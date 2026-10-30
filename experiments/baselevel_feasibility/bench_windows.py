"""Paired kernel timing of the engines on the base-level window panels (2026-10-20; pre-specified in REPORT.md).

Same protocol and statistics as the paper's 100-read benchmark (experiments/benchmark_100reads):
  - validate: every instance, engines per_haplotype, class (expanded), auto (guarded) and compact; values
    equal across engines; output signatures (values and pointers, or compact records) saved;
  - benchmark: one untimed warmup and --reps timed repetitions per instance in a fixed shuffled instance
    order, engine order shuffled per repetition, cold pair cache; every output's signature is checked after
    its timer stops; native kernel and Python packed-call times are recorded. A panel starts only when the
    host is idle (1-minute load below --max-load for --settle minutes) and no other project timing job runs;
  - summarize: per read, the median over repetitions of the summed two-instance time; median of paired
    per-read ratios with a 2,000-resample whole-read bootstrap, totals and wins, using summarize_reads of
    experiments/benchmark_100reads/summarize.py unchanged.

  python3.13 -B experiments/baselevel_feasibility/bench_windows.py --phase validate --out experiments/baselevel_feasibility/timing/validation
  python3.13 -B experiments/baselevel_feasibility/bench_windows.py --phase benchmark --validation .../validation --out .../timing/run1
  python3.13 -B experiments/baselevel_feasibility/bench_windows.py --phase summarize --validation .../validation --run .../run1 --out .../summary.json
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import random
import statistics as st
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'research/current/new'))
import gap_chaining as gc  # noqa: E402
import gap_native as gn  # noqa: E402

GAMMA, SCALE = 10000, 200
MODES = ('per_haplotype', 'class', 'auto', 'compact')
WINDOWS = ('w31250k_250kb', 'w31700k_250kb', 'w32375k_100kb')   # class I, class III, class II (REPORT.md)
SOURCES = [Path(__file__)] + [ROOT / 'research/current/new' / f for f in ('gap_chaining.py', 'gap_native.py', 'gap_native.cpp')]
OTHER_TIMING = ('benchmark_100reads/run.py', 'run_when_idle.py', 'minichain_integrated/run.py --phase timing',
                'minichain-engine/minichain')


def file_hash(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def digest(z, compact=False):
    fields = ('off', 'value', 'origin', 'kind', 'wj', 'wg', 'ww', 'masks', 'B', 'R') if compact else ('C', 'B', 'R', 'back')
    h = hashlib.sha256()
    for f in fields:
        a = getattr(z, f)
        h.update(len(a).to_bytes(8, 'little'))
        h.update(memoryview(a))
    return h.hexdigest()


def solve(G, batch, mode):
    return G.solve_compact(batch, warm=False) if mode == 'compact' else G.solve_packed(batch, mode=mode, warm=False)


def load_panel(window: str):
    panel = HERE / 'windows' / window / 'panel'
    graph = json.loads((panel / 'graph.json').read_text())
    cases = [json.loads(line) for line in (panel / 'instances.jsonl').open()]
    p = gc.Graph(graph['lengths'], [tuple(e) for e in graph['edges']], graph['paths']).prepare()
    hashes = {f: file_hash(panel / f) for f in ('graph.json', 'instances.jsonl')}
    return panel, gn.GapGraph(p), cases, hashes


def host():
    swap = subprocess.run(['sysctl', '-n', 'vm.swapusage'], capture_output=True, text=True).stdout.strip()
    return dict(time=datetime.now(timezone.utc).isoformat(), load=os.getloadavg(), swap=swap)


def other_timing_running() -> bool:
    out = subprocess.run(['ps', '-Ao', 'args'], capture_output=True, text=True).stdout
    return any(pat in line for line in out.splitlines() for pat in OTHER_TIMING if 'bench_windows.py' not in line)


def wait_quiet(max_load: float, settle: int):
    quiet = 0
    while quiet < settle:
        quiet = quiet + 1 if (os.getloadavg()[0] < max_load and not other_timing_running()) else 0
        if quiet < settle:
            time.sleep(60)


def validate(args):
    args.out.mkdir(parents=True)
    manifest = dict(phase='validate', started=datetime.now(timezone.utc).isoformat(), platform=platform.platform(),
                    source_hashes={str(p.relative_to(ROOT)): file_hash(p) for p in SOURCES},
                    library=gn.build()._name, gamma=GAMMA, weight_scale=SCALE, auto_factor=gn.AUTO_FACTOR, panels=[])
    for w in WINDOWS:
        panel, G, cases, hashes = load_panel(w)
        rows = []
        for c in cases:
            A = [gc.Anchor(*x[:5], x[5] * SCALE) for x in c['anchors']]
            b = G.batch(A, GAMMA)
            sig, stats, values = {}, {}, set()
            for m in MODES:
                z = solve(G, b, m)
                values.add((z.expand() if m == 'compact' else z).arrays_digest())
                sig[m] = digest(z, compact=m == 'compact')
                stats[m] = z.stats
            assert len(values) == 1, c['id']
            assert stats['auto']['auto_fallback_available'] and stats['auto']['auto_work'] <= stats['auto']['auto_budget']
            rows.append(dict(id=c['id'], read=c['read'], N=b.n, K=b.K, signatures=sig,
                             fallback=bool(stats['auto']['fallback']), compact_records=stats['compact']['compact_records'],
                             class_stats={k: v for k, v in stats['class'].items() if isinstance(v, (int, float))}))
        (args.out / f'{w}.jsonl').write_text(''.join(json.dumps(r) + '\n' for r in rows))
        manifest['panels'].append(dict(window=w, input_hashes=hashes, instances=len(rows), states=sum(r['K'] for r in rows)))
        print(w, 'validated', len(rows), flush=True)
    manifest['status'] = 'passed'
    (args.out / 'manifest.json').write_text(json.dumps(manifest, indent=1) + '\n')


def benchmark(args):
    valid = json.loads((args.validation / 'manifest.json').read_text())
    assert valid['status'] == 'passed'
    hashes = {str(p.relative_to(ROOT)): file_hash(p) for p in SOURCES}
    assert hashes == valid['source_hashes'], 'sources changed since validation'
    args.out.mkdir(parents=True)
    manifest = dict(phase='benchmark', started=datetime.now(timezone.utc).isoformat(), platform=platform.platform(),
                    source_hashes=hashes, library=gn.build()._name, repetitions=args.reps, pair_cache='cold',
                    max_load=args.max_load, settle_minutes=args.settle, panels=[], status='running')
    try:
        for w in WINDOWS:
            expected = {r['id']: r for r in map(json.loads, (args.validation / f'{w}.jsonl').read_text().splitlines())}
            panel, G, cases, hash_now = load_panel(w)
            vpanel = next(p for p in valid['panels'] if p['window'] == w)
            assert vpanel['input_hashes'] == hash_now
            wait_quiet(args.max_load, args.settle)
            before = host()
            rng = random.Random(20261020 + sum(w.encode()))
            order = list(cases)
            rng.shuffle(order)
            rows = []
            for c in order:
                A = [gc.Anchor(*x[:5], x[5] * SCALE) for x in c['anchors']]
                b = G.batch(A, GAMMA)
                ref = expected[c['id']]
                times = {m: [] for m in MODES}
                packed = {m: [] for m in MODES}
                for rep in range(-1, args.reps):
                    for m in rng.sample(list(MODES), len(MODES)):
                        z = solve(G, b, m)
                        assert digest(z, compact=m == 'compact') == ref['signatures'][m], (c['id'], m, rep)
                        if rep >= 0:
                            times[m].append(z.stats['kernel_seconds'])
                            packed[m].append(z.stats['packed_seconds'])
                        del z
                rows.append(dict(id=c['id'], read=c['read'], N=b.n, K=b.K, kernel_seconds=times, packed_seconds=packed,
                                 guard=dict(fallback=ref['fallback']), compact_records=ref['compact_records'],
                                 prior_primary=False))
            (args.out / f'{w}.jsonl').write_text(''.join(json.dumps(r) + '\n' for r in rows))
            manifest['panels'].append(dict(window=w, instances=len(rows), host_before=before, host_after=host()))
            (args.out / 'manifest.json').write_text(json.dumps(manifest, indent=1) + '\n')
            print(w, 'timed', flush=True)
        manifest['status'] = 'passed'
    finally:
        manifest['finished'] = datetime.now(timezone.utc).isoformat()
        (args.out / 'manifest.json').write_text(json.dumps(manifest, indent=1) + '\n')


def summarize(args):
    spec = importlib.util.spec_from_file_location('b100sum', ROOT / 'experiments/benchmark_100reads/summarize.py')
    b100 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(b100)
    man = json.loads((args.run / 'manifest.json').read_text())
    assert man['status'] == 'passed'
    out = dict(run=str(args.run), validation=str(args.validation), repetitions=man['repetitions'], windows=[])
    for p in man['panels']:
        w = p['window']
        rows = [json.loads(line) for line in (args.run / f'{w}.jsonl').read_text().splitlines()]
        groups = defaultdict(list)
        for r in rows:
            groups[r['read']].append(r)
        reads = []
        for name, comps in groups.items():
            r = dict(read=name, fallback=any(c['guard']['fallback'] for c in comps),
                     N=sum(c['N'] for c in comps), K=sum(c['K'] for c in comps),
                     records=sum(c['compact_records'] for c in comps))
            for boundary in ('kernel_seconds', 'packed_seconds'):
                r[boundary] = {m: st.median(sum(c[boundary][m][j] for c in comps) for j in range(man['repetitions']))
                               for m in MODES}
            reads.append(r)
        seed = 20261020 + sum(w.encode())
        out['windows'].append(dict(
            window=w, reads=len(reads), instances=len(rows), states=sum(r['K'] for r in reads),
            anchors=sum(r['N'] for r in reads), fallback_instances=sum(c['guard']['fallback'] for c in rows),
            fallback_reads=sum(r['fallback'] for r in reads),
            records_per_anchor=st.median(r['records'] / r['N'] for r in reads),
            states_per_anchor=st.median(r['K'] / r['N'] for r in reads),
            perhap_ms=st.median(r['kernel_seconds']['per_haplotype'] for r in reads) * 1e3,
            guarded_ms=st.median(r['kernel_seconds']['auto'] for r in reads) * 1e3,
            kernel=b100.summarize_reads(reads, 'kernel_seconds', seed),
            packed=b100.summarize_reads(reads, 'packed_seconds', seed),
            host_before=p['host_before'], host_after=p['host_after']))
    args.out.write_text(json.dumps(out, indent=1) + '\n')
    for x in out['windows']:
        a = x['kernel']['auto']
        print(x['window'], f"guarded {a['median_speedup']:.2f} {a['median_ci95']} total {a['total_speedup']:.2f} "
              f"wins {a['wins']}/{x['reads']} compact {x['kernel']['compact']['median_speedup']:.2f} fallbacks {x['fallback_instances']}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--phase', choices=('validate', 'benchmark', 'summarize'), required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--validation', type=Path)
    ap.add_argument('--run', type=Path)
    ap.add_argument('--reps', type=int, default=5)
    ap.add_argument('--max-load', type=float, default=1.5)
    ap.add_argument('--settle', type=int, default=5)
    args = ap.parse_args()
    {'validate': validate, 'benchmark': benchmark, 'summarize': summarize}[args.phase](args)


if __name__ == '__main__':
    main()
