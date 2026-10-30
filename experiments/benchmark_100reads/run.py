"""Expand the saved panel benchmark to all 100 reads, without changing solvers.

Run validation before timing. Each phase requires a fresh output directory.
Timings use a reused prepared graph/batch and cold pair caches. The native
kernel and Python packed-call boundaries are recorded separately. Compact
output is unguarded and does not include expansion in either timing boundary.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import statistics
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(ROOT/'research/current/new'), str(ROOT/'experiments/gap_fallback')]
import auto_real as data
import gap_chaining as gc
import gap_native as gn
import gap_validation as gv

MODES = ('per_haplotype', 'class', 'auto', 'compact')
SOURCES = [Path(__file__), ROOT/'experiments/gap_fallback/auto_real.py'] + [
    ROOT/'research/current/new'/f for f in (
        'gap_chaining.py', 'gap_native.py', 'gap_native.cpp',
        'gap_validation.py', 'gap_validation.cpp')]


def file_hash(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def digest(z, compact=False, pointers=False):
    fields = ('off', 'value', 'origin', 'kind', 'wj', 'wg', 'ww', 'masks', 'B', 'R') if compact else (
        ('C', 'B', 'R', 'back') if pointers else ('C', 'B', 'R'))
    h = hashlib.sha256()
    for field in fields:
        a = getattr(z, field)
        h.update(len(a).to_bytes(8, 'little'))
        h.update(memoryview(a))
    return h.hexdigest()


def save(p, value):
    p.write_text(json.dumps(value, indent=2)+'\n')


def read_id(c):
    return c['id'].rsplit(':', 1)[0]


def solve(G, batch, mode):
    return G.solve_compact(batch, warm=False) if mode == 'compact' else G.solve_packed(batch, mode=mode, warm=False)


def assert_guard(z):
    assert z.stats['auto_fallback_available']
    assert z.stats['auto_work'] <= z.stats['auto_budget']


def direct_checks(p, G, validator, cases, old_ids):
    # These are drawn only from the previously untimed 80 reads. All their
    # full components are independently validated below as well.
    fresh = sorted((c for c in cases if c['id'] not in old_ids and c['anchors']),
                   key=lambda c: len(c['anchors']))
    selected = [fresh[round(i*(len(fresh)-1)/3)] for i in range(4)]
    rows = []
    for c in selected:
        full = [gc.Anchor(*x[:5], x[5]*data.SCALE) for x in c['anchors']]
        n = min(32, len(full))
        indices = sorted({round(i*(len(full)-1)/(n-1)) for i in range(n)}) if n > 1 else [0]
        A = [full[i] for i in indices]
        for gamma in (0, 1, 10000, 20000):
            expected = gc.direct_oracle(p, A, gamma)
            v = validator.batch(A, gamma)
            z = G.solve_packed(A, gamma, mode='auto')
            v.check(z); v.close()
            C, ptr, B, R = z.as_dicts(A)
            assert (C, B, R) == expected[:3], c['id']
            gc.verify_pointers(p, A, gamma, C, B, ptr)
            rows.append(dict(id=c['id'], anchor_indices=indices, gamma=gamma, states=len(z.C)))
    return rows


def run_panel(args, locus, name, directory):
    started = time.perf_counter()
    graph, cases = data.load(directory)
    counts = Counter(read_id(c) for c in cases)
    assert len(cases) == 200 and len(counts) == 100 and set(counts.values()) == {2}
    old_ids = {c['id'] for c in data.primary(cases, locus)}
    assert len(old_ids) == 40
    input_hashes = {f: file_hash(directory/f) for f in ('graph.json', 'instances.jsonl')}
    t = time.perf_counter()
    p = gc.Graph(graph['lengths'], [tuple(e) for e in graph['edges']], graph['paths']).prepare()
    python_prepare_seconds = time.perf_counter()-t
    G = gn.GapGraph(p)
    assert G.warm_cache_info()['pairs'] == 0
    panel = dict(locus=locus, panel=name, H=p.H, V=len(graph['lengths']),
                 input_directory=str(directory.relative_to(ROOT)), input_hashes=input_hashes,
                 reads=100, components=200, python_prepare_seconds=python_prepare_seconds,
                 native_index_seconds=G.native_index_seconds, adapter_index_seconds=G.index_seconds,
                 load_before=os.getloadavg(), status='running')
    key = locus+'-'+name
    rows_path = args.out/(key+'.jsonl')
    validator = gv.Validator(graph) if args.phase == 'validate' else None
    expected = None
    if args.phase == 'benchmark':
        validated_panel = json.loads((args.validation/(key+'.json')).read_text())
        assert validated_panel['status'] == 'passed' and validated_panel['input_hashes'] == input_hashes
        expected = {r['id']: r for r in map(json.loads, (args.validation/(key+'.jsonl')).read_text().splitlines())}
        assert set(expected) == {c['id'] for c in cases}
    # Fixed shuffled component order; identical inputs and mode pairs within each repetition.
    rng = random.Random(20261017+sum(key.encode()))
    order = list(cases)
    rng.shuffle(order)
    rows = []
    with rows_path.open('x') as f:
        for k, c in enumerate(order):
            A = [gc.Anchor(*x[:5], x[5]*data.SCALE) for x in c['anchors']]
            t = time.perf_counter(); batch = G.batch(A, data.GAMMA)
            packing_seconds = time.perf_counter()-t
            row = dict(id=c['id'], read=read_id(c), prior_primary=c['id'] in old_ids,
                       N=batch.n, K=batch.K, anchor_packing_seconds=packing_seconds)
            if args.phase == 'validate':
                independent = validator.batch(A, data.GAMMA)
                values, signatures, stats = set(), {}, {}
                for mode in MODES:
                    z = solve(G, batch, mode)
                    expanded = z.expand() if mode == 'compact' else z
                    independent.check(expanded)
                    values.add(digest(expanded))
                    signatures[mode] = digest(z, compact=mode == 'compact', pointers=True)
                    stats[mode] = z.stats
                    if mode == 'auto': assert_guard(z)
                    del expanded, z
                independent.close()
                assert len(values) == 1, c['id']
                row.update(digest=values.pop(), signatures=signatures, stats=stats)
            else:
                ref = expected[c['id']]
                assert (row['N'], row['K']) == (ref['N'], ref['K'])
                times = {m: [] for m in MODES}
                packed = {m: [] for m in MODES}
                # Untimed warmup plus five paired measurements. Pair-cache reuse
                # is disabled; every output is checked after its timer stops.
                for rep in range(-1, args.reps):
                    for mode in rng.sample(list(MODES), len(MODES)):
                        z = solve(G, batch, mode)
                        assert digest(z, compact=mode == 'compact', pointers=True) == ref['signatures'][mode], (c['id'], mode, rep)
                        if mode == 'auto': assert_guard(z)
                        if rep >= 0:
                            times[mode].append(z.stats['kernel_seconds'])
                            packed[mode].append(z.stats['packed_seconds'])
                        del z
                row.update(kernel_seconds=times, packed_seconds=packed,
                           guard=ref['stats']['auto'], class_stats=ref['stats']['class'],
                           compact_records=ref['stats']['compact']['compact_records'])
            assert G.warm_cache_info()['pairs'] == 0
            f.write(json.dumps(row)+'\n'); f.flush(); rows.append(row)
            if (k+1) % 25 == 0:
                print(f'{args.phase} {key}: {k+1}/200 components', flush=True)
    if validator:
        panel['direct_checks'] = direct_checks(p, G, validator, cases, old_ids)
        validator.close()
    panel.update(status='passed', states=sum(r['K'] for r in rows), anchors=sum(r['N'] for r in rows),
                 seconds=time.perf_counter()-started, load_after=os.getloadavg())
    save(args.out/(key+'.json'), panel)
    print(json.dumps(panel if not validator else {k:v for k,v in panel.items() if k != 'direct_checks'}), flush=True)
    return panel


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--phase', choices=('validate', 'benchmark'), required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--validation', type=Path)
    ap.add_argument('--reps', type=int, default=5)
    args = ap.parse_args()
    assert args.reps >= 3
    args.out.mkdir()  # Refuse to overwrite any prior run, including a partial one.
    hashes = {str(p.relative_to(ROOT)): file_hash(p) for p in SOURCES}
    if args.phase == 'benchmark':
        assert args.validation
        valid = json.loads((args.validation/'manifest.json').read_text())
        assert valid['status'] == 'passed' and valid['source_hashes'] == hashes
    started = time.perf_counter()
    manifest = dict(status='running', phase=args.phase, command=sys.argv,
                    started_utc=datetime.now(timezone.utc).isoformat(),
                    platform=platform.platform(), python=sys.version, assertions=__debug__,
                    cpu=platform.processor(), threads=1, source_hashes=hashes,
                    compiler=subprocess.run(['clang++', '--version'], capture_output=True, text=True, check=True).stdout,
                    library=gn.build()._name, gamma=data.GAMMA, weight_scale=data.SCALE,
                    auto_factor=gn.AUTO_FACTOR, repetitions=args.reps, pair_cache='cold',
                    timing_boundary='kernel: native setup/init/sweep/export; packed: adapter allocations and call; graph preparation, anchor packing, verification and compact expansion excluded',
                    panels=[])
    try:
        for locus, name, directory in data.panels():
            manifest['panels'].append(run_panel(args, locus, name, directory))
            save(args.out/'manifest.json', manifest)
        assert len(manifest['panels']) == 14
        assert hashes == {str(p.relative_to(ROOT)): file_hash(p) for p in SOURCES}
        manifest['status'] = 'passed'
    except BaseException as e:
        manifest['status'] = 'failed'; manifest['error'] = repr(e)
        raise
    finally:
        manifest['seconds'] = time.perf_counter()-started
        save(args.out/'manifest.json', manifest)


if __name__ == '__main__':
    main()
