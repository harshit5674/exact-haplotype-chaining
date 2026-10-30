"""Minichain with the project's engines as a drop-in chainer: correctness and end-to-end timing (2026-10-20).

Kernel timings are not mapper timings, so this also times the whole mapper.
vendor/minichain-engine is a copy of experiments/minichain_matched/vendor/minichain-dp (pinned Minichain
d21e694 with the seed-export and DP-dump instrumentation) plus a drop-in chainer (graphUtils.cpp,
namespace recomb_engine; see build.log for the build). With RECOMB_ENGINE unset the binary runs Minichain's
own DP. With RECOMB_ENGINE=ph|class|auto the chaining DP is computed by gap_ph_solve, gap_class_packed or
gap_auto_packed from the production engine library, given Minichain's own graph, topological order,
last2reach and Distance tables; Minichain's backtracking and output then run unchanged. Anchors are passed
in reverse order (RECOMB_ENGINE_REVERSE=1), which aligns the engines' tie-breaking among equal keys with
Minichain's range query; it does not change any value.

Phase 'correctness': every variant on every dataset, with RECOMB_DP_DUMP. Checks: (1) the engine-off run's
GAF equals the pinned binary's saved GAF; (2) every engine reproduces the engine-off DP hashes (all state
values, B = C_m, max(w, R) = C_p, best score, valid states) for every (read, component); (3) GAF agreement
per read, with every difference described. Phase 'timing': shuffled paired end-to-end repetitions, run only
when the host is idle (load and swap recorded), without dumps; RECOMB_DP_TIME records each read's chaining
time and Minichain's log its total wall time and peak RSS.

  python3.13 -B experiments/minichain_integrated/run.py --phase correctness --out experiments/minichain_integrated/results/correctness
  python3.13 -B experiments/minichain_integrated/run.py --phase timing --reps 3 --out experiments/minichain_integrated/results/timing
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import re
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'research/current/new'))
import gap_native as gn  # noqa: E402

BINARY = HERE / 'vendor/minichain-engine/minichain'
MATCHED = ROOT / 'experiments/minichain_matched/data'
DATASETS = {   # name: (graph, reads, saved GAF of the pinned binary on the same input)
    'mhc2104': ('experiments/mhc_pilot/data/MHC-CHM13.0.gfa.gz', 'experiments/mhc_pilot/data/MHC_CHM13_PacBio_filt.fq.gz',
                MATCHED / 'dp_R10000/minichain.gaf'),
    **{f'hprc{h}': (f'experiments/hprc_mhc_panel/data/hprc{h}/graph.gfa.gz', 'experiments/mhc_pilot/data/sample100.fq',
                    MATCHED / f'time_hprc{h}_R10000_rep1/minichain.gaf') for h in (60, 120, 240, 463)},
}
VARIANTS = {'off': {}, 'ph': {'RECOMB_ENGINE': 'ph'}, 'class': {'RECOMB_ENGINE': 'class'}, 'auto': {'RECOMB_ENGINE': 'auto'}}


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def host_state():
    swap = subprocess.run(['sysctl', '-n', 'vm.swapusage'], capture_output=True, text=True).stdout.strip()
    return dict(load=os.getloadavg(), swap=swap)


def run_one(dataset: str, variant: str, out: Path, lib: Path, dump: bool, time_dump: bool) -> dict:
    graph, reads, _ = DATASETS[dataset]
    out.mkdir(parents=True)
    env = dict(os.environ, OMP_NUM_THREADS='1')
    for k in [k for k in env if k.startswith('RECOMB_')]:
        del env[k]
    env.update(VARIANTS[variant])
    if variant != 'off':
        env.update(RECOMB_ENGINE_LIB=str(lib), RECOMB_ENGINE_REVERSE='1')
    if dump:
        env['RECOMB_DP_DUMP'] = str(out / 'dump.tsv')
    if time_dump:
        env['RECOMB_DP_TIME'] = str(out / 'time.tsv')
    cmd = [str(BINARY), '-x', 'lr', '-t', '1', '-R10000', graph, reads]
    before = host_state()
    t0 = time.perf_counter()
    with (out / 'minichain.gaf').open('w') as gaf, (out / 'minichain.log').open('w') as log:
        rc = subprocess.run(cmd, cwd=ROOT, env=env, stdout=gaf, stderr=log).returncode
    wall = time.perf_counter() - t0
    log = (out / 'minichain.log').read_text()
    assert rc == 0, log[-2000:]
    m = re.search(r'Real time: ([0-9.]+) sec; CPU: ([0-9.]+) sec; Peak RSS: ([0-9.]+) GB', log)
    e = re.search(r'built (\d+) component graphs in ([0-9.]+) s', log)
    f = re.search(r'(\d+) engine calls, (\d+) fallbacks, (\d+) states', log)
    row = dict(dataset=dataset, variant=variant, command=cmd, switches={k: v for k, v in env.items() if k.startswith('RECOMB_')},
               wall_seconds=wall, minichain_real_seconds=float(m.group(1)), cpu_seconds=float(m.group(2)),
               peak_rss_gb=float(m.group(3)), engine_build_seconds=float(e.group(2)) if e else None,
               engine_calls=int(f.group(1)) if f else None, fallbacks=int(f.group(2)) if f else None,
               engine_states=int(f.group(3)) if f else None, host_before=before, host_after=host_state(),
               finished=datetime.now(timezone.utc).isoformat())
    (out / 'run.json').write_text(json.dumps(row, indent=1) + '\n')
    return row


def read_dump(p: Path):
    out, key = {}, None
    for line in p.open():
        a = line.rstrip('\n').split('\t')
        if a[0] == 'Q':
            key = (a[1], int(a[2]))
            assert key not in out
            out[key] = dict(N=int(a[3]), K=int(a[4]))
        elif a[0] == 'H':
            out[key].update(ha=a[1], hc=a[2], hb=a[3], hp=a[4], best=int(a[5]), nvalid=int(a[6]))
    return out


def gaf_by_read(p: Path):
    rows = {}
    for line in p.open():
        a = line.rstrip('\n').split('\t')
        rows.setdefault(a[0], []).append(a)
    return rows


def compare(base: Path, other: Path) -> dict:
    d0, d1 = read_dump(base / 'dump.tsv'), read_dump(other / 'dump.tsv')
    assert set(d0) == set(d1), 'different (read, component) sets'
    dp_mismatch = [k for k in d0 if d0[k] != d1[k]]
    g0, g1 = gaf_by_read(base / 'minichain.gaf'), gaf_by_read(other / 'minichain.gaf')
    assert set(g0) == set(g1)
    same = [r for r in g0 if g0[r] == g1[r]]
    diff = [r for r in g0 if g0[r] != g1[r]]
    target_same = [r for r in diff if [x[5:9] for x in g0[r]] == [x[5:9] for x in g1[r]]]
    fields = {}
    for r in diff:
        for x0, x1 in zip(g0[r], g1[r]):
            for i, (u, v) in enumerate(zip(x0, x1)):
                if u != v:
                    name = (x0[i].split(':')[0] if ':' in x0[i] and i >= 12 else f'col{i + 1}')
                    fields[name] = fields.get(name, 0) + 1
    return dict(components=len(d0), dp_hash_mismatches=len(dp_mismatch), dp_mismatch_examples=dp_mismatch[:5],
                reads=len(g0), gaf_identical_reads=len(same), gaf_different_reads=len(diff),
                different_reads_with_same_path_and_target_interval=len(target_same), differing_fields=fields)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--phase', choices=('correctness', 'timing'), required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--reps', type=int, default=3)
    ap.add_argument('--datasets', default=','.join(DATASETS))
    ap.add_argument('--max-load', type=float, default=1.5, help='timing: wait until the 1-minute load is below this')
    args = ap.parse_args()
    args.out.mkdir(parents=True)
    lib = Path(gn.build()._name)
    datasets = args.datasets.split(',')
    manifest = dict(phase=args.phase, started=datetime.now(timezone.utc).isoformat(), platform=platform.platform(),
                    binary=str(BINARY.relative_to(ROOT)), binary_sha256=sha(BINARY),
                    engine_library=str(lib.relative_to(ROOT)), engine_library_sha256=sha(lib),
                    engine_source_sha256=sha(ROOT / 'research/current/new/gap_native.cpp'),
                    patch_source_sha256=sha(HERE / 'vendor/minichain-engine/graphUtils.cpp'),
                    datasets={d: dict(graph=DATASETS[d][0], reads=DATASETS[d][1], graph_sha256=sha(ROOT / DATASETS[d][0]),
                                      reads_sha256=sha(ROOT / DATASETS[d][1])) for d in datasets}, runs=[], status='running')
    path = args.out / 'manifest.json'
    try:
        if args.phase == 'correctness':
            manifest['comparisons'] = {}
            for d in datasets:
                for v in VARIANTS:
                    manifest['runs'].append(run_one(d, v, args.out / d / v, lib, dump=True, time_dump=False))
                    path.write_text(json.dumps(manifest, indent=1) + '\n')
                pinned = DATASETS[d][2]
                manifest['comparisons'][d] = dict(
                    off_gaf_equals_pinned=(args.out / d / 'off/minichain.gaf').read_bytes() == pinned.read_bytes(),
                    pinned_gaf=str(pinned.relative_to(ROOT)),
                    **{v: compare(args.out / d / 'off', args.out / d / v) for v in VARIANTS if v != 'off'})
                path.write_text(json.dumps(manifest, indent=1) + '\n')
                print(d, json.dumps(manifest['comparisons'][d]), flush=True)
        else:
            rng = random.Random(20261020)
            for rep in range(args.reps):
                for d in datasets:
                    order = list(VARIANTS)
                    rng.shuffle(order)
                    for v in order:
                        while os.getloadavg()[0] > args.max_load:   # run only on an idle host
                            time.sleep(60)
                        manifest['runs'].append(dict(rep=rep, **run_one(d, v, args.out / f'rep{rep}' / d / v, lib,
                                                                       dump=False, time_dump=True)))
                        path.write_text(json.dumps(manifest, indent=1) + '\n')
        manifest['status'] = 'passed'
    except BaseException as e:
        manifest.update(status='failed', error=repr(e))
        raise
    finally:
        manifest['finished'] = datetime.now(timezone.utc).isoformat()
        path.write_text(json.dumps(manifest, indent=1) + '\n')


if __name__ == '__main__':
    main()
