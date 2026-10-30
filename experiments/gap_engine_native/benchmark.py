"""Paired native benchmark: origin-class engine versus per-haplotype Fenwick baseline.

Objective: the gap-sensitive contract (docs/MATHEMATICAL_MODEL.md), Minichain weight
scaling (anchor weight x 200) and gamma 10000 or 20000. Inputs: the 40 primary MHC
reads on panels of H haplotypes (the 60 real paths plus seed-6060 mosaics).

Modes, all from research/current/new/gap_native.cpp:
  per_haplotype        query-major per-haplotype Fenwick trees, packed per-state output
  class_expanded       class engine with the per-state packed output written in the call
  class_compact        class engine, per-anchor records only; expand() is timed separately
  class_compact_warm   secondary: the same with the pair-range cache kept across reads
                       (primed by one untimed pass over the level's reads)

Rules: one process, single-threaded kernels, sequential runs; graph indexing,
compilation and verification are outside the timings; gc.collect() before every run.
Every read runs `--repeat` repetitions of all modes in a freshly shuffled order (paired,
randomized). kernel_seconds is the native call (validation, sorting, sweep and output
export); packed_seconds adds Python input packing and output allocation and copying.
Exactness on every read and level: the class compact records equal the Python
specification's records (gap_chaining.solve) and every mode's C, B, R arrays are equal.
Pointers are checked separately by verify_native.py.
"""
from __future__ import annotations
import argparse
import gc
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
sys.path.insert(0, str(HERE))
import common as cm  # noqa: E402
import gap_chaining as gcs  # noqa: E402
import gap_native as gn  # noqa: E402

MODES = ('per_haplotype', 'class_expanded', 'class_compact', 'class_compact_warm')
STAGES = ('setup_seconds', 'init_seconds', 'sweep_seconds', 'export_seconds')
COUNTERS = ('merges', 'contested_merges', 'pairs', 'pairs_same_vertex', 'pairs_loose', 'pairs_range',
            'pairs_sliced', 'sliced_thresholds', 'sliced_minima', 'pair_ranges', 'range_labels',
            'route_planes', 'cut_records', 'prefix_records', 'remote_groups', 'remote_candidates',
            'remote_sliced', 'published_records', 'output_records', 'max_node_records',
            'ph_label_steps', 'ph_remote_queries')


def run_mode(G, anchors, gamma, mode):
    """One timed call; anchors are packed inside it, so packed_seconds includes packing."""
    gc.collect()
    if mode == 'per_haplotype':
        res = G.solve_packed(anchors, gamma, mode='per_haplotype')
        return res, dict(kernel=res.stats['kernel_seconds'], packed=res.stats['packed_seconds'], expansion=None)
    if mode == 'class_expanded':
        res = G.solve_packed(anchors, gamma, mode='class')
        return res, dict(kernel=res.stats['kernel_seconds'], packed=res.stats['packed_seconds'], expansion=None)
    comp = G.solve_compact(anchors, gamma, warm=(mode == 'class_compact_warm'))
    res = comp.expand()
    return (comp, res), dict(kernel=comp.stats['kernel_seconds'], packed=comp.stats['packed_seconds'],
                             expansion=res.stats['expansion_seconds'])


def level(G, p, H, gamma, cases, args, rng, out):
    batches = [G.batch(batch_anchors(c), gamma) for c in cases]
    G.warm_cache_reset()
    t = time.perf_counter()
    for b in batches:  # prime the warm cache (untimed)
        G.solve_compact(b, warm=True)
    warm_prime = time.perf_counter() - t
    warm_info = G.warm_cache_info()
    order = list(range(len(cases)))
    rng.shuffle(order)
    for number, k in enumerate(order):
        case, batch = cases[k], batches[k]
        anchors = batch_anchors(case)
        spec = gcs.solve(p, anchors, gamma)
        times = {m: [] for m in MODES}
        digests, meta = {}, {}
        for rep in range(args.repeat):
            modes = list(MODES)
            rng.shuffle(modes)
            for mode in modes:
                res, tm = run_mode(G, anchors, gamma, mode)
                times[mode].append(tm)
                if rep == 0:
                    if isinstance(res, tuple):
                        comp, res = res
                        assert comp.records() == spec.records, (H, gamma, case['id'], mode)
                        assert list(comp.B) == spec.B and comp.R_values() == spec.R, (H, gamma, case['id'], mode)
                        meta[mode] = {k2: v for k2, v in comp.stats.items() if isinstance(v, (int, float))}
                    else:
                        meta[mode] = {k2: v for k2, v in res.stats.items() if isinstance(v, (int, float))}
                    digests[mode] = res.arrays_digest()
                del res
        assert len(set(digests.values())) == 1, (H, gamma, case['id'])
        row = dict(H=H, gamma=gamma, id=case['id'], N=batch.n, K=batch.K, W=G.W, digest=digests['per_haplotype'],
                   times=times, stats=meta, warm_cache_after_prime=warm_info, warm_prime_seconds=warm_prime)
        for field in ('kernel', 'packed', 'expansion'):
            row[field + '_median'] = {m: statistics.median(x[field] for x in times[m]) for m in MODES
                                      if times[m][0][field] is not None}
        out.write(json.dumps(row) + '\n')
        out.flush()
        km = row['kernel_median']
        print(f'H={H} gamma={gamma} {number + 1}/{len(cases)} {case["id"]} N={batch.n}',
              {m: round(v * 1000, 3) for m, v in km.items()},
              'speedup compact %.2f expanded %.2f' % (km['per_haplotype'] / km['class_compact'],
                                                      km['per_haplotype'] / km['class_expanded']), flush=True)


_ANCHORS = {}


def batch_anchors(case):
    if case['id'] not in _ANCHORS:
        _ANCHORS[case['id']] = cm.anchors(case)
    return _ANCHORS[case['id']]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, required=True, help='JSONL rows; a manifest is written beside it')
    ap.add_argument('--H', type=lambda s: [int(x) for x in s.split(',')], default=list(cm.H_LEVELS))
    ap.add_argument('--gamma', type=lambda s: [int(x) for x in s.split(',')], default=[10000, 20000])
    ap.add_argument('--repeat', type=int, default=5)
    ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--seed', type=int, default=8080)
    args = ap.parse_args()
    meta_path = args.out.with_suffix('.meta.json')
    if args.out.exists() or meta_path.exists():
        raise SystemExit('refusing to overwrite')
    lib_path = gn.build()._name
    cases = cm.primary_cases()
    if args.limit:
        cases = cases[:args.limit]
    rng = random.Random(args.seed)
    load_before = os.getloadavg()
    started = time.time()
    with args.out.open('w') as out:
        for H in args.H:
            t = time.perf_counter()
            p = cm.prepared(H)
            G = gn.GapGraph(p)
            prep = dict(python_prepared_seconds=time.perf_counter() - t, native_index_seconds=G.native_index_seconds,
                        packing_index_seconds=G.index_seconds, remote_groups=G.remote_groups,
                        remote_entries=G.remote_entries, W=G.W, b=G.b)
            print('H', H, prep, flush=True)
            for gamma in args.gamma:
                level(G, p, H, gamma, cases, args, rng, out)
            del G
    compiler = subprocess.run(['clang++', '--version'], capture_output=True, text=True).stdout.splitlines()[0]
    meta = dict(script='experiments/gap_engine_native/benchmark.py', argv=sys.argv[1:], python=sys.version,
                platform=platform.platform(), machine=platform.machine(), compiler=compiler,
                compile_flags=['-std=c++20', '-O3', '-shared', '-fPIC'], library=lib_path,
                load_average_before=load_before, load_average_after=os.getloadavg(),
                wall_seconds=round(time.time() - started, 1), modes=MODES, repeat=args.repeat, seed=args.seed)
    meta_path.write_text(json.dumps(meta, indent=1) + '\n')
    print(json.dumps(meta), flush=True)


if __name__ == '__main__':
    main()
