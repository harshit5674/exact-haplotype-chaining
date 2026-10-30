"""Matched timing: the project's native engines on exactly the anchors Minichain chains.

For every read of a seed-export run (Minichain's own seed hits, split into Minichain's
components in its own order), this times the native engines on the compatibility objective
(Algorithm 2 with Minichain's Distance table, distance='minichain') or on Algorithm 2 as
transcribed (distance='true'). It verifies, in the same run, that the per-state table equals
Minichain's (FNV hash from a DP-dump run of the same reads and penalty) whenever --dump is given.

Modes, all from research/current/new/gap_native.cpp:
  per_haplotype    per-haplotype Fenwick trees, packed per-state output (every state + pointer)
  class_expanded   class engine, packed per-state output written in the call
  class_compact    class engine, per-anchor records (every state recoverable in O(1))
kernel_seconds is the native call (validation, sorting, sweep, output); packing of the
anchors into the native batch and graph indexing are excluded for every engine, as in
experiments/gap_engine_native/benchmark.py. Minichain's own times come from a separate,
sequential run of the instrumented binary with RECOMB_DP_TIME (run_minichain.py --time).
Run on an otherwise idle machine; one process; single-threaded kernels.
"""
from __future__ import annotations

import argparse
import gc
import json
import os
import platform
import random
import statistics
import subprocess
import sys
import time
from array import array
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mcdp  # noqa: E402
import gap_native as gn  # noqa: E402

MODES = ('per_haplotype', 'class_expanded', 'class_compact')


def call(G, batch, mode):
    gc.collect()
    if mode == 'class_compact':
        res = G.solve_compact(batch)
    else:
        res = G.solve_packed(batch, mode='per_haplotype' if mode == 'per_haplotype' else 'class')
    return res, res.stats['kernel_seconds']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--seeds', required=True, type=Path, help='seed-export run (seeds.tsv, index.tsv)')
    ap.add_argument('--dump', type=Path, help='DP-dump run of the same reads and penalty, for identity checks')
    ap.add_argument('--R', type=int, default=10000)
    ap.add_argument('--distance', choices=('minichain', 'true'), default='minichain')
    ap.add_argument('--repeat', type=int, default=3)
    ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--seed', type=int, default=20261006)
    ap.add_argument('--out', required=True, type=Path)
    args = ap.parse_args()
    meta_path = args.out.with_suffix('.meta.json')
    if args.out.exists() or meta_path.exists():
        raise SystemExit('refusing to overwrite')

    comps, header = mcdp.load_index(args.seeds / 'index.tsv')
    t = time.perf_counter()
    graphs = {c.cid: gn.GapGraph(mcdp.prepare(c, args.distance)) for c in comps}
    index_seconds = time.perf_counter() - t
    seeds, order = mcdp.read_seed_export(args.seeds / 'seeds.tsv')
    dump = mcdp.read_dp_dump(args.dump / 'dump.tsv') if args.dump else None
    if dump:
        assert all(m['recomb'] == args.R for m in dump.values()), 'DP dump penalty differs'
    split = mcdp.Splitter(comps)
    if args.limit:
        order = order[:args.limit]
    rng = random.Random(args.seed)
    order = list(order)
    rng.shuffle(order)
    sf = 200
    load_before = os.getloadavg()
    started = time.time()
    checked = 0
    with open(args.out, 'w') as out:
        for number, name in enumerate(order):
            parts = split.split(seeds[name]['hits'])
            row = dict(read=name, R=args.R, distance=args.distance, components=[])
            totals = {m: 0.0 for m in MODES}
            for c in comps:
                a0 = parts[c.cid]
                if not a0:
                    continue
                wt = (a0[0][4] - a0[0][3] + 1) * sf
                anchors = mcdp.to_gc_anchors(c, a0, wt)
                G = graphs[c.cid]
                batch = G.batch(anchors, args.R)
                times = {m: [] for m in MODES}
                digests = {}
                for rep in range(args.repeat):
                    modes = list(MODES)
                    rng.shuffle(modes)
                    for mode in modes:
                        res, k = call(G, batch, mode)
                        times[mode].append(k)
                        if rep == 0:
                            packed = res.expand() if mode == 'class_compact' else res
                            digests[mode] = packed.arrays_digest()
                            if mode == 'per_haplotype' and dump is not None:
                                m = dump[(name, c.cid)]
                                assert m['N'] == len(a0) and mcdp.input_hash(a0) == m['ha']
                                assert mcdp.fnv(packed.C) == m['hc'] and mcdp.fnv(packed.B) == m['hb'], (name, c.cid)
                                checked += 1
                            del packed
                        del res
                assert len(set(digests.values())) == 1, (name, c.cid, digests)
                med = {m: statistics.median(times[m]) for m in MODES}
                for m in MODES:
                    totals[m] += med[m]
                row['components'].append(dict(cid=c.cid, N=batch.n, K=batch.K, kernel_median=med, kernel_all=times))
            row['kernel_total'] = totals
            out.write(json.dumps(row) + '\n')
            if number % 100 == 0:
                print(number, name, {m: round(v * 1000, 3) for m, v in totals.items()}, flush=True)
    compiler = subprocess.run(['clang++', '--version'], capture_output=True, text=True).stdout.splitlines()[0]
    meta = dict(script='experiments/minichain_matched/bench_matched.py', argv=sys.argv[1:], python=sys.version,
                platform=platform.platform(), machine=platform.machine(), compiler=compiler,
                library=gn.build()._name, index_seconds=index_seconds, load_average_before=load_before,
                load_average_after=os.getloadavg(), wall_seconds=round(time.time() - started, 1), modes=MODES,
                repeat=args.repeat, seed=args.seed, identity_checked_pairs=checked)
    meta_path.write_text(json.dumps(meta, indent=1) + '\n')
    print(json.dumps(meta), flush=True)


if __name__ == '__main__':
    main()
