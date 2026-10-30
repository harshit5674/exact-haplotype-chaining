"""Calibrate the per-haplotype baseline's dense/sparse crossover (`dense_cut`).

An anchor with d labels uses dense all-H lanes iff 256*d >= dense_cut*H, and otherwise loops
over its own labels. The cutoff is tuned for the BASELINE's benefit, on instances that are
not among the 40 benchmarked primary ones (every 8th of the other 160 by size), on the
463-haplotype panel: per-haplotype kernel time for each cutoff, interleaved per instance,
median of repetitions. The cutoff with the smallest total time becomes gap_native.DENSE_CUT.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import random
import statistics
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panel_common as pc  # noqa: E402
import gap_native as gn  # noqa: E402

CUTS = (0, 16, 32, 64, 96, 128, 192, 257)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--panel', default='hprc463')
    ap.add_argument('--repeat', type=int, default=5)
    args = ap.parse_args()
    if args.out.exists():
        raise SystemExit('refusing to overwrite')
    p = pc.prepared(args.panel)
    G = gn.GapGraph(p)
    _, cases = pc.load_panel(args.panel)
    primary = set(pc.primary_ids())
    others = sorted((c for i, c in cases.items() if i not in primary and c['anchors']), key=lambda c: len(c['anchors']))
    picks = others[::8]
    rng = random.Random(1717)
    rows = []
    for case in picks:
        b = G.batch(pc.anchors(case), pc.GAMMA)
        times = {c: [] for c in CUTS}
        digests = set()
        for _ in range(args.repeat):
            order = list(CUTS)
            rng.shuffle(order)
            for cut in order:
                res = G.solve_packed(b, mode='per_haplotype', dense_cut=cut)
                times[cut].append(res.stats['kernel_seconds'])
                digests.add(res.arrays_digest())
        assert len(digests) == 1
        med = {c: statistics.median(v) for c, v in times.items()}
        sparse = G.solve_packed(b, mode='per_haplotype', dense_cut=gn.DENSE_CUT).stats['ph_sparse_anchors']
        rows.append(dict(id=case['id'], N=b.n, K=b.K, density=b.K / (b.n * p.H), median_seconds=med,
                         sparse_anchors_at_default=sparse))
        print(case['id'], b.n, {c: round(v * 1e3, 3) for c, v in med.items()}, flush=True)
    totals = {c: sum(r['median_seconds'][c] for r in rows) for c in CUTS}
    best = min(totals, key=totals.get)
    out = dict(script='experiments/hprc_mhc_panel/calibrate_dense.py', panel=args.panel, H=p.H, instances=len(rows),
               cuts=CUTS, total_seconds=totals, best_cut=best,
               speed_vs_always_dense={c: totals[0] / totals[c] for c in CUTS}, rows=rows)
    args.out.write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps({k: v for k, v in out.items() if k != 'rows'}, indent=1))


if __name__ == '__main__':
    main()
