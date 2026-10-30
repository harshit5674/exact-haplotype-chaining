"""Seeded random haplotype subsets of hprc463 for Phase 4 (panel sampling).

Giraffe can index at most 64 path IDs, and its own evaluations sample panels to 16 haplotypes
(+ reference) with vg haplotypes, which is k-mer-guided and builds synthetic recombinants. Here the
subsets are uniform random (a simplification, stated in REPORT.md), always without CHM13 unless drawn.

  python3.13 experiments/exactness_payoff/make_subsets.py --out experiments/exactness_payoff/phase4_subsets.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import random

ROOT = Path(__file__).resolve().parents[2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--seed', type=int, default=4630)
    ap.add_argument('--sizes', default='16,64')
    ap.add_argument('--replicates', type=int, default=3)
    args = ap.parse_args()
    out = Path(args.out)
    if out.exists():
        raise SystemExit(f'{out} exists')
    g = json.loads((ROOT / 'experiments/hprc_mhc_panel/data/hprc463/graph.json').read_text())
    H = len(g['paths'])
    rng = random.Random(args.seed)
    subsets = []
    for size in (int(x) for x in args.sizes.split(',')):
        for rep in range(args.replicates):
            haps = sorted(rng.sample(range(H), size))
            subsets.append(dict(name=f's{size}_r{rep}', size=size, haplotypes=haps,
                                names=[g['haplotypes'][h] for h in haps]))
    out.write_text(json.dumps(dict(seed=args.seed, H=H, subsets=subsets), indent=1) + '\n')
    print(json.dumps([(s['name'], s['haplotypes'][:5]) for s in subsets]))


if __name__ == '__main__':
    main()
