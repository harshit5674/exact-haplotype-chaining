"""Giraffe's single-record rule against the exact optimum on the base-level window panels (2026-10-20).

Pre-specified in REPORT.md: the exactness half of the base-level experiment. The unchanged Phase 2 runner of
experiments/exactness_payoff (faithful emulation of vg v1.77.0's DP, consistency bonus 0 and 12, against the
exact per-haplotype DP of the same objective; every reported chain re-scored independently in Python) runs on
the windows' real reads. Giraffe indexes at most 64 path identifiers, so, as in that study's Phase 4, the panel
is a seeded uniform random 64-haplotype subset without CHM13 (the reads' source), pruned to the subset's
vertices and edges; anchors on dropped vertices are removed. Three replicates per window; rho = 2 (Giraffe's
default), 10, 28 and 56. This module only registers the window panels with run_phase2; nothing in it is new
scoring code.

  python3.13 -B experiments/baselevel_feasibility/exactness_windows.py --out experiments/baselevel_feasibility/exactness
  python3.13 -B experiments/baselevel_feasibility/exactness_windows.py --summarize --out experiments/baselevel_feasibility/exactness
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'experiments/exactness_payoff'))
import run_phase2 as p2  # noqa: E402

WINDOWS = ('w31250k_250kb', 'w31700k_250kb', 'w32375k_100kb')
SIZE, REPLICATES, SEED = 64, 3, 20261020
RHOS = '2,10,28,56'


def names():
    return [f'{w}_s{SIZE}_r{r}' for w in WINDOWS for r in range(REPLICATES)]


def register():
    """Make the window panels and their subsets known to run_phase2 (also in its spawned workers)."""
    for w in WINDOWS:
        p2.DATASETS[w] = HERE / 'windows' / w / 'panel'
        for r in range(REPLICATES):
            p2.SUBSET16[f'{w}_s{SIZE}_r{r}'] = (w, SIZE, SEED + 1000 * WINDOWS.index(w) + r)


register()


def summarize(out: Path):
    res = {}
    for name in names():
        rows = [json.loads(line) for line in (out / f'{name}.jsonl').open()]
        meta = json.loads((out / f'{name}.jsonl.meta.json').read_text())
        by_rho = defaultdict(list)
        for r in rows:
            by_rho[r['rho']].append(r)
        res[name] = dict(H=meta['H'], V=meta['V'], instances=meta['instances'], dropped_anchors=meta['dropped_anchors'],
                         rhos={})
        for rho, rs in sorted(by_rho.items()):
            entry = dict(instances=len(rs), exact_with_switch=sum(1 for r in rs if r['exact_switches'] > 0))
            for b in ('h0', 'h12'):
                gaps = [r for r in rs if r[b]['gap'] > 0]
                entry[b] = dict(instances_below=len(gaps),
                                max_rel_gap=max((r[b]['rel_gap'] for r in gaps), default=0.0),
                                max_gap=max((r[b]['gap'] for r in gaps), default=0))
            res[name]['rhos'][str(rho)] = entry
    (out / 'summary.json').write_text(json.dumps(res, indent=1) + '\n')
    for name, v in res.items():
        print(name, v['instances'], {rho: (e['h0']['instances_below'], e['h12']['instances_below'],
                                         round(100 * max(e['h0']['max_rel_gap'], e['h12']['max_rel_gap']), 3),
                                         e['exact_with_switch']) for rho, e in v['rhos'].items()})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--workers', type=int, default=3)
    ap.add_argument('--summarize', action='store_true')
    ap.add_argument('--only', default='', help='comma-separated subset names (default: all)')
    ap.add_argument('--limit', type=int, default=0, help='first n instances only (smoke tests)')
    args = ap.parse_args()
    if args.summarize:
        summarize(args.out)
        return
    args.out.mkdir(parents=True, exist_ok=True)
    for name in [n for n in names() if not args.only or n in args.only.split(',')]:
        target = args.out / f'{name}.jsonl'
        if target.exists():
            print(name, 'exists; skipped', flush=True)
            continue
        sys.argv = ['run_phase2.py', '--dataset', name, '--out', str(target), '--rhos', RHOS,
                    '--workers', str(args.workers)] + (['--limit', str(args.limit)] if args.limit else [])
        p2.main()


if __name__ == '__main__':
    main()
