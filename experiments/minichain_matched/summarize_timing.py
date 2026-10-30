"""Join Minichain's own per-read timers (RECOMB_DP_TIME runs) with bench_matched.py rows.

Per read, Minichain's numbers are the steady-clock time of graphUtils::Chaining() and of
its DP section (tuple construction, sort and sweep), summed over its components; ours are
the median native kernel times summed over the same components. Several Minichain timing
runs of the same reads are combined by taking each read's median. Ratios are paired per
read (Minichain time / our time). Totals are sums over reads.
"""
from __future__ import annotations

import argparse
import json
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mcdp  # noqa: E402


def quantiles(xs):
    xs = sorted(xs)
    if not xs:
        return None
    q = lambda p: xs[min(len(xs) - 1, max(0, round(p * (len(xs) - 1))))]  # noqa: E731
    return dict(min=xs[0], p10=q(0.10), median=statistics.median(xs), p90=q(0.90), max=xs[-1], n=len(xs))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--minichain', nargs='+', required=True, type=Path, help='timing run directories (time.tsv)')
    ap.add_argument('--ours', required=True, type=Path, help='bench_matched.py JSONL')
    ap.add_argument('--label', required=True)
    ap.add_argument('--out', required=True, type=Path)
    args = ap.parse_args()
    runs = [mcdp.read_time_dump(d / 'time.tsv') for d in args.minichain]
    reads = set(runs[0])
    for r in runs[1:]:
        assert set(r) == reads, 'timing runs cover different reads'
    mc = {name: dict(total=statistics.median(r[name]['total_ns'] for r in runs) / 1e9,
                     dp=statistics.median(r[name]['dp_ns'] for r in runs) / 1e9,
                     n=runs[0][name]['n']) for name in reads}
    ours = {}
    with open(args.ours) as f:
        for line in f:
            row = json.loads(line)
            ours[row['read']] = row['kernel_total']
    common = sorted(set(mc) & set(ours))
    missing = sorted(set(mc) ^ set(ours))
    modes = list(next(iter(ours.values())).keys())
    out = dict(label=args.label, reads=len(common), unmatched_reads=len(missing),
               minichain_runs=[str(d) for d in args.minichain],
               totals_seconds=dict(minichain_chaining=sum(mc[r]['total'] for r in common),
                                   minichain_dp=sum(mc[r]['dp'] for r in common),
                                   **{m: sum(ours[r][m] for r in common) for m in modes}),
               per_read_ms=dict(minichain_chaining=quantiles([mc[r]['total'] * 1e3 for r in common]),
                                minichain_dp=quantiles([mc[r]['dp'] * 1e3 for r in common]),
                                **{m: quantiles([ours[r][m] * 1e3 for r in common]) for m in modes}),
               paired_ratio=dict())
    for m in modes:
        out['paired_ratio'][f'minichain_chaining/{m}'] = quantiles([mc[r]['total'] / ours[r][m] for r in common])
        out['paired_ratio'][f'minichain_dp/{m}'] = quantiles([mc[r]['dp'] / ours[r][m] for r in common])
    for m in modes:
        if m != 'per_haplotype':
            out['paired_ratio'][f'per_haplotype/{m}'] = quantiles([ours[r]['per_haplotype'] / ours[r][m]
                                                                  for r in common])
    args.out.write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(dict(label=out['label'], reads=out['reads'], totals_seconds=out['totals_seconds'],
                          median_ratios={k: round(v['median'], 2) for k, v in out['paired_ratio'].items()}),
                     indent=1))


if __name__ == '__main__':
    main()
