"""Summarize the integrated Minichain runs (2026-10-20): correctness comparisons and end-to-end timing.

Correctness: copies the per-dataset comparisons from the correctness manifest and totals them.
Timing: for every dataset and variant, the median over repetitions of Minichain's own reported real time
(whole process: graph loading and indexing, seeding, chaining, backtracking, output), CPU time and peak RSS,
and of the summed per-read chaining time (RECOMB_DP_TIME: Chaining() and its DP section). Speedups are
ratios of medians (off over variant) and medians of per-repetition paired ratios.

  python3.13 -B experiments/minichain_integrated/summarize.py --correctness experiments/minichain_integrated/results/correctness \
      [--timing experiments/minichain_integrated/results/timing] --out experiments/minichain_integrated/results/summary.json
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import statistics

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def time_tsv(p: Path):
    total = dp = 0
    reads = 0
    for line in p.open():
        a = line.rstrip('\n').split('\t')
        assert a[0] == 'T'
        total += int(a[2])
        dp += int(a[3])
        reads += 1
    return dict(reads=reads, chaining_seconds=total / 1e9, dp_seconds=dp / 1e9)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--correctness', type=Path, required=True)
    ap.add_argument('--timing', type=Path)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    cm = json.loads((args.correctness / 'manifest.json').read_text())
    assert cm['status'] == 'passed', cm['status']
    out = dict(correctness=cm['comparisons'], totals={})
    for v in ('ph', 'class', 'auto'):
        t = defaultdict(int)
        for d, c in cm['comparisons'].items():
            assert c['off_gaf_equals_pinned'], d
            for k in ('components', 'dp_hash_mismatches', 'reads', 'gaf_identical_reads', 'gaf_different_reads',
                      'different_reads_with_same_path_and_target_interval'):
                t[k] += c[v][k]
        out['totals'][v] = dict(t)
    if args.timing:
        tm = json.loads((args.timing / 'manifest.json').read_text())
        assert tm['status'] == 'passed', tm['status']
        groups = defaultdict(list)
        for r in tm['runs']:
            run_dir = args.timing / f"rep{r['rep']}" / r['dataset'] / r['variant']
            groups[(r['dataset'], r['variant'])].append(dict(r, **time_tsv(run_dir / 'time.tsv')))
        timing = {}
        for (d, v), rs in sorted(groups.items()):
            timing.setdefault(d, {})[v] = dict(
                reps=len(rs), real_seconds=statistics.median(x['minichain_real_seconds'] for x in rs),
                cpu_seconds=statistics.median(x['cpu_seconds'] for x in rs),
                peak_rss_gb=statistics.median(x['peak_rss_gb'] for x in rs),
                chaining_seconds=statistics.median(x['chaining_seconds'] for x in rs),
                dp_seconds=statistics.median(x['dp_seconds'] for x in rs),
                engine_build_seconds=statistics.median(x['engine_build_seconds'] for x in rs) if rs[0]['engine_build_seconds'] is not None else None,
                load_before=[x['host_before']['load'][0] for x in rs])
        for d, vs in timing.items():
            off = vs['off']
            for v, s in vs.items():
                if v == 'off':
                    continue
                pairs = [(a['minichain_real_seconds'], b['minichain_real_seconds'])
                         for a in groups[(d, 'off')] for b in groups[(d, v)] if a['rep'] == b['rep']]
                s['end_to_end_speedup'] = off['real_seconds'] / s['real_seconds']
                s['end_to_end_speedup_paired_median'] = statistics.median(x / y for x, y in pairs)
                s['chaining_speedup'] = off['chaining_seconds'] / s['chaining_seconds']
        out['timing'] = timing
    args.out.write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out['totals'], indent=1))
    if 'timing' in out:
        for d, vs in out['timing'].items():
            print(d, {v: (round(s['real_seconds'], 2), s.get('end_to_end_speedup')) for v, s in vs.items()})


if __name__ == '__main__':
    main()
