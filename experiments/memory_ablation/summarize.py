"""Summarize experiments/memory_ablation/run.py output per panel (2026-10-20).

Per read, the two instances (components) run one after the other, so a read's native peak is the larger
of its two instance peaks; its working total adds the returned arrays of that instance. Per-anchor counts
sum both instances, then divide by the read's anchors, as the paper's counters do. Medians are over the
100 reads of a panel; maxima are over reads. Minichain's whole-process peak RSS on the HPRC MHC panels is
read from the logs of the matched timing runs (experiments/minichain_matched/data/time_*), for context.

  python3.13 -B experiments/memory_ablation/summarize.py --run experiments/memory_ablation/results/run1 \
      --out experiments/memory_ablation/results/summary.json
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import re
import statistics

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MODES = ('per_haplotype', 'class', 'auto', 'compact')
COUNTS = ('states', 'keys', 'scores', 'score_values', 'records', 'published')
MB = 1 << 20


def per_read(rows):
    reads = defaultdict(list)
    for r in rows:
        reads[r['read']].append(r)
    assert len(reads) == 100 and all(len(v) == 2 for v in reads.values()), len(reads)
    out = []
    for name, (a, b) in sorted(reads.items()):
        rec = dict(read=name, N=a['N'] + b['N'], K=a['K'] + b['K'])
        for m in MODES:
            rec[m + '.native_peak'] = max(a[m]['native_peak'], b[m]['native_peak'])
            rec[m + '.output'] = max(a[m]['output_bytes'], b[m]['output_bytes'])
            rec[m + '.total'] = max(a[m]['native_peak'] + a[m]['output_bytes'] + a['input_bytes'],
                                    b[m]['native_peak'] + b[m]['output_bytes'] + b['input_bytes'])
        rec['fallback'] = a['auto']['fallback'] or b['auto']['fallback']
        for k in COUNTS:
            rec['per_anchor.' + k] = (a['counts'][k] + b['counts'][k]) / rec['N']
        out.append(rec)
    return out


def minichain_rss():
    out = {}
    for d in sorted((ROOT / 'experiments/minichain_matched/data').glob('time_hprc*_R10000_rep*')):
        m = re.search(r'Peak RSS: ([0-9.]+) GB', (d / 'minichain.log').read_text())
        h = int(re.search(r'hprc(\d+)', d.name).group(1))
        out.setdefault(h, []).append(float(m.group(1)))
    return {h: dict(reps=len(v), median_gb=statistics.median(v), reads=100) for h, v in sorted(out.items())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    manifest = json.loads((args.run / 'manifest.json').read_text())
    assert manifest['status'] == 'passed', manifest['status']
    panels = []
    for p in manifest['panels']:
        key = f"{p['locus']}-{p['panel']}"
        rows = [json.loads(line) for line in open(args.run / f'{key}.jsonl')]
        assert len(rows) == 200
        reads = per_read(rows)
        s = dict(locus=p['locus'], panel=p['panel'], H=p['H'], V=p['V'], W=p['W'],
                 index_native_mb=p['index_native_bytes'] / MB, index_native_peak_mb=p['index_native_peak'] / MB,
                 fallback_reads=sum(r['fallback'] for r in reads))
        for m in MODES:
            for part in ('native_peak', 'output', 'total'):
                xs = [r[f'{m}.{part}'] / MB for r in reads]
                s[f'{m}.{part}.median_mb'] = statistics.median(xs)
                s[f'{m}.{part}.max_mb'] = max(xs)
        for k in COUNTS:
            s['per_anchor.' + k] = statistics.median(r['per_anchor.' + k] for r in reads)
        for k in ('keys', 'scores'):
            s[f'ratio.{k}_over_records'] = statistics.median(r['per_anchor.' + k] / r['per_anchor.records'] for r in reads)
        s['ratio.per_haplotype_over_class_native'] = statistics.median(
            r['per_haplotype.native_peak'] / r['class.native_peak'] for r in reads)
        s['ratio.expanded_over_compact_output'] = statistics.median(
            r['per_haplotype.output'] / r['compact.output'] for r in reads)
        panels.append(s)
    out = dict(run=str(args.run.resolve().relative_to(ROOT)), panels=panels, minichain_peak_rss=minichain_rss())
    args.out.write_text(json.dumps(out, indent=1) + '\n')
    hdr = (f"{'panel':16} {'H':>4} {'idxMB':>6} | native peak MB (median): {'ph':>6} {'class':>6} {'auto':>6} "
           f"{'cmpct':>6} | output MB: {'exp':>6} {'cmpct':>6} | per anchor: {'states':>7} {'keys':>7} "
           f"{'scores':>7} {'recs':>5}")
    print(hdr)
    for s in panels:
        print(f"{s['locus'] + '-' + s['panel']:16} {s['H']:>4} {s['index_native_mb']:6.1f} | "
              f"{'':24}{s['per_haplotype.native_peak.median_mb']:6.2f} {s['class.native_peak.median_mb']:6.2f} "
              f"{s['auto.native_peak.median_mb']:6.2f} {s['compact.native_peak.median_mb']:6.2f} | {'':11}"
              f"{s['per_haplotype.output.median_mb']:6.2f} {s['compact.output.median_mb']:6.2f} | {'':12}"
              f"{s['per_anchor.states']:7.1f} {s['per_anchor.keys']:7.2f} {s['per_anchor.scores']:7.2f} "
              f"{s['per_anchor.records']:5.2f}")
    print('Minichain peak RSS (GB, 100 reads):', {h: v['median_gb'] for h, v in out['minichain_peak_rss'].items()})


if __name__ == '__main__':
    main()
