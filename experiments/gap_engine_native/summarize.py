"""Summarize benchmark.jsonl: paired speedups, absolute times, stages and work counters.

Ratios are paired per read: per-haplotype time divided by class time on the same read,
each the median of its repetitions. Summaries give the median over reads with the
[minimum, maximum] range and the number of reads where the class mode was faster.
"compact+expand" adds the separately timed expansion to the compact kernel time.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import statistics


def rng3(xs):
    return [round(min(xs), 3), round(statistics.median(xs), 3), round(max(xs), 3)]


def summarize(rows):
    out = {}
    for H in sorted({r['H'] for r in rows}):
        for gamma in sorted({r['gamma'] for r in rows if r['H'] == H}):
            rs = [r for r in rows if r['H'] == H and r['gamma'] == gamma]
            km = lambda r, m: r['kernel_median'][m]  # noqa: E731
            pm = lambda r, m: r['packed_median'][m]  # noqa: E731
            ex = lambda r, m: r['expansion_median'][m]  # noqa: E731
            ph = [km(r, 'per_haplotype') for r in rs]
            variants = {
                'class_compact': [km(r, 'class_compact') for r in rs],
                'class_expanded': [km(r, 'class_expanded') for r in rs],
                'class_compact_plus_expansion': [km(r, 'class_compact') + ex(r, 'class_compact') for r in rs],
                'class_compact_warm': [km(r, 'class_compact_warm') for r in rs],
            }
            item = dict(reads=len(rs), anchors=sum(r['N'] for r in rs), states=sum(r['K'] for r in rs), W=rs[0]['W'],
                        kernel_ms_median={'per_haplotype': round(statistics.median(ph) * 1e3, 3),
                                          **{k: round(statistics.median(v) * 1e3, 3) for k, v in variants.items()}},
                        kernel_speedup={}, kernel_wins={},
                        packed_ms_median={m: round(statistics.median(pm(r, m) for r in rs) * 1e3, 3)
                                          for m in ('per_haplotype', 'class_expanded', 'class_compact')},
                        packed_speedup={}, expansion_ms_median=round(statistics.median(ex(r, 'class_compact') for r in rs) * 1e3, 3))
            for k, v in variants.items():
                ratios = [a / b for a, b in zip(ph, v)]
                item['kernel_speedup'][k] = rng3(ratios)
                item['kernel_wins'][k] = sum(x > 1 for x in ratios)
            for k in ('class_expanded', 'class_compact'):
                ratios = [pm(r, 'per_haplotype') / pm(r, k) for r in rs]
                item['packed_speedup'][k] = rng3(ratios)
            item['total_kernel_speedup'] = {k: round(sum(ph) / sum(v), 3) for k, v in variants.items()}
            # Stage breakdowns (median over reads, milliseconds) from the first repetition's stats.
            item['stages_ms'] = {m: {s: round(statistics.median(r['stats'][m][s] for r in rs) * 1e3, 3)
                                     for s in ('setup_seconds', 'init_seconds', 'sweep_seconds', 'export_seconds')}
                                 for m in ('per_haplotype', 'class_expanded', 'class_compact')}
            # Work counters per anchor (class compact; per-haplotype steps and queries).
            per = lambda key, m='class_compact': [r['stats'][m][key] / r['N'] for r in rs]  # noqa: E731
            pairs = lambda key: [r['stats']['class_compact'][key] / max(1, r['stats']['class_compact']['pairs']) for r in rs]  # noqa: E731
            item['counters_per_anchor'] = {k: rng3(per(k)) for k in (
                'merges', 'contested_merges', 'pairs', 'pairs_sliced', 'sliced_thresholds', 'sliced_minima',
                'range_labels', 'route_planes', 'cut_records', 'prefix_records', 'remote_groups',
                'remote_candidates', 'remote_sliced', 'published_records', 'output_records')}
            item['pair_decision_shares'] = {k: rng3(pairs(k)) for k in ('pairs_same_vertex', 'pairs_loose', 'pairs_range', 'pairs_sliced')}
            item['pair_ranges_per_read'] = rng3([r['stats']['class_compact']['pair_ranges'] for r in rs])
            item['max_node_records'] = rng3([r['stats']['class_compact']['max_node_records'] for r in rs])
            item['per_haplotype_label_steps_per_anchor'] = rng3(per('ph_label_steps', 'per_haplotype'))
            item['per_haplotype_remote_queries_per_anchor'] = rng3(per('ph_remote_queries', 'per_haplotype'))
            item['warm_cache_after_prime'] = rs[0]['warm_cache_after_prime']
            item['warm_prime_seconds'] = round(rs[0]['warm_prime_seconds'], 3)
            out[f'H{H}_gamma{gamma}'] = item
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('rows', type=Path)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    if args.out.exists():
        raise SystemExit('refusing to overwrite')
    rows = [json.loads(s) for s in args.rows.read_text().splitlines()]
    summary = summarize(rows)
    args.out.write_text(json.dumps(dict(source=str(args.rows), levels=summary), indent=1) + '\n')
    print(f"{'level':>16} {'reads':>5} {'ph ms':>8} {'cmp ms':>7} {'exp ms':>7} {'cmp x':>20} {'exp x':>20} {'cmp+exp x':>20} {'warm x':>20} wins(c/e)")
    for key, it in summary.items():
        km, ks = it['kernel_ms_median'], it['kernel_speedup']
        print(f"{key:>16} {it['reads']:>5} {km['per_haplotype']:>8} {km['class_compact']:>7} {km['class_expanded']:>7} "
              f"{str(ks['class_compact']):>20} {str(ks['class_expanded']):>20} {str(ks['class_compact_plus_expansion']):>20} "
              f"{str(ks['class_compact_warm']):>20} {it['kernel_wins']['class_compact']}/{it['kernel_wins']['class_expanded']}")


if __name__ == '__main__':
    main()
