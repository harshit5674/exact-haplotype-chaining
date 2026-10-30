"""Per-anchor work counters of the native engines on every panel (no timing is recorded).

For each panel and instance this runs, once each, the native class engine (solve_compact,
cold pair cache) and the native per-haplotype baseline (solve_packed mode='per_haplotype'),
records their deterministic work counters (gap_native.STAT_NAMES; every *_seconds field is
dropped), and checks that the expanded class output and the per-haplotype output have
identical C, B, R arrays (an exactness cross-check on every instance; the full checks are in
verify_locus.py). Sources:
  kir, chr20  data/<locus>/<panel>/ with the simulated reads (primary = first 20 reads x 2)
  mhc_sim     data/mhc_sim/<panel>/: simulated CHM13 reads on the MHC panels (realism check)
  mhc         experiments/hprc_mhc_panel/data/<panel>/ with the real CHM13 PacBio reads
              (primary = the MHC pilot's 40 instances): recomputes the MHC report's counters
Output: counters_<source>.jsonl (one row per panel, gamma and instance) and, with
--summarize, counters_summary.json: medians [min, median, max] over the 40 primary instances
(and over all instances) computed exactly as experiments/gap_engine_native/summarize.py does.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import statistics
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import locus_common as lc  # noqa: E402
import locus_panels as lp  # noqa: E402

lp.native_library_ready()
import gap_chaining as gc  # noqa: E402
import gap_native as gn  # noqa: E402

COUNTERS = ('merges', 'contested_merges', 'pairs', 'pairs_same_vertex', 'pairs_loose', 'pairs_range',
            'pairs_sliced', 'sliced_thresholds', 'sliced_minima', 'pair_ranges', 'range_labels',
            'route_planes', 'cut_records', 'prefix_records', 'remote_groups', 'remote_candidates',
            'remote_sliced', 'published_records', 'output_records', 'max_node_records')
PH_COUNTERS = ('ph_label_steps', 'ph_remote_queries', 'ph_sparse_anchors')
PER_ANCHOR = ('merges', 'contested_merges', 'pairs', 'pairs_sliced', 'sliced_thresholds', 'sliced_minima',
              'range_labels', 'route_planes', 'cut_records', 'prefix_records', 'remote_groups',
              'remote_candidates', 'remote_sliced', 'published_records', 'output_records')
SHARES = ('pairs_same_vertex', 'pairs_loose', 'pairs_range', 'pairs_sliced')
MHC_PANELS = ('hprc60', 'hprc120', 'hprc240', 'hprc463', 'hprc462_noCHM13')


def source_panels(source):
    if source == 'mhc':
        return list(MHC_PANELS)
    if source == 'mhc_sim':
        return ['hprc60', 'hprc463']
    info = lp.panels(source)
    return [info['all'], *info['subsets'], info['leaveout']]


def source_data(source, name):
    d = lc.MHC / 'data' / name if source == 'mhc' else lc.DATA / source / name
    graph = json.loads((d / 'graph.json').read_text())
    with (d / 'instances.jsonl').open() as f:
        cases = [json.loads(s) for s in f]
    if source == 'mhc':
        rows = (lc.ROOT / 'experiments/mhc_pilot/final_hybrid.jsonl').read_text().splitlines()
        primary = {json.loads(s)['id'] for s in rows}
        groups = {}
    else:
        primary = set(json.loads((lc.DATA / source / 'simulation.json').read_text())['primary_ids'])
        groups = lp.read_groups(source)
    assert len(primary) == 40 and primary <= {c['id'] for c in cases}
    return graph, cases, primary, groups


def run(source, name, gamma, out, scope):
    graph, cases, primary, groups = source_data(source, name)
    p = gc.Graph(graph['lengths'], [tuple(e) for e in graph['edges']], graph['paths']).prepare()
    G = gn.GapGraph(p)
    bad = 0
    for case in cases:
        if scope == 'primary' and case['id'] not in primary:
            continue
        A = [gc.Anchor(*x[:5], x[5] * lc.SCALE) for x in case['anchors']]
        b = G.batch(A, gamma)
        comp = G.solve_compact(b)
        exp = comp.expand()
        ph = G.solve_packed(b, mode='per_haplotype')
        same = exp.arrays_digest() == ph.arrays_digest()
        bad += not same
        row = dict(source=source, panel=name, H=p.H, V=p.V, W=G.W, b=G.b, gamma=gamma, id=case['id'],
                   primary=case['id'] in primary, group=groups.get(case['read'], {}).get('group', 'present'),
                   qlen=case['qlen'], N=b.n, K=b.K, digest=exp.arrays_digest(), class_equals_per_haplotype=same,
                   compact_records=comp.stats['compact_records'],
                   **{k: int(comp.stats[k]) for k in COUNTERS}, **{k: int(ph.stats[k]) for k in PH_COUNTERS})
        out.write(json.dumps(row) + '\n')
        del comp, exp, ph
    out.flush()
    print(source, name, gamma, 'mismatches', bad, flush=True)
    return bad


def rng3(xs):
    return [round(min(xs), 3), round(statistics.median(xs), 3), round(max(xs), 3)]


def summarize_rows(rs):
    per = lambda key: [r[key] / r['N'] for r in rs]  # noqa: E731
    item = dict(instances=len(rs), anchors=sum(r['N'] for r in rs), states=sum(r['K'] for r in rs),
                states_per_anchor=round(sum(r['K'] for r in rs) / sum(r['N'] for r in rs), 1),
                anchors_per_instance=rng3([r['N'] for r in rs]),
                counters_per_anchor={k: rng3(per(k)) for k in PER_ANCHOR},
                pair_decision_shares={k: rng3([r[k] / max(1, r['pairs']) for r in rs]) for k in SHARES},
                pooled_pair_decision_shares={k: round(sum(r[k] for r in rs) / max(1, sum(r['pairs'] for r in rs)), 4)
                                             for k in SHARES},
                pooled_per_anchor={k: round(sum(r[k] for r in rs) / sum(r['N'] for r in rs), 4)
                                   for k in ('pairs_sliced', 'sliced_thresholds', 'sliced_minima', 'output_records',
                                             'published_records')},
                max_node_records=rng3([r['max_node_records'] for r in rs]),
                pair_ranges_per_read=rng3([r['pair_ranges'] for r in rs]),
                per_haplotype_label_steps_per_anchor=rng3(per('ph_label_steps')),
                per_haplotype_remote_queries_per_anchor=rng3(per('ph_remote_queries')),
                class_equals_per_haplotype=sum(r['class_equals_per_haplotype'] for r in rs))
    return item


def summarize(paths, out_path):
    rows = [json.loads(s) for p in paths for s in p.read_text().splitlines()]
    out = {}
    for key in sorted({(r['source'], r['panel'], r['gamma']) for r in rows}):
        rs = [r for r in rows if (r['source'], r['panel'], r['gamma']) == key]
        entry = dict(source=key[0], panel=key[1], gamma=key[2], H=rs[0]['H'], V=rs[0]['V'], W=rs[0]['W'], b=rs[0]['b'],
                     primary=summarize_rows([r for r in rs if r['primary']]))
        if len(rs) > 40:
            entry['all'] = summarize_rows(rs)
            for g in sorted({r['group'] for r in rs}):
                entry[f'all_{g}'] = summarize_rows([r for r in rs if r['group'] == g])
        out[f'{key[0]}/{key[1]}/g{key[2]}'] = entry
    lc.dump(out_path, dict(script='experiments/second_locus/counters.py', sources=[str(p.name) for p in paths], levels=out))
    print(f"{'level':>34} {'inst':>4} {'st/anc':>6} {'same':>6} {'loose':>6} {'range':>6} {'slice':>6} "
          f"{'pairs_sl/a':>10} {'thr/a':>6} {'min/a':>6} {'out/a':>6} {'maxnode':>7}")
    for key, e in out.items():
        s = e['primary']
        sh, ca = s['pair_decision_shares'], s['counters_per_anchor']
        print(f"{key:>34} {s['instances']:>4} {s['states_per_anchor']:>6} {sh['pairs_same_vertex'][1]:>6} {sh['pairs_loose'][1]:>6} "
              f"{sh['pairs_range'][1]:>6} {sh['pairs_sliced'][1]:>6} {ca['pairs_sliced'][1]:>10} {ca['sliced_thresholds'][1]:>6} "
              f"{ca['sliced_minima'][1]:>6} {ca['output_records'][1]:>6} {s['max_node_records'][1]:>7}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sources', default='kir,chr20,mhc_sim,mhc')
    ap.add_argument('--gamma', type=lambda s: [int(x) for x in s.split(',')], default=[10000, 20000])
    ap.add_argument('--scope', choices=('all', 'primary'), default='all', help="'all' 200 instances or the 40 primary")
    ap.add_argument('--summarize', type=Path, help='only summarize existing counters_*.jsonl into this file')
    args = ap.parse_args()
    if args.summarize:
        if args.summarize.exists():
            raise SystemExit('refusing to overwrite')
        summarize([HERE / f'counters_{s}.jsonl' for s in args.sources.split(',')], args.summarize)
        return 0
    gn.build()
    bad = 0
    for source in args.sources.split(','):
        path = HERE / f'counters_{source}.jsonl'
        if path.exists():
            raise SystemExit(f'refusing to overwrite {path}')
        with path.open('w') as out:
            for name in source_panels(source):
                for gamma in args.gamma:
                    bad += run(source, name, gamma, out, args.scope)
    print('total class/per-haplotype mismatches', bad)
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
