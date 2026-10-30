"""Stage A cost model: certified fast paths versus bit-slicing, and predicted native work.

For every primary MHC read (Minichain scaling: weights x 200, gamma 10000) and every
panel size H in {60, 120, 240, 480, 960} (the real 60 paths plus seed-6060 mosaics),
this runs the specification (gap_chaining.OriginSolver) and two plans of the fast-path
solver (gap_fastpath.FastPathSolver), each with a cold pair cache per read, and
requires identical records from both:

  prompt  the Stage A plan: certified merge tests (i) same origin vertex,
          (ii) loose per-vertex bounds, (iii) exact cached pair range; range fast
          paths for floors, switches and minima; one remote point query per
          missing label, scanning short node lists;
  full    the same plus two further certified steps: a grouped remote query per
          last2reach vertex, and the range floor cut in prefix queries.

Per read it records how the specification's bit-sliced merges, all contested merges
and all record pairs are decided; floor, switch and minimum decisions; bit-sliced
columns (the specification's merge columns versus every column still needed, by
source); distinct origin-vertex pairs; and gap_fastpath's unit-cost prediction:
class word operations (compact, and expanded with 2K writes; cold pair ranges, and
warm with ranges precomputed for the graph) versus per-haplotype Fenwick steps.
Ablations (H in --ablate) rerun the prompt plan with tests 'none', 'i' and 'ii'.
Only operation counts are reported; Python timings predict nothing native.
"""
from __future__ import annotations
import argparse
import json
from multiprocessing import get_context
from pathlib import Path
import statistics
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import common as cm  # noqa: E402
import gap_chaining as gc  # noqa: E402
import gap_fastpath as fp  # noqa: E402

_PREP = {}
PLANS = dict(prompt=dict(tests='iii', remote='labels', cut='none'),
             full=dict(tests='iii', remote='grouped', cut='range'))


def graph(H):
    if H not in _PREP:
        p = cm.prepared(H)
        _PREP[H] = (p, fp.VertexBounds(p), fp.RemoteGroups(p))
    return _PREP[H]


def plan_row(rf, W, K):
    c = rf.cost
    pred = fp.predicted_word_ops(c, W, K)
    warm = pred['class_compact'] - c['labels']  # pair ranges precomputed for the graph
    row = dict(cost=dict(sorted(c.items())), predicted=pred,
               gain_compact=pred['per_haplotype'] / pred['class_compact'],
               gain_expanded=pred['per_haplotype'] / pred['class_expanded'],
               gain_compact_warm=pred['per_haplotype'] / warm,
               gain_expanded_warm=pred['per_haplotype'] / (warm + 2 * K))
    if W == 15:  # native templates round 15 words up to 16
        alt = 16 * (c['mask'] + c['col']) + c['scalar'] + c['labels']
        row['gain_compact_W16'] = pred['per_haplotype'] / alt
        row['gain_expanded_W16'] = pred['per_haplotype'] / (alt + 2 * K)
    return row


def one(job):
    H, idx, gamma, ablate = job
    p, bounds, groups = graph(H)
    case = cm.primary_cases()[idx]
    A = cm.anchors(case)
    n = len(A)
    rs = gc.OriginSolver(p, A, gamma).run()
    K, W = rs.stats['valid_states'], p.W
    s = rs.stats
    row = dict(H=H, id=case['id'], N=n, K=K, W=W, b=p.b, gamma=gamma,
               missing_incidences=s['missing_incidences'],
               spec=dict(merges=s['merges'], scalar_fast_merges=s['scalar_fast_merges'],
                         coded_merges=s['direct_coded_merges'] + s['clustered_merges'],
                         clustered_merges=s['clustered_merges'], coded_columns=s['coded_columns'],
                         route_calls=s['route_calls'], route_columns=s['route_columns'],
                         threshold_calls=s['threshold_calls'], minimum_calls=s['minimum_calls'],
                         lambda_records=s['lambda_records'], output_classes=s['output_classes'],
                         published_classes=s['published_classes'], max_node_classes=s['max_node_classes']),
               plans={})
    for name, cfg in PLANS.items():
        fast = fp.FastPathSolver(p, A, gamma, bounds=bounds, groups=groups, **cfg)
        rf = fast.run()
        assert rf.records == rs.records and rf.B == rs.B and rf.R == rs.R, (H, case['id'], name)
        row['plans'][name] = plan_row(rf, W, K)
        row['plans'][name]['pair_keys'] = sorted(fast.cache.ranges)
    if ablate:
        row['ablation'] = {}
        for tests in ('none', 'i', 'ii'):
            rf = fp.FastPathSolver(p, A, gamma, bounds=bounds, tests=tests).run()
            assert rf.records == rs.records and rf.B == rs.B and rf.R == rs.R, (H, case['id'], tests)
            pr = fp.predicted_word_ops(rf.cost, W, K)
            row['ablation'][tests] = dict(col=rf.cost['col'], mask=rf.cost['mask'], scalar=rf.cost['scalar'],
                                          labels=rf.cost['labels'], predicted=pr,
                                          gain_compact=pr['per_haplotype'] / pr['class_compact'],
                                          gain_expanded=pr['per_haplotype'] / pr['class_expanded'])
    return row


def med(xs):
    return statistics.median(xs) if xs else None


def rng3(xs):
    return [round(min(xs), 4), round(med(xs), 4), round(max(xs), 4)] if xs else None


def summarize(rows, union_labels):
    out = {}
    for H in sorted({r['H'] for r in rows}):
        for gamma in sorted({r['gamma'] for r in rows if r['H'] == H}):
            rs = [r for r in rows if r['H'] == H and r['gamma'] == gamma]
            per = lambda f: [f(r) / r['N'] for r in rs]  # noqa: E731
            pc = lambda r, k: r['plans']['prompt']['cost'].get(k, 0)  # noqa: E731
            fc = lambda r, k: r['plans']['full']['cost'].get(k, 0)  # noqa: E731
            share = lambda num, den, get=pc: [get(r, num) / max(1, get(r, den)) for r in rs]  # noqa: E731
            coded = lambda lv: [pc(r, 'coded_level%d' % lv) / max(1, r['spec']['coded_merges']) for r in rs]  # noqa: E731
            item = dict(
                reads=len(rs), anchors=sum(r['N'] for r in rs), states=sum(r['K'] for r in rs),
                W=rs[0]['W'], b=rs[0]['b'],
                states_per_anchor=rng3(per(lambda r: r['K'])),
                missing_per_anchor=rng3(per(lambda r: r['missing_incidences'])),
                spec=dict(
                    coded_merge_share=rng3([r['spec']['coded_merges'] / r['spec']['merges'] for r in rs]),
                    coded_columns_per_read=rng3([r['spec']['coded_columns'] for r in rs]),
                    coded_columns_per_anchor=rng3(per(lambda r: r['spec']['coded_columns'])),
                    output_classes_per_anchor=rng3(per(lambda r: r['spec']['output_classes'])),
                    max_node_classes=rng3([r['spec']['max_node_classes'] for r in rs])),
                prompt_plan=dict(
                    coded_merges_decided=dict(
                        same_vertex=rng3(coded(1)), loose=rng3(coded(2)), pair_range=rng3(coded(3)),
                        undecided=rng3(coded(4)),
                        by_loose_or_better=rng3([(pc(r, 'coded_level1') + pc(r, 'coded_level2'))
                                                 / max(1, r['spec']['coded_merges']) for r in rs])),
                    contested_merges_per_anchor=rng3(per(lambda r: pc(r, 'merges_contested'))),
                    contested_merges_decided=dict(
                        same_vertex=rng3(share('merges_level1', 'merges_contested')),
                        loose=rng3(share('merges_level2', 'merges_contested')),
                        pair_range=rng3(share('merges_level3', 'merges_contested')),
                        undecided=rng3(share('merges_level4', 'merges_contested'))),
                    pairs_per_anchor=rng3(per(lambda r: pc(r, 'pairs'))),
                    pairs_decided=dict(
                        same_vertex=rng3(share('pairs_level1', 'pairs')),
                        loose=rng3(share('pairs_level2', 'pairs')),
                        pair_range=rng3(share('pairs_level3', 'pairs')),
                        undecided=rng3(share('pairs_level4', 'pairs'))),
                    floor_sliced_share=rng3(share('floor_sliced', 'floor_tests')),
                    switch_sliced_share=rng3(share('switch_sliced', 'switch_tests')),
                    min_sliced_share=rng3(share('min_sliced', 'min_tests')),
                    min_pruned_share=rng3(share('min_pruned', 'min_tests')),
                    columns_per_read=rng3([pc(r, 'col') for r in rs]),
                    columns_per_anchor=rng3(per(lambda r: pc(r, 'col'))),
                    columns_by_source_total={k: sum(pc(r, k) for r in rs)
                                             for k in ('col_merge', 'col_floor', 'col_switch', 'col_min')},
                    column_reduction_vs_spec_merges=rng3([r['spec']['coded_columns'] / max(1, pc(r, 'col'))
                                                          for r in rs]),
                    distinct_origin_vertex_pairs_in_merges_per_read=rng3([pc(r, 'distinct_merge_pairs') for r in rs]),
                    distinct_pairs_per_read=rng3([pc(r, 'cached_ranges') for r in rs]),
                    remote_queries_per_anchor=rng3(per(lambda r: pc(r, 'remote_queries'))),
                    **work_block(rs, 'prompt')),
                full_plan=dict(
                    remote_groups_per_anchor=rng3(per(lambda r: fc(r, 'remote_groups'))),
                    remote_candidates_single_share=rng3(share('remote_single', 'remote_groups', fc)),
                    remote_sliced_per_anchor=rng3(per(lambda r: fc(r, 'remote_sliced'))),
                    prefix_records_cut_share=rng3(share('cut_records', 'prefix_records', fc)),
                    contested_merges_per_anchor=rng3(per(lambda r: fc(r, 'merges_contested'))),
                    columns_per_anchor=rng3(per(lambda r: fc(r, 'col'))),
                    columns_by_source_total={k: sum(fc(r, k) for r in rs)
                                             for k in ('col_merge', 'col_floor', 'col_switch', 'col_min', 'col_remote')},
                    distinct_pairs_per_read=rng3([fc(r, 'cached_ranges') for r in rs]),
                    pair_ranges_union_over_reads=union_labels[(H, gamma)][0],
                    pair_range_labels_union_over_reads=union_labels[(H, gamma)][1],
                    pair_range_labels_sum_over_reads=sum(fc(r, 'labels') for r in rs),
                    **work_block(rs, 'full')))
            abl = [r for r in rs if 'ablation' in r]
            if abl:
                item['ablation'] = {t: dict(gain_compact=rng3([r['ablation'][t]['gain_compact'] for r in abl]),
                                            gain_expanded=rng3([r['ablation'][t]['gain_expanded'] for r in abl]),
                                            columns_per_anchor=rng3([r['ablation'][t]['col'] / r['N'] for r in abl]))
                                    for t in ('none', 'i', 'ii')}
            out[f'H{H}_gamma{gamma}'] = item
    return out


def work_block(rs, plan):
    per = lambda f: [f(r) / r['N'] for r in rs]  # noqa: E731
    c = lambda r, k: r['plans'][plan]['cost'].get(k, 0)  # noqa: E731
    g = lambda k: [r['plans'][plan][k] for r in rs if k in r['plans'][plan]]  # noqa: E731
    tot = lambda k: sum(r['plans'][plan]['predicted'][k] for r in rs)  # noqa: E731
    phases = ('remote', 'local', 'floor', 'output', 'update')
    block = dict(
        mask_ops_per_anchor=rng3(per(lambda r: c(r, 'mask'))),
        mask_ops_by_phase_per_anchor={ph: rng3(per(lambda r, ph=ph: c(r, 'mask_' + ph))) for ph in phases},
        scalar_ops_per_anchor=rng3(per(lambda r: c(r, 'scalar'))),
        scalar_ops_by_phase_per_anchor={ph: rng3(per(lambda r, ph=ph: c(r, 'scalar_' + ph))) for ph in phases},
        range_labels_per_anchor=rng3(per(lambda r: c(r, 'labels'))),
        predicted_per_anchor={k: rng3(per(lambda r, k=k: r['plans'][plan]['predicted'][k]))
                              for k in ('class_compact', 'class_expanded', 'per_haplotype')},
        per_haplotype_steps_per_anchor={k: rng3(per(lambda r, k=k: c(r, 'ph_' + k))) for k in ('local', 'update', 'remote')},
        predicted_gain_compact=rng3(g('gain_compact')),
        predicted_gain_expanded=rng3(g('gain_expanded')),
        predicted_gain_compact_warm=rng3(g('gain_compact_warm')),
        predicted_gain_expanded_warm=rng3(g('gain_expanded_warm')),
        predicted_gain_total_compact=tot('per_haplotype') / tot('class_compact'),
        predicted_gain_total_expanded=tot('per_haplotype') / tot('class_expanded'))
    if g('gain_compact_W16'):
        block['predicted_gain_compact_W16'] = rng3(g('gain_compact_W16'))
        block['predicted_gain_expanded_W16'] = rng3(g('gain_expanded_W16'))
    return block


def union_pair_labels(rows):
    """Distinct pair ranges over all reads of a level, and the labels to build them once."""
    out = {}
    for H in sorted({r['H'] for r in rows}):
        p = graph(H)[0]
        for gamma in sorted({r['gamma'] for r in rows if r['H'] == H}):
            keys = set()
            for r in rows:
                if r['H'] == H and r['gamma'] == gamma:
                    keys.update(map(tuple, r['plans']['full']['pair_keys']))
            out[(H, gamma)] = (len(keys), sum((p.A[u] & p.A[w]).bit_count() for u, w in keys))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, required=True, help='JSONL rows; a .summary.json is written beside it')
    ap.add_argument('--H', type=lambda s: [int(x) for x in s.split(',')], default=list(cm.H_LEVELS))
    ap.add_argument('--gamma', type=lambda s: [int(x) for x in s.split(',')], default=[cm.GAMMA])
    ap.add_argument('--ablate', type=lambda s: [int(x) for x in s.split(',')], default=[60, 960],
                    help='H levels that also run the tests=none/i/ii ablations')
    ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--workers', type=int, default=4)
    args = ap.parse_args()
    summary_path = args.out.with_suffix('.summary.json')
    if args.out.exists() or summary_path.exists():
        raise SystemExit('refusing to overwrite')
    start = time.perf_counter()
    n = len(cm.primary_cases()) if not args.limit else args.limit
    rows = []
    with args.out.open('w') as f, get_context('spawn').Pool(args.workers) as pool:
        for H in args.H:
            for gamma in args.gamma:
                jobs = [(H, k, gamma, H in args.ablate) for k in range(n)]
                for row in pool.imap(one, jobs):
                    rows.append(row)
                    f.write(json.dumps(row) + '\n')
                    f.flush()
                print(f'H={H} gamma={gamma} done after {time.perf_counter() - start:.0f}s', flush=True)
    summary = dict(script='experiments/gap_engine_native/stageA_cost.py', argv=sys.argv[1:], python=sys.version,
                   plans=PLANS, unit_model=fp.__doc__.split('COST MODEL.')[1].strip(),
                   levels=summarize(rows, union_pair_labels(rows)), seconds=round(time.perf_counter() - start, 1))
    summary_path.write_text(json.dumps(summary, indent=1) + '\n')
    for key, item in summary['levels'].items():
        for plan in ('prompt_plan', 'full_plan'):
            print(key, plan, 'gain compact', item[plan]['predicted_gain_compact'],
                  'expanded', item[plan]['predicted_gain_expanded'], 'cols/anchor', item[plan]['columns_per_anchor'])


if __name__ == '__main__':
    main()
