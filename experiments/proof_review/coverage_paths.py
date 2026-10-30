"""Which code paths do the test families exercise? (REPORT.md, coverage table)

For a deterministic sample of every family it runs the specification, the fast path (prompt
and full plans) and the native class engine, and sums their branch counters:
  spec    scalar_fast_merges, direct_coded_merges, clustered_merges, threshold_calls,
          minimum_calls, route_calls, missing_incidences;
  fast    pairs_level1..4 (decided by same vertex / loose / exact range / sliced),
          floor_/switch_ sliced|all|none|self, min_single|min_sliced|min_pruned,
          remote_single|remote_sliced|remote_pruned, cut_records;
  native  pairs_sliced, sliced_thresholds, sliced_minima, remote_sliced, cut_records;
  witness kinds of the specification's states, and labels with undefined last2reach.
It does not check correctness (run_families.py and hostile.py do). Writes
results/coverage.json. Single process; usage: python3.13 -B coverage_paths.py [--per 4000]
"""
from __future__ import annotations
import argparse
import collections
import itertools
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import harness as H  # noqa: E402
import families as fam  # noqa: E402
import hostile as ho  # noqa: E402
gc, fp, gn = H.gc, H.fp, H.gn

SPEC_KEYS = ('scalar_fast_merges', 'direct_coded_merges', 'clustered_merges', 'threshold_calls',
             'minimum_calls', 'route_calls', 'missing_incidences')
FP_KEYS = ('pairs_level1', 'pairs_level2', 'pairs_level3', 'pairs_level4', 'floor_sliced', 'floor_all',
           'floor_none', 'switch_sliced', 'switch_all', 'switch_none', 'min_single', 'min_sliced',
           'min_pruned', 'remote_single', 'remote_sliced', 'remote_pruned', 'cut_records')
NAT_KEYS = ('pairs_sliced', 'sliced_thresholds', 'sliced_minima', 'remote_sliced', 'cut_records', 'pairs')


def samples(per):
    yield 'EXA', itertools.islice(itertools.chain.from_iterable(
        fam.exa_instances(V, 3, fam.ANCHOR_SETS_A, 0, 9) for V in (1, 2, 3)), per)
    yield 'EXA4', itertools.islice(fam.exa_instances(4, 2, fam.ANCHOR_SETS_A4, 0, 61), per)
    yield 'EXB2', itertools.islice(itertools.chain.from_iterable(
        fam.exb2_instances(n, 0, 97) for n in fam.HOSTILE_GRAPHS), per)
    yield 'EXB3', itertools.islice(itertools.chain.from_iterable(
        fam.exb3_instances(n, 0, 37) for n in fam.HOSTILE_GRAPHS_3), per)
    yield 'RS', (fam.random_stress_instance(s) for s in range(per))
    for g, f in ho.GENERATORS.items():
        yield 'H_' + g, (f(s) for s in range(min(per, 600)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--per', type=int, default=3000)
    args = ap.parse_args()
    out = {}
    for family, it in samples(args.per):
        c = collections.Counter()
        cache = {}
        for name, (L, E, P), A, gamma in it:
            key = (tuple(L), tuple(map(tuple, E)), tuple(map(tuple, P)))
            if key not in cache:
                cache.clear()
                cache[key] = H.GraphCache(L, E, P)
            g = cache[key]
            anchors = [gc.Anchor(*a) for a in A]
            p = g.p
            c['instances'] += 1
            for z in anchors:
                for h in range(p.H):
                    if not p.on[z.v][h]:
                        c['missing_labels'] += 1
                        c['missing_labels_without_last2reach'] += p.last[z.v][h] < 0
            rs = gc.solve(p, anchors, gamma)
            for k in SPEC_KEYS:
                c['spec_' + k] += rs.stats.get(k, 0)
            for recs in rs.records:
                for r in recs:
                    c['spec_witness_' + r.witness.kind] += r.mask.bit_count()
            for plan, cfg in (('prompt', H.FP_CONFIGS['fp_prompt']), ('full', H.FP_CONFIGS['fp_full'])):
                rf = fp.FastPathSolver(p, anchors, gamma, **cfg).run()
                for k in FP_KEYS:
                    c['fast_%s_%s' % (plan, k)] += rf.cost.get(k, 0)
            if g.native is not None:
                try:
                    comp = g.native.solve_compact(anchors, gamma)
                    for k in NAT_KEYS:
                        c['native_' + k] += int(comp.stats.get(k, 0))
                except (gn.NativeRangeError, ValueError):
                    c['native_rejected'] += 1
        out[family] = dict(sorted(c.items()))
        print(family, {k: v for k, v in out[family].items() if 'sliced' in k or 'cluster' in k or 'without' in k
                       or k == 'instances'})
    (HERE / 'results').mkdir(exist_ok=True)
    (HERE / 'results/coverage.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
