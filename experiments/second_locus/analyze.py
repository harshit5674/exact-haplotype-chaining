"""Assemble the report tables from this experiment's outputs (no solver runs).

Reads survey.json, data/<locus>/panels.json, simulation.json, every panel's conversion.json
and instances.jsonl, verify_<locus>.json, counters_summary.json and realism_check.json, and
writes analysis.json plus markdown tables on stdout (pasted into REPORT.md).
"""
from __future__ import annotations
import argparse
import collections
import json
from pathlib import Path
import statistics
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import locus_common as lc  # noqa: E402

LOCI = ('kir', 'chr20')


def load(path):
    return json.loads(path.read_text()) if path.exists() else None


def survey_table(out):
    s = load(HERE / 'survey.json')
    rows = []
    print('\n### Survey\n')
    print('| Candidate | Window | CHM13 window | Bubbles | Bubbles/Mb | Kept / 464 | Excluded (reason: haplotypes) |')
    print('| --- | --- | --- | ---: | ---: | ---: | --- |')
    for r in s['results']:
        exc = ', '.join(f'{k}: {v}' for k, v in r.get('excluded', {}).items()) or 'none'
        row = dict(name=r['name'], window=r['window'], chm13=r['chm13'], bubbles=r['bubbles'],
                   bubbles_per_mb=r.get('bubbles_per_mb'), kept=r.get('kept'), excluded=r.get('excluded'),
                   with_reverse_alleles=r.get('with_reverse_alleles'), multi_contig=r.get('multi_contig'))
        rows.append(row)
        print(f"| {r['name']} | {r['window']} | {r['chrom']}:{r['chm13'][0]:,}-{r['chm13'][1]:,} | {r['bubbles']} | "
              f"{r.get('bubbles_per_mb')} | {r.get('kept')} | {exc} |")
    out['survey'] = rows


def panel_table(out):
    print('\n### Panels\n')
    print('| Locus | Panel | H | Segments | Edges | Segment bases | Haplotype length (min-max) | Anchors (100 reads) | '
          'Primary N (min / median / max) | States/anchor (primary) | Primary states |')
    print('| --- | --- | ---: | ---: | ---: | ---: | --- | ---: | --- | ---: | ---: |')
    out['panels'] = {}
    for locus in LOCI:
        meta = load(lc.DATA / locus / 'panels.json')
        # N and K of the 40 primary instances, from the counter rows (gamma 10000), not re-parsing instances
        crow = [json.loads(s) for s in (HERE / f'counters_{locus}.jsonl').read_text().splitlines()]
        for p in meta['panels']:
            d = lc.DATA / locus / p['panel']
            conv = load(d / 'conversion.json')
            rs = [r for r in crow if r['panel'] == p['panel'] and r['gamma'] == lc.GAMMA and r['primary']]
            assert len(rs) == 40
            Ns, Ks = [r['N'] for r in rs], [r['K'] for r in rs]
            item = dict(H=p['H'], segments=p['segments'], edges=p['edges'], segment_bases=p['total_segment_length'],
                        haplotype_length=[p['haplotype_length_min'], p['haplotype_length_max']], anchors=conv['anchors'],
                        instances=conv['orientation_instances'], primary_N=[min(Ns), statistics.median(Ns), max(Ns)],
                        primary_states=sum(Ks), states_per_anchor=round(sum(Ks) / sum(Ns), 1))
            out['panels'][f'{locus}/{p["panel"]}'] = item
            print(f"| {locus} | {p['panel']} | {p['H']} | {p['segments']:,} | {p['edges']:,} | {p['total_segment_length']:,} | "
                  f"{p['haplotype_length_min']:,}-{p['haplotype_length_max']:,} | {conv['anchors']:,} | "
                  f"{min(Ns):,} / {statistics.median(Ns):,.0f} / {max(Ns):,} | {item['states_per_anchor']} | {sum(Ks):,} |")


def verify_table(out):
    print('\n### Verification\n')
    print('| Locus | Panel | H | gamma | Instances | Passed | Anchors | States | Pointer kinds (restart / local / remote / same) |')
    print('| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |')
    out['verification'] = {}
    for locus in LOCI:
        v = load(HERE / f'verify_{locus}.json')
        if v is None:
            print(f'| {locus} | (not run) | | | | | | | |')
            continue
        out['verification'][locus] = dict(status=v['status'], failures=v['failures'], instances=v['instances'],
                                          states=v['states'], seconds=v['seconds'], levels=[])
        for lv in v['levels']:
            kinds = collections.Counter()
            groups = collections.defaultdict(collections.Counter)
            for r in lv['rows']:
                pc = r.get('pointer_counts', {}).get('class_compact', {})
                kinds.update(pc)
                groups[r['group']].update(pc)
            out['verification'][locus]['levels'].append(dict(panel=lv['panel'], H=lv['H'], gamma=lv['gamma'],
                                                             instances=lv['reads'], passed=lv['passed'], states=lv['states'],
                                                             anchors=lv['anchors'], pointer_kinds=dict(kinds),
                                                             pointer_kinds_by_group={g: dict(c) for g, c in groups.items()},
                                                             failures=lv['failures']))
            print(f"| {locus} | {lv['panel']} | {lv['H']} | {lv['gamma']} | {lv['reads']} | {lv['passed']} | {lv['anchors']:,} | "
                  f"{lv['states']:,} | {kinds['restart']:,} / {kinds['local']:,} / {kinds['remote']:,} / {kinds['same']:,} |")
        print(f"| {locus} | **total** | | | {v['instances']} | {v['instances'] - v['failures']} | | {v['states']:,} | {v['status']} |")


def counters_table(out):
    c = load(HERE / 'counters_summary.json')
    if c is None:
        return
    out['counters'] = c['levels']
    for gamma in (10000, 20000):
        print(f'\n### Counters, 40 primary instances, gamma {gamma} (median over instances [min, max] where shown)\n')
        print('| Source | Panel | H | States/anchor | Same vertex | Loose | Range | Sliced | Sliced pairs/anchor | '
              'Sliced thresholds/anchor | Sliced minima/anchor | Published/anchor | Output records/anchor | Longest node list |')
        print('| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | --- | ---: | ---: | --- |')
        for key, e in c['levels'].items():
            if e['gamma'] != gamma:
                continue
            s = e['primary']
            sh, ca = s['pair_decision_shares'], s['counters_per_anchor']
            f3 = lambda x: f'{x[1]} [{x[0]}, {x[2]}]'  # noqa: E731
            print(f"| {e['source']} | {e['panel']} | {e['H']} | {s['states_per_anchor']} | {sh['pairs_same_vertex'][1]:.3f} | "
                  f"{sh['pairs_loose'][1]:.3f} | {sh['pairs_range'][1]:.3f} | {sh['pairs_sliced'][1]:.3f} | {f3(ca['pairs_sliced'])} | "
                  f"{f3(ca['sliced_thresholds'])} | {f3(ca['sliced_minima'])} | {ca['published_records'][1]} | "
                  f"{ca['output_records'][1]} | {f3(s['max_node_records'])} |")
    print('\n### Counters by read source group (all 200 instances, gamma 10000; medians)\n')
    print('| Source | Panel | H | Group | Instances | Sliced pairs/anchor | Sliced thresholds/anchor | Sliced minima/anchor | '
          'Output records/anchor | Pooled sliced pairs/anchor | Class = per-haplotype |')
    print('| --- | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |')
    for key, e in c['levels'].items():
        if e['gamma'] != 10000 or 'all' not in e:
            continue
        for g in ('absent', 'present'):
            s = e.get(f'all_{g}')
            if not s:
                continue
            ca = s['counters_per_anchor']
            print(f"| {e['source']} | {e['panel']} | {e['H']} | {g} | {s['instances']} | {ca['pairs_sliced'][1]} | "
                  f"{ca['sliced_thresholds'][1]} | {ca['sliced_minima'][1]} | {ca['output_records'][1]} | "
                  f"{s['pooled_per_anchor']['pairs_sliced']} | {s['class_equals_per_haplotype']}/{s['instances']} |")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, default=HERE / 'analysis.json')
    args = ap.parse_args()
    if args.out.exists():
        raise SystemExit('refusing to overwrite ' + str(args.out))
    out = {}
    survey_table(out)
    panel_table(out)
    verify_table(out)
    counters_table(out)
    r = load(HERE / 'realism_check.json')
    if r:
        out['realism'] = r['panels']
    lc.dump(args.out, dict(script='experiments/second_locus/analyze.py', **out))


if __name__ == '__main__':
    main()
