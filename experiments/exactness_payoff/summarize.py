"""Summaries of the Phase 2 (optimality gap) and Phase 3/4 (simulated reads) outputs.

  python3.13 experiments/exactness_payoff/summarize.py phase2 phase2_*.jsonl --out phase2_summary.json
  python3.13 experiments/exactness_payoff/summarize.py phase3 phase3_hprc60.jsonl --out phase3_summary.json
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import statistics as st


def quantiles(xs):
    if not xs:
        return None
    xs = sorted(xs)
    q = lambda f: xs[min(len(xs) - 1, int(f * (len(xs) - 1) + 0.5))]
    return dict(n=len(xs), min=xs[0], median=q(0.5), p90=q(0.9), max=xs[-1], mean=sum(xs) / len(xs))


# ------------------------------------------------------------ phase 2

def phase2(paths):
    out = {}
    for path in paths:
        rows = [json.loads(s) for s in Path(path).read_text().splitlines()]
        ds = Path(path).stem.replace('phase2_', '')
        by_rho = defaultdict(list)
        for r in rows:
            by_rho[r['rho']].append(r)
        res = {}
        for rho, rs in sorted(by_rho.items()):
            reads = defaultdict(list)
            for r in rs:
                reads[r['read']].append(r)
            entry = dict(instances=len(rs), reads=len(reads),
                         exact_with_switch=sum(1 for r in rs if r['exact_switches'] > 0),
                         anchors=sum(r['N'] for r in rs), transitions=sum(r['transitions'] for r in rs))
            for b in ('h0', 'h12'):
                gaps = [r[b]['gap'] for r in rs]
                pos = [r for r in rs if r[b]['gap'] > 0]
                primary = [max(v, key=lambda x: x['opt']) for v in reads.values()]
                entry[b] = dict(
                    instances_with_gap=len(pos), share_instances=len(pos) / len(rs),
                    reads_primary_with_gap=sum(1 for r in primary if r[b]['gap'] > 0),
                    reads_any_with_gap=sum(1 for v in reads.values() if any(r[b]['gap'] > 0 for r in v)),
                    gap_abs=quantiles([r[b]['gap'] for r in pos]), gap_rel=quantiles([r[b]['rel_gap'] for r in pos]),
                    total_gap=sum(gaps), heuristic_with_switch=sum(1 for r in rs if r[b]['switches'] > 0),
                    anchors_below_exact=sum(r[b]['anchors_below_exact'] for r in rs),
                    instances_with_anchor_below_exact=sum(1 for r in rs if r[b]['anchors_below_exact'] > 0))
            tie = [r for r in rs if 'tie_reverse' in r]
            if tie:
                entry['tie_reverse_changed'] = {b: sum(1 for r in tie if r['tie_reverse'][b] != r[b]['best'])
                                                for b in ('h0', 'h12')}
            res[str(rho)] = entry
        out[ds] = res
    return out


def phase2_table(summary):
    lines = ['| Dataset | ρ | Instances | Exact optimum uses a switch | Gap > 0, bonus 0 | Gap > 0, bonus 12 | '
             'Reads (primary) with gap, bonus 0 / 12 | Anchors with stored < exact, bonus 0 / 12 |',
             '| --- | ---: | ---: | ---: | ---: | ---: | --- | --- |']
    for ds, res in summary.items():
        for rho, e in res.items():
            lines.append(f"| {ds} | {rho} | {e['instances']} | {e['exact_with_switch']} | "
                         f"{e['h0']['instances_with_gap']} | {e['h12']['instances_with_gap']} | "
                         f"{e['h0']['reads_primary_with_gap']} / {e['h12']['reads_primary_with_gap']} | "
                         f"{e['h0']['anchors_below_exact']} / {e['h12']['anchors_below_exact']} of {e['anchors']} |")
    return '\n'.join(lines)


# ------------------------------------------------------------ phase 3

def match(pred, true, tol):
    pairs = sorted((abs(p - t), i, j) for i, p in enumerate(pred) for j, t in enumerate(true) if abs(p - t) <= tol)
    used_p, used_t, m = set(), set(), 0
    for _, i, j in pairs:
        if i not in used_p and j not in used_t:
            used_p.add(i)
            used_t.add(j)
            m += 1
    return m, used_t


def phase3(path, classes=None):
    rows = [json.loads(s) for s in Path(path).read_text().splitlines()]
    methods = sorted({m for r in rows for m in r['methods']})
    out = {}
    for k in ('all', 1, 2, 3):
        rs = [r for r in rows if k == 'all' or r['k'] == k]
        if not rs:
            continue
        cls = {}
        for m in methods:
            ev = [(r, r['methods'][m]) for r in rs if m in r['methods']]
            if not ev:
                continue
            e = dict(reads=len(ev))
            e['span_cov'] = quantiles([x['span_cov'] for _, x in ev])
            e['placed_frac'] = quantiles([x['placed_frac'] for _, x in ev])
            e['chain_len'] = quantiles([x['chain_len'] for _, x in ev])
            e['switches'] = quantiles([x['switches'] for _, x in ev])
            e['switches_equal_truth'] = sum(1 for r, x in ev if x['switches'] == len(r['true_breakpoints']))
            e['switches_equal_truth_min'] = sum(1 for r, x in ev if x['switches'] == r['truth_min_switches'])
            hc = sum(x['hap_correct'] for _, x in ev)
            ht = sum(x['hap_total'] for _, x in ev)
            e['hap_accuracy_pooled'] = hc / ht if ht else None
            e['hap_accuracy_reads'] = quantiles([x['hap_correct'] / x['hap_total'] for _, x in ev if x['hap_total']])
            e['reads_hap_all_correct'] = sum(1 for _, x in ev if x['hap_total'] and x['hap_correct'] == x['hap_total'])
            if all(x.get('walk_correct') is not None for _, x in ev):
                wc = sum(x['walk_correct'] for _, x in ev)
                e['walk_accuracy_pooled'] = wc / ht if ht else None
            e['mean_set_size'] = quantiles([x['mean_set_size'] for _, x in ev if x['mean_set_size'] is not None])
            if 'gap' in ev[0][1]:
                e['gap_reads'] = sum(1 for _, x in ev if x['gap'] > 0)
                e['gap_abs'] = quantiles([x['gap'] for _, x in ev if x['gap'] > 0])
                e['gap_rel'] = quantiles([x['gap'] / (x['score'] + x['gap']) for _, x in ev if x['gap'] > 0])
            for tol in (1000, 5000):
                tp = npred = ntrue = tpo = ntrue_o = 0
                for r, x in ev:
                    mm, used = match(x['breakpoints'], r['true_breakpoints'], tol)
                    tp += mm
                    npred += len(x['breakpoints'])
                    ntrue += len(r['true_breakpoints'])
                    for j, o in enumerate(r['observable_5kb']):
                        if o:
                            ntrue_o += 1
                            tpo += j in used
                e[f'bp_{tol}'] = dict(true=ntrue, predicted=npred, matched=tp,
                                      recall=tp / ntrue if ntrue else None, precision=tp / npred if npred else None,
                                      observable_true=ntrue_o, observable_matched=tpo,
                                      observable_recall=tpo / ntrue_o if ntrue_o else None)
            cls[m] = e
        out[str(k)] = dict(reads=len(rs), true_breakpoints=sum(len(r['true_breakpoints']) for r in rs),
                           observable_breakpoints=sum(sum(r['observable_5kb']) for r in rs),
                           truth_needs_switch=sum(1 for r in rs if r['truth_min_switches'] > 0), methods=cls)
    return out


def phase3_paired(path):
    """Exact vs Giraffe per read at the same rho: suboptimal reads vs equal-score reads."""
    rows = [json.loads(s) for s in Path(path).read_text().splitlines()]
    out = {}
    rhos = sorted({int(m.split('rho')[1]) for r in rows for m in r['methods'] if m.startswith('exact_rho')})
    for rho in rhos:
        for b in (0, 12):
            ex, gi = f'exact_rho{rho}', f'giraffe_b{b}_rho{rho}'
            pairs = [(r, r['methods'][ex], r['methods'][gi]) for r in rows if ex in r['methods'] and gi in r['methods']]
            if not pairs:
                continue
            groups = dict(suboptimal=[p for p in pairs if p[2]['gap'] > 0],
                          equal_score=[p for p in pairs if p[2]['gap'] == 0])
            e = {}
            for gname, ps in groups.items():
                acc = lambda x: x['hap_correct'] / x['hap_total'] if x['hap_total'] else None
                ae = [acc(x) for _, x, _ in ps if acc(x) is not None]
                ag = [acc(y) for _, _, y in ps if acc(y) is not None]
                diff = [acc(x) - acc(y) for _, x, y in ps if acc(x) is not None and acc(y) is not None]
                bp = {}
                for tol in (1000, 5000):
                    me = sum(match(x['breakpoints'], r['true_breakpoints'], tol)[0] for r, x, _ in ps)
                    mg = sum(match(y['breakpoints'], r['true_breakpoints'], tol)[0] for r, _, y in ps)
                    bp[tol] = dict(exact_matched=me, giraffe_matched=mg,
                                   exact_predicted=sum(len(x['breakpoints']) for _, x, _ in ps),
                                   giraffe_predicted=sum(len(y['breakpoints']) for _, _, y in ps),
                                   true=sum(len(r['true_breakpoints']) for r, _, _ in ps))
                e[gname] = dict(reads=len(ps), exact_hap_mean=sum(ae) / len(ae) if ae else None,
                                giraffe_hap_mean=sum(ag) / len(ag) if ag else None,
                                reads_exact_better=sum(1 for d in diff if d > 1e-12),
                                reads_giraffe_better=sum(1 for d in diff if d < -1e-12),
                                acc_diff=quantiles([d for d in diff if abs(d) > 1e-12]),
                                gap=quantiles([y['gap'] for _, _, y in ps]) if gname == 'suboptimal' else None,
                                switches_exact=sum(x['switches'] for _, x, _ in ps),
                                switches_giraffe=sum(y['switches'] for _, _, y in ps), breakpoints=bp)
            out[f'rho{rho}_b{b}'] = e
    return out


def phase4(full_path, subset_paths, subsets_json):
    """Pair each subset run with the full-panel run on the same reads."""
    full = {json.loads(s)['name']: json.loads(s) for s in Path(full_path).read_text().splitlines()}
    subsets = {x['name']: x for x in json.loads(Path(subsets_json).read_text())['subsets']}
    out = dict(full=phase3(full_path))
    for path in subset_paths:
        name = Path(path).stem.replace('phase4_', '')
        rows = [json.loads(s) for s in Path(path).read_text().splitlines()]
        e = dict(size=subsets[name]['size'], reads=len(rows), summary=phase3(path))
        e['reads_all_sources_in_panel'] = sum(1 for r in rows if all(r['true_haps_in_panel']))
        e['reads_no_source_in_panel'] = sum(1 for r in rows if not any(r['true_haps_in_panel']))
        e['dropped_anchors'] = quantiles([r['dropped'] / (r['N'] + r['dropped']) for r in rows if r['N'] + r['dropped']])
        e['truth_min_switches_increase'] = quantiles([r['truth_min_switches'] - full[r['name']]['truth_min_switches']
                                                      for r in rows])
        paired = {}
        for m in sorted({m for r in rows for m in r['methods']}):
            if not m.startswith(('exact', 'cgj')):
                continue
            drop, rel, acc_full, acc_sub, cov_full, cov_sub, wf, ws = [], [], [], [], [], [], [], []
            for r in rows:
                f = full[r['name']]['methods'].get(m)
                s = r['methods'].get(m)
                if not f or not s:
                    continue
                drop.append(f['score'] - s['score'])
                rel.append((f['score'] - s['score']) / f['score'] if f['score'] else 0.0)
                if f['hap_total'] and s['hap_total']:
                    acc_full.append(f['hap_correct'] / f['hap_total'])
                    acc_sub.append(s['hap_correct'] / s['hap_total'])
                    if f.get('walk_correct') is not None and s.get('walk_correct') is not None:
                        wf.append(f['walk_correct'] / f['hap_total'])
                        ws.append(s['walk_correct'] / s['hap_total'])
                cov_full.append(f['span_cov'])
                cov_sub.append(s['span_cov'])
            paired[m] = dict(score_drop=quantiles(drop), rel_score_drop=quantiles(rel),
                             reads_with_drop=sum(1 for x in drop if x > 0),
                             hap_acc_full_mean=sum(acc_full) / len(acc_full) if acc_full else None,
                             hap_acc_subset_mean=sum(acc_sub) / len(acc_sub) if acc_sub else None,
                             walk_acc_full_mean=sum(wf) / len(wf) if wf else None,
                             walk_acc_subset_mean=sum(ws) / len(ws) if ws else None,
                             reads_walk_worse=sum(1 for a, b in zip(wf, ws) if b < a - 1e-12),
                             reads_walk_better=sum(1 for a, b in zip(wf, ws) if b > a + 1e-12),
                             span_cov_full_mean=sum(cov_full) / len(cov_full) if cov_full else None,
                             span_cov_subset_mean=sum(cov_sub) / len(cov_sub) if cov_sub else None)
        e['paired_vs_full'] = paired
        out[name] = e
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('kind', choices=('phase2', 'phase3', 'phase4'))
    ap.add_argument('--subsets', help='phase4: phase4_subsets.json')
    ap.add_argument('inputs', nargs='+')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    if args.kind == 'phase2':
        s = phase2(args.inputs)
        Path(args.out).write_text(json.dumps(s, indent=1) + '\n')
        print(phase2_table(s))
    elif args.kind == 'phase4':
        s = phase4(args.inputs[0], args.inputs[1:], args.subsets)
        Path(args.out).write_text(json.dumps(s, indent=1) + '\n')
        for name, e in s.items():
            if name == 'full':
                continue
            print(name, e['reads_all_sources_in_panel'], e['reads_no_source_in_panel'],
                  json.dumps({m: (v['reads_with_drop'], v['hap_acc_full_mean'], v['hap_acc_subset_mean'])
                              for m, v in e['paired_vs_full'].items()}))
    else:
        s = {Path(p).stem: dict(phase3(p), paired=phase3_paired(p)) for p in args.inputs}
        Path(args.out).write_text(json.dumps(s, indent=1) + '\n')
        for name, res in s.items():
            print(name, json.dumps({k: {m: (e['hap_accuracy_pooled'], e['bp_5000']['recall'], e['bp_5000']['precision'])
                                        for m, e in v['methods'].items()} for k, v in res.items() if k != 'paired'},
                                   indent=1))


if __name__ == '__main__':
    main()
