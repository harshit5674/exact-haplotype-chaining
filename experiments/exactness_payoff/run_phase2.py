"""Phase 2: optimality gap of Giraffe's single-record chaining DP on real anchors.

For every (read, strand) instance of a 60-haplotype panel and each penalty rho, run the faithful
emulation of vg v1.77.0's DP (consistency bonus 0 and 12) and the exact per-haplotype DP over the
same transitions, items, gap terms and masks (SPEC_GIRAFFE.md). Every reported chain is re-scored
independently in Python (own shortest distances, own transition rule, own greedy and label DP).

  python3.13 experiments/exactness_payoff/run_phase2.py --dataset hprc60 --out .../phase2_hprc60.jsonl
  python3.13 experiments/exactness_payoff/run_phase2.py --counterexample --out .../phase2_counterexample.json
"""
from __future__ import annotations

import argparse
import json
from multiprocessing import Pool
import os
from pathlib import Path
import platform
import random
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import emul  # noqa: E402

ROOT = HERE.parents[1]
DATASETS = {
    'hprc60': ROOT / 'experiments/hprc_mhc_panel/data/hprc60',
    'mhc60': ROOT / 'experiments/mhc_pilot/data',
}
# Leave-one-out panels: the same graph minus the reads' own haplotype (CHM13), pruned to the
# vertices and edges of the remaining 59 haplotypes; anchors on dropped vertices are removed.
LEAVE_OUT = {'hprc60_noCHM13': ('hprc60', 'CHM13#0|1|MHC'), 'mhc60_noCHM13': ('mhc60', 'CHM13#0|1|MHC')}
# Random 16-haplotype subsets without CHM13 (vg's evaluations sample panels to 16 haplotypes; vg's
# sampler is k-mer guided, these are uniform). Seed 1600 + replicate.
SUBSET16 = {f'{base}_s16_r{r}': (base, 16, 1600 + r) for base in ('hprc60', 'mhc60') for r in range(3)}
RHOS = (1, 2, 4, 10, 28)
BONUSES = (0, 12)

_PANEL = None


def make_panel(dataset):
    """Panel for a dataset name and the old->new vertex map (identity for full panels)."""
    if dataset in SUBSET16:
        base, size, seed = SUBSET16[dataset]
        panel = emul.load_panel(DATASETS[base] / 'graph.json', base)
        pool = [h for h, name in enumerate(panel.hap_names) if not name.startswith('CHM13')]
        keep = sorted(random.Random(seed).sample(pool, size))
        return emul.subset_panel(panel, keep, dataset)
    base, drop = LEAVE_OUT.get(dataset, (dataset, None))
    panel = emul.load_panel(DATASETS[base] / 'graph.json', base)
    if drop is None:
        return panel, list(range(panel.V))
    keep = [h for h, name in enumerate(panel.hap_names) if name != drop]
    assert len(keep) == panel.H - 1, f'{drop} not found exactly once'
    return emul.subset_panel(panel, keep, dataset)


def map_case(case, vmap):
    anchors = [[vmap[a[0]], *a[1:]] for a in case['anchors'] if vmap[a[0]] >= 0]
    return dict(case, anchors=anchors, dropped_anchors=len(case['anchors']) - len(anchors))


def _init(dataset):
    global _PANEL
    _PANEL, _ = make_panel(dataset)
    _PANEL.handle  # build the native index once per worker


def verify_heuristic_chain(panel, r, b, params):
    """Independent re-scoring of the traced heuristic chain."""
    chain = emul.heuristic_chain(r, b)
    end = chain[-1]
    A = r.anchors
    j0s = []
    for s, a in zip(chain, chain[1:]):
        j0 = emul.py_transition(panel, A[s], A[a], params['graph_lookback'], params['indel_limit'])
        assert j0 is not None, 'heuristic chain uses an illegal transition'
        j0s.append(j0)
    masks = [panel.A[A[k][0]] for k in chain]
    resets, _ = emul.greedy_resets(masks)
    minsw, _ = emul.min_switch_labels(masks, panel.H)
    assert resets == minsw, 'greedy resets differ from the minimum label changes'
    assert resets == r.hrec[b * r.n + end]
    value = sum(r.item[k] for k in chain) + sum(j0s) - r.rho * resets
    assert value == r.hs[b * r.n + end], (value, r.hs[b * r.n + end])
    return chain, resets


def verify_exact_chain(panel, r, params):
    ch = emul.exact_chain(r)
    A = r.anchors
    j0s = []
    for (s, _), (a, _) in zip(ch, ch[1:]):
        j0 = emul.py_transition(panel, A[s], A[a], params['graph_lookback'], params['indel_limit'])
        assert j0 is not None, 'exact chain uses an illegal transition'
        j0s.append(j0)
    for i, h in ch:
        assert panel.A[A[i][0]] >> h & 1, 'label not on the anchor vertex'
    switches = sum(1 for (_, g), (_, h) in zip(ch, ch[1:]) if g != h)
    value = sum(r.item[i] for i, _ in ch) + sum(j0s) - r.rho * switches
    opt = max(r.B)
    assert value == opt, (value, opt)
    masks = [panel.A[A[i][0]] for i, _ in ch]
    minsw, _ = emul.min_switch_labels(masks, panel.H)
    assert minsw == switches, 'exact chain labeling is not switch-minimal for its chain'
    return ch, switches


def run_case(case):
    panel = _PANEL
    rows = []
    for rho in case['_rhos']:
        t = time.perf_counter()
        r = emul.solve(panel, case['anchors'], case['qlen'], rho=rho, bonuses=BONUSES, tie='graph',
                       prune=case['_prune'])
        dt = time.perf_counter() - t
        opt = max(r.B)
        n = r.n
        row = dict(id=case['id'], read=case['read'], strand=case['strand'], qlen=case['qlen'], rho=rho, N=n,
                   dropped_anchors=case.get('dropped_anchors', 0), transitions=r.stats['transitions'],
                   seconds=round(dt, 4), opt=opt)
        ech, esw = verify_exact_chain(panel, r, r.params)
        row.update(exact_chain_len=len(ech), exact_switches=esw,
                   exact_span=[r.anchors[ech[0][0]][3], r.anchors[ech[-1][0]][4]])
        for b, cb in enumerate(BONUSES):
            S = r.hs[b * n:(b + 1) * n]
            below = sum(1 for i in range(n) if S[i] < r.B[i])
            assert all(S[i] <= r.B[i] for i in range(n)), 'heuristic exceeds exact on some anchor'
            chain, resets = verify_heuristic_chain(panel, r, b, r.params)
            best = S[chain[-1]]
            row[f'h{cb}'] = dict(best=best, gap=opt - best, rel_gap=(opt - best) / opt if opt else 0.0,
                                 chain_len=len(chain), switches=resets, anchors_below_exact=below,
                                 span=[r.anchors[chain[0]][3], r.anchors[chain[-1]][4]],
                                 end_B=r.B[chain[-1]])
        if rho == 2:
            # Tie-order sensitivity: std::sort leaves identical read intervals unordered.
            rr = emul.solve(panel, case['anchors'], case['qlen'], rho=rho, bonuses=BONUSES, tie='reverse',
                            do_exact=False)
            row['tie_reverse'] = {f'h{cb}': max(rr.hs[b * n:(b + 1) * n]) for b, cb in enumerate(BONUSES)}
        rows.append(row)
    return rows


def counterexample(out: Path):
    """The novelty audit's four anchors, executed through the same native DP core."""
    names = ['A', 'B', 'C', 'D']
    items = [10, 11, 10, 10]
    masks = [0b01, 0b10, 0b11, 0b01]  # h1 = bit 0, h2 = bit 1
    trans = [(0, 2, 0), (1, 2, 0), (2, 3, 0)]
    res = {}
    for rho in (5, 2):
        r = emul.solve_explicit(4, 2, items, masks, trans, rho=rho, bonuses=(0, 12))
        entry = dict(exact_B=list(r.B), exact_opt=max(r.B),
                     exact_chain=[[names[i], f'h{h + 1}'] for i, h in emul.exact_chain(r)])
        for b, cb in enumerate((0, 12)):
            ch = emul.heuristic_chain(r, b)
            entry[f'bonus{cb}'] = dict(scores=list(r.hs[b * 4:(b + 1) * 4]),
                                       sources=[names[s] if s >= 0 else None for s in r.hsrc[b * 4:(b + 1) * 4]],
                                       best=r.hs[b * 4 + ch[-1]], chain=[names[i] for i in ch],
                                       recombinations=r.hrec[b * 4 + ch[-1]])
        res[f'rho{rho}'] = entry
    doc = dict(description='Novelty-audit counterexample: A{h1} 10, B{h2} 11, C{h1,h2} 10, D{h1} 10; '
                           'transitions A->C, B->C, C->D with zero gap cost', results=res,
               library=emul.build().identity)
    out.write_text(json.dumps(doc, indent=1) + '\n')
    print(json.dumps(res, indent=1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dataset', choices=sorted(DATASETS) + sorted(LEAVE_OUT) + sorted(SUBSET16))
    ap.add_argument('--out', required=True)
    ap.add_argument('--rhos', default=','.join(map(str, RHOS)))
    ap.add_argument('--workers', type=int, default=3)
    ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--counterexample', action='store_true')
    ap.add_argument('--no-prune', action='store_true',
                    help='exact DP without source pruning (identical values; the engine of the first runs)')
    args = ap.parse_args()
    out = Path(args.out)
    if out.exists():
        raise SystemExit(f'{out} exists; choose a fresh output path')
    if args.counterexample:
        counterexample(out)
        return
    rhos = [int(x) for x in args.rhos.split(',')]
    base = SUBSET16[args.dataset][0] if args.dataset in SUBSET16 else LEAVE_OUT.get(args.dataset, (args.dataset, None))[0]
    lines = (DATASETS[base] / 'instances.jsonl').read_text().splitlines()
    panel, vmap = make_panel(args.dataset)
    cases = [map_case(json.loads(s), vmap) for s in lines]
    cases = [c for c in cases if c['anchors']]
    if args.limit:
        cases = cases[:args.limit]
    for c in cases:
        c['_rhos'] = rhos
        c['_prune'] = not args.no_prune
    # Largest first, so the pool's tail is short.
    cases.sort(key=lambda c: -len(c['anchors']))
    t0 = time.perf_counter()
    workers = min(args.workers, 3)
    with Pool(workers, initializer=_init, initargs=(args.dataset,)) as pool, out.open('x') as f:
        for k, rows in enumerate(pool.imap_unordered(run_case, cases, chunksize=1)):
            for row in rows:
                f.write(json.dumps(row) + '\n')
            f.flush()
            if (k + 1) % 20 == 0:
                print(f'{k + 1}/{len(cases)} instances, {time.perf_counter() - t0:.0f} s', flush=True)
    meta = dict(dataset=args.dataset, source=str(DATASETS[base].relative_to(ROOT)), H=panel.H, V=panel.V,
                haplotypes=panel.hap_names,
                edges=len(panel.edges), instances=len(cases), dropped_anchors=sum(c['dropped_anchors'] for c in cases),
                rhos=rhos, bonuses=list(BONUSES), prune=not args.no_prune, workers=workers,
                seconds=round(time.perf_counter() - t0, 1),
                library=emul.build().identity, python=platform.python_version(), host=platform.machine(),
                argv=sys.argv, assertions=__debug__)
    Path(str(out) + '.meta.json').write_text(json.dumps(meta, indent=1) + '\n')
    print(json.dumps(meta))


if __name__ == '__main__':
    main()
