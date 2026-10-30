"""Haplotype switches on one optimal chain per instance: does the leave-out design bite?

For each primary instance (gamma 10000, weights x 200) on every panel of a locus, this runs the
native class engine (solve_packed mode='class'), takes the first anchor attaining max B and its
first state attaining B, and follows the class pointers to a restart. It counts the pointer
kinds on that chain (local; remote and same-anchor are haplotype switches) and, for simulated
reads, whether the read's source haplotype alone attains the optimum (the best state value on
the source label equals max B; tie-independent, unlike the labels on the traced chain, which
are the lowest-index tied labels). Ties: other optimal chains may switch differently; the
switch count is for one legal optimal chain, not a property of all of them.
The MHC source (real CHM13 reads, earlier experiment's panels) is included for reference.
Output: chain_switches.jsonl (one row per instance and panel) and chain_switches_summary.json.
"""
from __future__ import annotations
import argparse
import bisect
import collections
import json
from pathlib import Path
import statistics
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import locus_common as lc  # noqa: E402
import locus_panels as lp  # noqa: E402

lp.native_library_ready()
import counters as ct  # noqa: E402  (source_panels, source_data)
import gap_chaining as gc  # noqa: E402
import gap_native as gn  # noqa: E402


def trace(res, A, p, labels_of):
    """Pointer kinds and labels along one optimal chain of a packed class result."""
    B, off, back = res.B, res.batch.offsets, res.back
    best = max(B)
    i = next(k for k in range(len(B)) if B[k] == best)
    st = next(s for s in range(off[i], off[i + 1]) if res.C[s] == best)
    kinds, chain_labels, anchors = collections.Counter(), [], set()
    while True:
        a = bisect.bisect_right(off, st) - 1
        h = labels_of(A[a].v)[st - off[a]]
        chain_labels.append(h)
        anchors.add(a)
        nxt = back[st]
        if nxt < 0:
            kinds['restart'] += 1
            break
        j = bisect.bisect_right(off, nxt) - 1
        g = labels_of(A[j].v)[nxt - off[j]]
        kinds['same' if j == a else 'local' if g == h else 'remote'] += 1
        st = nxt
    return best, kinds, chain_labels, len(anchors)


def run(source, name, gamma, out):
    graph, cases, primary, groups = ct.source_data(source, name)
    p = gc.Graph(graph['lengths'], [tuple(e) for e in graph['edges']], graph['paths']).prepare()
    G = gn.GapGraph(p)
    names = [h.split('|')[0] for h in graph['haplotypes']]
    cache = {}

    def labels_of(v):
        if v not in cache:
            cache[v] = list(gc.labels(p.A[v]))
        return cache[v]
    for case in cases:
        if case['id'] not in primary:
            continue
        A = [gc.Anchor(*x[:5], x[5] * lc.SCALE) for x in case['anchors']]
        res = G.solve_packed(G.batch(A, gamma), mode='class')
        best, kinds, chain_labels, n_anchors = trace(res, A, p, labels_of)
        info = groups.get(case['read'], dict(group='present', source='CHM13#0'))
        src = names.index(info['source']) if info['source'] in names else None
        source_best = None
        if src is not None:
            off = res.batch.offsets
            vals = [res.C[off[i] + labels_of(z.v).index(src)] for i, z in enumerate(A) if p.A[z.v] >> src & 1]
            source_best = max(vals) if vals else None
        row = dict(source=source, panel=name, H=p.H, gamma=gamma, id=case['id'], group=info['group'],
                   read_source=info['source'], source_in_panel=src is not None, N=len(A), best=best,
                   chain_anchors=n_anchors, chain_states=len(chain_labels), kinds=dict(kinds),
                   switches=kinds['remote'] + kinds['same'], source_best=source_best,
                   source_attains_optimum=None if source_best is None else source_best == best)
        out.write(json.dumps(row) + '\n')
    out.flush()


def summarize(rows):
    out = {}
    for key in sorted({(r['source'], r['panel'], r['H']) for r in rows}, key=lambda k: (k[0], k[2])):
        for g in ('absent', 'present'):
            rs = [r for r in rows if (r['source'], r['panel'], r['H']) == key and r['group'] == g]
            if not rs:
                continue
            item = dict(instances=len(rs), source_in_panel=sum(r['source_in_panel'] for r in rs),
                        with_switch=sum(r['switches'] > 0 for r in rs),
                        switches_median=statistics.median(r['switches'] for r in rs),
                        switches_max=max(r['switches'] for r in rs),
                        chain_anchors_median=statistics.median(r['chain_anchors'] for r in rs))
            att = [r['source_attains_optimum'] for r in rs if r['source_attains_optimum'] is not None]
            if att:
                item['source_attains_optimum'] = sum(att)
                gaps = [r['best'] - r['source_best'] for r in rs if r['source_best'] is not None and r['source_best'] < r['best']]
                item['source_shortfall_median_when_below'] = statistics.median(gaps) if gaps else 0
            out[f'{key[0]}/{key[1]}/{g}'] = item
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sources', default='kir,chr20,mhc')
    ap.add_argument('--gamma', type=int, default=lc.GAMMA)
    args = ap.parse_args()
    path = HERE / 'chain_switches.jsonl'
    if path.exists():
        raise SystemExit('refusing to overwrite')
    gn.build()
    with path.open('w') as out:
        for source in args.sources.split(','):
            for name in ct.source_panels(source):
                run(source, name, args.gamma, out)
                print(source, name, flush=True)
    rows = [json.loads(s) for s in path.read_text().splitlines()]
    summary = summarize(rows)
    lc.dump(HERE / 'chain_switches_summary.json', dict(script='experiments/second_locus/chain_switches.py',
                                                       gamma=args.gamma, levels=summary))
    for k, v in summary.items():
        print(k, v)


if __name__ == '__main__':
    main()
