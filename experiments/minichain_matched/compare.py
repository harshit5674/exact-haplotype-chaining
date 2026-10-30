"""Score identity between Minichain's own chaining DP and the project's native engines.

For every (read, component) that Minichain chained in a DP-dump run, rebuild the exact anchor
list Minichain chained (seed export of the same reads, split by component in Minichain's
order, input hash checked), run the native engines with the chosen distance table, and compare
the FNV-1a hashes of the full state table C (every valid (anchor, path) state, in Minichain's
(anchor, path) order), of B (Minichain's C_m), of max(w, R) (Minichain's C_p), the number of
valid states and the best score.

distance='minichain' is the compatibility mode (Algorithm 2 with Minichain's order-dependent
Distance table); distance='true' is Algorithm 2 as transcribed (shortest distances).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from array import array
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mcdp  # noqa: E402
import gap_native as gn  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', required=True, type=Path, help='DP-dump run directory (dump.tsv)')
    ap.add_argument('--seeds', required=True, type=Path, help='seed-export run directory (seeds.tsv, index.tsv)')
    ap.add_argument('--distance', choices=('minichain', 'true'), required=True)
    ap.add_argument('--engines', default='per_haplotype,class')
    ap.add_argument('--out', required=True, type=Path, help='per-pair JSONL (must not exist)')
    ap.add_argument('--limit', type=int, default=0)
    args = ap.parse_args()
    if args.out.exists():
        raise SystemExit(f'{args.out} exists')
    engines = args.engines.split(',')

    index = args.run / 'index.tsv' if (args.run / 'index.tsv').exists() else args.seeds / 'index.tsv'
    comps, header = mcdp.load_index(index)
    seed_comps, _ = mcdp.load_index(args.seeds / 'index.tsv')
    for a, b in zip(comps, seed_comps):   # same graph and component numbering in both runs
        assert (a.glob, a.inn, a.top, a.paths, a.l2r, a.d2b) == (b.glob, b.inn, b.top, b.paths, b.l2r, b.d2b)
    graphs = {c.cid: gn.GapGraph(mcdp.prepare(c, args.distance)) for c in comps}
    seeds, order = mcdp.read_seed_export(args.seeds / 'seeds.tsv')
    dump = mcdp.read_dp_dump(args.run / 'dump.tsv')
    split = mcdp.Splitter(comps)

    totals = dict(pairs=0, states=0, anchors=0, reads=0)
    match = {e: dict(C=0, B=0, Cp=0, best=0, nvalid=0, all=0) for e in engines}
    mismatches = []
    t0 = time.perf_counter()
    with open(args.out, 'w') as fo:
        for n_read, name in enumerate(order):
            if args.limit and n_read >= args.limit:
                break
            parts = split.split(seeds[name]['hits'])
            totals['reads'] += 1
            for c in comps:
                a0 = parts[c.cid]
                key = (name, c.cid)
                if not a0:
                    assert key not in dump
                    continue
                m = dump[key]
                assert m['N'] == len(a0) and mcdp.input_hash(a0) == m['ha'], ('input mismatch', key)
                wt = (a0[0][4] - a0[0][3] + 1) * m['sf']
                assert wt == m['wt'], ('weight', key, wt, m['wt'])
                anchors = mcdp.to_gc_anchors(c, a0, wt)
                row = dict(read=name, cid=c.cid, N=len(a0), recomb=m['recomb'], gate=m['gate'],
                           mc_best=m['best'], mc_nvalid=m['nvalid'])
                for e in engines:
                    res = graphs[c.cid].solve_packed(anchors, m['recomb'], mode=e)
                    hc, hb = mcdp.fnv(res.C), mcdp.fnv(res.B)
                    hp = mcdp.fnv(array('q', [max(wt, r) for r in res.R]))
                    ok = dict(C=hc == m['hc'], B=hb == m['hb'], Cp=hp == m['hp'], best=max(res.C) == m['best'],
                              nvalid=len(res.C) == m['nvalid'])
                    ok['all'] = all(ok.values())
                    for k, v in ok.items():
                        match[e][k] += v
                    row[e] = dict(ok, best=max(res.C))
                    if not ok['all']:
                        mismatches.append(dict(read=name, cid=c.cid, engine=e, **{k: v for k, v in ok.items()}))
                totals['pairs'] += 1
                totals['states'] += m['nvalid']
                totals['anchors'] += len(a0)
                fo.write(json.dumps(row) + '\n')
    summary = dict(run=str(args.run), seeds=str(args.seeds), distance=args.distance, engines=engines,
                   index_header=header, totals=totals, match=match, mismatches=mismatches[:200],
                   mismatch_count=len(mismatches), seconds=time.perf_counter() - t0)
    args.out.with_suffix('.summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(dict(totals=totals, match=match, mismatch_count=len(mismatches)), indent=1))


if __name__ == '__main__':
    main()
