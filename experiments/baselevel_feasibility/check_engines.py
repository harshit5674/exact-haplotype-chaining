"""Exactness checks and work counters of the native engines on a base-level window panel (2026-10-20, spike).

For every instance of windows/<name>/panel (real reads, pinned Minichain seeds, weights x200, gamma 10000):
  - the per-haplotype, class (expanded), guarded and compact engines must give identical C, B, R
    (compact after expansion), and every pointer of every engine is checked by gap_chaining.verify_pointers
    on a bounded subsample;
  - bounded direct-recurrence checks (gap_chaining.direct_oracle) on evenly spaced anchor subsets of
    instances spread over the anchor-count range, at gamma 0, 1, 10000 and 20000, as the pilot did;
  - the class engine's counters (records, pair decisions, bit-sliced decisions) and guard fallbacks.
The full independent graph-major validator stores all-pairs distances (V^2 words), too large for this
window on this host; that limit is reported, not hidden. Nothing here is timed.

  python3.13 -B experiments/baselevel_feasibility/check_engines.py --window windows/w31250k_250kb
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import statistics
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'research/current/new'))
import gap_chaining as gc  # noqa: E402
import gap_native as gn  # noqa: E402

GAMMA, SCALE = 10000, 200
MODES = ('per_haplotype', 'class', 'auto', 'compact')


def solve(G, batch, mode):
    return G.solve_compact(batch, warm=False) if mode == 'compact' else G.solve_packed(batch, mode=mode, warm=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--window', type=Path, required=True)
    ap.add_argument('--direct', type=int, default=8, help='instances for bounded direct-recurrence checks')
    args = ap.parse_args()
    win = args.window if args.window.is_absolute() else HERE / args.window
    panel = win / 'panel'
    graph = json.loads((panel / 'graph.json').read_text())
    cases = [json.loads(line) for line in (panel / 'instances.jsonl').open()]
    t0 = time.perf_counter()
    p = gc.Graph(graph['lengths'], [tuple(e) for e in graph['edges']], graph['paths']).prepare()
    G = gn.GapGraph(p)
    prep_seconds = time.perf_counter() - t0
    rows = []
    for c in cases:
        A = [gc.Anchor(*x[:5], x[5] * SCALE) for x in c['anchors']]
        batch = G.batch(A, GAMMA)
        out = {m: solve(G, batch, m) for m in MODES}
        ref = out['per_haplotype'].arrays_digest()
        assert out['class'].arrays_digest() == ref and out['auto'].arrays_digest() == ref, c['id']
        assert out['compact'].expand().arrays_digest() == ref, c['id']
        cs, au = out['class'].stats, out['auto'].stats
        assert au['auto_fallback_available'] and au['auto_work'] <= au['auto_budget']
        rows.append(dict(id=c['id'], read=c['read'], N=batch.n, K=batch.K,
                         records=out['compact'].stats['compact_records'], fallback=bool(au['fallback']),
                         pairs=cs['contested_pairs'] if 'contested_pairs' in cs else None,
                         stats={k: cs[k] for k in cs if isinstance(cs[k], (int, float))}))
    # Bounded direct-recurrence and pointer checks, spread over the anchor-count range.
    by_size = sorted((c for c in cases if c['anchors']), key=lambda c: len(c['anchors']))
    picks = [by_size[round(i * (len(by_size) - 1) / max(1, args.direct - 1))] for i in range(args.direct)]
    direct = []
    for c in picks:
        full = [gc.Anchor(*x[:5], x[5] * SCALE) for x in c['anchors']]
        n = min(32, len(full))
        idx = sorted({round(i * (len(full) - 1) / (n - 1)) for i in range(n)}) if n > 1 else [0]
        A = [full[i] for i in idx]
        for gamma in (0, 1, 10000, 20000):
            expected = gc.direct_oracle(p, A, gamma)
            for mode in ('per_haplotype', 'class', 'auto'):
                z = G.solve_packed(A, gamma, mode=mode)
                C, ptr, B, R = z.as_dicts(A)
                assert (C, B, R) == expected[:3], (c['id'], gamma, mode)
                gc.verify_pointers(p, A, gamma, C, B, ptr)
            direct.append(dict(id=c['id'], anchors=len(A), gamma=gamma, states=sum(len(r) for r in expected[0])))
    reads = defaultdict(list)
    for r in rows:
        reads[r['read']].append(r)
    def tot(v, key):
        return sum(x['stats'][key] for x in v)
    per_read = [dict(read=k, N=sum(x['N'] for x in v), K=sum(x['K'] for x in v), records=sum(x['records'] for x in v),
                     fallback=any(x['fallback'] for x in v), pairs=tot(v, 'pairs'), sliced=tot(v, 'pairs_sliced'),
                     same=tot(v, 'pairs_same_vertex'), loose=tot(v, 'pairs_loose'), rng=tot(v, 'pairs_range'))
                for k, v in reads.items()]
    shares = {k: statistics.median(r[k] / r['pairs'] for r in per_read if r['pairs'])
              for k in ('same', 'loose', 'rng', 'sliced')}
    summary = dict(panel=str(panel.relative_to(ROOT)), H=p.H, V=p.V, W=G.W, instances=len(rows), reads=len(per_read),
                   anchors=sum(r['N'] for r in rows), states=sum(r['K'] for r in rows), prep_seconds=prep_seconds,
                   states_per_anchor_median=statistics.median(r['K'] / r['N'] for r in per_read),
                   records_per_anchor_median=statistics.median(r['records'] / r['N'] for r in per_read),
                   sliced_per_anchor_median=statistics.median(r['sliced'] / r['N'] for r in per_read),
                   decision_shares_median=shares,
                   fallback_reads=sum(r['fallback'] for r in per_read),
                   fallback_instances=sum(r['fallback'] for r in rows),
                   direct_checks=len(direct), direct_states=sum(d['states'] for d in direct),
                   stat_keys=sorted(rows[0]['stats']))
    (win / 'engine_check.json').write_text(json.dumps(dict(summary=summary, rows=rows, direct=direct), indent=1) + '\n')
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
