"""CGJ's prose Problem 5 versus their Algorithm 2 (the project's contract), on a two-anchor instance.

CGJ (Genome Research 2024, Europe PMC PMC11529843) state Problem 5 as maximizing
    sum weight(M[s_j.i]) - sum_{j>=2} gap(M[s_{j-1}.i], M[s_j.i]) - gamma * r(S),
"assuming the same definition of function gap(M[i], M[j]) as used by Chandra and Jain (2023);
see Problem 2a". Problem 2a's gap is a function of the two anchors only: the minimum, over cover
paths i through the predecessor's vertex, of the route along P_i to last2reach(target, i) followed
by a shortest path. Read literally, a chain may therefore stay on haplotype h1 while being charged
another haplotype's shorter route, with no recombination penalty. Algorithm 2 (and Minichain, and
docs/MATHEMATICAL_MODEL.md's contract) instead charges a step along the haplotype the chain
follows, and a route through another haplotype only as a remote switch costing gamma.

Instance: h1 runs v0 -> X (length 100) -> v4; h2 runs v0 -> Y -> Z -> v5, and an edge Z -> v4 lies
on no haplotype. Anchor j sits on v0 (labels {h1, h2}), anchor i on v4 (label {h1}); weights 200,
gamma 5.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'research/current/new'))
import gap_chaining as gc  # noqa: E402

LENGTHS = [1, 100, 1, 1, 1, 1]
EDGES = [(0, 1), (1, 4), (0, 2), (2, 3), (3, 5), (3, 4)]
PATHS = [[0, 1, 4], [0, 2, 3, 5]]
ANCHORS = [gc.Anchor(0, 1, 1, 1, 1, 200), gc.Anchor(4, 1, 1, 2, 2, 200)]
GAMMA = 5


def main():
    p = gc.Graph(LENGTHS, EDGES, PATHS).prepare()
    res = gc.solve(p, ANCHORS, GAMMA)
    oracle = gc.direct_oracle(p, ANCHORS, GAMMA)
    j, i = ANCHORS
    cands = []
    for g in range(len(PATHS)):
        if not p.on[j.v][g]:
            continue
        w, D = p.last[i.v][g], p.distance[i.v][g]
        gap_graph = p.phi[w][g] - p.phi[j.v][g] + D + i.x - j.y - 1
        gap_query = i.c - j.d - 1
        cands.append(dict(g=g, last2reach=w, D=D, gap=gap_graph + gap_query))
    prose = j.weight + i.weight - min(c['gap'] for c in cands)   # chain (j,h1) -> (i,h1), r(S) = 0
    out = dict(contract_B=res.B, contract_R=res.R, contract_records=[[(r.coef, r.origin, r.mask, r.witness.kind)
                                                                     for r in row] for row in res.records],
               direct_oracle=repr(oracle)[:400], problem2a_candidates=cands,
               prose_problem5_value_at_i_h1=prose,
               note='Algorithm 2 gives 393 (remote switch from h2, paying gamma); the literal prose reading gives '
                    '398 by staying on h1 while charging h2\'s route without a recombination penalty.')
    assert res.B[1] == 393 and prose == 398
    print(json.dumps(out, indent=1, default=str))
    (HERE / 'results' / 'problem5_prose_example.json').write_text(json.dumps(out, indent=1, default=str) + '\n')


if __name__ == '__main__':
    main()
