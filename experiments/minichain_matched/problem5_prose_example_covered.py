"""CGJ's prose Problem 5 versus their Algorithm 2, on an instance whose graph is fully covered by haplotypes.

Added 2026-10-20. The instance of problem5_prose_example.py (kept
unchanged with its saved result) has an edge Z -> v4 on no haplotype, and CGJ's graph model requires every
vertex and edge to lie on a haplotype. Here a third haplotype h3 = (Z, v4) covers that edge without passing
the first anchor's vertex v0, so every vertex and edge lies on a haplotype and the distinction survives in
CGJ's own domain.

Instance: h1 runs v0 -> X (length 100) -> v4; h2 runs v0 -> Y -> Z -> v5; h3 runs Z -> v4. Anchor j sits on
v0 (labels {h1, h2}), anchor i on v4 (labels {h1, h3}); weights 200, gamma 5. Algorithm 2 (the project's
contract) gives 393, a remote switch from h2 paying gamma; the literal prose reading gives 398 by staying on
h1 while being charged h2's shorter route without a recombination penalty. The prose value is computed as in
problem5_prose_example.py: Problem 2a's gap is the minimum, over haplotypes through the predecessor's
vertex, of the route along that haplotype to last2reach(target) followed by a shortest path.

    python3.13 -B experiments/minichain_matched/problem5_prose_example_covered.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'research/current/new'))
import gap_chaining as gc  # noqa: E402

LENGTHS = [1, 100, 1, 1, 1, 1]                      # v0, X, Y, Z, v4, v5
EDGES = [(0, 1), (1, 4), (0, 2), (2, 3), (3, 5), (3, 4)]
PATHS = [[0, 1, 4], [0, 2, 3, 5], [3, 4]]           # h1, h2, h3
ANCHORS = [gc.Anchor(0, 1, 1, 1, 1, 200), gc.Anchor(4, 1, 1, 2, 2, 200)]
GAMMA = 5
OUT = HERE / 'results' / 'problem5_prose_example_covered.json'


def main():
    covered_edges = {(p[k], p[k + 1]) for p in PATHS for k in range(len(p) - 1)}
    covered_vertices = {v for p in PATHS for v in p}
    uncovered_edges = sorted(set(EDGES) - covered_edges)
    uncovered_vertices = sorted(set(range(len(LENGTHS))) - covered_vertices)
    assert not uncovered_edges and not uncovered_vertices

    p = gc.Graph(LENGTHS, EDGES, PATHS).prepare()
    res = gc.solve(p, ANCHORS, GAMMA)
    C, B, R, back = gc.direct_oracle(p, ANCHORS, GAMMA)
    assert list(res.B) == list(B) and list(res.R) == list(R)
    j, i = ANCHORS
    assert not p.on[j.v][2], 'h3 must not pass the first anchor'
    cands = []
    for g in range(len(PATHS)):
        if not p.on[j.v][g]:
            continue
        w, D = p.last[i.v][g], p.distance[i.v][g]
        gap_graph = p.phi[w][g] - p.phi[j.v][g] + D + i.x - j.y - 1
        gap_query = i.c - j.d - 1
        cands.append(dict(g=g, last2reach=w, D=D, gap=gap_graph + gap_query))
    prose = j.weight + i.weight - min(c['gap'] for c in cands)   # chain (j,h1) -> (i,h1), r(S) = 0
    out = dict(instance=dict(lengths=LENGTHS, edges=EDGES, paths=PATHS, gamma=GAMMA,
                             anchors=[[a.v, a.x, a.y, a.c, a.d, a.weight] for a in ANCHORS]),
               uncovered_edges=uncovered_edges, uncovered_vertices=uncovered_vertices,
               contract_B=res.B, contract_R=res.R,
               contract_records=[[(r.coef, r.origin, r.mask, r.witness.kind) for r in row] for row in res.records],
               direct_oracle_C=[{str(h): v for h, v in row.items()} for row in C],
               direct_oracle_B=list(B), direct_oracle_R=list(R),
               problem2a_candidates=cands, prose_problem5_value_at_i_h1=prose,
               note='Every vertex and edge lies on a haplotype. Algorithm 2 gives 393 at both labels of the second '
                    'anchor (remote switch from h2, paying gamma); the literal prose reading gives 398.')
    assert res.B[1] == 393 and prose == 398
    assert all(v == 393 for v in C[1].values())
    print(json.dumps(out, indent=1, default=str))
    OUT.write_text(json.dumps(out, indent=1, default=str) + '\n')


if __name__ == '__main__':
    main()
