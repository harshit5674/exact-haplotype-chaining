"""Concrete instances for PROOF.md Remark 2.5(c)-(d): the contract is a specific restriction.

The contract (Algorithm 2 as transcribed) computes the optimum over the chains of PROOF.md
Definition 2.2 (verified here by the brute force and by gap_chaining.solve). Two natural
alternative readings of "recombination-aware gap" give strictly higher optima on these tiny
instances, so a paper must state which one it solves:
  (d) shortest-exit remote: a remote switch from (j, g) pays the true shortest graph route
      from the end of j to the start of i, i.e. D(v_j, v_i) - y_j + x_i - 1, instead of
      leaving P_g at last(v_i, g);
  (c) switch-anywhere: a switch from label g to label h between anchors j and i may use any
      route that follows g to some vertex, then any graph route, then h (here: the shortest
      route), also when g passes through v_i.
Writes results/semantics_examples.json.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import brute_chains as bc  # noqa: E402
import harness as H  # noqa: E402


def contract_values(L, E, P, A, gamma):
    facts = bc.Facts(L, E, P)
    C, Cnot, B, R, n = bc.enumerate_chains(facts, A, gamma)
    g = H.GraphCache(L, E, P)
    rs = H.gc.solve(g.p, [H.gc.Anchor(*a) for a in A], gamma)
    C2, _ = rs.expand()
    assert H.flat(C2) == C and rs.B == B and rs.R == R
    return facts, C, B, R


def main():
    out = []
    gamma = 5
    # (d) remote exit fixed at last2reach. P_g = 0 -> 1 -> 2 with len(1) = 100; vertex 3 is on h only;
    # edges 0 -> 3 and 2 -> 3. last(3, g) = 2, so the contract routes 0 -> 1 -> 2 -> 3.
    L = [1, 100, 1, 1]
    E = [(0, 1), (1, 2), (0, 3), (2, 3)]
    P = [[0, 1, 2], [3]]
    A = [(0, 1, 1, 1, 1, 200), (3, 1, 1, 2, 2, 200)]
    facts, C, B, R = contract_values(L, E, P, A, gamma)
    contract_gap = (facts.phi[0][2] + facts.D[2][3] + 1) - (facts.phi[0][0] + 1) + (2 - 1) - 2
    shortest_gap = (facts.D[0][3] - 1 + 1 - 1) + (2 - 1 - 1)
    out.append(dict(example='(d) remote exit fixed at last(v_i, g)', lengths=L, edges=E, paths=P, anchors=A,
                    gamma=gamma, last_3_g=facts.last[3][0], contract_C_of_anchor1=C[(1, 1)],
                    contract_remote_gap=contract_gap, shortest_route_gap=shortest_gap,
                    shortest_exit_value=200 + 200 - shortest_gap - gamma,
                    note='contract value = 200 + 200 - contract_gap - gamma; a shortest-route reading gives more'))
    # (c) no remote switch when g passes through v_i. j on vertex 0 (g only); i on vertex 3 (g and h).
    # g = 0 -> 1 -> 3 with len(1) = 100; h = 2 -> 3; edge 0 -> 2 on no haplotype.
    L = [1, 100, 1, 1]
    E = [(0, 1), (1, 3), (0, 2), (2, 3)]
    P = [[0, 1, 3], [2, 3]]
    A = [(0, 1, 1, 1, 1, 200), (3, 1, 1, 2, 2, 200)]
    facts, C, B, R = contract_values(L, E, P, A, gamma)
    local_gap = (facts.phi[0][3] + 1) - (facts.phi[0][0] + 1) + (2 - 1) - 2
    switch_route_gap = (facts.D[0][3] - 1 + 1 - 1) + (2 - 1 - 1)
    out.append(dict(example='(c) no switch between anchors when g passes v_i', lengths=L, edges=E, paths=P,
                    anchors=A, gamma=gamma, contract_C=dict((str(k), v) for k, v in C.items()), contract_R=R,
                    contract_best_via_g=200 + 200 - local_gap, switch_anywhere_value=200 + 200 - switch_route_gap - gamma,
                    note='contract: local along g (gap through the long vertex 1), optionally then switch to h; '
                         'a switch-anywhere reading may leave g at 0 and join h at 2'))
    (HERE / 'results').mkdir(exist_ok=True)
    (HERE / 'results/semantics_examples.json').write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
