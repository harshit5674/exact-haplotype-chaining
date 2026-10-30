"""Deliberately broken variants (mutants) of the specification, to measure test sensitivity.

Each mutant is a MODIFIED COPY of gap_chaining.OriginSolver logic (the package and the working
copy are not changed). A test family is sensitive to a mutant if some instance makes the
harness report a C/B/R mismatch with the brute force or a pointer failure.

Non-equivalent mutants (must be detected):
  ends_first      ENDs before STARTs at equal query coordinate (allows d_j = c_i)
  thr_nonstrict   same-vertex predecessor allowed when y_j = x_i
  remote_gate     remote value applied only to labels with a local predecessor (Minichain v1.3)
  no_same         same-anchor switch omitted
  remote_skip_w   remote query excludes anchors on the last2reach vertex w itself
  switch_strict   same-anchor retention strict (old state loses ties): pointer cycles at gamma=0
  minichain_dist  last2reach distance by the pinned Minichain relaxation (order-dependent)
  remote_any_g    remote switches also from labels g in A_i (not in the contract)
Equivalent mutants (values unchanged; pointers still legal; must NOT be reported):
  floor_ties_local  local candidate wins ties with the floor
  merge_ties_right  merges give ties to the right list
  best_largest      h* is the largest attaining label
"""
from __future__ import annotations
import sys
from bisect import bisect_left, bisect_right
from dataclasses import replace
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'research/current/new'))
import gap_chaining as gc  # noqa: E402
from gap_chaining import Record, Witness, coalesce, labels  # noqa: E402

NON_EQUIVALENT = ('ends_first', 'thr_nonstrict', 'remote_gate', 'no_same', 'remote_skip_w',
                  'switch_strict', 'minichain_dist', 'remote_any_g')
EQUIVALENT = ('floor_ties_local', 'merge_ties_right', 'best_largest')


class MinichainDistPrepared(gc.Prepared):
    """Prepared with the pinned Minichain distance update Distance[v] = len(u) + min(Distance[v],
    Distance[u]) over predecessors in increasing id order (graphUtils.cpp:792 as quoted in
    experiments/gap_solution_review/REVIEW.md)."""

    def __init__(self, graph):
        super().__init__(graph)
        for h in range(self.H):
            for v in range(self.V):
                if self.on[v][h] or self.last[v][h] < 0:
                    continue
                w = self.last[v][h]
                cur = None
                for u in self.pred[v]:
                    if self.last[u][h] != w:
                        continue
                    du = self.distance[u][h] if not self.on[u][h] else 0
                    cur = graph.lengths[u] + (du if cur is None else min(cur, du))
                self.distance[v][h] = cur


class Mutant(gc.OriginSolver):
    def __init__(self, prepared, anchors, gamma, mutant):
        super().__init__(prepared, anchors, gamma)
        self.mutant = mutant
        if mutant == 'thr_nonstrict':
            self.threshold = [bisect_right(self.ends, (z.v, z.x, self.n)) for z in anchors]

    def merge(self, left, right):
        if self.mutant != 'merge_ties_right' or not left or not right:
            return super().merge(left, right)
        # Ties to the right: evaluate per label explicitly (a reference merge, not bit-sliced).
        p = self.p
        out = []
        allm = 0
        for r in left + right:
            allm |= r.mask
        owner = {}
        for side, lst in ((0, left), (1, right)):
            for r in lst:
                for h in labels(r.mask):
                    val = r.coef + p.phi[self.a[r.origin].v][h]
                    if h not in owner or val >= owner[h][0]:
                        owner[h] = (val, side, id(r), r)
        keep = {}
        for h, (_, _, rid, r) in owner.items():
            keep.setdefault(rid, [r, 0])[1] |= 1 << h
        for r in left + right:
            if id(r) in keep:
                out.append(replace(r, mask=keep[id(r)][1]))
        return out

    def remote(self, i):
        if self.mutant not in ('remote_skip_w', 'remote_any_g'):
            return super().remote(i)
        z, p = self.a[i], self.p
        value, witness = None, None
        cand_labels = p.U & ~p.A[z.v] if self.mutant == 'remote_skip_w' else p.U
        for g in labels(cand_labels):
            w = p.last[z.v][g]
            if w < 0:
                continue
            if p.on[z.v][g]:
                w = z.v  # remote_any_g: treat a label through v_i as leaving at v_i itself
            if self.mutant == 'remote_skip_w':
                r = bisect_left(self.ends, (w, -1, -1))
            else:
                r = bisect_right(self.ends, (w, p.graph.lengths[w], self.n)) if w != z.v else self.threshold[i]
            best, bj = None, -1
            while r:
                sel = self.selectors[r]
                if sel is not None:
                    rec, _ = sel.get(g)
                    if rec is not None:
                        v = rec.coef + p.phi[self.a[rec.origin].v][g]
                        if best is None or v > best:
                            best, bj = v, rec.witness.j
                r -= r & -r
            if best is not None:
                dist = p.distance[z.v][g] if not p.on[z.v][g] else 0
                cand = best + z.weight - z.x - z.c + 2 - p.phi[w][g] - dist - self.gamma
                if value is None or cand > value:
                    value, witness = cand, Witness('remote', bj, g, w)
        return value, witness

    def run(self):
        mut = self.mutant
        states = [[] for _ in range(self.n)]
        Bs, Rs = [0] * self.n, [None] * self.n
        phase_start, phase_end = (1, 0) if mut == 'ends_first' else (0, 1)
        events = sorted([(z.c, phase_start, i) for i, z in enumerate(self.a)] +
                        [(z.d, phase_end, i) for i, z in enumerate(self.a)])
        for _, task, i in events:
            z = self.a[i]
            if task == phase_start:
                self.route_cache = {}
                remote, rw = self.remote(i)
                Rs[i] = remote
                floor, fw = z.weight, Witness('restart')
                if remote is not None and remote > floor:
                    floor, fw = remote, rw
                local = self.prefix(self.threshold[i], self.p.A[z.v])
                alpha = z.weight - z.x - z.c + 2
                local = [replace(r, coef=r.coef + alpha) for r in local]
                if mut == 'remote_gate':
                    # Labels without a local predecessor keep the restart value only.
                    pre = self.floor_states(local, i, floor, fw, keep_equal=False)
                    have = 0
                    for r in local:
                        have |= r.mask
                    fixed = []
                    for r in pre:
                        if r.origin == i and r.witness.kind == 'remote':
                            nolocal = r.mask & ~have
                            if nolocal:
                                fixed.append(Record(z.weight, i, nolocal, Witness('restart')))
                            if r.mask & have:
                                fixed.append(replace(r, mask=r.mask & have))
                        else:
                            fixed.append(r)
                    pre = fixed
                else:
                    pre = self.floor_states(local, i, floor, fw, keep_equal=(mut == 'floor_ties_local'))
                best, bh = None, -1
                for r in pre:
                    if r.origin == i:
                        mn, h = 0, ((r.mask & -r.mask).bit_length() - 1 if mut != 'best_largest'
                                    else r.mask.bit_length() - 1)
                    else:
                        mn, h = self.minimum(self.route(r.origin, i), r.mask)
                        if mut == 'best_largest':
                            route = self.route(r.origin, i)
                            hs = [g for g in labels(r.mask)
                                  if sum(((route[k] >> g) & 1) << k for k in range(len(route))) == mn]
                            h = max(hs)
                    val = r.coef - mn
                    if best is None or val > best or (val == best and (h < bh if mut != 'best_largest' else h > bh)):
                        best, bh = val, h
                Bs[i] = best
                if mut == 'no_same':
                    post = pre
                else:
                    post = self.floor_states(pre, i, best - self.gamma, Witness('same', i, bh),
                                             keep_equal=(mut != 'switch_strict'))
                states[i] = coalesce(post)
            else:
                incoming = coalesce([replace(r, coef=r.coef + z.y + z.d) for r in states[i]], publish=i)
                r = self.rank[i]
                while r <= self.n:
                    self.nodes[r] = self.merge(self.nodes[r], incoming)
                    self.selectors[r] = gc.Selector(self.nodes[r])
                    r += r & -r
        return gc.Result(self.a, self.p, states, Bs, Rs, self.stats)


def mutant_outputs(lengths, edges, paths, anchors_t, gamma, mutant):
    graph = gc.Graph(list(lengths), [tuple(e) for e in edges], [list(q) for q in paths])
    p = MinichainDistPrepared(graph) if mutant == 'minichain_dist' else graph.prepare()
    anchors = [gc.Anchor(*a) for a in anchors_t]
    res = Mutant(p, anchors, gamma, mutant).run()
    C, ptr = res.expand()
    return p, anchors, C, res.B, res.R, ptr
