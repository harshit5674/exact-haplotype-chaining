"""Certified scalar fast paths for the gap-sensitive origin-class solver (Stage A).

FastPathSolver subclasses gap_chaining.OriginSolver, the byte-identical working
copy of the external specification, which is not modified. It takes every decision
of the specification - merge winners, floor and switch thresholds, masked minima,
remote point queries - and returns the same records. It decides with scalar tests
whenever they are certified, and uses bit-sliced route planes only when they are not:

  (i)   two records whose origins lie on the same VERTEX: the coordinate cancels;
  (ii)  loose per-vertex bounds: phi_h(v) lies in [lo(v), hi(v)] over h in A_v;
  (iii) the exact range [min, max] of the route length l_h(u, w) = phi_h(w) - phi_h(u)
        over the haplotypes through both u and w (u <= w), computed once per vertex
        pair and cached; a range with min == max decides a masked minimum too.

Invariant (R): every record's labels pass through its origin vertex. Floor and switch
records created at anchor i have origin i and labels in A_i; a local candidate
restricts a published record's mask to A_i; publication and trimming only shrink or
keep masks. Hence for any two records the overlap of their masks lies in
A_u & A_w for their origin vertices u, w, and a floor, switch or minimum at anchor i
acts on a mask inside A_{v_o} & A_{v_i}. The range of (iii) bounds every compared
lane, so its decisions are exact. Topological numbering puts v_o <= v_i on every
path through both, so route lengths are nonnegative.

Undecided merge pairs are resolved by one bit-sliced threshold on the route planes of
their two origin vertices, restricted to the pair's overlap. Every label of an overlap
lies in exactly one left and one right record, so this equals the specification's
aggregate comparison. Ties go left, exactly as in the specification.

COST MODEL. Alongside the computation the solver counts the work of the planned native
engine (research/current/new/gap_native.cpp), in these units:
  mask    one W-word operation on a label mask (AND/OR/ANDNOT with a fused
          nonzero test, or a record copy), W = ceil(H/64);
  col     one bit-sliced column step on W words (route-plane subtraction,
          threshold comparison, masked-minimum step);
  scalar  one word-sized step: Fenwick node visit, pair decision, comparison,
          cache lookup, remote record scan, sort comparison;
  labels  labels scanned to build an exact pair range (one scalar each).
predicted class word operations = W*(mask + col) + scalar + labels, plus 2K if
states are expanded. For a per-haplotype Fenwick solver of the same recurrence it
counts Fenwick steps: |A_i| x (prefix nodes at thr(i) + update nodes at er(i)) per
anchor, plus prefix nodes at the last END on last(v_i, g) for each missing label g.
Python big-integer timings do not predict native cost; only these counters are used.
"""
from __future__ import annotations
from bisect import bisect_right
from collections import Counter
from dataclasses import replace
from typing import Optional

from gap_chaining import (Anchor, Graph, Prepared, Record, Result, Witness,
                          OriginSolver, coalesce, labels, solve)

__all__ = ['VertexBounds', 'PairCache', 'RemoteGroups', 'FastPathSolver', 'solve_fast',
           'per_haplotype_steps', 'predicted_word_ops']


class VertexBounds:
    """Reusable per-graph offset bounds lo(v) <= phi_h(v) <= hi(v) over h in A_v; O(HV)."""
    def __init__(self, p: Prepared):
        self.p = p
        self.lo, self.hi = [0] * p.V, [0] * p.V
        for v in range(p.V):
            offs = [p.phi[v][h] for h in labels(p.A[v])]
            if offs:
                self.lo[v], self.hi[v] = min(offs), max(offs)


class PairCache:
    """Exact route ranges and route bitplanes per vertex pair u < w, built lazily.

    range(u, w) is (min, max) of phi_h(w) - phi_h(u) over h in A_u & A_w.
    planes(u, w) are b bitplanes of that difference (valid on A_u & A_w only).
    A fresh cache per read is the conservative (cold) accounting.
    """
    def __init__(self, p: Prepared):
        self.p = p
        self.ranges: dict = {}
        self.planes_cache: dict = {}

    def range(self, u: int, w: int, cost: Counter) -> tuple[int, int]:
        assert u < w
        r = self.ranges.get((u, w))
        cost['scalar'] += 1
        if r is None:
            common = self.p.A[u] & self.p.A[w]
            cost['mask'] += 1
            vals = [self.p.phi[w][h] - self.p.phi[u][h] for h in labels(common)]
            cost['labels'] += len(vals)
            cost['pair_ranges'] += 1
            assert vals, 'range requested for vertices without a common haplotype'
            r = (min(vals), max(vals))
            assert 0 <= r[0] <= r[1] <= self.p.S
            self.ranges[(u, w)] = r
        return r

    def planes(self, u: int, w: int, cost: Counter) -> list[int]:
        assert u < w
        ans = self.planes_cache.get((u, w))
        cost['scalar'] += 1
        if ans is None:
            p, U = self.p, self.p.U
            left, right = p.phibits[w], p.phibits[u]
            borrow, ans = 0, []
            for x, y in zip(left, right):
                ans.append((x ^ y ^ borrow) & U)
                borrow = ((~x & (y | borrow)) | (y & borrow)) & U
            cost['col'] += p.b
            cost['route_planes'] += 1
            self.planes_cache[(u, w)] = ans
        return ans


class RemoteGroups:
    """Missing labels of every vertex grouped by their last2reach vertex; O(HV), reusable.

    For a label g not on v with w = last(v, g) defined, the remote transition of an
    anchor on v subtracts phi_g(w) + D(w, v), where D(w, v) = distance[v][g] is the
    shortest character distance from w to v and so depends on w alone (checked
    here). groups[v] lists (w, mask of such g, D(w, v)) by w.
    """
    def __init__(self, p: Prepared):
        self.p = p
        self.groups = []
        for v in range(p.V):
            by = {}
            for g in labels(p.U & ~p.A[v]):
                w = p.last[v][g]
                if w < 0:
                    continue
                d = p.distance[v][g]
                if w in by:
                    assert by[w][1] == d, 'last2reach distance must depend on w only'
                    by[w][0] |= 1 << g
                else:
                    by[w] = [1 << g, d]
            self.groups.append(sorted((w, m, d) for w, (m, d) in by.items()))


def per_haplotype_steps(n: int, r: int, up: bool) -> int:
    """Fenwick nodes visited by a prefix query at r (up=False) or an update at r."""
    steps = 0
    if up:
        while r <= n:
            steps += 1
            r += r & -r
    else:
        while r:
            steps += 1
            r -= r & -r
    return steps


def _sort_cost(m: int) -> int:
    return m * max(1, (m - 1).bit_length()) if m > 1 else 0


class FastPathSolver(OriginSolver):
    """Origin-class solver with certified scalar decisions; identical records to the spec."""

    def __init__(self, prepared: Prepared, anchors: list[Anchor], gamma: int,
                 bounds: Optional[VertexBounds] = None, cache: Optional[PairCache] = None,
                 tests: str = 'iii', remote: str = 'labels', cut: str = 'none',
                 groups: Optional[RemoteGroups] = None):
        """tests: 'none' | 'i' | 'ii' | 'iii' (cumulative certified tests).
        remote: 'labels' scans node lists once per missing label (the planned native
        point query); 'grouped' answers each last2reach group at once.
        cut: 'none' | 'zero' | 'range' drops prefix records that cannot beat the floor.
        """
        super().__init__(prepared, anchors, gamma)
        if tests not in ('none', 'i', 'ii', 'iii'):
            raise ValueError('tests must be none, i, ii or iii')
        if remote not in ('labels', 'grouped') or cut not in ('none', 'zero', 'range'):
            raise ValueError('unknown remote or cut option')
        if (remote == 'grouped' or cut == 'range') and tests != 'iii':
            raise ValueError('grouped remote and range cut use the exact pair ranges')
        self.tests, self.remote_mode, self.cut = tests, remote, cut
        self.groups = groups if groups is not None else (RemoteGroups(prepared) if remote == 'grouped' else None)
        assert self.groups is None or self.groups.p is prepared
        self._cut = None
        self._tau = {}
        self.bounds = bounds if bounds is not None else VertexBounds(prepared)
        assert self.bounds.p is prepared
        self.cache = cache if cache is not None else PairCache(prepared)
        assert self.cache.p is prepared
        self.cost = Counter()
        self.merge_pairs = set()   # distinct (u, w) origin-vertex pairs met in contested merges
        self.target_pairs = set()  # distinct (v_o, v_i) pairs met in floors, minima, switches
        self.vertex = [z.v for z in anchors]

    # ----- certified decisions -------------------------------------------------
    def _route_range(self, u: int, w: int) -> tuple[int, int]:
        """Range of phi_h(w) - phi_h(u) over A_u & A_w, for any vertex order."""
        if u == w:
            return 0, 0
        if u < w:
            return self.cache.range(u, w, self.cost)
        lo, hi = self.cache.range(w, u, self.cost)
        return -hi, -lo

    def _decide(self, ua: int, ub: int, t: int):
        """Right wins at h iff delta(h) = phi_h(ub) - phi_h(ua) > t.

        Returns (True | False | None, level): True = right wins the whole overlap,
        False = left wins it, None = undecided; level 1 same vertex, 2 loose bounds,
        3 exact pair range, 4 undecided.
        """
        c = self.cost
        c['scalar'] += 1
        if ua == ub and self.tests != 'none':
            return (0 > t), 1
        if self.tests in ('ii', 'iii'):
            c['scalar'] += 1
            b = self.bounds
            lo, hi = b.lo[ub] - b.hi[ua], b.hi[ub] - b.lo[ua]
            if lo > t:
                return True, 2
            if hi <= t:
                return False, 2
        if self.tests == 'iii':
            lo, hi = self._route_range(ua, ub)
            if lo > t:
                return True, 3
            if hi <= t:
                return False, 3
        return None, 4

    def _target_range(self, vo: int, vi: int) -> tuple[int, int]:
        """Certified bounds on l_h(vo, vi) >= 0 for the enabled tests."""
        if vo == vi and self.tests != 'none':
            return 0, 0
        if self.tests == 'iii':
            self.target_pairs.add((vo, vi))
            return self._route_range(vo, vi)
        if self.tests == 'ii':
            self.cost['scalar'] += 1
            b = self.bounds
            return max(0, b.lo[vi] - b.hi[vo]), b.hi[vi] - b.lo[vo]
        return 0, self.p.S + 1

    def _threshold(self, planes: list[int], theta: int, mask: int) -> int:
        """{h in mask : value(h) < theta} by a bit-sliced comparison (b columns)."""
        self.cost['sliced_thresholds'] += 1
        self.stats['threshold_calls'] += 1
        if theta <= 0:
            return 0
        if theta > self.p.S:
            return mask
        self.cost['col'] += self.b
        rhs = [self.U if (theta >> k) & 1 else 0 for k in range(self.b)]
        return self.gt(rhs, planes, mask)

    def _sliced_right_wins(self, ua: int, ub: int, t: int, inter: int) -> int:
        self.cost['sliced_pairs'] += 1
        if ua == ub:  # only reached with tests='none'
            planes = [0] * self.b
            self.cost['col'] += self.b
            return inter & ~self._threshold(planes, t + 1, inter)
        if ua < ub:   # delta = l(ua, ub) >= 0; right wins iff l > t, i.e. not l < t+1
            lt = self._threshold(self.cache.planes(ua, ub, self.cost), t + 1, inter)
            self.cost['mask'] += 1
            return inter & ~lt
        # delta = -l(ub, ua); right wins iff l(ub, ua) < -t
        return self._threshold(self.cache.planes(ub, ua, self.cost), -t, inter)

    # ----- merges ------------------------------------------------------------------
    def merge(self, left: list[Record], right: list[Record]) -> list[Record]:
        """Pointwise maximum; ties go to the left list; never enumerate labels."""
        c, st = self.cost, self.stats
        m = len(left) + len(right)
        st['merges'] += 1
        st['lambda_records'] += m
        c['merges'] += 1
        if not left or not right:
            c['merges_trivial'] += 1
            ans = list(left or right)
            c['mask'] += len(ans)  # copy into the destination list
            st['lambda_records'] += len(ans)
            return ans
        lm = rm = 0
        for r in left:
            lm |= r.mask
        for r in right:
            rm |= r.mask
        # Native lists keep their union masks; one op for the overlap.
        overlap = lm & rm
        c['mask'] += 1
        if not overlap:
            c['merges_disjoint'] += 1
            ans = left + right
            c['mask'] += len(right) + 1  # append right records, update the union
            st['lambda_records'] += len(ans)
            return ans
        origins = {r.origin for r in left} | {r.origin for r in right}
        spec_kind = 'scalar' if len(origins) == 1 else 'coded'
        c['merges_contested'] += 1
        c['merges_spec_' + spec_kind] += 1
        win = rm & ~lm
        c['mask'] += 1
        level = 0
        for a in left:
            ia = a.mask & overlap
            c['mask'] += 1
            if not ia:
                continue
            ua = self.vertex[a.origin]
            for b in right:
                inter = ia & b.mask
                c['mask'] += 1
                if not inter:
                    continue
                ub = self.vertex[b.origin]
                c['pairs'] += 1
                verdict, how = self._decide(ua, ub, a.coef - b.coef)
                level = max(level, how)
                c['pairs_level%d' % how] += 1
                if how >= 2:
                    self.merge_pairs.add((min(ua, ub), max(ua, ub)))
                if verdict is True:
                    win |= inter
                    c['mask'] += 1
                elif verdict is None:
                    before = c['col']
                    win |= self._sliced_right_wins(ua, ub, a.coef - b.coef, inter)
                    c['col_merge'] += c['col'] - before
                    c['mask'] += 1
        c['merges_level%d' % level] += 1
        if spec_kind == 'coded':
            c['coded_level%d' % level] += 1
        ans = []
        for r in left:
            if r.mask & overlap:
                c['mask'] += 1
                if (keep := r.mask & ~win):
                    ans.append(replace(r, mask=keep))
            else:
                ans.append(r)
        for r in right:
            c['mask'] += 1
            if (keep := r.mask & win):
                ans.append(replace(r, mask=keep))
        c['mask'] += 1  # union of the result
        st['lambda_records'] += len(ans)
        return ans

    def _floor_cut(self, z: Record, filtered: bool) -> bool:
        """True if record z cannot beat the floor on any label of the current START.

        Its local value at h is K + alpha_i - l_h(v_o, v_i) <= K + alpha_i - lo, and a
        local candidate wins only if strictly above the floor f. Dropping z leaves every
        floor_states output unchanged: wherever z would win the merge, the winner is
        at most f and is replaced by the floor record either way.
        """
        i, bound = self._cut  # bound = f - alpha_i
        self.cost['scalar'] += 1
        if not filtered:  # before the A_i filter: the free bound l >= 0
            lo = 0
        elif self.cut == 'range':  # after it: z shares a label with v_i, so the range exists
            vo, vi = self.vertex[z.origin], self.vertex[i]
            lo = 0 if vo == vi else self.cache.range(vo, vi, self.cost)[0]
        else:
            return False
        if z.coef - lo <= bound:
            self.cost['cut_records'] += 1
            return True
        return False

    def prefix(self, r: int, mask: int) -> list[Record]:
        ans = []
        while r:
            self.cost['scalar'] += 1  # node visit
            node = self.nodes[r]
            self.stats['lambda_records'] += len(node)
            filtered = []
            for z in node:
                self.cost['prefix_records'] += 1
                if self._cut is not None and self._floor_cut(z, False):
                    continue
                self.cost['mask'] += 1  # restrict to A_i (writes the filtered record)
                if (m := z.mask & mask):
                    if self.cut == 'range' and self._cut is not None and self._floor_cut(z, True):
                        continue
                    filtered.append(replace(z, mask=m))
            self.cost['mask'] += len(filtered)  # union of the filtered list
            ans = self.merge(ans, filtered) if filtered else ans
            r -= r & -r
        return ans

    # ----- remote point queries by scanning short node lists ---------------------
    def remote(self, i: int):
        z, p, c = self.a[i], self.p, self.cost
        value, witness = None, None
        missing = self.U & ~p.A[z.v]
        self.stats['missing_incidences'] += missing.bit_count()
        c['mask'] += 1
        for g in labels(missing):
            w = p.last[z.v][g]
            c['scalar'] += 1
            if w < 0:
                continue
            c['remote_queries'] += 1
            r = bisect_right(self.ends, (w, p.graph.lengths[w], self.n))
            c['ph_remote'] += per_haplotype_steps(self.n, r, False)
            best, bj = None, -1
            while r:
                c['scalar'] += 1
                c['remote_nodes'] += 1
                for rec in self.nodes[r]:
                    c['scalar'] += 1  # one mask word tested
                    c['remote_scans'] += 1
                    if rec.mask >> g & 1:
                        v = rec.coef + p.phi[self.vertex[rec.origin]][g]
                        c['scalar'] += 2
                        if best is None or v > best:
                            best, bj = v, rec.witness.j
                        break
                r -= r & -r
            if best is not None:
                cand = (best + z.weight - z.x - z.c + 2
                        - p.phi[w][g] - p.distance[z.v][g] - self.gamma)
                c['scalar'] += 3
                if value is None or cand > value:
                    value, witness = cand, Witness('remote', bj, g, w)
        return value, witness

    def tau(self, w: int) -> int:
        """END rank of the last anchor end on or before vertex w (cached per read)."""
        r = self._tau.get(w)
        if r is None:
            r = bisect_right(self.ends, (w, self.p.graph.lengths[w], self.n))
            self.cost['scalar'] += max(1, self.n.bit_length())
            self._tau[w] = r
        return r

    def remote_grouped(self, i: int):
        """R_i and its witness from one pass per last2reach group.

        For a group (w, G, D): every g in G has candidate
        max_records [K - l_g(v_o, w)] + alpha_i - D - gamma over records holding g in the
        prefix nodes of tau(w). So the group's maximum is the maximum over (node, record)
        of K - min_{g in mask & G} l_g(v_o, w), decided by the cached pair range when it
        is single-valued and otherwise by a masked minimum. Ties replicate the
        specification: the smallest label, then the earliest node in walk order.
        """
        z, p, c = self.a[i], self.p, self.cost
        self.stats['missing_incidences'] += (self.U & ~p.A[z.v]).bit_count()
        alpha = z.weight - z.x - z.c + 2
        best, witness = None, None  # best = (candidate, -label, -node order)
        for w, G, D in self.groups.groups[z.v]:
            c['scalar'] += 2
            c['remote_groups'] += 1
            r = self.tau(w)
            labels_in_group = G.bit_count()
            c['remote_queries'] += labels_in_group
            c['ph_remote'] += labels_in_group * per_haplotype_steps(self.n, r, False)
            const = alpha - D - self.gamma
            cands, order = [], 0
            while r:
                c['scalar'] += 1
                c['remote_nodes'] += 1
                node = self.nodes[r]
                union = 0
                for rec in node:
                    union |= rec.mask
                c['mask'] += 1  # maintained node union against G
                if union & G:
                    for rec in node:
                        inter = rec.mask & G
                        c['mask'] += 1
                        c['remote_scans'] += 1
                        if inter:
                            vo = self.vertex[rec.origin]
                            lo, hi = (0, 0) if vo == w else self.cache.range(vo, w, c)
                            cands.append((rec.coef - lo, order, rec, inter, lo, hi, vo))
                            c['scalar'] += 1
                order += 1
                r -= r & -r
            cands.sort(key=lambda t: -t[0])
            c['scalar'] += _sort_cost(len(cands))
            for ub, order, rec, inter, lo, hi, vo in cands:
                low = (inter & -inter).bit_length() - 1
                c['mask'] += 1
                c['scalar'] += 1
                if best is not None and (ub + const, -low, -order) < best:
                    c['remote_pruned'] += 1
                    continue
                if lo == hi:
                    val, g = ub, low
                    c['remote_single'] += 1
                else:
                    before = c['col']
                    mn, g = self._masked_minimum(self.cache.planes(vo, w, c), inter)
                    c['col_remote'] += c['col'] - before
                    val = rec.coef - mn
                    c['remote_sliced'] += 1
                key = (val + const, -g, -order)
                if best is None or key > best:
                    best, witness = key, Witness('remote', rec.witness.j, g, w)
        return (None, None) if best is None else (best[0], witness)

    # ----- floor, minimum and switch -----------------------------------------------
    def floor_states(self, records: list[Record], i: int, floor: int,
                     witness: Witness, keep_equal: bool) -> list[Record]:
        c = self.cost
        stage = 'switch' if keep_equal else 'floor'
        mask = self.p.A[self.a[i].v]
        vi = self.vertex[i]
        out, used = [], 0
        for r in records:
            c['scalar'] += 1
            c[stage + '_tests'] += 1
            if r.origin == i:
                ok = r.coef >= floor if keep_equal else r.coef > floor
                keep = r.mask if ok else 0
                c[stage + '_self'] += 1
            else:
                vo = self.vertex[r.origin]
                theta = r.coef - floor + int(keep_equal)  # keep h iff l_h(vo, vi) < theta
                lo, hi = self._target_range(vo, vi)
                if hi < theta:
                    keep = r.mask
                    c[stage + '_all'] += 1
                elif lo >= theta:
                    keep = 0
                    c[stage + '_none'] += 1
                else:
                    c[stage + '_sliced'] += 1
                    before = c['col']
                    planes = [0] * self.b if vo == vi else self.cache.planes(vo, vi, c)
                    keep = self._threshold(planes, theta, r.mask)
                    c['col_' + stage] += c['col'] - before
            if keep:
                c['mask'] += 2  # write the kept record, accumulate used
                out.append(replace(r, mask=keep))
                used |= keep
        c['mask'] += 1
        if (rest := mask & ~used):
            out.append(Record(floor, i, rest, witness))
            c[stage + '_new_record'] += 1
        self.stats['lambda_records'] += len(records) + len(out)
        return out

    def _masked_minimum(self, planes: list[int], mask: int) -> tuple[int, int]:
        self.cost['col'] += self.b
        self.cost['sliced_minima'] += 1
        return self.minimum(planes, mask)

    def best_state(self, records: list[Record], i: int) -> tuple[int, int]:
        """B_i and the smallest label attaining it; exact best-first with range bounds."""
        c = self.cost
        vi = self.vertex[i]
        cand = []
        for r in records:
            c['scalar'] += 1
            c['min_tests'] += 1
            if r.origin == i:
                lo = hi = 0
            else:
                lo, hi = self._target_range(self.vertex[r.origin], vi)
            cand.append((r.coef - lo, r, lo, hi))
        prune = self.tests != 'none'
        if prune:
            cand.sort(key=lambda z: -z[0])
            c['scalar'] += _sort_cost(len(cand))
        best, bh = None, -1
        for ub, r, lo, hi in cand:
            low = (r.mask & -r.mask).bit_length() - 1
            c['scalar'] += 1
            if prune and best is not None and (ub < best or (ub == best and low > bh)):
                c['min_pruned'] += 1
                continue
            if lo == hi:
                val, h = ub, low
                c['min_self' if r.origin == i else 'min_single'] += 1
            else:
                vo = self.vertex[r.origin]
                before = c['col']
                planes = [0] * self.b if vo == vi else self.cache.planes(vo, vi, c)
                mn, h = self._masked_minimum(planes, r.mask)
                c['col_min'] += c['col'] - before
                val = r.coef - mn
                c['min_sliced'] += 1
            if best is None or val > best or (val == best and h < bh):
                best, bh = val, h
        assert best is not None
        return best, bh

    # ----- the sweep --------------------------------------------------------------
    def _phase(self, name: Optional[str]):
        """Attribute mask, scalar, column and range-label counts to sweep phases."""
        c = self.cost
        now = (c['mask'], c['scalar'], c['col'], c['labels'])
        if getattr(self, '_current', None) is not None:
            for unit, a, b in zip(('mask', 'scalar', 'col', 'labels'), self._mark, now):
                c['%s_%s' % (unit, self._current)] += b - a
        self._current, self._mark = name, now

    def run(self) -> Result:
        c, n = self.cost, self.n
        states = [[] for _ in range(n)]
        Bs, Rs = [0] * n, [None] * n
        events = sorted([(z.c, 0, i) for i, z in enumerate(self.a)] +
                        [(z.d, 1, i) for i, z in enumerate(self.a)])
        c['scalar'] += _sort_cost(2 * n)
        for _, task, i in events:
            z = self.a[i]
            if task == 0:
                k = self.p.A[z.v].bit_count()
                c['ph_local'] += k * per_haplotype_steps(n, self.threshold[i], False)
                c['ph_update'] += k * per_haplotype_steps(n, self.rank[i], True)
                c['ph_states'] += k
                self._phase('remote')
                remote, rw = self.remote_grouped(i) if self.remote_mode == 'grouped' else self.remote(i)
                Rs[i] = remote
                floor, fw = z.weight, Witness('restart')
                c['scalar'] += 2
                if remote is not None and remote > floor:
                    floor, fw = remote, rw
                self._phase('local')
                alpha = z.weight - z.x - z.c + 2
                self._cut = None if self.cut == 'none' else (i, floor - alpha)
                local = self.prefix(self.threshold[i], self.p.A[z.v])
                self._cut = None
                local = [replace(r, coef=r.coef + alpha) for r in local]
                c['scalar'] += len(local)
                self._phase('floor')
                pre = self.floor_states(local, i, floor, fw, keep_equal=False)
                best, bh = self.best_state(pre, i)
                Bs[i] = best
                post = self.floor_states(pre, i, best - self.gamma,
                                         Witness('same', i, bh), keep_equal=True)
                self._phase('output')
                states[i] = coalesce(post)
                c['scalar'] += _sort_cost(len(post))
                c['mask'] += (len(post) - len(states[i])) + len(states[i])  # OR merged, write output
                self.stats['output_classes'] += len(states[i])
                self.stats['origin_value_classes'] += len({(r.coef, r.origin) for r in states[i]})
                self.stats['lambda_records'] += len(post) + len(states[i])
            else:
                self._phase('update')
                incoming = coalesce([replace(r, coef=r.coef + z.y + z.d) for r in states[i]], publish=i)
                c['scalar'] += len(states[i]) + _sort_cost(len(states[i]))
                c['mask'] += (len(states[i]) - len(incoming)) + len(incoming)  # OR merged, union
                self.stats['published_classes'] += len(incoming)
                self.stats['lambda_records'] += len(states[i]) + len(incoming)
                r = self.rank[i]
                while r <= n:
                    c['scalar'] += 1
                    self.nodes[r] = self.merge(self.nodes[r], incoming)
                    self.stats['max_node_classes'] = max(self.stats['max_node_classes'], len(self.nodes[r]))
                    r += r & -r
        self._phase(None)
        self.stats['coordinate_bits'] = self.b
        self.stats['mask_words'] = self.p.W
        self.stats['valid_states'] = sum(self.p.A[z.v].bit_count() for z in self.a)
        c['distinct_merge_pairs'] = len(self.merge_pairs)
        c['distinct_target_pairs'] = len(self.target_pairs)
        c['distinct_pairs'] = len(self.merge_pairs | self.target_pairs)
        c['cached_ranges'] = len(self.cache.ranges)
        c['cached_planes'] = len(self.cache.planes_cache)
        res = Result(self.a, self.p, states, Bs, Rs, self.stats)
        res.cost = self.cost
        return res


def solve_fast(prepared: Prepared, anchors: list[Anchor], gamma: int, **kwargs) -> Result:
    return FastPathSolver(prepared, anchors, gamma, **kwargs).run()


def predicted_word_ops(cost: Counter, W: int, K: int) -> dict:
    """Predicted native word operations under the unit model in the module docstring."""
    compact = W * (cost['mask'] + cost['col']) + cost['scalar'] + cost['labels']
    per_hap = cost['ph_local'] + cost['ph_update'] + cost['ph_remote']
    return dict(class_compact=compact, class_expanded=compact + 2 * K, per_haplotype=per_hap)
