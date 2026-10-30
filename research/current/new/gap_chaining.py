"""Exact origin-class solver for PROBLEM.md, including remote and same-anchor switches.

Pure-Python executable specification: Python integers implement H-bit masks, not
unit-cost H-lane operations. In the analysis each mask operation costs ceil(H/64).
No external dependencies. Vertex ids MUST already be a topological ordering.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
from bisect import bisect_left, bisect_right
from typing import Optional, Iterable


def labels(mask: int) -> Iterable[int]:
    while mask:
        b = mask & -mask
        yield b.bit_length() - 1
        mask ^= b


@dataclass(frozen=True)
class Anchor:
    v: int
    x: int
    y: int
    c: int
    d: int
    weight: int


@dataclass(frozen=True)
class Witness:
    kind: str  # restart, local, remote, same
    j: int = -1
    g: int = -1  # local: -1 means the current label
    w: int = -1  # last2reach vertex for remote switches


@dataclass(frozen=True)
class Record:
    coef: int
    origin: int  # an anchor id, not a vertex id
    mask: int
    witness: Witness


@dataclass
class Graph:
    lengths: list[int]
    edges: list[tuple[int, int]]
    paths: list[list[int]]

    def prepare(self) -> 'Prepared':
        return Prepared(self)


class Prepared:
    """Reusable O(H(V+E)) graph index, plus transposed coordinate bitplanes.

    The Python transpose below uses O(H V b) elementary operations. A native
    word-matrix transpose reduces this term; it is entirely reusable graph prep.
    """
    def __init__(self, graph: Graph):
        self.graph = graph
        self.V, self.H = len(graph.lengths), len(graph.paths)
        if not self.V or not self.H or any(type(z) is not int or z <= 0 for z in graph.lengths):
            raise ValueError('Nonempty graph/panel and positive integer vertex lengths required')
        eset = set(graph.edges)
        if len(eset) != len(graph.edges):
            raise ValueError('Duplicate edge')
        if any(not (0 <= u < v < self.V) for u, v in eset):
            raise ValueError('Edges must follow the supplied topological vertex numbering')
        self.pred = [[] for _ in range(self.V)]
        for u, v in graph.edges:
            self.pred[v].append(u)
        for a in self.pred:
            a.sort()
        self.A = [0] * self.V
        self.phi = [[0] * self.H for _ in range(self.V)]
        self.on = [[False] * self.H for _ in range(self.V)]
        for h, path in enumerate(graph.paths):
            if not path or any(not 0 <= v < self.V for v in path):
                raise ValueError('Invalid haplotype path')
            if any((u, v) not in eset for u, v in zip(path, path[1:])):
                raise ValueError('Consecutive haplotype vertices must be graph edges')
            off = 0
            for v in path:
                self.phi[v][h] = off
                self.on[v][h] = True
                self.A[v] |= 1 << h
                off += graph.lengths[v]
        self.U = (1 << self.H) - 1
        self.W = (self.H + 63) // 64
        self.S = sum(graph.lengths)  # safe upper bound on every absolute coordinate and route length
        self.b = max(1, self.S.bit_length())
        self.phibits = [[0] * self.b for _ in range(self.V)]
        for v in range(self.V):
            for h in labels(self.A[v]):
                off = self.phi[v][h]
                for k in range(self.b):
                    if (off >> k) & 1:
                        self.phibits[v][k] |= 1 << h
        self.last = [[-1] * self.H for _ in range(self.V)]
        self.distance = [[None] * self.H for _ in range(self.V)]
        self.dparent = [[-1] * self.H for _ in range(self.V)]
        for h in range(self.H):
            for v in range(self.V):
                if self.on[v][h]:
                    self.last[v][h], self.distance[v][h] = v, 0
                    continue
                w = max((self.last[u][h] for u in self.pred[v]), default=-1)
                if w < 0:
                    continue
                choices = [(self.distance[u][h] + graph.lengths[u], u)
                           for u in self.pred[v] if self.last[u][h] == w]
                dist, par = min(choices)
                self.last[v][h], self.distance[v][h], self.dparent[v][h] = w, dist, par

    def check_anchors(self, anchors: list[Anchor], gamma: int):
        if type(gamma) is not int or gamma < 0:
            raise ValueError('gamma must be a nonnegative integer')
        for a in anchors:
            if any(type(t) is not int for t in (a.v, a.x, a.y, a.c, a.d, a.weight)):
                raise ValueError('Integer anchors required')
            if not (0 <= a.v < self.V and 1 <= a.x <= a.y <= self.graph.lengths[a.v]
                    and 1 <= a.c <= a.d and a.weight >= 0 and self.A[a.v]):
                raise ValueError(f'Invalid or uncovered anchor: {a}')


class Selector:
    """Balanced tree over DISJOINT record masks. O(log(number of records)) label lookup.

    No dense per-label owner table is constructed.
    """
    def __init__(self, records: list[Record]):
        self.records = records
        self.n = len(records)
        self.size = 1 << max(0, (self.n - 1).bit_length())
        self.tree = [0] * (2 * self.size)
        for i, r in enumerate(records):
            self.tree[self.size + i] = r.mask
        for i in range(self.size - 1, 0, -1):
            self.tree[i] = self.tree[2*i] | self.tree[2*i+1]

    def get(self, h: int) -> tuple[Optional[Record], int]:
        bit = 1 << h
        if not self.n or not self.tree[1] & bit:
            return None, 1
        p, steps = 1, 1
        while p < self.size:
            p = 2*p if self.tree[2*p] & bit else 2*p+1
            steps += 1
        return self.records[p-self.size], steps


def coalesce(records: list[Record], publish: Optional[int] = None) -> list[Record]:
    """Deterministic grouping; never combines different immediate witnesses."""
    items = []
    for r in records:
        if not r.mask:
            continue
        w = r.witness if publish is None else Witness('local', publish)
        key = (r.coef, r.origin, w.kind, w.j, w.g, w.w)
        items.append((key, r.mask, w))
    items.sort(key=lambda z: z[0])
    ans = []
    for key, mask, w in items:
        if ans and (ans[-1].coef, ans[-1].origin, ans[-1].witness) == (key[0], key[1], w):
            ans[-1] = replace(ans[-1], mask=ans[-1].mask | mask)
        else:
            ans.append(Record(key[0], key[1], mask, w))
    return ans


@dataclass
class Result:
    anchors: list[Anchor]
    prepared: Prepared
    records: list[list[Record]]
    B: list[int]
    R: list[Optional[int]]
    stats: dict

    def expand(self):
        """One enumeration per output mask; O(K) writes after mask traversal.

        Mask-word traversal costs O(W D_out) and is already charged to output
        record processing in the theorem; no score recomputation DP is performed.
        """
        C, pointers = [], []
        p = self.prepared
        for i, rs in enumerate(self.records):
            row, ptr = {}, {}
            for r in rs:
                ov = self.anchors[r.origin].v
                for h in labels(r.mask):
                    if h in row:
                        raise AssertionError('Overlapping output classes')
                    row[h] = r.coef - (p.phi[self.anchors[i].v][h] - p.phi[ov][h])
                    ptr[h] = r.witness
            if sum(1 << h for h in row) != p.A[self.anchors[i].v]:
                raise AssertionError('Incomplete output partition')
            C.append(row)
            pointers.append(ptr)
        return C, pointers


class OriginSolver:
    """Origin-class Fenwick tree with bit-sliced comparisons and route minima.

    Relative coding is used for short score ranges. Otherwise coefficients
    separated by more than S are encoded in ordered clusters. Coordinate b and
    log(number of merge records), not absolute integer score magnitude, bound
    the bit-sliced comparison width. Scalar big-integer costs remain explicit.
    """
    def __init__(self, prepared: Prepared, anchors: list[Anchor], gamma: int):
        prepared.check_anchors(anchors, gamma)
        self.p, self.a, self.gamma = prepared, anchors, gamma
        self.n = len(anchors)
        self.U, self.b = prepared.U, prepared.b
        self.stats = dict(merges=0, lambda_records=0, coded_records=0,
                          coded_columns=0, route_calls=0, route_columns=0,
                          selector_steps=0, remote_nodes=0, published_classes=0,
                          output_classes=0, max_node_classes=0, scalar_fast_merges=0,
                          missing_incidences=0, origin_value_classes=0,
                          threshold_calls=0, minimum_calls=0,
                          direct_coded_merges=0, clustered_merges=0)
        ends = sorted((z.v, z.y, i) for i, z in enumerate(anchors))
        self.ends = ends
        self.rank = [0] * self.n
        for r, (_, _, i) in enumerate(ends, 1):
            self.rank[i] = r
        self.threshold = [bisect_left(ends, (z.v, z.x, -1)) for z in anchors]
        self.nodes = [[] for _ in range(self.n + 1)]
        self.selectors = [None] * (self.n + 1)
        self.route_cache = {}

    def add_const(self, bits: list[int], k: int, width: int) -> list[int]:
        if k < 0:
            raise ValueError('Internal add_const expects nonnegative coding constants')
        out, carry = [], 0
        for t in range(width):
            x = bits[t] if t < len(bits) else 0
            y = self.U if (k >> t) & 1 else 0
            xy = x ^ y
            out.append(xy ^ carry)
            carry = (x & y) | (carry & xy)
        if carry:
            raise AssertionError('Coding overflow')
        return out

    def gt(self, left: list[int], right: list[int], mask: int) -> int:
        eq, win = mask, 0
        for x, y in zip(reversed(left), reversed(right)):
            win |= eq & x & ~y
            eq &= ~(x ^ y)
            if not eq:
                break
        return win

    def route(self, origin: int, target: int) -> list[int]:
        key = (origin, target)
        if key in self.route_cache:
            return self.route_cache[key]
        ov, tv = self.a[origin].v, self.a[target].v
        self.stats['route_calls'] += 1
        self.stats['route_columns'] += self.b
        if ov == tv:
            ans = [0] * self.b
        else:
            left, right = self.p.phibits[tv], self.p.phibits[ov]
            borrow, ans = 0, []
            for x, y in zip(left, right):
                ans.append((x ^ y ^ borrow) & self.U)
                borrow = ((~x & (y | borrow)) | (y & borrow)) & self.U
        # Values outside the caller's legal mask are immaterial and may wrap.
        self.route_cache[key] = ans
        return ans

    def less_than(self, bits: list[int], threshold: int, mask: int) -> int:
        self.stats['threshold_calls'] += 1
        if threshold <= 0:
            return 0
        if threshold > self.p.S:
            return mask
        rhs = [self.U if (threshold >> k) & 1 else 0 for k in range(self.b)]
        return self.gt(rhs, bits, mask)

    def minimum(self, bits: list[int], mask: int) -> tuple[int, int]:
        self.stats['minimum_calls'] += 1
        if not mask:
            raise ValueError('Minimum over empty mask')
        val, live = 0, mask
        for k in range(len(bits)-1, -1, -1):
            zero = live & ~bits[k]
            if zero:
                live = zero
            else:
                val |= 1 << k
                live &= bits[k]
        return val, (live & -live).bit_length()-1

    def merge(self, left: list[Record], right: list[Record]) -> list[Record]:
        """Pointwise maximum; ties go to the left list. Never enumerate labels."""
        m = len(left) + len(right)
        self.stats['merges'] += 1
        self.stats['lambda_records'] += m
        if not left or not right:
            ans = list(left or right)
            self.stats['lambda_records'] += len(ans)
            return ans
        lm, rm = 0, 0
        for r in left:
            lm |= r.mask
        for r in right:
            rm |= r.mask
        overlap = lm & rm
        all_records = left + right
        if not overlap:
            ans = all_records
        elif len({r.origin for r in all_records}) == 1:
            # Coordinate term is common: ordinary colored maxima.
            self.stats['scalar_fast_merges'] += 1
            used, ans = 0, []
            for r in sorted(all_records, key=lambda r: -r.coef):
                keep = r.mask & ~used
                if keep:
                    ans.append(replace(r, mask=keep))
                    used |= keep
        else:
            lo = min(r.coef for r in all_records)
            hi = max(r.coef for r in all_records)
            direct_width = max(1, (hi-lo+self.p.S).bit_length())
            if direct_width <= self.b + 2*(m+1).bit_length() + 2:
                offsets = [r.coef-lo for r in all_records]
                width = direct_width
                self.stats['direct_coded_merges'] += 1
            else:
                order = sorted(range(m), key=lambda i: all_records[i].coef)
                offsets = [0] * m
                cluster, base, prev = -1, 0, None
                shift = max(1, (m*self.p.S + 1).bit_length())
                for idx in order:
                    k = all_records[idx].coef
                    if prev is None or k-prev > self.p.S:
                        cluster += 1
                        base = k
                    offsets[idx] = (cluster << shift) + k - base
                    prev = k
                width = max(1, (max(offsets) + self.p.S).bit_length())
                self.stats['clustered_merges'] += 1
            lp, rp = [0]*width, [0]*width
            for idx, rec in enumerate(all_records):
                pb = self.p.phibits[self.a[rec.origin].v]
                eb = self.add_const(pb, offsets[idx], width)
                dest = lp if idx < len(left) else rp
                for t in range(width):
                    dest[t] |= eb[t] & rec.mask
            right_wins = (rm & ~lm) | self.gt(rp, lp, overlap)
            ans = [replace(r, mask=keep) for r in left
                   if (keep := r.mask & ~right_wins)]
            ans += [replace(r, mask=keep) for r in right
                    if (keep := r.mask & right_wins)]
            self.stats['coded_records'] += m
            self.stats['coded_columns'] += width*m
        self.stats['lambda_records'] += len(ans)
        return ans

    def prefix(self, r: int, mask: int) -> list[Record]:
        ans = []
        while r:
            self.stats['lambda_records'] += len(self.nodes[r])
            filtered = [replace(z, mask=m) for z in self.nodes[r] if (m := z.mask & mask)]
            ans = self.merge(ans, filtered)
            r -= r & -r
        return ans

    def remote(self, i: int) -> tuple[Optional[int], Optional[Witness]]:
        z, p = self.a[i], self.p
        value, witness = None, None
        missing = self.U & ~p.A[z.v]
        self.stats['missing_incidences'] += missing.bit_count()
        for g in labels(missing):
            w = p.last[z.v][g]
            if w < 0:
                continue
            r = bisect_right(self.ends, (w, p.graph.lengths[w], self.n))
            best, bj = None, -1
            while r:
                self.stats['remote_nodes'] += 1
                sel = self.selectors[r]
                if sel is not None:
                    rec, steps = sel.get(g)
                    self.stats['selector_steps'] += steps
                    if rec is not None:
                        v = rec.coef + p.phi[self.a[rec.origin].v][g]
                        if best is None or v > best:
                            best, bj = v, rec.witness.j
                r -= r & -r
            if best is not None:
                cand = (best + z.weight - z.x - z.c + 2
                        - p.phi[w][g] - p.distance[z.v][g] - self.gamma)
                if value is None or cand > value:
                    value, witness = cand, Witness('remote', bj, g, w)
        return value, witness

    def floor_states(self, records: list[Record], i: int, floor: int,
                     witness: Witness, keep_equal: bool) -> list[Record]:
        mask = self.p.A[self.a[i].v]
        out, used = [], 0
        for r in records:
            if r.origin == i:
                ok = r.coef >= floor if keep_equal else r.coef > floor
                keep = r.mask if ok else 0
            else:
                route = self.route(r.origin, i)
                keep = self.less_than(route, r.coef-floor+int(keep_equal), r.mask)
            if keep:
                out.append(replace(r, mask=keep))
                used |= keep
        if (rest := mask & ~used):
            out.append(Record(floor, i, rest, witness))
        self.stats['lambda_records'] += len(records) + len(out)
        return out

    def run(self) -> Result:
        states = [[] for _ in range(self.n)]
        Bs, Rs = [0]*self.n, [None]*self.n
        events = sorted([(z.c, 0, i) for i, z in enumerate(self.a)] +
                        [(z.d, 1, i) for i, z in enumerate(self.a)])
        for _, task, i in events:
            z = self.a[i]
            if task == 0:
                self.route_cache = {}  # no pairwise-origin preprocessing table
                remote, rw = self.remote(i)
                Rs[i] = remote
                floor, fw = z.weight, Witness('restart')
                if remote is not None and remote > floor:
                    floor, fw = remote, rw
                local = self.prefix(self.threshold[i], self.p.A[z.v])
                alpha = z.weight-z.x-z.c+2
                local = [replace(r, coef=r.coef+alpha) for r in local]
                pre = self.floor_states(local, i, floor, fw, keep_equal=False)
                best, bh = None, -1
                self.stats['lambda_records'] += len(pre)
                for r in pre:
                    if r.origin == i:
                        mn, h = 0, (r.mask & -r.mask).bit_length()-1
                    else:
                        mn, h = self.minimum(self.route(r.origin, i), r.mask)
                    val = r.coef-mn
                    if best is None or val > best or (val == best and h < bh):
                        best, bh = val, h
                assert best is not None
                Bs[i] = best
                post = self.floor_states(pre, i, best-self.gamma,
                                         Witness('same', i, bh), keep_equal=True)
                states[i] = coalesce(post)
                self.stats['output_classes'] += len(states[i])
                self.stats['origin_value_classes'] += len({(r.coef, r.origin) for r in states[i]})
                self.stats['lambda_records'] += len(post) + len(states[i])
            else:
                incoming = coalesce([replace(r, coef=r.coef+z.y+z.d) for r in states[i]], publish=i)
                self.stats['published_classes'] += len(incoming)
                self.stats['lambda_records'] += len(states[i]) + len(incoming)
                r = self.rank[i]
                while r <= self.n:
                    self.nodes[r] = self.merge(self.nodes[r], incoming)
                    self.selectors[r] = Selector(self.nodes[r])
                    self.stats['lambda_records'] += len(self.nodes[r])  # selector construction
                    self.stats['max_node_classes'] = max(self.stats['max_node_classes'], len(self.nodes[r]))
                    r += r & -r
        self.stats['coordinate_bits'] = self.b
        self.stats['mask_words'] = self.p.W
        self.stats['xi_calls'] = (self.stats['route_calls'] + self.stats['threshold_calls']
                                  + self.stats['minimum_calls'])
        self.stats['valid_states'] = sum(self.p.A[z.v].bit_count() for z in self.a)
        return Result(self.a, self.p, states, Bs, Rs, self.stats)


def solve(prepared: Prepared, anchors: list[Anchor], gamma: int) -> Result:
    return OriginSolver(prepared, anchors, gamma).run()


def direct_oracle(p: Prepared, a: list[Anchor], gamma: int):
    """Independent explicit predecessor-state recurrence; O(N^2 H)."""
    p.check_anchors(a, gamma)
    n = len(a)
    C, pointers = [dict() for _ in a], [dict() for _ in a]
    B, R = [0]*n, [None]*n
    for i in sorted(range(n), key=lambda i: (a[i].c, i)):
        z = a[i]
        rv, rw = None, None
        for g in range(p.H):
            if p.on[z.v][g]:
                continue
            w = p.last[z.v][g]
            if w < 0:
                continue
            for j, u in enumerate(a):
                if u.d >= z.c or not p.on[u.v][g] or u.v > w:
                    continue
                gap = ((p.phi[w][g] + p.distance[z.v][g] + z.x)
                       -(p.phi[u.v][g]+u.y)+(z.c-u.d)-2)
                v = C[j][g]+z.weight-gap-gamma
                if rv is None or v > rv:
                    rv, rw = v, Witness('remote', j, g, w)
        R[i] = rv
        pre, pp = {}, {}
        for h in range(p.H):
            if not p.on[z.v][h]:
                continue
            best, ptr = z.weight, Witness('restart')
            if rv is not None and rv > best:
                best, ptr = rv, rw
            for j, u in enumerate(a):
                if not (u.d < z.c and p.on[u.v][h]
                        and (u.v < z.v or (u.v == z.v and u.y < z.x))):
                    continue
                gap = ((p.phi[z.v][h]+z.x)-(p.phi[u.v][h]+u.y)+(z.c-u.d)-2)
                v = C[j][h]+z.weight-gap
                if v > best:
                    best, ptr = v, Witness('local', j)
            pre[h], pp[h] = best, ptr
        B[i] = max(pre.values())
        bh = min(h for h in pre if pre[h] == B[i])
        for h in pre:
            if B[i]-gamma > pre[h]:
                C[i][h], pointers[i][h] = B[i]-gamma, Witness('same', i, bh)
            else:
                C[i][h], pointers[i][h] = pre[h], pp[h]
    return C, B, R, pointers


def scalar_fenwick(p: Prepared, a: list[Anchor], gamma: int):
    """Per-haplotype O((K+M) log N) reference, independent of origin normalization."""
    p.check_anchors(a, gamma)
    n = len(a)
    ends = sorted((z.v, z.y, i) for i, z in enumerate(a))
    rank = [0]*n
    for r, (_, _, i) in enumerate(ends, 1):
        rank[i] = r
    trees = [dict() for _ in range(p.H)]
    C, ptrs, B, R = [dict() for _ in a], [dict() for _ in a], [0]*n, [None]*n
    def query(h, r):
        val, j = None, -1
        while r:
            item = trees[h].get(r)
            if item is not None and (val is None or item[0] > val):
                val, j = item
            r -= r & -r
        return val, j
    for _, task, i in sorted([(z.c, 0, i) for i, z in enumerate(a)] +
                             [(z.d, 1, i) for i, z in enumerate(a)]):
        z = a[i]
        if task:
            for h in labels(p.A[z.v]):
                value = C[i][h]+z.y+z.d+p.phi[z.v][h]
                r = rank[i]
                while r <= n:
                    old = trees[h].get(r)
                    if old is None or value > old[0]:
                        trees[h][r] = (value, i)
                    r += r & -r
            continue
        rv, rw = None, None
        for g in labels(p.U & ~p.A[z.v]):
            w = p.last[z.v][g]
            if w < 0:
                continue
            f, j = query(g, bisect_right(ends, (w, p.graph.lengths[w], n)))
            if f is not None:
                cand = f+z.weight-z.x-z.c+2-p.phi[w][g]-p.distance[z.v][g]-gamma
                if rv is None or cand > rv:
                    rv, rw = cand, Witness('remote', j, g, w)
        R[i] = rv
        floor, fw = z.weight, Witness('restart')
        if rv is not None and rv > floor:
            floor, fw = rv, rw
        for h in labels(p.A[z.v]):
            f, j = query(h, bisect_left(ends, (z.v, z.x, -1)))
            local = None if f is None else f+z.weight-z.x-z.c+2-p.phi[z.v][h]
            if local is not None and local > floor:
                C[i][h], ptrs[i][h] = local, Witness('local', j)
            else:
                C[i][h], ptrs[i][h] = floor, fw
        B[i] = max(C[i].values())
        bh = min(h for h in C[i] if C[i][h] == B[i])
        for h in C[i]:
            if B[i]-gamma > C[i][h]:
                C[i][h], ptrs[i][h] = B[i]-gamma, Witness('same', i, bh)
    return C, B, R, ptrs


def verify_pointers(p: Prepared, a: list[Anchor], gamma: int, C, B, pointers) -> dict:
    """Checks legality, exact score, and full traceback acyclicity for EVERY state."""
    counts = dict(restart=0, local=0, remote=0, same=0)
    parent = {}
    for i, z in enumerate(a):
        assert B[i] == max(C[i].values())
        for h in labels(p.A[z.v]):
            wit = pointers[i][h]
            counts[wit.kind] += 1
            key = (i, h)
            if wit.kind == 'restart':
                assert C[i][h] == z.weight
                parent[key] = None
            elif wit.kind == 'same':
                assert wit.j == i and wit.g != h and p.on[z.v][wit.g]
                assert C[i][h] == C[i][wit.g]-gamma
                parent[key] = (i, wit.g)
            else:
                j = wit.j
                u = a[j]
                g = h if wit.kind == 'local' else wit.g
                assert 0 <= j < len(a) and p.on[u.v][g] and u.d < z.c
                if wit.kind == 'local':
                    assert p.on[z.v][h] and (u.v < z.v or (u.v == z.v and u.y < z.x))
                    gap = p.phi[z.v][h]+z.x-p.phi[u.v][h]-u.y+z.c-u.d-2
                    value = C[j][h]+z.weight-gap
                else:
                    assert wit.kind == 'remote' and not p.on[z.v][g]
                    w = p.last[z.v][g]
                    assert wit.w == w and w >= 0 and u.v <= w
                    gap = p.phi[w][g]+p.distance[z.v][g]+z.x-p.phi[u.v][g]-u.y+z.c-u.d-2
                    value = C[j][g]+z.weight-gap-gamma
                assert value == C[i][h], (i, h, wit, value, C[i][h])
                parent[key] = (j, g)
    done = set()
    for key in parent:
        path, seen, cur = [], set(), key
        while cur is not None and cur not in done:
            assert cur not in seen, ('cycle', cur, path)
            seen.add(cur)
            path.append(cur)
            cur = parent[cur]
        done.update(path)
    return counts
