"""Independent brute-force chain semantics for the gap-sensitive contract.

Written from the chain-level definition in PROOF.md section 2 alone. It imports nothing
from the repository: no gap_chaining, gap_fastpath, gap_native, review oracle or index.
It does not assume topological vertex numbering, does not use END ranks, Fenwick trees,
last2reach dynamic programming or the recurrence's algebraic form.

Instance: lengths (positive ints), edges [(u, v)], paths [[v, ...]], anchors
[(v, x, y, c, d, w)], gamma. Anchor ids are list positions.

Definitions (PROOF.md Definitions 2.1-2.6):
  reach(u, v)  some directed route u -> ... -> v exists (u reaches itself);
  D(w, v)      minimum over routes w = u_0 -> ... -> u_k = v of len(u_0)+...+len(u_{k-1});
  phi_h(v)     characters before v on path h;
  A(v)         labels whose path contains v;
  last(v, g)   the vertex of path g with the largest position on g that reaches v.
  Steps of a chain (sigma_{t-1} -> sigma_t):
    local  (j,h) -> (i,h): h in A(v_j) & A(v_i), d_j < c_i, reach(v_j, v_i),
           (v_j != v_i or y_j < x_i); gain w_i - g_h(j,i),
           g_h(j,i) = (phi_h(v_i)+x_i) - (phi_h(v_j)+y_j) + (c_i-d_j) - 2;
    remote (j,g) -> (i,h): g not in A(v_i), w = last(v_i, g) exists, pos_g(v_j) <= pos_g(w),
           d_j < c_i, h in A(v_i); gain w_i - G_g(j,i) - gamma,
           G_g(j,i) = (phi_g(w)+D(w,v_i)+x_i) - (phi_g(v_j)+y_j) + (c_i-d_j) - 2;
    same   (i,h') -> (i,h): h' != h, both in A(v_i); gain -gamma; never two in a row.
  Score: w of the first anchor plus the gains.

Two evaluators:
  enumerate_chains  explicit depth-first enumeration of every chain (exponential; tiny inputs);
  chain_dag_optimum longest path in the explicit chain-step DAG, whose source paths are exactly
                    the chains (layer 0: entered by start/local/remote, layer 1: entered by a
                    same-anchor step). Used to cross-check enumeration and for larger inputs.
Both return C[(i,h)], Cnot[(i,h)] (best chain whose last step is not 'same'), B[i], R[i].
"""
from __future__ import annotations

NEG = None  # "no chain": represented as None, never as a sentinel number


class Facts:
    """Graph facts computed by explicit route enumeration (tiny graphs only)."""

    def __init__(self, lengths, edges, paths):
        V = len(lengths)
        self.V, self.lengths, self.paths = V, list(lengths), [list(p) for p in paths]
        if any((type(z) is not int) or z <= 0 for z in lengths):
            raise ValueError('positive integer vertex lengths required')
        succ = [[] for _ in range(V)]
        for u, v in edges:
            if not (0 <= u < V and 0 <= v < V) or u == v:
                raise ValueError('bad edge')
            succ[u].append(v)
        self.succ = succ
        # D[w][v] by enumerating every route from w (DFS); detect cycles explicitly.
        INF = None
        self.D = [[INF] * V for _ in range(V)]
        for w in range(V):
            stack = [(w, 0, (w,))]
            while stack:
                u, dist, route = stack.pop()
                if self.D[w][u] is None or dist < self.D[w][u]:
                    self.D[w][u] = dist
                for x in succ[u]:
                    if x in route:
                        raise ValueError('graph has a directed cycle')
                    stack.append((x, dist + lengths[u], route + (x,)))
        self.reach = [[self.D[u][v] is not None for v in range(V)] for u in range(V)]
        H = len(paths)
        self.H = H
        self.phi = [dict() for _ in range(H)]
        self.pos = [dict() for _ in range(H)]
        eset = set(edges)
        for h, p in enumerate(self.paths):
            if not p:
                raise ValueError('empty path')
            off = 0
            for k, v in enumerate(p):
                if v in self.pos[h]:
                    raise ValueError('path revisits a vertex')
                if k and (p[k - 1], v) not in eset:
                    raise ValueError('path step is not an edge')
                self.phi[h][v] = off
                self.pos[h][v] = k
                off += lengths[v]
        self.A = [frozenset(h for h in range(H) if v in self.pos[h]) for v in range(V)]
        self.last = [[None] * H for _ in range(V)]
        for v in range(V):
            for g in range(H):
                cands = [w for w in self.paths[g] if self.reach[w][v]]
                if cands:
                    self.last[v][g] = max(cands, key=lambda w: self.pos[g][w])


def check_instance(facts, anchors, gamma):
    if type(gamma) is not int:
        raise ValueError('integer gamma')
    for (v, x, y, c, d, w) in anchors:
        if not (0 <= v < facts.V) or not facts.A[v]:
            raise ValueError('anchor vertex must lie on a haplotype')
        if not (1 <= x <= y <= facts.lengths[v]) or not (1 <= c <= d):
            raise ValueError('bad anchor coordinates')


def step_targets(facts, anchors, gamma, j, h):
    """All (target_state, gain, kind, extra) chain steps leaving state (j, h) to a later anchor.

    kind is 'local' or 'remote'; extra is w (the last2reach vertex) for remote steps.
    """
    vj, xj, yj, cj, dj, wj = anchors[j]
    out = []
    for i, (vi, xi, yi, ci, di, wi) in enumerate(anchors):
        if not dj < ci:
            continue
        if h in facts.A[vi]:
            if facts.reach[vj][vi] and (vj != vi or yj < xi):
                gap = (facts.phi[h][vi] + xi) - (facts.phi[h][vj] + yj) + (ci - dj) - 2
                out.append(((i, h), wi - gap, 'local', None))
        else:
            w = facts.last[vi][h]
            if w is not None and facts.pos[h][vj] <= facts.pos[h][w]:
                gap = (facts.phi[h][w] + facts.D[w][vi] + xi) - (facts.phi[h][vj] + yj) + (ci - dj) - 2
                for h2 in sorted(facts.A[vi]):
                    out.append(((i, h2), wi - gap - gamma, 'remote', w))
    return out


def step_table(facts, anchors, gamma):
    """step_targets for every state, computed once per instance (performance only)."""
    return {(j, h): step_targets(facts, anchors, gamma, j, h)
            for j in range(len(anchors)) for h in sorted(facts.A[anchors[j][0]])}


def _better(a, b):
    return b is None or a > b


def enumerate_chains(facts, anchors, gamma, cap=None):
    """Explicitly enumerate every chain. Returns (C, Cnot, B, R, number_of_chains).

    cap: stop and return None if more than cap chains would be enumerated.
    """
    check_instance(facts, anchors, gamma)
    n = len(anchors)
    C, Cnot, R = {}, {}, [None] * n
    succ_cache = {}
    count = 0
    stack = []
    for i in range(n):
        for h in sorted(facts.A[anchors[i][0]]):
            stack.append(((i, h), anchors[i][5], 'start'))
    while stack:
        (i, h), score, kind = stack.pop()
        count += 1
        if cap is not None and count > cap:
            return None
        key = (i, h)
        if _better(score, C.get(key)):
            C[key] = score
        if kind != 'same' and _better(score, Cnot.get(key)):
            Cnot[key] = score
        if kind == 'remote' and _better(score, R[i]):
            R[i] = score
        if kind != 'same':
            for h2 in sorted(facts.A[anchors[i][0]]):
                if h2 != h:
                    stack.append(((i, h2), score - gamma, 'same'))
        if key not in succ_cache:
            succ_cache[key] = step_targets(facts, anchors, gamma, i, h)
        for tgt, gain, k2, _w in succ_cache[key]:
            stack.append((tgt, score + gain, k2))
    B = [max(C[(i, h)] for h in facts.A[anchors[i][0]]) for i in range(n)]
    return C, Cnot, B, R, count


def chain_dag_optimum(facts, anchors, gamma, table=None):
    """Longest path in the explicit chain-step DAG (same semantics, polynomial time).

    Nodes: (i, h, layer); layer 0 entered by start/local/remote, layer 1 by a same step.
    Edges follow step_targets from both layers and same steps from layer 0 only.
    A chain is exactly a path from the virtual source (start edges into layer 0).
    """
    check_instance(facts, anchors, gamma)
    n = len(anchors)
    best = {}
    R = [None] * n
    for i in range(n):
        for h in facts.A[anchors[i][0]]:
            best[(i, h, 0)] = anchors[i][5]
    # Every local/remote step strictly increases c (c_j <= d_j < c_i): process anchors by c,
    # and within one anchor layer 0 before layer 1. Incoming steps into anchor i come only from
    # anchors with smaller c, all finalized before i is processed.
    order = sorted(range(n), key=lambda t: (anchors[t][3], t))
    incoming = {i: [] for i in range(n)}
    for j in range(n):
        for h in facts.A[anchors[j][0]]:
            succ = table[(j, h)] if table is not None else step_targets(facts, anchors, gamma, j, h)
            for (i, h2), gain, kind, _w in succ:
                incoming[i].append((j, h, h2, gain, kind))
    for i in order:
        for (j, h, h2, gain, kind) in incoming[i]:
            src = [best[(j, h, L)] for L in (0, 1) if (j, h, L) in best]
            for s in src:
                val = s + gain
                if _better(val, best.get((i, h2, 0))):
                    best[(i, h2, 0)] = val
                if kind == 'remote' and _better(val, R[i]):
                    R[i] = val
        labs = sorted(facts.A[anchors[i][0]])
        for h in labs:
            for h2 in labs:
                if h2 != h:
                    val = best[(i, h, 0)] - gamma
                    if _better(val, best.get((i, h2, 1))):
                        best[(i, h2, 1)] = val
    C, Cnot = {}, {}
    for i in range(n):
        for h in facts.A[anchors[i][0]]:
            Cnot[(i, h)] = best[(i, h, 0)]
            C[(i, h)] = max(v for v in (best.get((i, h, 0)), best.get((i, h, 1))) if v is not None)
    B = [max(C[(i, h)] for h in facts.A[anchors[i][0]]) for i in range(n)]
    return C, Cnot, B, R


def trace_and_score(facts, anchors, gamma, C, pointers, table=None):
    """Independently check one solver's pointers under the chain semantics.

    C: {(i, h): value}; pointers: {(i, h): (kind, j, g, w)} with kind in
    restart/local/remote/same (local: g is ignored, the label is the current one; remote: w
    must equal last(v_i, g) as computed here). For every state, follow pointers back to a
    restart; every step must be a legal chain step, no two same-anchor steps may be
    consecutive, there must be no cycle, and the chain's score, recomputed from the chain
    definition, must equal C[(i, h)]. Returns a list of problems (empty if none).
    """
    problems = []
    memo = {}
    n = len(anchors)

    def gain_of(src, dst, kind, w):
        succ = table[src] if table is not None else step_targets(facts, anchors, gamma, src[0], src[1])
        for tgt, gain, k2, wv in succ:
            if tgt == dst and k2 == kind and (kind == 'local' or wv == w):
                return gain
        return None

    def score(s, visiting):
        if s in memo:
            return memo[s]
        if s in visiting:
            problems.append(('cycle', s))
            return None
        visiting.add(s)
        i, h = s
        kind, j, g, w = pointers[s]
        val = None
        if kind == 'restart':
            val = anchors[i][5]
        elif kind == 'same':
            if j != i or g == h or g not in facts.A[anchors[i][0]]:
                problems.append(('illegal same', s, pointers[s]))
            elif pointers[(i, g)][0] == 'same':
                problems.append(('consecutive same', s, pointers[s]))
            else:
                base = score((i, g), visiting)
                val = None if base is None else base - gamma
        elif kind in ('local', 'remote') and 0 <= j < n:
            src = (j, h) if kind == 'local' else (j, g)
            if src not in C:
                problems.append(('source state invalid', s, pointers[s]))
            else:
                gain = gain_of(src, s, kind, w)
                if gain is None:
                    problems.append(('illegal ' + kind, s, pointers[s]))
                else:
                    base = score(src, visiting)
                    val = None if base is None else base + gain
        else:
            problems.append(('bad pointer', s, pointers[s]))
        visiting.discard(s)
        memo[s] = val
        return val

    for s in sorted(C):
        v = score(s, set())
        if v is not None and v != C[s]:
            problems.append(('score', s, v, C[s]))
    return problems
