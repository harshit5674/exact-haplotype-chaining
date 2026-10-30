"""Python side of the Giraffe emulation (giraffe_emul.cpp) and independent checks.

Standard library only. The native library is built into experiments/exactness_payoff/build/.
Anchors are the project's (v, x, y, c, d, w) tuples, 1-based closed, inside one vertex.

Giraffe hifi parameters (vg v1.77.0; SPEC_GIRAFFE.md section 3):
  item points = match(1) * span + item_bonus(2); gap_scale 0.27579; rho 2; consistency bonus 12;
  graph lookback max(20000, 0.10501002120802233 * read length); no read lookback;
  indel limit max(5000, 2.45 * read length); L = mean seed length (17 here).
"""
from __future__ import annotations

from array import array
import ctypes as ct
import hashlib
import heapq
import json
import math
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
NEG = -(1 << 29)
SIZE_MAX = (1 << 64) - 1

HIFI = dict(item_bonus=2, match=1, gap_scale=0.27579, rho=2, consistency_bonus=12,
            graph_lookback=20000, graph_lookback_per_base=0.10501002120802233,
            read_lookback=None, indel=5000, indel_per_base=2.45)


class Params(ct.Structure):
    _fields_ = [('graph_lookback', ct.c_int64), ('read_lookback', ct.c_int64), ('indel_limit', ct.c_int64),
                ('base_seed_length', ct.c_int32), ('rho', ct.c_int32), ('gap_scale', ct.c_double),
                ('n_bonus', ct.c_int32), ('bonus', ct.c_int32 * 4), ('do_exact', ct.c_int32),
                ('do_heuristic', ct.c_int32), ('prune', ct.c_int32)]


_LIB = None


def build():
    global _LIB
    if _LIB is not None:
        return _LIB
    src = HERE / 'giraffe_emul.cpp'
    flags = ['-std=c++20', '-O3', '-shared', '-fPIC']
    ident = hashlib.sha256(src.read_bytes() + repr(flags).encode()).hexdigest()[:16]
    target = HERE / 'build' / f'giraffe_emul-{ident}.dylib'
    target.parent.mkdir(exist_ok=True)
    if not target.exists():
        tmp = target.with_suffix('.tmp')
        subprocess.run(['clang++', *flags, str(src), '-o', str(tmp)], check=True)
        tmp.replace(target)
    lib = ct.CDLL(str(target))
    P = ct.POINTER
    i32, i64, u64, i16, dbl = P(ct.c_int32), P(ct.c_int64), P(ct.c_uint64), P(ct.c_int16), P(ct.c_double)
    lib.ge_error.restype = ct.c_char_p
    lib.ge_graph_new.argtypes = [ct.c_int, ct.c_int, ct.c_int, i64, ct.c_int, i32, i32, u64]
    lib.ge_graph_new.restype = ct.c_void_p
    lib.ge_graph_free.argtypes = [ct.c_void_p]
    lib.ge_graph_dist.argtypes = [ct.c_void_p, ct.c_int, ct.c_int]
    lib.ge_graph_dist.restype = ct.c_int32
    lib.ge_solve.argtypes = [ct.c_void_p, ct.c_int, i64, i32, P(Params), i32, i32, u64, i32, i32, i32, i16, i32, dbl]
    lib.ge_solve.restype = ct.c_int
    lib.ge_solve_explicit.argtypes = [ct.c_int, ct.c_int, ct.c_int, i32, u64, i64, i32, i32, P(Params),
                                      i32, i32, u64, i32, i32, i32, i16, i32]
    lib.ge_solve_explicit.restype = ct.c_int
    lib.ge_transitions.argtypes = [ct.c_void_p, ct.c_int, i64, P(Params), ct.c_int64, i32, i32, i32, i32]
    lib.ge_transitions.restype = ct.c_int64
    lib.identity = target.name
    _LIB = lib
    return lib


def buf(a):
    typ = {'i': ct.c_int32, 'q': ct.c_int64, 'Q': ct.c_uint64, 'h': ct.c_int16, 'd': ct.c_double}[a.typecode]
    assert a.itemsize == ct.sizeof(typ)
    return (typ * len(a)).from_buffer(a) if len(a) else None


def words(mask: int, W: int):
    return [(mask >> (64 * q)) & 0xFFFFFFFFFFFFFFFF for q in range(W)]


# ---------------------------------------------------------------- graph

class Panel:
    """Graph JSON (lengths, edges, paths in topological vertex order) with masks and the native index."""

    def __init__(self, graph: dict, name: str = ''):
        self.name = name
        self.lengths = list(graph['lengths'])
        self.edges = [tuple(e) for e in graph['edges']]
        self.paths = [list(p) for p in graph['paths']]
        self.hap_names = list(graph.get('haplotypes', [str(h) for h in range(len(self.paths))]))
        self.V, self.H = len(self.lengths), len(self.paths)
        self.W = max(1, (self.H + 63) // 64)
        self.A = [0] * self.V
        for h, p in enumerate(self.paths):
            for v in p:
                self.A[v] |= 1 << h
        assert all(self.A), 'every vertex must lie on a haplotype'
        self.out = [[] for _ in range(self.V)]
        for u, v in self.edges:
            assert 0 <= u < v < self.V, 'vertex IDs must be topological'
            self.out[u].append(v)
        self._handle = None
        self._dist_cache = {}

    @property
    def handle(self):
        if self._handle is None:
            lib = build()
            lengths = array('q', self.lengths)
            eu = array('i', [u for u, _ in self.edges])
            ev = array('i', [v for _, v in self.edges])
            vm = array('Q', [w for v in range(self.V) for w in words(self.A[v], self.W)])
            h = lib.ge_graph_new(self.V, self.H, self.W, buf(lengths), len(self.edges), buf(eu), buf(ev), buf(vm))
            if not h:
                raise ValueError(lib.ge_error().decode())
            self._handle = h
        return self._handle

    def __del__(self):
        if self._handle is not None and _LIB is not None:
            _LIB.ge_graph_free(self._handle)
            self._handle = None

    def native_dist(self, u: int, v: int) -> int | None:
        d = build().ge_graph_dist(self.handle, u, v)
        return None if d < 0 else d

    def py_dist_from(self, u: int) -> dict:
        """Independent check: Dijkstra over the DAG; D(u,v) = start-of-u to start-of-v."""
        if u in self._dist_cache:
            return self._dist_cache[u]
        dist = {u: 0}
        heap = [(0, u)]
        while heap:
            d, x = heapq.heappop(heap)
            if d > dist.get(x, math.inf):
                continue
            nd = d + self.lengths[x]
            for y in self.out[x]:
                if nd < dist.get(y, math.inf):
                    dist[y] = nd
                    heapq.heappush(heap, (nd, y))
        self._dist_cache[u] = dist
        return dist

    def py_dist(self, u: int, v: int) -> int | None:
        return self.py_dist_from(u).get(v)


def load_panel(path: Path, name: str = '') -> Panel:
    return Panel(json.loads(Path(path).read_text()), name or Path(path).parent.name)


def subset_panel(panel: Panel, haps: list[int], name: str = '') -> tuple[Panel, list[int]]:
    """Panel induced by a subset of haplotypes: their vertices and edges, renumbered topologically.

    Returns the new panel and old->new vertex map (-1 for dropped vertices)."""
    keep = sorted({v for h in haps for v in panel.paths[h]})
    new_id = [-1] * panel.V
    for k, v in enumerate(keep):
        new_id[v] = k
    edges = sorted({(new_id[u], new_id[v]) for h in haps for u, v in zip(panel.paths[h], panel.paths[h][1:])})
    graph = dict(lengths=[panel.lengths[v] for v in keep], edges=[list(e) for e in edges],
                 paths=[[new_id[v] for v in panel.paths[h]] for h in haps],
                 haplotypes=[panel.hap_names[h] for h in haps])
    return Panel(graph, name), new_id


# ---------------------------------------------------------------- parameters and instances

def hifi_params(qlen: int, rho: int = HIFI['rho'], bonuses=(0, HIFI['consistency_bonus']), L: int = 17,
                do_exact=True, do_heuristic=True, prune=True) -> Params:
    p = Params()
    p.graph_lookback = max(HIFI['graph_lookback'], int(HIFI['graph_lookback_per_base'] * qlen))
    p.read_lookback = (1 << 62)  # SIZE_MAX in vg: no read limit
    p.indel_limit = max(HIFI['indel'], int(HIFI['indel_per_base'] * qlen))
    p.base_seed_length = L
    p.rho = rho
    p.gap_scale = HIFI['gap_scale']
    p.n_bonus = len(bonuses)
    for k, b in enumerate(bonuses):
        p.bonus[k] = b
    p.do_exact = int(do_exact)
    p.do_heuristic = int(do_heuristic)
    p.prune = int(prune)
    return p


def sort_order(anchors: list, tie: str = 'graph') -> list[int]:
    """Giraffe's view order: read start asc, read end desc (CI 96-106); std::sort leaves ties unspecified.

    tie='graph': ties by (v, x), then input index; tie='reverse': the reverse of that."""
    if tie == 'graph':
        key = lambda k: (anchors[k][3], -anchors[k][4], anchors[k][0], anchors[k][1], k)
    elif tie == 'reverse':
        key = lambda k: (anchors[k][3], -anchors[k][4], -anchors[k][0], -anchors[k][1], -k)
    else:
        raise ValueError(tie)
    return sorted(range(len(anchors)), key=key)


def item_points(a, match=HIFI['match'], item_bonus=HIFI['item_bonus']) -> int:
    span = a[2] - a[1] + 1
    assert span == a[4] - a[3] + 1, 'equal graph and read spans expected (exact matches)'
    return match * span + item_bonus


class Result:
    pass


def solve(panel: Panel, anchors: list, qlen: int, rho: int, bonuses=(0, 12), tie='graph', do_exact=True,
          do_heuristic=True, L: int = 17, item_bonus: int = HIFI['item_bonus'], prune: bool = True) -> Result:
    lib = build()
    order = sort_order(anchors, tie)
    sa = [anchors[k] for k in order]
    n, H, W = len(sa), panel.H, panel.W
    data = array('q', [t for a in sa for t in a[:6]])
    item = array('i', [item_points(a, item_bonus=item_bonus) for a in sa])
    p = hifi_params(qlen, rho, bonuses, L, do_exact, do_heuristic and W == 1, prune)
    nb = len(bonuses) if p.do_heuristic else 0
    hs, hsrc, hrec = (array('i', bytes(4 * nb * n)) for _ in range(3))
    hmask = array('Q', bytes(8 * nb * n))
    K = n * H if do_exact else 0
    C, Pj = array('i', bytes(4 * K)), array('i', bytes(4 * K))
    Pg = array('h', bytes(2 * K))
    B = array('i', bytes(4 * (n if do_exact else 0)))
    stats = array('d', bytes(8 * 10))
    err = lib.ge_solve(panel.handle, n, buf(data), buf(item), ct.byref(p), buf(hs), buf(hsrc), buf(hmask), buf(hrec),
                       buf(C), buf(Pj), buf(Pg), buf(B), buf(stats))
    if err:
        raise ValueError(lib.ge_error().decode())
    r = Result()
    r.order, r.anchors, r.item, r.n, r.H, r.rho, r.bonuses = order, sa, item, n, H, rho, list(bonuses)
    r.hs, r.hsrc, r.hmask, r.hrec, r.nb = hs, hsrc, hmask, hrec, nb
    r.C, r.Pj, r.Pg, r.B = C, Pj, Pg, B
    r.stats = dict(zip(('transitions', 'scanned', 'rej_lookback', 'rej_graph_gap', 'rej_read_lookback', 'rej_indel',
                        'rej_unreachable', 'max_in_degree', 'exact_source_updates'), [int(x) for x in stats]))
    r.params = dict(graph_lookback=p.graph_lookback, indel_limit=p.indel_limit, rho=rho, L=L,
                    gap_scale=p.gap_scale, item_bonus=item_bonus, bonuses=list(bonuses))
    return r


def solve_explicit(n, H, items, masks, trans, rho, bonuses=(0, 12), do_exact=True, do_heuristic=True,
                   prune=True) -> Result:
    """trans: list of (src, dst, jump0) with src < dst. masks: Python ints (H <= 64 for the heuristic)."""
    lib = build()
    W = max(1, (H + 63) // 64)
    by_dst = [[] for _ in range(n)]
    for s, d, j0 in trans:
        assert 0 <= s < d < n
        by_dst[d].append((s, j0))
    off, tsrc, tj0 = [0], [], []
    for d in range(n):
        for s, j0 in sorted(by_dst[d]):
            tsrc.append(s)
            tj0.append(j0)
        off.append(len(tsrc))
    p = Params()
    p.rho, p.n_bonus, p.do_exact, p.do_heuristic = rho, len(bonuses), int(do_exact), int(do_heuristic and W == 1)
    p.prune = int(prune)
    for k, b in enumerate(bonuses):
        p.bonus[k] = b
    nb = len(bonuses) if p.do_heuristic else 0
    item = array('i', items)
    am = array('Q', [w for m in masks for w in words(m, W)])
    hs, hsrc, hrec = (array('i', bytes(4 * nb * n)) for _ in range(3))
    hmask = array('Q', bytes(8 * nb * n))
    K = n * H if do_exact else 0
    C, Pj = array('i', bytes(4 * K)), array('i', bytes(4 * K))
    Pg = array('h', bytes(2 * K))
    B = array('i', bytes(4 * (n if do_exact else 0)))
    err = lib.ge_solve_explicit(n, H, W, buf(item), buf(am), buf(array('q', off)), buf(array('i', tsrc)),
                                buf(array('i', tj0)), ct.byref(p), buf(hs), buf(hsrc), buf(hmask), buf(hrec),
                                buf(C), buf(Pj), buf(Pg), buf(B))
    if err:
        raise ValueError(lib.ge_error().decode())
    r = Result()
    r.n, r.H, r.rho, r.bonuses, r.item = n, H, rho, list(bonuses), item
    r.hs, r.hsrc, r.hmask, r.hrec, r.nb = hs, hsrc, hmask, hrec, nb
    r.C, r.Pj, r.Pg, r.B = C, Pj, Pg, B
    return r


# ---------------------------------------------------------------- tracebacks

def heuristic_best(r: Result, b: int):
    """Giraffe's first traceback start: max (score, source), nowhere = SIZE_MAX (CI 655-664)."""
    n = r.n
    S, SR = r.hs[b * n:(b + 1) * n], r.hsrc[b * n:(b + 1) * n]
    key = lambda i: (S[i], SIZE_MAX if SR[i] < 0 else SR[i])
    return max(range(n), key=key)


def heuristic_chain(r: Result, b: int, end: int | None = None) -> list[int]:
    n = r.n
    SR = r.hsrc[b * n:(b + 1) * n]
    i = heuristic_best(r, b) if end is None else end
    chain = []
    while i >= 0:
        chain.append(i)
        i = SR[i]
    return chain[::-1]


def exact_best(r: Result) -> tuple[int, int]:
    i = max(range(r.n), key=lambda k: (r.B[k], -k))
    H = r.H
    row = r.C[i * H:(i + 1) * H]
    h = max(range(H), key=lambda g: (row[g], -g))
    assert row[h] == r.B[i]
    return i, h


def exact_chain(r: Result, start: tuple[int, int] | None = None) -> list[tuple[int, int]]:
    i, h = exact_best(r) if start is None else start
    H = r.H
    out = []
    while True:
        out.append((i, h))
        j, g = r.Pj[i * H + h], r.Pg[i * H + h]
        if j < 0:
            break
        i, h = j, g
    return out[::-1]


# ---------------------------------------------------------------- independent recomputation (Python)

def score_chain_gap(d: int, L: int) -> int:
    if d == 0:
        return 0
    return int(0.01 * L * d + 0.5 * math.log2(d))


def jump0(indel: int, L: int = 17, gap_scale: float = HIFI['gap_scale']) -> int:
    return int(-score_chain_gap(indel, L) * gap_scale)


def py_transition(panel: Panel, s, a, graph_lookback: int, indel_limit: int, L: int = 17):
    """Independent re-derivation of one transition s -> a; returns jump0 or None if illegal."""
    if not s[4] < a[3]:
        return None
    D = panel.py_dist(s[0], a[0])
    if D is None:
        return None
    Z = D + a[1] - s[1]
    if Z <= 0 or Z > graph_lookback:
        return None
    gg = Z - (s[2] - s[1] + 1)
    if gg < 0:
        return None
    rg = a[3] - 1 - s[4]
    indel = abs(rg - gg)
    if indel > indel_limit:
        return None
    return jump0(indel, L)


def greedy_resets(masks: list[int]) -> tuple[int, list[int]]:
    """Intersect-or-reset scan along a chain; returns (#resets, positions of the anchors that reset)."""
    cur, resets = masks[0], []
    for k, m in enumerate(masks[1:], 1):
        if cur & m:
            cur &= m
        else:
            cur = m
            resets.append(k)
    return len(resets), resets


def min_switch_labels(masks: list[int], H: int) -> tuple[int, list[set]]:
    """Min label changes along a chain and, per position, the labels used by some optimal labeling.

    Forward-backward DP over (position, label); independent of the greedy scan."""
    k = len(masks)
    INF = 1 << 30
    lab = [[h for h in range(H) if m >> h & 1] for m in masks]
    f = [dict.fromkeys(lab[0], 0)]
    for t in range(1, k):
        best_prev = min(f[-1].values())
        f.append({h: min(f[-1].get(h, INF), best_prev + 1) for h in lab[t]})
    g = [None] * k
    g[k - 1] = dict.fromkeys(lab[k - 1], 0)
    for t in range(k - 2, -1, -1):
        best_next = min(g[t + 1].values())
        g[t] = {h: min(g[t + 1].get(h, INF), best_next + 1) for h in lab[t]}
    opt = min(f[-1].values())
    feas = [{h for h in lab[t] if f[t][h] + g[t][h] == opt} for t in range(k)]
    return opt, feas
