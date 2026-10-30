"""Rebuild Minichain v1.3's chaining inputs exactly, for score-identity checks and matched timing.

Minichain (pinned d21e694) splits seed hits by connected component of its bidirected vertex
set, where vertex = 2*segment + strand. On the MHC graph these are the forward component
and its mirror, the reverse component, whose walks are the reversed haplotype walks. Each
component is chained on its own, in 0-based coordinates, with every anchor weighted by the
first anchor's span times the scale factor.

This module reads the component index dumped by the instrumented build
(vendor/minichain-dp, RECOMB_DP_INDEX). It rebuilds a component as a gap_chaining.Graph
whose vertices are renumbered by Minichain's own topological order. It converts Minichain's
seed hits to 1-based closed anchors. That is a uniform +1 shift: the gap-sensitive
recurrence uses only coordinate differences and strict comparisons, so the shift changes
nothing (see REPORT.md, "Coordinate shift"). Minichain's dependence on absolute
coordinates (its RMQ default of (0,0)) is handled separately in minichain_ref.py.
"""
from __future__ import annotations

import ctypes as ct
import subprocess
import sys
from array import array
from dataclasses import dataclass, field
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'research/current/new'))
import gap_chaining as gc  # noqa: E402

INT_MAX = 2**31 - 1
FNV0 = 1469598103934665603

_FNV = None


def fnv_lib():
    """Compile fnvhash.c once into build/ and return the ctypes handle."""
    global _FNV
    if _FNV is None:
        out = HERE / 'build' / 'libfnvhash.dylib'
        src = HERE / 'fnvhash.c'
        if not out.exists() or out.stat().st_mtime < src.stat().st_mtime:
            out.parent.mkdir(exist_ok=True)
            subprocess.run(['clang', '-O2', '-shared', '-o', str(out), str(src)], check=True)
        lib = ct.CDLL(str(out))
        lib.fnv_int64.argtypes = [ct.c_void_p, ct.c_long, ct.c_uint64]
        lib.fnv_int64.restype = ct.c_uint64
        _FNV = lib
    return _FNV


def fnv(values, h: int = FNV0) -> int:
    """FNV-1a over int64 little-endian bytes, identical to recomb_dp::fnv in the instrumented build."""
    a = values if isinstance(values, array) and values.typecode == 'q' else array('q', values)
    if not len(a):
        return h
    buf = (ct.c_int64 * len(a)).from_buffer(a)
    return fnv_lib().fnv_int64(buf, len(a), h)


@dataclass
class Component:
    cid: int
    V: int
    K: int
    glob: list = field(default_factory=list)      # local vertex -> global vertex (2*seg + strand)
    seg: list = field(default_factory=list)       # local vertex -> segment name
    strand: list = field(default_factory=list)
    length: list = field(default_factory=list)
    inn: list = field(default_factory=list)       # local vertex -> in-neighbours, Minichain's order
    top: list = field(default_factory=list)       # Minichain's topological order (local ids)
    paths: list = field(default_factory=list)     # path k -> local vertices
    names: list = field(default_factory=list)     # path k -> walk sample name
    l2r: list = field(default_factory=list)       # [k][v] last2reach index on path k, or -1
    dist: list = field(default_factory=list)      # [k][v] Minichain's Distance as dumped
    d2b: list = field(default_factory=list)       # [k][v] dist2begin

    def __post_init__(self):
        self._rank = None

    @property
    def rank(self):
        if self._rank is None:
            self._rank = [0] * self.V
            for r, v in enumerate(self.top):
                self._rank[v] = r
        return self._rank


def load_index(path: Path):
    """Parse RECOMB_DP_INDEX output into Component objects (index by cid) and its header."""
    comps, cur, header = [], None, {}
    with open(path) as f:
        for line in f:
            a = line.rstrip('\n').split('\t')
            tag = a[0]
            if tag == 'I':
                header = dict(x.split('=') for x in a[1:])
            elif tag == 'CID':
                cur = Component(int(a[1]), int(a[2]), int(a[3]))
                comps.append(cur)
            elif tag == 'V':
                v = int(a[1])
                assert v == len(cur.glob)
                cur.glob.append(int(a[2]))
                cur.seg.append(a[3])
                cur.strand.append(int(a[4]))
                cur.length.append(int(a[5]))
                ins = a[6]
                assert ins.startswith('in')
                cur.inn.append([int(u) for u in ins[2:].split(',') if u])
            elif tag == 'TOP':
                cur.top = [int(v) for v in a[1:]]
            elif tag == 'P':
                assert int(a[1]) == len(cur.paths)
                cur.names.append(a[2])
                cur.paths.append([int(v) for v in a[3:]])
            elif tag in ('L2R', 'DIST', 'D2B'):
                row = [int(x) for x in a[2:]]
                assert len(row) == cur.V
                {'L2R': cur.l2r, 'DIST': cur.dist, 'D2B': cur.d2b}[tag].append(row)
    for c in comps:
        assert sorted(c.top) == list(range(c.V)), 'topological order must be a permutation'
        assert len(c.paths) == len(c.l2r) == len(c.dist) == len(c.d2b) == c.K
    return comps, header


def minichain_distance(c: Component, truedist: bool = False):
    """Replicate graphUtils.cpp MPC_index's Distance computation, including its in-neighbour
    order dependence (line 792 of the pinned source); truedist=True uses min(D[v], D[u]+len(u))."""
    out = []
    for k in range(c.K):
        D = [INT_MAX // 2] * c.V
        for v in c.paths[k]:
            D[v] = 0
        l2r = c.l2r[k]
        for v in c.top:
            for u in c.inn[v]:
                if l2r[v] == l2r[u] and l2r[v] != -1:
                    if truedist:
                        D[v] = min(D[v], c.length[u] + D[u])
                    else:
                        D[v] = c.length[u] + min(D[v], D[u])
                    assert -INT_MAX - 1 <= D[v] <= INT_MAX, 'int overflow in Minichain distance'
        for v in range(c.V):
            if l2r[v] == -1:
                D[v] = 0
        out.append(D)
    return out


def prepare(c: Component, distance: str = 'minichain', check: bool = True) -> gc.Prepared:
    """gap_chaining.Prepared for component c, vertices renumbered by Minichain's topological order.

    distance='true' keeps gap_chaining's shortest distances (Algorithm 2 as transcribed);
    'minichain' substitutes Minichain's order-dependent Distance (MPC_index, line 792).
    """
    r = c.rank
    edges = sorted({(r[u], r[v]) for v in range(c.V) for u in c.inn[v]})
    lengths = [c.length[v] for v in c.top]
    paths = [[r[v] for v in p] for p in c.paths]
    p = gc.Graph(lengths, edges, paths).prepare()
    if check:
        for k in range(c.K):
            on = set(c.paths[k])
            for v in range(c.V):
                rv = r[v]
                if v in on:
                    assert p.on[rv][k] and p.phi[rv][k] == c.d2b[k][v], (k, v)
                l = c.l2r[k][v]
                want = -1 if l < 0 else r[c.paths[k][l]]
                assert p.last[rv][k] == want, ('last2reach mismatch', k, v, p.last[rv][k], want)
    if distance == 'minichain':
        D = minichain_distance(c)
        for k in range(c.K):
            for v in range(c.V):
                rv = r[v]
                if p.last[rv][k] >= 0 and not p.on[rv][k]:
                    p.distance[rv][k] = D[k][v]
    elif distance != 'true':
        raise ValueError(distance)
    return p


def distance_report(c: Component):
    """Compare Minichain's dumped Distance, its replication, and true shortest distances."""
    mc = minichain_distance(c)
    tr = minichain_distance(c, truedist=True)
    assert mc == c.dist or tr == c.dist, 'dumped Distance matches neither relaxation'
    pairs = diff = 0
    worst = 0
    for k in range(c.K):
        on = set(c.paths[k])
        for v in range(c.V):
            if c.l2r[k][v] >= 0 and v not in on:
                pairs += 1
                if mc[k][v] != tr[k][v]:
                    diff += 1
                    assert mc[k][v] > tr[k][v]
                    worst = max(worst, mc[k][v] - tr[k][v])
    return dict(remote_pairs=pairs, differing_pairs=diff, max_overestimate=worst,
                dumped_is=('minichain' if mc == c.dist else 'true'))


def read_seed_export(path: Path):
    """Seed hits per read from the RECOMB_DUMP_ANCHORS export, in Minichain's array order."""
    reads, cur = {}, None
    order = []
    with open(path) as f:
        for line in f:
            a = line.rstrip('\n').split('\t')
            if a[0] == 'R':
                assert a[1] not in reads, 'duplicate read in export'
                cur = reads[a[1]] = dict(qlen=int(a[2]), n=int(a[3]), hits=[])
                order.append(a[1])
            else:
                cur['hits'].append((a[1], int(a[2]), int(a[3]), int(a[4]), int(a[5])))
    for name in order:
        assert len(reads[name]['hits']) == reads[name]['n'], name
    return reads, order


class Splitter:
    """Split a read's seed hits into Minichain's per-component anchor lists (M[cid])."""
    def __init__(self, comps):
        self.where = {}
        for c in comps:
            for v, (s, st) in enumerate(zip(c.seg, c.strand)):
                self.where[(s, st)] = (c.cid, v)
        self.ncomp = len(comps)

    def split(self, hits):
        out = [[] for _ in range(self.ncomp)]
        for seg, strand, pos, qpos, span in hits:
            cid, lv = self.where[(seg, strand)]
            out[cid].append((lv, pos - span + 1, pos, qpos - span + 1, qpos, span))
        return out


def input_hash(anchors0) -> int:
    """Hash of (local vertex, x, y, c, d) in Minichain's 0-based coordinates, as recomb_dp writes it."""
    return fnv([t for z in anchors0 for t in z[:5]])


def to_gc_anchors(c: Component, anchors0, weight: int):
    r = c.rank
    return [gc.Anchor(r[lv], x + 1, y + 1, cc + 1, d + 1, weight) for lv, x, y, cc, d, _ in anchors0]


def read_dp_dump(path: Path):
    """Parse RECOMB_DP_DUMP: {(read, cid): dict(N, K, wt, recomb, gate, sf, ha, hc, hb, hp, best, nvalid, full)}."""
    out, cur = {}, None
    with open(path) as f:
        for line in f:
            a = line.rstrip('\n').split('\t')
            if a[0] == 'Q':
                key = (a[1], int(a[2]))
                assert key not in out, ('duplicate (read, component) in DP dump', key)
                cur = out[key] = dict(N=int(a[3]), K=int(a[4]), wt=int(a[5]), recomb=int(a[6]),
                                      gate=int(a[7]), sf=int(a[8]), full=None)
            elif a[0] == 'H':
                cur.update(ha=int(a[1], 16), hc=int(a[2], 16), hb=int(a[3], 16), hp=int(a[4], 16),
                           best=int(a[5]), nvalid=int(a[6]))
            elif a[0] == 'A':
                if cur['full'] is None:
                    cur['full'] = []
                row = dict(a=int(a[1]), lv=int(a[2]), x=int(a[3]), y=int(a[4]), c=int(a[5]), d=int(a[6]),
                           Cm=int(a[7]), Cp=int(a[8]),
                           C=[None if s == '.' else tuple(int(t) for t in s.split(':')) for s in a[12:]])
                cur['full'].append(row)
    return out


def read_time_dump(path: Path):
    out = {}
    with open(path) as f:
        for line in f:
            a = line.rstrip('\n').split('\t')
            assert a[0] == 'T'
            assert a[1] not in out, ('duplicate read in time dump', a[1])
            out[a[1]] = dict(total_ns=int(a[2]), dp_ns=int(a[3]), tuples=int(a[4]), n=int(a[5]),
                             per_cid=[int(x) for x in a[6:]])
    return out
