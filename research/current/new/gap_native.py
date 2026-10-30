"""Native exact gap-sensitive chaining: origin-class engine and per-haplotype baseline.

Kernels: gap_native.cpp (new; native_chaining.* and class_chaining.* are unchanged).
The objective is the gap-sensitive contract of docs/MATHEMATICAL_MODEL.md, as computed
by gap_chaining.py.

GapGraph(prepared) packs a validated gap_chaining.Prepared index once (label masks,
offsets phi, last2reach vertices and distances); the kernel derives offset bounds,
coordinate bitplanes and remote groups. Per read:

  solve_compact(anchors, gamma)  class engine; every state as a per-anchor record
                                 (value N, origin, witness rule, label mask), no
                                 per-state output. CompactGapResult.expand() writes the
                                 packed arrays; records() gives gap_chaining Records.
  solve_packed(anchors, gamma, mode='class' | 'per_haplotype' | 'auto')
                                 packed per-state arrays from one native call. 'auto'
                                 runs the class engine under a work budget (auto_factor
                                 times the per-haplotype engine's lane operations) and
                                 falls back to the per-haplotype engine if it is exceeded.

The class engine decides merges, floors, switches, minima and remote switches with the
certified fast paths of gap_fastpath.py ("full" plan) and is record-identical to the
specification on every tested input. The per-haplotype engine runs query-major Fenwick
trees over global END ranks in node-major lanes: all H lanes (vectorized) for anchors with
many labels, the anchor's own labels otherwise (`dense_cut`); remote queries are grouped
by last2reach vertex.

Packed output: anchor i owns states offsets[i]:offsets[i+1] in increasing label order.
C holds exact integer values. back holds the packed index of the predecessor state, or
-1 for a restart: the same anchor means a same-anchor switch (i, h*), the same label a
local transition, another label g a remote switch with w = last(v_i, g). R uses
NATIVE_NEG for "no remote predecessor". Pointers are legal tie-dependent witnesses.

Arithmetic is int64 with an explicit bound: with S the total graph length,
  bound = sum(w) + 2S + max x + max y + max c + max d + max w + gamma + 16 < 2^60,
which keeps every intermediate within +-2^62; the per-haplotype lanes also need
published keys sum(w) + S + max y + max d < 2^38 and fewer than 2^24 anchors.
Inputs outside these domains raise NativeRangeError; nothing is rounded.
"""
from __future__ import annotations
from array import array
from itertools import accumulate
import ctypes as ct
import hashlib
import os
from pathlib import Path
import subprocess
import sys
import time

import gap_chaining as gc

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
NATIVE_NEG = -(1 << 61)
LIMIT = 1 << 60
STAT_NAMES = ('kernel_seconds', 'setup_seconds', 'init_seconds', 'sweep_seconds', 'export_seconds', 'K', 'N',
              'reserved7', 'merges', 'contested_merges', 'pairs', 'pairs_same_vertex', 'pairs_loose',
              'pairs_range', 'pairs_sliced', 'sliced_thresholds', 'sliced_minima', 'pair_ranges',
              'range_labels', 'route_planes', 'cut_records', 'prefix_records', 'remote_groups',
              'remote_candidates', 'remote_sliced', 'published_records', 'output_records', 'max_node_records',
              'ph_label_steps', 'ph_remote_queries', 'ph_sparse_anchors')
KINDS = ('restart', 'local', 'remote', 'same')
AUTO_FACTOR = 1.0  # auto mode: class work budget in units of the per-haplotype engine's lane operations
DENSE_CUT = 128  # per-haplotype: dense lanes iff an anchor has >= H/2 labels (experiments/hprc_mhc_panel/calibration_dense.json)
_LIB = None


class NativeRangeError(ValueError):
    """Input outside the explicitly checked native integer or index domain."""


def build():
    global _LIB
    if _LIB is not None:
        return _LIB
    source = HERE / 'gap_native.cpp'
    compiler = os.environ.get('CXX', 'clang++')
    flags = ['-std=c++20', '-O3', '-shared', '-fPIC']
    identity = source.read_bytes() + repr((compiler, flags)).encode()
    target = ROOT / 'results/native_build' / ('gap-' + hashlib.sha256(identity).hexdigest()[:16] + ('.dylib' if sys.platform == 'darwin' else '.so'))
    target.parent.mkdir(parents=True, exist_ok=True)
    if not target.exists():
        temp = target.with_name(target.name + '.' + str(os.getpid()) + '.tmp')
        subprocess.run([compiler, *flags, str(source), '-o', str(temp)], check=True)
        temp.replace(target)
    lib = ct.CDLL(str(target))
    ip, qp, up, dp, vp = (ct.POINTER(ct.c_int), ct.POINTER(ct.c_int64), ct.POINTER(ct.c_uint64),
                          ct.POINTER(ct.c_double), ct.c_void_p)
    lib.gap_graph_new.argtypes = [ct.c_int, ct.c_int, ct.c_int, qp, up, qp, ip, qp, dp]
    lib.gap_graph_new.restype = vp
    lib.gap_graph_free.argtypes = [vp]
    lib.gap_graph_free.restype = None
    lib.gap_class_compact.argtypes = [vp, ct.c_int, qp, ct.c_int64, qp, qp, dp, ct.c_int]
    lib.gap_class_compact.restype = vp
    lib.gap_class_packed.argtypes = [vp, ct.c_int, qp, ct.c_int64, ip, qp, qp, qp, qp, dp, ct.c_int]
    lib.gap_warm_info.argtypes = [vp, dp]
    lib.gap_warm_info.restype = None
    lib.gap_warm_reset.argtypes = [vp]
    lib.gap_warm_reset.restype = None
    lib.gap_class_packed.restype = ct.c_int
    lib.gap_ph_solve.argtypes = [vp, ct.c_int, qp, ct.c_int64, ip, qp, qp, qp, qp, dp, ct.c_int]
    lib.gap_ph_solve.restype = ct.c_int
    lib.gap_auto_packed.argtypes = [vp, ct.c_int, qp, ct.c_int64, ip, qp, qp, qp, qp, dp, ct.c_int, ct.c_double, ct.c_int]
    lib.gap_auto_packed.restype = ct.c_int
    lib.gap_compact_size.argtypes = [vp]
    lib.gap_compact_size.restype = ct.c_longlong
    lib.gap_compact_copy.argtypes = [vp, ip, qp, ip, ip, ip, ip, ip, up]
    lib.gap_compact_copy.restype = None
    lib.gap_compact_free.argtypes = [vp]
    lib.gap_compact_free.restype = None
    lib.gap_expand.argtypes = [vp, ct.c_int, ip, ip, ip, qp, ip, ip, ip, ip, up, qp, qp]
    lib.gap_expand.restype = ct.c_int
    lib.gap_error.restype = ct.c_char_p
    _LIB = lib
    return lib


def _buf(a):
    typ = {'i': ct.c_int, 'q': ct.c_int64, 'Q': ct.c_uint64}[a.typecode]
    assert a.itemsize == ct.sizeof(typ)
    return (typ * len(a)).from_buffer(a) if len(a) else None


def _words(mask: int, W: int):
    return [(mask >> (64 * q)) & 0xFFFFFFFFFFFFFFFF for q in range(W)]


class AnchorBatch:
    """Packed, range-checked anchors of one read; reusable across engines and repetitions."""
    def __init__(self, graph: 'GapGraph', anchors, gamma: int):
        if type(gamma) is not int or gamma < 0:
            raise ValueError('gamma must be a nonnegative integer')
        n = len(anchors)
        flat = [t for z in anchors for t in (z.v, z.x, z.y, z.c, z.d, z.weight)]
        if any(type(t) is not int for t in flat):
            raise ValueError('integer anchors required')
        try:
            data = array('q', flat)
        except OverflowError as e:
            raise NativeRangeError('anchor coordinates and weights must fit in int64') from e
        if n:
            vs, xs, ys, cs, ds, ws = (data[k::6] for k in range(6))
            if min(vs) < 0 or max(vs) >= graph.V:
                raise ValueError('anchor vertex out of bounds')
            bound = sum(ws) + 2 * graph.S + max(xs) + max(ys) + max(cs) + max(ds) + max(ws) + gamma + 16
            self.max_key = sum(ws) + graph.S + max(ys) + max(ds)
        else:
            vs, bound, self.max_key = array('q'), gamma + 16, 0
        if bound >= LIMIT or n >= 1 << 27:
            raise NativeRangeError('input exceeds the checked native integer domain')
        self.gamma, self.n, self.data = gamma, n, data
        self.vertex = array('i', vs)
        degree = graph.degree
        self.offsets = array('i', accumulate((degree[v] for v in vs), initial=0))
        if self.offsets[-1] >= 1 << 31:
            raise NativeRangeError('too many states for 32-bit offsets')
        self.K = self.offsets[-1]


class PackedGapResult:
    def __init__(self, graph, batch, C, back, B, R, stats):
        self.graph, self.batch, self.C, self.back, self.B, self.R, self.stats = graph, batch, C, back, B, R, stats

    def arrays_digest(self) -> str:
        """SHA-256 of the exact C, B, R arrays (pointers excluded: they may differ on ties)."""
        h = hashlib.sha256()
        for a in (self.C, self.B, self.R):
            h.update(len(a).to_bytes(8, 'little'))
            h.update(a.tobytes())
        return h.hexdigest()

    def as_dicts(self, anchors):
        """C, pointers, B, R in gap_chaining's dictionary form (for oracles/verify_pointers)."""
        p, n = self.graph.prepared, self.batch.n
        owner = [(i, h) for i in range(n) for h in gc.labels(p.A[anchors[i].v])]
        C, ptrs = [], []
        off = self.batch.offsets
        for i in range(n):
            v = anchors[i].v
            row, prow = {}, {}
            for st in range(off[i], off[i + 1]):
                h = owner[st][1]
                row[h] = self.C[st]
                b = self.back[st]
                if b < 0:
                    prow[h] = gc.Witness('restart')
                else:
                    j, g = owner[b]
                    if j == i:
                        prow[h] = gc.Witness('same', i, g)
                    elif g == h:
                        prow[h] = gc.Witness('local', j)
                    else:
                        prow[h] = gc.Witness('remote', j, g, p.last[v][g])
            C.append(row)
            ptrs.append(prow)
        R = [None if r == NATIVE_NEG else r for r in self.R]
        return C, ptrs, list(self.B), R


class CompactGapResult:
    """Every state as a per-anchor record: C(i,h) = value - (phi_h(v_i) - phi_h(v_origin))."""
    def __init__(self, graph, batch, off, value, origin, kind, wj, wg, ww, masks, B, R, stats):
        self.graph, self.batch = graph, batch
        self.off, self.value, self.origin, self.kind = off, value, origin, kind
        self.wj, self.wg, self.ww, self.masks = wj, wg, ww, masks
        self.B, self.R, self.stats = B, R, stats

    def expand(self) -> PackedGapResult:
        t = time.perf_counter()
        b, g = self.batch, self.graph
        C, back = array('q', bytes(8 * b.K)), array('q', bytes(8 * b.K))
        if b.n and b.K:
            err = build().gap_expand(g.handle, b.n, _buf(b.vertex), _buf(b.offsets), _buf(self.off),
                                     _buf(self.value), _buf(self.origin), _buf(self.kind), _buf(self.wj),
                                     _buf(self.wg), _buf(self.masks), _buf(C), _buf(back))
            if err:
                raise ValueError(build().gap_error().decode())
        stats = dict(self.stats, expansion_seconds=time.perf_counter() - t)
        return PackedGapResult(g, b, C, back, self.B, self.R, stats)

    def records(self):
        """Per-anchor gap_chaining.Record lists, in gap_chaining.coalesce order."""
        W, out = self.graph.W, []
        for i in range(self.batch.n):
            row = []
            for z in range(self.off[i], self.off[i + 1]):
                mask = sum(self.masks[z * W + q] << (64 * q) for q in range(W))
                k = KINDS[self.kind[z]]
                w = (gc.Witness('restart') if k == 'restart' else gc.Witness('local', self.wj[z]) if k == 'local'
                     else gc.Witness('remote', self.wj[z], self.wg[z], self.ww[z]) if k == 'remote'
                     else gc.Witness('same', self.wj[z], self.wg[z]))
                row.append(gc.Record(self.value[z], self.origin[z], mask, w))
            row.sort(key=lambda r: (r.coef, r.origin, r.witness.kind, r.witness.j, r.witness.g, r.witness.w))
            out.append(row)
        return out

    def R_values(self):
        return [None if r == NATIVE_NEG else r for r in self.R]


class GapGraph:
    """Packed reusable graph index for both native engines (one panel, any reads)."""
    def __init__(self, prepared: 'gc.Prepared'):
        t = time.perf_counter()
        p = self.prepared = prepared
        V, H = p.V, p.H
        W = 1
        while 64 * W < H:
            W *= 2
        if W > 16:
            raise NativeRangeError('native engines support at most 1024 haplotypes')
        self.V, self.H, self.W, self.S = V, H, W, p.S
        if p.S >= 1 << 40:
            raise NativeRangeError('total graph length exceeds 2^40')
        self.degree = [p.A[v].bit_count() for v in range(V)]
        lengths = array('q', p.graph.lengths)
        vmask = array('Q', [w for v in range(V) for w in _words(p.A[v], W)])
        phi = array('q', [x for row in p.phi for x in row])
        last = array('i', [x for row in p.last for x in row])
        dist = array('q', [0 if x is None else x for row in p.distance for x in row])
        stats = (ct.c_double * 5)()
        lib = build()
        handle = lib.gap_graph_new(V, H, W, _buf(lengths), _buf(vmask), _buf(phi), _buf(last), _buf(dist), stats)
        if not handle:
            raise ValueError(lib.gap_error().decode())
        self.handle = handle
        self.native_index_seconds, self.b = stats[0], int(stats[2])
        self.remote_groups, self.remote_entries = int(stats[3]), int(stats[4])
        self.index_seconds = time.perf_counter() - t

    def __del__(self):
        handle = getattr(self, 'handle', None)
        if handle and _LIB is not None:
            _LIB.gap_graph_free(handle)
            self.handle = None

    def batch(self, anchors, gamma: int) -> AnchorBatch:
        return AnchorBatch(self, anchors, gamma)

    def _stats(self, raw, mode, packing, started):
        meta = dict(zip(STAT_NAMES, raw))
        meta.update(mode=mode, W=self.W, H=self.H, packing_seconds=packing,
                    packed_seconds=time.perf_counter() - started)
        return meta

    def warm_cache_info(self) -> dict:
        out = (ct.c_double * 2)()
        build().gap_warm_info(self.handle, out)
        return dict(pairs=int(out[0]), plane_words=int(out[1]))

    def warm_cache_reset(self):
        build().gap_warm_reset(self.handle)

    def solve_compact(self, anchors, gamma=None, warm=False) -> CompactGapResult:
        """warm=True keeps exact pair ranges in a graph-persistent cache across reads
        (lazily built reusable graph data); the default rebuilds them for every read."""
        lib = build()
        t = time.perf_counter()
        b = anchors if isinstance(anchors, AnchorBatch) else AnchorBatch(self, anchors, gamma)
        n = b.n
        B, R = array('q', bytes(8 * n)), array('q', bytes(8 * n))
        raw = (ct.c_double * 40)()
        packing = time.perf_counter() - t
        handle = lib.gap_class_compact(self.handle, n, _buf(b.data), b.gamma, _buf(B), _buf(R), raw, int(warm))
        if not handle:
            raise ValueError(lib.gap_error().decode())
        try:
            total = lib.gap_compact_size(handle)
            off = array('i', bytes(4 * (n + 1)))
            value = array('q', bytes(8 * total))
            origin, kind, wj, wg, ww = (array('i', bytes(4 * total)) for _ in range(5))
            masks = array('Q', bytes(8 * total * self.W))
            lib.gap_compact_copy(handle, _buf(off), _buf(value), _buf(origin), _buf(kind), _buf(wj), _buf(wg),
                                 _buf(ww), _buf(masks))
        finally:
            lib.gap_compact_free(handle)
        meta = self._stats(raw, 'class_compact_warm' if warm else 'class_compact', packing, t)
        meta['compact_records'] = total
        return CompactGapResult(self, b, off, value, origin, kind, wj, wg, ww, masks, B, R, meta)

    def solve_packed(self, anchors, gamma=None, mode='class', warm=False, dense_cut=DENSE_CUT,
                     auto_factor=AUTO_FACTOR) -> PackedGapResult:
        """dense_cut (per_haplotype and auto, units of 1/256): an anchor with d labels uses dense
        all-H lanes iff 256*d >= dense_cut*H, else its own labels; 0 means always dense.

        mode='auto' runs the class engine under a work budget of auto_factor times the
        per-haplotype engine's lane operations on this read (plus 64*W words per anchor) and
        re-solves the read with the per-haplotype engine if the budget is exceeded; stats
        report 'fallback', 'auto_work', 'auto_estimate' and 'auto_budget'. Work is charged
        before execution, including cache probes, rebuilding and output traversal. These
        are accounting units that bound word work up to a constant, not elapsed time.
        The guarantee applies only when 'auto_fallback_available' is true (keys < 2^38,
        n < 2^24); otherwise the budget is infinite and the class engine runs unguarded.
        auto_factor must be finite and nonnegative; a finite budget must be <= 2^53 so
        the counter represents every admitted integer unit exactly."""
        if mode not in ('class', 'per_haplotype', 'auto'):
            raise ValueError('mode must be class, per_haplotype or auto')
        if warm and mode == 'per_haplotype':
            raise ValueError('warm applies to the class engine only')
        lib = build()
        t = time.perf_counter()
        b = anchors if isinstance(anchors, AnchorBatch) else AnchorBatch(self, anchors, gamma)
        if mode == 'per_haplotype' and (b.max_key >= 1 << 38 or b.n >= 1 << 24):
            raise NativeRangeError('per-haplotype lanes need keys below 2^38 and fewer than 2^24 anchors')
        n, K = b.n, b.K
        C, back = array('q', bytes(8 * K)), array('q', bytes(8 * K))
        B, R = array('q', bytes(8 * n)), array('q', bytes(8 * n))
        raw = (ct.c_double * 40)()
        packing = time.perf_counter() - t
        args = (self.handle, n, _buf(b.data), b.gamma, _buf(b.offsets), _buf(C), _buf(back), _buf(B), _buf(R), raw)
        if mode == 'class':
            err = lib.gap_class_packed(*args, int(warm))
        elif mode == 'per_haplotype':
            err = lib.gap_ph_solve(*args, int(dense_cut))
        else:
            err = lib.gap_auto_packed(*args, int(warm), float(auto_factor), int(dense_cut))
        if err:
            raise ValueError(lib.gap_error().decode())
        name = 'class_expanded' + ('_warm' if warm else '') if mode == 'class' else mode
        meta = self._stats(raw, name, packing, t)
        if mode == 'auto':
            meta.update(fallback=bool(raw[7]), auto_work=raw[36], auto_estimate=raw[37], auto_budget=raw[38],
                        auto_fallback_available=bool(raw[39]), auto_factor=float(auto_factor),
                        auto_cache_probes=raw[31], auto_cache_rehashes=raw[32],
                        auto_cache_rehash_aborts=raw[33], auto_denied_work=raw[34])
        return PackedGapResult(self, b, C, back, B, R, meta)
