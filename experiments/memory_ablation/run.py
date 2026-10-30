"""Native memory per engine and record counts under alternative groupings (2026-10-20).

Earlier drafts of the paper reported no
memory, and no ablation separating origin normalization from simpler sharing. One pass over every panel
of the 100-read benchmark (experiments/benchmark_100reads, same inputs, weights x200, gamma 10000):

Memory. The engines run from a separately built copy of research/current/new/gap_native.cpp linked with
alloc_counter.cpp, which counts this library's own heap (every operator new/delete, by malloc_size).
The production library and adapter are unchanged; the adapter's ctypes declarations are reused. For each
instance and engine (per_haplotype, class expanded, guarded auto, compact) we record the peak native heap
above the level before the call, and the bytes of the arrays the adapter returns (Python-allocated, so
not in the native peak). The native index of each panel is measured once. These are heap bytes, not RSS.

Ablation. From the exact all-state output (the per-haplotype engine's C and pointers) and the compact
records, per anchor i with labels A_i:
  states        |A_i|                                   (per-haplotype lanes)
  keys          distinct key_h(i) = C(i,h)+y_i+d_i+phi_h(v_i) over h in A_i
                (what one publication contributes to a range-maximum structure that shares equal keys)
  scores        distinct (C(i,h), witness rule) over h  (state records when sharing equal current scores,
                each normalized at its own anchor; witness rules as in the compact records)
  score_values  distinct C(i,h)                         (the same without witnesses)
  records       origin-normalized compact records       (ours; witness-bearing)
  published     distinct (value, origin) among them     (records entering the tree after coalescing)
Every engine's output is checked against the per-haplotype engine (C, B, R digests), and the compact
output's expansion likewise. Nothing is timed.

  python3.13 -B experiments/memory_ablation/run.py --out experiments/memory_ablation/results/run1
"""
from __future__ import annotations

import argparse
import ctypes as ct
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(ROOT / 'research/current/new'), str(ROOT / 'experiments/gap_fallback')]
import auto_real as data  # noqa: E402
import gap_chaining as gc  # noqa: E402
import gap_native as gn  # noqa: E402

MODES = ('per_haplotype', 'class', 'auto', 'compact')
SOURCES = [Path(__file__), HERE / 'alloc_counter.cpp', ROOT / 'experiments/gap_fallback/auto_real.py'] + [
    ROOT / 'research/current/new' / f for f in ('gap_chaining.py', 'gap_native.py', 'gap_native.cpp')]


def file_hash(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def instrumented_library() -> Path:
    source, counter = ROOT / 'research/current/new/gap_native.cpp', HERE / 'alloc_counter.cpp'
    compiler = os.environ.get('CXX', 'clang++')
    flags = ['-std=c++20', '-O3', '-shared', '-fPIC']   # the adapter's own flags
    ident = source.read_bytes() + counter.read_bytes() + repr((compiler, flags)).encode()
    target = HERE / 'build' / ('gap-mem-' + hashlib.sha256(ident).hexdigest()[:16] + '.dylib')
    target.parent.mkdir(exist_ok=True)
    if not target.exists():
        temp = target.with_name(f'{target.name}.{os.getpid()}.tmp')
        subprocess.run([compiler, *flags, str(source), str(counter), '-o', str(temp)], check=True)
        temp.replace(target)
    return target


def install():
    """Load the counting library in place of the production one, reusing gap_native.build()'s declarations."""
    assert gn._LIB is None
    path = instrumented_library()
    real = ct.CDLL
    ct.CDLL = lambda _ignored, *a, **k: real(str(path), *a, **k)
    try:
        lib = gn.build()
    finally:
        ct.CDLL = real
    assert lib._name == str(path)
    for f in ('mem_current', 'mem_peak', 'mem_calls'):
        getattr(lib, f).restype, getattr(lib, f).argtypes = ct.c_longlong, []
    lib.mem_reset_peak.restype, lib.mem_reset_peak.argtypes = None, []
    return lib, path


def solve(G, batch, mode):
    return G.solve_compact(batch, warm=False) if mode == 'compact' else G.solve_packed(batch, mode=mode, warm=False)


def output_bytes(z) -> int:
    arrays = ((z.off, z.value, z.origin, z.kind, z.wj, z.wg, z.ww, z.masks, z.B, z.R)
              if isinstance(z, gn.CompactGapResult) else (z.C, z.back, z.B, z.R))
    return sum(len(a) * a.itemsize for a in arrays)


def ablation_counts(p, anchors, batch, ph, comp) -> dict:
    """Per-instance totals over anchors (see the module docstring)."""
    off, C, back = batch.offsets, ph.C, ph.back
    n, K = batch.n, batch.K
    labels = [list(gc.labels(p.A[a.v])) for a in anchors]   # gc.labels is a one-shot iterator
    state_anchor, state_label = [0] * K, [0] * K
    for i in range(n):
        for st, h in zip(range(off[i], off[i + 1]), labels[i]):
            state_anchor[st], state_label[st] = i, h
    tot = dict(states=K, keys=0, scores=0, score_values=0, records=0, published=0)
    for i in range(n):
        assert len(labels[i]) == off[i + 1] - off[i]
        phi_row = p.phi[anchors[i].v]
        keys, scores, values = set(), set(), set()
        for st, h in zip(range(off[i], off[i + 1]), labels[i]):
            c, b = C[st], back[st]
            keys.add(c + phi_row[h])
            values.add(c)
            if b < 0:
                rule = 0
            else:
                j, g = state_anchor[b], state_label[b]
                rule = (1, g) if j == i else (2, j) if g == h else (3, j, g)
            scores.add((c, rule))
        tot['keys'] += len(keys)
        tot['scores'] += len(scores)
        tot['score_values'] += len(values)
        lo, hi = comp.off[i], comp.off[i + 1]
        tot['records'] += hi - lo
        tot['published'] += len({(comp.value[z], comp.origin[z]) for z in range(lo, hi)})
    return tot


def run_panel(lib, locus, name, directory, out: Path) -> dict:
    started = time.perf_counter()
    graph, cases = data.load(directory)
    assert len(cases) == 200
    input_hashes = {f: file_hash(directory / f) for f in ('graph.json', 'instances.jsonl')}
    p = gc.Graph(graph['lengths'], [tuple(e) for e in graph['edges']], graph['paths']).prepare()
    before = lib.mem_current()
    lib.mem_reset_peak()
    G = gn.GapGraph(p)
    index = dict(bytes=lib.mem_current() - before, peak=lib.mem_peak() - before)
    panel = dict(locus=locus, panel=name, H=p.H, V=p.V, W=G.W, b=G.b, input_hashes=input_hashes,
                 index_native_bytes=index['bytes'], index_native_peak=index['peak'],
                 load_before=os.getloadavg(), status='running')
    rows_path = out / f'{locus}-{name}.jsonl'
    with rows_path.open('x') as fh:
        for k, c in enumerate(cases):
            A = [gc.Anchor(*x[:5], x[5] * data.SCALE) for x in c['anchors']]
            batch = G.batch(A, data.GAMMA)
            row = dict(id=c['id'], read=c['id'].rsplit(':', 1)[0], N=batch.n, K=batch.K,
                       input_bytes=sum(len(x) * x.itemsize for x in (batch.data, batch.offsets, batch.vertex)))
            res = {}
            for mode in MODES:
                base = lib.mem_current()
                lib.mem_reset_peak()
                z = solve(G, batch, mode)
                row[mode] = dict(native_peak=lib.mem_peak() - base, native_retained=lib.mem_current() - base,
                                 output_bytes=output_bytes(z))
                if mode == 'auto':
                    row[mode]['fallback'] = bool(z.stats['fallback'])
                    assert z.stats['auto_fallback_available'] and z.stats['auto_work'] <= z.stats['auto_budget']
                if mode == 'compact':
                    row[mode]['records'] = z.stats['compact_records']
                res[mode] = z
            ref = res['per_haplotype'].arrays_digest()
            assert res['class'].arrays_digest() == ref and res['auto'].arrays_digest() == ref, c['id']
            assert res['compact'].expand().arrays_digest() == ref, c['id']
            assert all(row[m]['native_retained'] == 0 for m in MODES), (c['id'], [row[m]['native_retained'] for m in MODES])
            row['counts'] = ablation_counts(p, A, batch, res['per_haplotype'], res['compact'])
            assert row['counts']['records'] == row['compact']['records']
            fh.write(json.dumps(row) + '\n')
            fh.flush()
            del res
            if (k + 1) % 50 == 0:
                print(f'{locus}-{name}: {k + 1}/200', flush=True)
    panel.update(status='passed', seconds=time.perf_counter() - started, load_after=os.getloadavg(),
                 graph_retained_after=lib.mem_current() - before)
    (out / f'{locus}-{name}.json').write_text(json.dumps(panel, indent=1) + '\n')
    print(json.dumps(panel), flush=True)
    return panel


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--panels', default='', help='comma-separated locus-panel keys (default: all 14)')
    args = ap.parse_args()
    args.out.mkdir(parents=True)        # refuse to reuse a directory
    lib, path = install()
    want = set(filter(None, args.panels.split(',')))
    manifest = dict(status='running', command=sys.argv, started_utc=datetime.now(timezone.utc).isoformat(),
                    platform=platform.platform(), python=sys.version, assertions=__debug__,
                    source_hashes={str(q.relative_to(ROOT)): file_hash(q) for q in SOURCES},
                    instrumented_library=str(path.relative_to(ROOT)), gamma=data.GAMMA, weight_scale=data.SCALE,
                    auto_factor=gn.AUTO_FACTOR, pair_cache='cold', panels=[])
    try:
        for locus, name, directory in data.panels():
            if want and f'{locus}-{name}' not in want:
                continue
            manifest['panels'].append(run_panel(lib, locus, name, directory, args.out))
            (args.out / 'manifest.json').write_text(json.dumps(manifest, indent=1) + '\n')
        manifest['status'] = 'passed'
    except BaseException as e:
        manifest.update(status='failed', error=repr(e))
        raise
    finally:
        manifest['finished_utc'] = datetime.now(timezone.utc).isoformat()
        (args.out / 'manifest.json').write_text(json.dumps(manifest, indent=1) + '\n')


if __name__ == '__main__':
    main()
