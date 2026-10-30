"""Phase 3 (and Phase 4 restricted panels): downstream accuracy on simulated recombinant reads.

For each simulated read, on its true-strand instance (anchors exported by pinned Minichain):
  exact_rho<r>       exact chain under Giraffe's objective (basic objective, Giraffe transitions/gaps)
  giraffe_b<b>_rho<r> Giraffe's single-record chain (consistency bonus b), same transitions/gaps
  cgj                exact chain under the CGJ gap-sensitive objective (Problem 5 as transcribed),
                     native class engine (research/current/new/gap_native.py), weights x200, gamma 10000.

Truth: every anchor is classed as truly placed (its graph position is where the read base came from)
or not, and carries the true haplotype of its read base. Evaluation (definitions in REPORT.md):
  - a chain's segments: basic-objective chains use the intersect-or-reset scan along the chain (min
    switches for that chain; the same rule for exact and Giraffe chains); CGJ uses its traced labels;
  - an anchor's assignment is correct iff the true haplotype passes through every anchored vertex of the
    anchor's segment (the segment's label set contains the truth, or a haplotype identical to it on all
    of the segment's anchored vertices when the truth is not in the panel);
  - breakpoints: midpoint of each switch's ambiguity interval (basic objective) or of the gap where the
    traced label changes (CGJ), against the true cuts in read coordinates.

  python3.13 experiments/exactness_payoff/run_phase3.py --prefix experiments/exactness_payoff/sim/hprc60 \
      --panel hprc60 --out experiments/exactness_payoff/phase3_hprc60.jsonl
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from multiprocessing import Pool
import os
from pathlib import Path
import platform
import random
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / 'research/current/new'))
import emul  # noqa: E402
import simulate_reads as sim  # noqa: E402

SCALE, GAMMA = 200, 10000
PANELS = {name: ROOT / 'experiments/hprc_mhc_panel/data' / name for name in ('hprc60', 'hprc463')}


def gap_native_guard():
    """Use the shared native build only if it already exists (never build outside this folder)."""
    src = ROOT / 'research/current/new/gap_native.cpp'
    compiler = os.environ.get('CXX', 'clang++')
    flags = ['-std=c++20', '-O3', '-shared', '-fPIC']
    ident = hashlib.sha256(src.read_bytes() + repr((compiler, flags)).encode()).hexdigest()[:16]
    target = ROOT / 'results/native_build' / f'gap-{ident}.dylib'
    if not target.exists():
        raise SystemExit(f'{target} missing; refusing to build the shared gap_native library here')
    return target.name


# ------------------------------------------------------------ worker state

W = {}


def _init(panel_name, subset, rhos, do_cgj):
    full = emul.load_panel(PANELS[panel_name] / 'graph.json', panel_name)
    hs = sim.Haplotypes(PANELS[panel_name])
    assert hs.paths == full.paths
    if subset:
        chain_panel, vmap = emul.subset_panel(full, subset, f'{panel_name}-subset{len(subset)}')
    else:
        chain_panel, vmap = full, list(range(full.V))
    chain_panel.handle
    W.update(full=full, hs=hs, panel=chain_panel, vmap=vmap, rhos=rhos, subset=subset)
    if do_cgj:
        import gap_chaining as gc
        import gap_native as gn
        p = chain_panel
        W['gc'] = gc
        W['gap'] = gn.GapGraph(gc.Graph(p.lengths, p.edges, p.paths).prepare())


# ------------------------------------------------------------ truth

def true_walk(read):
    """Vertices (>= 17 bp, full-panel ids) of the read's true walk with their forward-read intervals."""
    hs = W['hs']
    rmap = sim.read_to_mosaic(read['error_events'], read['mosaic_length'])
    inv = [-1] * read['mosaic_length']
    for r, m in enumerate(rmap):
        if m >= 0 and inv[m] < 0:
            inv[m] = r
    last = 0
    for m in range(read['mosaic_length']):  # deleted bases take the next read position
        if inv[m] < 0:
            inv[m] = last
        last = inv[m]
    out = []
    for seg in read['segments']:
        h, a = seg['hap'], seg['hap_from']
        L = seg['read_to'] - seg['read_from']
        k0, k1 = hs.locate(h, a)[0], hs.locate(h, a + L - 1)[0]
        for k in range(k0, k1 + 1):
            v = hs.paths[h][k]
            if len(hs.seq[v]) < 17:
                continue
            v0 = max(hs.starts[h][k], a)
            v1 = min(hs.starts[h][k] + len(hs.seq[v]), a + L)
            m0, m1 = seg['read_from'] + v0 - a, seg['read_from'] + v1 - a
            out.append((inv[m0], inv[m1 - 1], v))
    return out


def anchor_truth(read, anchors):
    """Per anchor (full-panel vertex ids): (truly placed, true haplotype or None)."""
    hs = W['hs']
    rmap = sim.read_to_mosaic(read['error_events'], read['mosaic_length'])
    segs = read['segments']
    out = []
    for a in anchors:
        m = rmap[a[3] - 1] if a[3] - 1 < len(rmap) else -1
        if m < 0:
            out.append((False, None))
            continue
        seg = next(s for s in segs if s['read_from'] <= m < s['read_to'])
        q = seg['hap_from'] + (m - seg['read_from'])
        _, v, off = hs.locate(seg['hap'], q)
        out.append((v == a[0] and off + 1 == a[1], seg['hap']))
    return out


# ------------------------------------------------------------ chain evaluation

def basic_segments(masks):
    """Intersect-or-reset scan: [(first, last, intersection)] and backward left boundaries."""
    segs, cur, start = [], masks[0], 0
    for k in range(1, len(masks)):
        if cur & masks[k]:
            cur &= masks[k]
        else:
            segs.append((start, k - 1, cur))
            start, cur = k, masks[k]
    segs.append((start, len(masks) - 1, cur))
    lefts, cur = [], masks[-1]
    for k in range(len(masks) - 2, -1, -1):
        if cur & masks[k]:
            cur &= masks[k]
        else:
            lefts.append(k)
            cur = masks[k]
    lefts.reverse()
    assert len(lefts) == len(segs) - 1
    return segs, lefts


def label_segments(labels, masks):
    segs, start = [], 0
    for k in range(1, len(labels) + 1):
        if k == len(labels) or labels[k] != labels[start]:
            cur = masks[start]
            for t in range(start + 1, k):
                cur &= masks[t]
            assert cur >> labels[start] & 1, 'traced label not on every anchor of its run'
            segs.append((start, k - 1, cur))
            start = k
    return segs


def evaluate(chain, labels, anchors_chain_ids, anchors_full, truth, qlen, walk=None):
    """chain: indices into the instance anchors (read order); labels: None (basic) or chain-panel labels.

    walk_correct: the segment's label set contains a chaining-panel haplotype that traverses every
    vertex (>= 17 bp) of the read's true walk inside the segment's read span (a truth allele missing
    from the panel makes the segment wrong)."""
    panel, full, vmap = W['panel'], W['full'], W['vmap']
    A = [anchors_chain_ids[i] for i in chain]
    masks = [panel.A[a[0]] for a in A]
    full_masks = [full.A[anchors_full[i][0]] for i in chain]
    if labels is None:
        segs, lefts = basic_segments(masks)
        intervals = []
        for k, (first, _, _) in enumerate(segs[1:]):
            lft = lefts[k]
            assert lft < first, 'backward boundary after forward reset'
            intervals.append((A[lft][4], A[first][3] - 1))  # 0-based read positions of the cut's range
    else:
        segs = label_segments(labels, masks)
        intervals = [(A[f - 1][4], A[f][3] - 1) for f, _, _ in segs[1:]]
    bps = [(lo + hi) / 2 for lo, hi in intervals]
    correct = total = walk_correct = 0
    set_sizes = []
    for first, last, inter in segs:
        fm = full_masks[first]
        for t in range(first + 1, last + 1):
            fm &= full_masks[t]
        need = inter
        if walk is not None:
            s0, s1 = A[first][3] - 1, A[last][4] - 1
            for r0, r1, v in walk:
                if r1 >= s0 and r0 <= s1:
                    nv = vmap[v]
                    need &= panel.A[nv] if nv >= 0 else 0
        for t in range(first, last + 1):
            placed, hap = truth[chain[t]]
            if placed:
                total += 1
                correct += fm >> hap & 1
                walk_correct += need != 0
                set_sizes.append(inter.bit_count())
    placed = sum(1 for i in chain if truth[i][0])
    return dict(chain_len=len(chain), placed_frac=placed / len(chain),
                span_cov=(A[-1][4] - A[0][3] + 1) / qlen, switches=len(segs) - 1,
                breakpoints=bps, intervals=intervals, hap_correct=correct, hap_total=total,
                walk_correct=walk_correct if walk is not None else None,
                mean_set_size=sum(set_sizes) / len(set_sizes) if set_sizes else None)


def run_read(item):
    read, inst = item
    panel, vmap = W['panel'], W['vmap']
    anchors_full = inst['anchors']
    keep = [k for k, a in enumerate(anchors_full) if vmap[a[0]] >= 0]
    anchors_full = [anchors_full[k] for k in keep]
    anchors = [[vmap[a[0]], *a[1:]] for a in anchors_full]
    truth = anchor_truth(read, anchors_full)
    walk = true_walk(read)
    qlen = inst['qlen']
    row = dict(name=read['name'], k=read['k'], strand=read['strand'], qlen=qlen, N=len(anchors),
               dropped=len(inst['anchors']) - len(anchors),
               true_haps=[s['hap'] for s in read['segments']], true_breakpoints=read['breakpoints_read'],
               placed_anchors=sum(1 for t in truth if t[0]))
    # Observability of each true cut: no single panel haplotype passes through every truly placed
    # anchored vertex within +-5 kb (chaining panel masks).
    obs = []
    for b in read['breakpoints_read']:
        m = (1 << panel.H) - 1
        for a, t in zip(anchors, truth):
            if t[0] and abs(a[3] - 1 - b) <= 5000:
                m &= panel.A[a[0]]
        obs.append(m == 0)
    row['observable_5kb'] = obs
    pl = [panel.A[a[0]] for a, t in zip(anchors, truth) if t[0]]
    row['truth_min_switches'] = emul.greedy_resets(pl)[0] if pl else 0
    if W['subset'] is not None:
        row['true_haps_in_panel'] = [h in W['subset'] for h in row['true_haps']]
    methods = {}
    if anchors:
        for rho in W['rhos']:
            t = time.perf_counter()
            r = emul.solve(panel, anchors, qlen, rho=rho, bonuses=(0, 12),
                           do_heuristic=panel.H <= 64)
            secs = time.perf_counter() - t
            order = r.order
            opt = max(r.B)
            ech = emul.exact_chain(r)
            ev = evaluate([order[i] for i, _ in ech], None, anchors, anchors_full, truth, qlen, walk)
            ev.update(score=opt, label_switches=sum(1 for (_, g), (_, h) in zip(ech, ech[1:]) if g != h),
                      seconds=round(secs, 3))
            methods[f'exact_rho{rho}'] = ev
            for b, cb in (enumerate((0, 12)) if r.nb else []):
                ch = emul.heuristic_chain(r, b)
                ev = evaluate([order[i] for i in ch], None, anchors, anchors_full, truth, qlen, walk)
                ev.update(score=r.hs[b * r.n + ch[-1]], gap=opt - r.hs[b * r.n + ch[-1]])
                methods[f'giraffe_b{cb}_rho{rho}'] = ev
        if 'gap' in W:
            gc = W['gc']
            ga = [gc.Anchor(*a[:5], a[5] * SCALE) for a in anchors]
            t = time.perf_counter()
            res = W['gap'].solve_packed(ga, GAMMA, mode='class')
            secs = time.perf_counter() - t
            prep = W['gap'].prepared
            owner = [(i, h) for i in range(len(ga)) for h in gc.labels(prep.A[ga[i].v])]
            best = max(res.B)
            st = next(s for s in range(len(res.C)) if res.C[s] == best)
            states = []
            while st >= 0:
                states.append(owner[st])
                st = res.back[st]
            states.reverse()
            # A same-anchor switch repeats the anchor; keep its outgoing label.
            dedup = []
            for i, h in states:
                if dedup and dedup[-1][0] == i:
                    dedup[-1] = (i, h)
                else:
                    dedup.append((i, h))
            chain = [i for i, _ in dedup]
            assert all(anchors[x][3] < anchors[y][3] for x, y in zip(chain, chain[1:]))
            ev = evaluate(chain, [h for _, h in dedup], anchors, anchors_full, truth, qlen, walk)
            ev.update(score=best, label_switches=sum(1 for (_, g), (_, h) in zip(states, states[1:]) if g != h),
                      seconds=round(secs, 3))
            methods['cgj'] = ev
    row['methods'] = methods
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--prefix', required=True)
    ap.add_argument('--panel', choices=sorted(PANELS), required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--rhos', default='2,10,56')
    ap.add_argument('--subset', default='', help='comma-separated haplotype indices: chain on this sub-panel')
    ap.add_argument('--no-cgj', action='store_true')
    ap.add_argument('--workers', type=int, default=3)
    ap.add_argument('--limit', type=int, default=0)
    args = ap.parse_args()
    out = Path(args.out)
    if out.exists():
        raise SystemExit(f'{out} exists')
    lib = None if args.no_cgj else gap_native_guard()
    truth = json.loads(Path(str(args.prefix) + '.truth.json').read_text())
    inst = {}
    plain = Path(str(args.prefix) + '.instances.jsonl')
    text = plain.read_text() if plain.exists() else gzip.decompress(Path(str(plain) + '.gz').read_bytes()).decode()
    for line in text.splitlines():
        c = json.loads(line)
        inst[c['id']] = c
    items = []
    for read in truth['reads']:
        c = inst.get(read['name'] + ':' + read['strand'])
        if c is None:
            c = dict(id=read['name'] + ':' + read['strand'], qlen=read['read_length'], anchors=[])
        items.append((read, c))
    if args.limit:
        items = items[:args.limit]
    items.sort(key=lambda x: -len(x[1]['anchors']))
    subset = [int(x) for x in args.subset.split(',')] if args.subset else None
    rhos = [int(x) for x in args.rhos.split(',')]
    t0 = time.perf_counter()
    with Pool(min(args.workers, 3), initializer=_init, initargs=(args.panel, subset, rhos, not args.no_cgj)) as pool, \
            out.open('x') as f:
        for k, row in enumerate(pool.imap_unordered(run_read, items, chunksize=1)):
            f.write(json.dumps(row) + '\n')
            f.flush()
            if (k + 1) % 25 == 0:
                print(f'{k + 1}/{len(items)} reads, {time.perf_counter() - t0:.0f} s', flush=True)
    meta = dict(prefix=str(args.prefix), panel=args.panel, subset=subset, rhos=rhos, reads=len(items),
                gap_native_library=lib, emul_library=emul.build().identity, scale=SCALE, gamma=GAMMA,
                seconds=round(time.perf_counter() - t0, 1), python=platform.python_version(), argv=sys.argv,
                assertions=__debug__)
    Path(str(out) + '.meta.json').write_text(json.dumps(meta, indent=1) + '\n')
    print(json.dumps(meta))


if __name__ == '__main__':
    main()
