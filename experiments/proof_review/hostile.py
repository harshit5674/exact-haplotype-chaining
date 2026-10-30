"""Targeted hostile instances aimed at the proof's weak points (deterministic generators).

Each generator yields (name, (lengths, edges, paths), anchors, gamma). Run with
  python3.13 -B experiments/proof_review/hostile.py --workers 3 [--count 1500]
which checks every instance with harness.check_instance (all solvers, brute force, both
pointer checkers, record identity) and writes results/hostile.json.

H1  tie_storm         equal and unequal allele lengths, equal weights, narrow coordinates:
                      merge, floor, switch and best-state ties at once
H2  zero_weights      all weights 0: every state value ties (gaps are nonnegative)
H3  gamma0_wide       gamma = 0, 5..130 labels on few vertices (word boundaries 63/64/65,
                      127/128/129): same-anchor switches with ties, pointer acyclicity
H4  shared_last       many missing labels share last(v_i, g) = w with different phi_g(w),
                      ties between groups and labels
H5  no_last           haplotypes that start after or end before the anchor vertices, so
                      last2reach is undefined for some missing labels
H6  route_variants    star bubbles: 3..8 alleles of different lengths between two shared
                      vertices; thresholds inside the exact route range (bit-sliced paths)
H7  duplicates        identical anchors repeated 2..6 times
H8  boundaries        x = y = len(v), x = 1, c = d, y_j + 1 = x_i, d_j + 1 = c_i, d_j = c_i
H9  huge_spread       weights spread over 2^13..2^41 so the specification's merges take the
                      cluster-coded branch (never reached with small weights); native domain kept
H10 medium            V 5..12, H 4..16 with mid-graph haplotypes, N 20..45 (deeper Fenwick trees,
                      longer lists); checked against the chain-step DAG (enumeration is too large)
"""
from __future__ import annotations
import argparse
import collections
import json
import multiprocessing as mp
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def _anchor(rng, lengths, verts, cmax, weights):
    v = rng.choice(verts)
    L = lengths[v]
    x = rng.choice([1, L, rng.randint(1, L)])
    y = rng.choice([x, L, rng.randint(x, L)])
    c = rng.randint(1, cmax)
    d = c + rng.choice([0, 0, 1])
    return (v, x, y, c, d, rng.choice(weights))


def _verts(paths):
    return sorted({v for p in paths for v in p})


def tie_storm(seed):
    rng = random.Random(100_000 + seed)
    a1, a2 = rng.choice([(1, 1), (2, 2), (1, 2), (1, 3), (2, 3)])
    lengths = [1, a1, a2, rng.choice([1, 2]), 1]
    edges = [(0, 1), (0, 2), (1, 3), (2, 3), (3, 4), (0, 3)]
    base = [[0, 1, 3, 4], [0, 2, 3, 4], [0, 3, 4], [1, 3], [0, 2, 3]]
    H = rng.randint(2, 5)
    paths = [list(base[rng.randrange(len(base))]) for _ in range(H)]
    verts = _verts(paths)
    w = rng.choice([1, 2, 3])
    A = [_anchor(rng, lengths, verts, 4, (w,)) for _ in range(rng.randint(3, 6))]
    return ('tie_storm', seed), (lengths, edges, paths), A, rng.choice([0, 0, 1, 2])


def zero_weights(seed):
    rng = random.Random(200_000 + seed)
    name, graph, A, gamma = tie_storm(seed)
    A = [(v, x, y, c, d, 0) for (v, x, y, c, d, _) in A]
    return ('zero_weights', seed), graph, A, rng.choice([0, 1])


def gamma0_wide(seed):
    rng = random.Random(300_000 + seed)
    H = rng.choice([5, 8, 63, 64, 65, 127, 128, 129, 130])
    lengths = [1, rng.randint(1, 3), rng.randint(1, 3), 1]
    edges = [(0, 1), (0, 2), (1, 3), (2, 3), (0, 3)]
    fams = [[0, 1, 3], [0, 2, 3], [0, 3], [1, 3], [2, 3], [0, 1], [3]]
    k = rng.randint(1, 4)
    pool = rng.sample(fams, k)
    paths = [list(rng.choice(pool)) for _ in range(H)]
    verts = _verts(paths)
    A = [_anchor(rng, lengths, verts, 5, (0, 1, 2)) for _ in range(rng.randint(2, 6))]
    return ('gamma0_wide', seed), (lengths, edges, paths), A, 0


def shared_last(seed):
    rng = random.Random(400_000 + seed)
    # Vertex 1 = w shared by several missing labels ending there (different offsets by
    # their start), long route 1 -> 2 -> 4 versus shortcut 1 -> 4; target vertex 4.
    lengths = [rng.randint(1, 3), rng.randint(1, 3), rng.randint(1, 4), rng.randint(1, 3), 1]
    edges = [(0, 1), (1, 2), (2, 4), (0, 3), (3, 4), (1, 4)]
    if rng.random() < .5:
        edges.append((0, 4))
    edges = sorted(set(edges))
    paths = [[0, 3, 4]]
    for _ in range(rng.randint(2, 5)):
        paths.append(rng.choice([[0, 1], [1], [0, 1, 2], [1, 2], [0, 1, 2, 4], [0, 3]]))
    rng.shuffle(paths)
    verts = _verts(paths)
    A = [_anchor(rng, lengths, verts, 6, (0, 1, 3, 6)) for _ in range(rng.randint(3, 6))]
    return ('shared_last', seed), (lengths, edges, paths), A, rng.choice([0, 1, 2, 5])


def no_last(seed):
    rng = random.Random(500_000 + seed)
    V = rng.randint(3, 6)
    lengths = [rng.randint(1, 3) for _ in range(V)]
    edges = {(v, v + 1) for v in range(V - 1)}
    for u in range(V):
        for v in range(u + 2, V):
            if rng.random() < .2:
                edges.add((u, v))
    edges = sorted(edges)
    paths = []
    for _ in range(rng.randint(2, 5)):
        s = rng.randrange(V)
        e = rng.randrange(s, V)
        paths.append(list(range(s, e + 1)))  # consecutive vertices: edges (v, v+1) exist
    verts = _verts(paths)
    A = [_anchor(rng, lengths, verts, 6, (0, 1, 2, 5)) for _ in range(rng.randint(3, 7))]
    return ('no_last', seed), (lengths, edges, paths), A, rng.choice([0, 1, 3])


def route_variants(seed):
    rng = random.Random(600_000 + seed)
    k = rng.randint(3, 8)
    # 0 -> allele_1..allele_k -> k+1; alleles have lengths 1..6 (some equal).
    lengths = [rng.randint(1, 2)] + [rng.randint(1, 6) for _ in range(k)] + [rng.randint(1, 2)]
    t = k + 1
    edges = [(0, a) for a in range(1, k + 1)] + [(a, t) for a in range(1, k + 1)]
    H = rng.randint(k, k + 3)
    paths = [[0, 1 + (h % k), t] for h in range(H)]
    rng.shuffle(paths)
    verts = _verts(paths)
    A = []
    for _ in range(rng.randint(3, 6)):
        v = rng.choice([0, t, 0, t] + verts)
        L = lengths[v]
        x = rng.randint(1, L)
        y = rng.randint(x, L)
        c = rng.randint(1, 12)
        A.append((v, x, y, c, c + rng.choice([0, 1]), rng.choice([2, 5, 9, 14])))
    return ('route_variants', seed), (lengths, edges, paths), A, rng.choice([0, 1, 3, 10])


def duplicates(seed):
    rng = random.Random(700_000 + seed)
    name, graph, A, gamma = tie_storm(seed)
    base = A[: rng.randint(1, 3)]
    out = []
    for a in base:
        out += [a] * rng.randint(2, 4)
    rng.shuffle(out)
    return ('duplicates', seed), graph, out[:7], gamma


def boundaries(seed):
    rng = random.Random(800_000 + seed)
    lengths = [rng.randint(2, 4), rng.randint(1, 3), rng.randint(2, 4)]
    edges = [(0, 1), (1, 2), (0, 2)]
    paths = rng.choice([[[0, 1, 2], [0, 2]], [[0, 1, 2], [0, 2], [1, 2]], [[0, 2], [0, 1]], [[0, 1, 2], [2]]])
    A = []
    c = 1
    for _ in range(rng.randint(3, 6)):
        v = rng.choice(_verts(paths))
        L = lengths[v]
        kind = rng.randrange(4)
        if kind == 0:
            x = y = L
        elif kind == 1:
            x = y = 1
        elif kind == 2:
            x, y = 1, L
        else:
            x = rng.randint(1, L)
            y = x
        d = c + rng.choice([0, 1])
        A.append((v, x, y, c, d, rng.choice([1, 2, 4])))
        c = d + rng.choice([0, 1, 1])  # next start equal to or right after this end
    # Same-vertex consecutive pair: y_j + 1 = x_i or y_j = x_i.
    v = rng.choice(_verts(paths))
    if lengths[v] >= 2:
        A.append((v, 1, 1, c, c, 3))
        A.append((v, rng.choice([1, 2]), 2, c + 1, c + 1, 3))
    return ('boundaries', seed), (lengths, edges, paths), A, rng.choice([0, 1, 2])


def huge_spread(seed):
    rng = random.Random(900_000 + seed)
    name, graph, A, gamma = (tie_storm if seed % 2 else route_variants)(seed)
    big = [0, 3, 1 << 13, (1 << 13) + 2, 1 << 20, (1 << 20) + 1, 1 << 41]
    A = [(v, x, y, c, d, rng.choice(big)) for (v, x, y, c, d, _) in A]
    return ('huge_spread', seed), graph, A, rng.choice([0, 1, 5, 1 << 14, 1 << 30])


def medium(seed):
    rng = random.Random(1_000_000 + seed)
    V = rng.randint(5, 12)
    lengths = [rng.choice([1, 1, 2, 3, 4, 8]) for _ in range(V)]
    H = rng.randint(4, 16)
    fams = []
    for _ in range(rng.randint(2, 5)):
        s0 = rng.choice([0, 0, rng.randrange(V)])
        e0 = rng.choice([V - 1, V - 1, rng.randrange(s0, V)])
        keep = rng.uniform(.3, .9)
        fams.append([s0] + [v for v in range(s0 + 1, e0) if rng.random() < keep] + ([e0] if e0 > s0 else []))
    paths = []
    for _ in range(H):
        p = list(rng.choice(fams))
        if rng.random() < .2 and len(p) > 2:  # private indel
            k = rng.randrange(1, len(p) - 1)
            p = p[:k] + p[k + 1:]
        paths.append(p)
    edges = {(u, v) for p in paths for u, v in zip(p, p[1:])}
    for u in range(V):
        for v in range(u + 1, V):
            if rng.random() < .05:
                edges.add((u, v))
    verts = _verts(paths)
    n = rng.randint(20, 45)
    A = [_anchor(rng, lengths, verts, rng.choice([15, 30, 60]), (0, 1, 2, 5, 9, 20, 40)) for _ in range(n)]
    return ('medium', seed), (lengths, sorted(edges), paths), A, rng.choice([0, 1, 3, 10, 100])


GENERATORS = dict(tie_storm=tie_storm, zero_weights=zero_weights, gamma0_wide=gamma0_wide,
                  shared_last=shared_last, no_last=no_last, route_variants=route_variants,
                  duplicates=duplicates, boundaries=boundaries, huge_spread=huge_spread, medium=medium)


# The independent traceback check is O(N^2 H^2) per output; cap the two expensive generators.
CAPS = dict(medium=300, gamma0_wide=400)


def worker(job):
    import harness as H
    gen, start, stop = job
    stats = collections.Counter()
    failures, nfail = [], 0
    for seed in range(start, stop):
        name, (L, E, P), A, gamma = GENERATORS[gen](seed)
        try:
            gcache = H.GraphCache(L, E, P)
            f = H.check_instance(gcache, A, gamma, stats)
        except Exception:
            import traceback
            f = [dict(solver='harness', kind='exception', detail=traceback.format_exc()[-1500:])]
        if f:
            nfail += 1
            if len(failures) < 50:
                failures.append(H.failure_record(L, E, P, A, gamma, f, list(map(str, name))))
    return gen, dict(stats), nfail, failures


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--count', type=int, default=1500)
    ap.add_argument('--workers', type=int, default=3)
    args = ap.parse_args()
    assert args.workers <= 3 and sys.flags.optimize == 0
    jobs = []
    for gen in GENERATORS:
        count = min(args.count, CAPS.get(gen, args.count))
        step = max(1, count // 4)
        for s in range(0, count, step):
            jobs.append((gen, s, min(count, s + step)))
    with mp.get_context('spawn').Pool(args.workers) as pool:
        parts = pool.map(worker, jobs, chunksize=1)
    out = {}
    for gen, st, nf, fl in parts:
        o = out.setdefault(gen, dict(stats=collections.Counter(), failing_instances=0, failures=[]))
        o['stats'].update(st)
        o['failing_instances'] += nf
        o['failures'] += fl
    report = dict(count_per_generator=args.count, caps=CAPS,
                  generators={g: dict(failing_instances=o['failing_instances'],
                                      stats=dict(sorted(o['stats'].items())), failures=o['failures'][:50])
                              for g, o in out.items()})
    (HERE / 'results').mkdir(exist_ok=True)
    (HERE / 'results/hostile.json').write_text(json.dumps(report, indent=1, default=str) + '\n')
    for g, o in report['generators'].items():
        s = o['stats']
        print(g, 'instances', s.get('instances'), 'states', s.get('states'), 'failing', o['failing_instances'],
              'enumerated', s.get('enumerated', 0), 'dag_only', s.get('dag_only', 0),
              'fp_pairs_level4', s.get('fp_pairs_level4', 0), 'ptr_same', s.get('ptr_same', 0),
              'ptr_remote', s.get('ptr_remote', 0))
    sys.exit(1 if any(o['failing_instances'] for o in out.values()) else 0)


if __name__ == '__main__':
    main()
