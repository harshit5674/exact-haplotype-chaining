"""Instance families for the adversarial tests (deterministic; every instance is reproducible).

An instance is (lengths, edges, paths, anchors, gamma) with anchors as (v, x, y, c, d, w).
Vertex ids are topological (edges u < v) because the solvers require it; the brute force
does not use that fact.

EXA   graph-exhaustive: every DAG on V <= 3 vertices (all edge subsets u < v), every vertex
      length vector in {1,2}^V, every ordered H-tuple (H <= 3, repetition allowed) of
      directed paths (single vertices included; paths may start and end anywhere).
      For each graph, ANCHOR_SETS_A deterministic pseudo-random anchor sets (N in 2..4,
      small coordinates, ties likely) with gamma cycling through GAMMAS.
EXA4  the same for V = 4 and H <= 2, ANCHOR_SETS_A4 anchor sets per graph.
EXB2  anchor-exhaustive, N = 2, on each hostile graph of HOSTILE_GRAPHS: every ordered pair
      of anchors from the alphabet (v covered; x, y in {1, len(v)} with x <= y; all
      1 <= c <= d <= 3; w in {0, 1, 3}) and every gamma in {0, 1, 2, 1000}.
EXB3  anchor-exhaustive, N = 3, on HOSTILE_GRAPHS_3: every multiset of three anchors from
      the reduced alphabet (v covered; (x, y) in {(1,1), (len,len)}; (c,d) in
      {(1,1),(2,2),(3,3),(1,2)}; w in {0, 2}) and gamma in {0, 1}.
RS    random stress: V in 2..7, H in 1..6, N in 2..9, mid-graph haplotype starts and ends,
      duplicate haplotypes, extra non-haplotype edges, uncovered vertices, duplicate
      anchors, boundary coordinates, zero weights, gamma in {0, 1, 2, 7, 10**6}.
"""
from __future__ import annotations
import itertools
import random

GAMMAS = (0, 1, 2, 10 ** 6)
ANCHOR_SETS_A = 24
ANCHOR_SETS_A4 = 3


def dags(V):
    pairs = [(u, v) for u in range(V) for v in range(u + 1, V)]
    for k in range(1 << len(pairs)):
        yield [pairs[t] for t in range(len(pairs)) if k >> t & 1]


def all_paths(V, edges):
    succ = {u: sorted(v for (a, v) in edges if a == u) for u in range(V)}
    out = []

    def ext(p):
        out.append(tuple(p))
        for v in succ[p[-1]]:
            ext(p + [v])
    for u in range(V):
        ext([u])
    return out


def graphs_exhaustive(V, Hmax, lengths_set=(1, 2)):
    """Yield (gid, lengths, edges, paths) for every graph of the family, in a fixed order."""
    gid = 0
    for edges in dags(V):
        P = all_paths(V, edges)
        for lengths in itertools.product(lengths_set, repeat=V):
            for H in range(1, Hmax + 1):
                for tup in itertools.product(P, repeat=H):
                    yield gid, list(lengths), edges, [list(p) for p in tup]
                    gid += 1


def covered(V, paths):
    return sorted({v for p in paths for v in p})


def random_anchor(rng, lengths, verts, cmax=5, weights=(0, 1, 2, 3, 5)):
    v = rng.choice(verts)
    L = lengths[v]
    r = rng.random()
    if r < .3:
        x = y = rng.choice([1, L])
    elif r < .5:
        x, y = 1, L
    else:
        x = rng.randint(1, L)
        y = rng.randint(x, L)
    c = rng.randint(1, cmax)
    d = c + rng.choice([0, 0, 0, 1, 2])
    return (v, x, y, c, d, rng.choice(weights))


def anchors_for_graph(seed, lengths, paths, nmin=2, nmax=4):
    rng = random.Random(seed)
    verts = covered(len(lengths), paths)
    n = rng.randint(nmin, nmax)
    A = [random_anchor(rng, lengths, verts) for _ in range(n)]
    if rng.random() < .15 and A:  # duplicate an anchor exactly
        A.append(A[rng.randrange(len(A))])
    return A


def exa_instances(V, Hmax, sets_per_graph, shard=0, nshards=1):
    """Graph-exhaustive family; deterministic anchors seeded by (V, gid, k)."""
    for gid, lengths, edges, paths in graphs_exhaustive(V, Hmax):
        if gid % nshards != shard:
            continue
        for k in range(sets_per_graph):
            seed = (V * 1_000_003 + gid) * 101 + k
            A = anchors_for_graph(seed, lengths, paths)
            gamma = GAMMAS[(gid + k) % len(GAMMAS)]
            yield ('EXA%d' % V, (gid, k)), (lengths, edges, paths), A, gamma


# ---------------------------------------------------------------------- hostile graphs
HOSTILE_GRAPHS = {
    # Length-variant bubble (route 0->3 is 1 via vertex 1, 3 via vertex 2); h2 duplicates h0.
    'bubble_dup': ([1, 1, 3, 1], [(0, 1), (0, 2), (1, 3), (2, 3)], [[0, 1, 3], [0, 2, 3], [0, 1, 3]]),
    # SOURCE_AUDIT.md remote-only graph: every edge on a haplotype.
    'audit_remote': ([1] * 5, [(0, 1), (1, 2), (2, 4), (0, 3), (3, 4), (0, 2), (2, 3)],
                     [[0, 1, 2, 4], [0, 3, 4], [0, 2, 3, 4]]),
    # Mid-graph starts and ends: last(0, h1) and last(0, h3) do not exist.
    'mid_paths': ([1, 2, 1, 2], [(0, 1), (1, 2), (2, 3), (0, 2), (1, 3)],
                  [[0, 1, 2, 3], [2, 3], [0, 1], [1, 3]]),
    # Two missing labels share last(4, .) = 1 with different offsets; shortcut 1->4 on no haplotype.
    'shared_last': ([1, 2, 1, 3, 1], [(0, 1), (1, 2), (2, 4), (0, 3), (3, 4), (1, 4)],
                    [[0, 1], [1], [0, 3, 4], [0, 1, 2, 4]]),
    # Three-allele bubble plus a path through two alleles: route lengths 1, 2, 3, 3.
    'triple_bubble': ([1, 1, 2, 3, 1], [(0, 1), (0, 2), (0, 3), (1, 4), (2, 4), (3, 4), (1, 2)],
                      [[0, 1, 4], [0, 2, 4], [0, 3, 4], [0, 1, 2, 4]]),
    # Uncovered vertex 2 carries the short route from 1 to 4: D(1,4) = len(1)+len(2) = 2, not
    # 5 via vertex 3; label 2 = [0,1] has last(4, 2) = 1. Edge 0->4 is on label 1 only.
    'uncovered_mid': ([1, 1, 1, 4, 1], [(0, 1), (1, 2), (2, 4), (1, 3), (3, 4), (0, 4)],
                      [[0, 1, 3, 4], [0, 4], [0, 1]]),
}
HOSTILE_GRAPHS_3 = ('bubble_dup', 'mid_paths', 'shared_last', 'uncovered_mid')


def alphabet_b2(lengths, paths):
    out = []
    for v in covered(len(lengths), paths):
        L = lengths[v]
        xys = sorted({(1, 1), (1, L), (L, L)})
        for (x, y) in xys:
            for c in range(1, 4):
                for d in range(c, 4):
                    for w in (0, 1, 3):
                        out.append((v, x, y, c, d, w))
    return out


def alphabet_b3(lengths, paths):
    out = []
    for v in covered(len(lengths), paths):
        L = lengths[v]
        for (x, y) in sorted({(1, 1), (L, L)}):
            for (c, d) in ((1, 1), (2, 2), (3, 3), (1, 2)):
                for w in (0, 2):
                    out.append((v, x, y, c, d, w))
    return out


def exb2_instances(name, shard=0, nshards=1):
    lengths, edges, paths = HOSTILE_GRAPHS[name]
    alpha = alphabet_b2(lengths, paths)
    k = 0
    for a in alpha:
        for b in alpha:
            for gamma in (0, 1, 2, 1000):
                if k % nshards == shard:
                    yield ('EXB2', name, k), (lengths, edges, paths), [a, b], gamma
                k += 1


def exb3_instances(name, shard=0, nshards=1):
    lengths, edges, paths = HOSTILE_GRAPHS[name]
    alpha = alphabet_b3(lengths, paths)
    k = 0
    for combo in itertools.combinations_with_replacement(range(len(alpha)), 3):
        for gamma in (0, 1):
            if k % nshards == shard:
                yield ('EXB3', name, k), (lengths, edges, paths), [alpha[t] for t in combo], gamma
            k += 1


def random_stress_instance(seed):
    rng = random.Random(9_000_000 + seed)
    V = rng.randint(2, 7)
    lengths = [rng.choice([1, 1, 2, 3, 5]) for _ in range(V)]
    H = rng.randint(1, 6)
    paths = []
    for _ in range(H):
        if paths and rng.random() < .15:
            paths.append(list(rng.choice(paths)))  # duplicate haplotype
            continue
        s = rng.randrange(V)
        e = rng.randrange(s, V)
        keep = rng.uniform(.2, 1)
        p = [s] + [v for v in range(s + 1, e) if rng.random() < keep] + ([e] if e > s else [])
        paths.append(p)
    edges = {(u, v) for p in paths for u, v in zip(p, p[1:])}
    for u in range(V):
        for v in range(u + 1, V):
            if rng.random() < .12:
                edges.add((u, v))
    edges = sorted(edges)
    verts = covered(V, paths)
    n = rng.randint(2, 9)
    cmax = rng.choice([3, 5, 8, 20])
    A = []
    for _ in range(n):
        if A and rng.random() < .12:
            A.append(A[rng.randrange(len(A))])
        else:
            A.append(random_anchor(rng, lengths, verts, cmax=cmax,
                                   weights=(0, 0, 1, 2, 3, 5, 8, 13, 40)))
    gamma = rng.choice([0, 0, 1, 2, 7, 10 ** 6])
    return ('RS', seed), (lengths, edges, paths), A, gamma
