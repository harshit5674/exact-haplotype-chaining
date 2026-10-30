#!/usr/bin/env python3
"""Differential tests of the native gap-sensitive engines in gap_native.cpp.

For every generated input inside the native integer domain it requires:
  * class compact records identical to the specification gap_chaining.solve
    (value, origin, label mask and witness of every record), and identical B, R;
  * expanded class, packed class (cold and warm pair cache) and per-haplotype states (dense,
    label-sparse and hybrid lanes) identical to gap_chaining.direct_oracle (C, B, R) and to
    each other;
  * gap_chaining.verify_pointers on the class and per-haplotype pointers.
Inputs outside the domain (for example the generator's 2^250 weights) must raise
NativeRangeError; they are counted, never rounded. Families: test_gap_chaining's
random_case (wide label counts every 4th seed) and fixed cases, 96-1024 labels around
word boundaries, inputs whose keys exceed the per-haplotype lane bound only, haplotypes that
start and end mid-graph, and the star panel on which the class engine's merge loop is
quadratic.

The auto mode (class engine under a work budget, falling back to per-haplotype) is checked on
every case at the default factor and with factor 0, which forces the fallback whenever the
class work exceeds 64*W words per anchor; on star panels the default must fall back.
"""
import argparse
import hashlib
import json
import random
import time
from pathlib import Path

import gap_chaining as gc
import gap_native as gn
from test_gap_chaining import random_case, fixed_cases


def wide_case(seed):
    rng = random.Random(20_000 + seed)
    H = rng.choice([96, 127, 128, 129, 255, 256, 257, 511, 512, 513, 700, 960, 1023, 1024])
    V = rng.randint(2, 9)
    lengths = [rng.randint(1, 9) for _ in range(V)]
    families = [[0] + [v for v in range(1, V - 1) if rng.random() < rng.uniform(.2, .9)] + [V - 1]
                for _ in range(rng.randint(1, 6))]
    paths = []
    for _ in range(H):
        base = rng.choice(families)
        if rng.random() < .15:
            base = [0] + [v for v in range(1, V - 1) if rng.random() < .5] + [V - 1]
        paths.append(list(base))
    edges = {(u, v) for path in paths for u, v in zip(path, path[1:])}
    edges |= {(u, v) for u in range(V) for v in range(u + 1, V) if rng.random() < .1}
    vertices = sorted({v for path in paths for v in path})
    anchors = []
    for _ in range(rng.randint(0, 20)):
        v = rng.choice(vertices)
        x = rng.randint(1, lengths[v]); y = rng.randint(x, lengths[v])
        c = rng.randint(1, 40); d = c + rng.randint(0, 5)
        anchors.append(gc.Anchor(v, x, y, c, d, rng.choice([0, 1, 7, 200, 3400]) * rng.randint(1, 4)))
    return gc.Graph(lengths, sorted(edges), paths), anchors, rng.choice([0, 1, 5, 100, 10000])


def midgraph_case(seed):
    """random_case with each haplotype cut to a random contiguous stretch (start/end mid-graph)."""
    rng = random.Random(30_000 + seed)
    g, a, gamma = random_case(seed, wide=seed % 4 == 0)
    paths = []
    for path in g.paths:
        i = rng.randrange(len(path)); j = rng.randrange(i, len(path))
        paths.append(path[i:j + 1] if rng.random() < .7 else list(path))
    covered = {v for path in paths for v in path}
    edges = sorted(set(g.edges) | {(u, v) for path in paths for u, v in zip(path, path[1:])})
    return gc.Graph(g.lengths, edges, paths), [z for z in a if z.v in covered], gamma


def star_case(m, k):
    """Star panel (experiments/proof_review/pair_cost.py): label h runs 0 -> h+1 -> sink; one private
    anchor per label with distinct weights, then k anchors on the sink. Contested merges have |L| = |R| = H."""
    sink = m + 1
    lengths = [1] + [1] * m + [4]
    edges = [(0, h + 1) for h in range(m)] + [(h + 1, sink) for h in range(m)]
    paths = [[0, h + 1, sink] for h in range(m)]
    anchors = [gc.Anchor(h + 1, 1, 1, h + 1, h + 1, 10 * m + h) for h in range(m)]
    for t in range(k):
        anchors.append(gc.Anchor(sink, 1 + t % 4, 1 + t % 4, m + 1 + t, m + 1 + t, 1))
    return gc.Graph(lengths, edges, paths), anchors, 10 ** 6


def big_key_case(seed):
    """Weights near 2^45: inside the class engine's domain, outside the lane packing."""
    g, a, gamma = random_case(seed)
    return g, [gc.Anchor(z.v, z.x, z.y, z.c, z.d, z.weight + (1 << 45)) for z in a], gamma


def run_case(name, case, report, expect_fallback=False):
    g, a, gamma = case
    p = g.prepare()
    try:
        G = gn.GapGraph(p)
        b = G.batch(a, gamma)
    except gn.NativeRangeError:
        report['rejected'] += 1
        return
    spec = gc.solve(p, a, gamma)
    comp = G.solve_compact(b)
    assert comp.records() == spec.records, (name, 'records')
    assert list(comp.B) == spec.B and comp.R_values() == spec.R, (name, 'B/R')
    eC, eB, eR, _ = gc.direct_oracle(p, a, gamma)
    outputs = [comp.expand(), G.solve_packed(b, mode='class'), G.solve_packed(b, mode='class', warm=True)]
    try:
        # Default hybrid, always dense lanes, always label-sparse.
        for cut in (gn.DENSE_CUT, 0, 257):
            outputs.append(G.solve_packed(b, mode='per_haplotype', dense_cut=cut))
    except gn.NativeRangeError:
        report['per_haplotype_rejected'] += 1
    auto = G.solve_packed(b, mode='auto')
    forced = G.solve_packed(b, mode='auto', auto_factor=0.0)
    outputs += [auto, forced]
    # Warm-cache reuse after interruption must also preserve scores and pointers.
    outputs.append(G.solve_packed(b, mode='auto', warm=True, auto_factor=0.25))
    for res in outputs:
        if res.stats['mode'] == 'auto' and res.stats['auto_fallback_available']:
            assert res.stats['auto_work'] <= res.stats['auto_budget'], (name, 'budget overshoot')
    report['auto_fallbacks'] += auto.stats['fallback']
    report['forced_fallbacks'] += forced.stats['fallback']
    assert auto.stats['fallback'] or not expect_fallback, (name, 'auto did not fall back')
    digests = set()
    for res in outputs:
        C, ptr, B, R = res.as_dicts(a)
        assert (C, B, R) == (eC, eB, eR), (name, res.stats['mode'])
        counts = gc.verify_pointers(p, a, gamma, C, B, ptr)
        if res.stats['mode'] in ('class_compact', 'per_haplotype'):
            for k, v in counts.items():
                report['pointer_counts'][res.stats['mode']][k] += v
        digests.add(res.arrays_digest())
    assert len(digests) == 1, name
    report['cases'] += 1
    report['states'] += b.K
    report['max_H'] = max(report['max_H'], p.H)
    report['case_hashes'].append([name, digests.pop()])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cases', type=int, default=2000)
    ap.add_argument('--wide', type=int, default=200)
    ap.add_argument('--bigkey', type=int, default=100)
    ap.add_argument('--midgraph', type=int, default=300)
    ap.add_argument('--out', default='verification_native.json')
    args = ap.parse_args()
    start = time.perf_counter()
    gn.build()
    report = dict(cases=0, rejected=0, per_haplotype_rejected=0, states=0, max_H=0,
                  pointer_counts={m: dict(restart=0, local=0, remote=0, same=0)
                                  for m in ('class_compact', 'per_haplotype')},
                  auto_fallbacks=0, forced_fallbacks=0, midgraph_cases=0, star_cases=0, case_hashes=[])
    for name, case in fixed_cases():
        run_case(name, case, report)
    for seed in range(args.cases):
        run_case(f'random_{seed}', random_case(seed, wide=seed % 4 == 0), report)
    for seed in range(args.wide):
        run_case(f'wide_{seed}', wide_case(seed), report)
    for seed in range(args.bigkey):
        run_case(f'bigkey_{seed}', big_key_case(seed), report)
    before = report['cases']
    for seed in range(args.midgraph):
        run_case(f'midgraph_{seed}', midgraph_case(seed), report)
    report['midgraph_cases'] = report['cases'] - before
    for m in (8, 32, 96):
        run_case(f'star_{m}', star_case(m, 8), report, expect_fallback=True)
        report['star_cases'] += 1
    # Out-of-domain inputs must be rejected, not rounded.
    g = gc.Graph([5], [], [[0], [0]])
    G = gn.GapGraph(g.prepare())
    for bad in ([gc.Anchor(0, 1, 1, 1, 1, 1 << 61)], [gc.Anchor(0, 1, 1, 1, 1, 1)]):
        gamma = 1 << 61 if bad[0].weight == 1 else 0
        try:
            G.solve_compact(bad, gamma)
        except gn.NativeRangeError:
            report['rejected_explicit'] = report.get('rejected_explicit', 0) + 1
        else:
            raise AssertionError('out-of-domain input was not rejected')
    assert report['rejected_explicit'] == 2
    assert report['per_haplotype_rejected'] >= args.bigkey // 2, 'lane bound not exercised'
    from test_gap_budget import run_budget_regressions
    report['budget_regressions'] = run_budget_regressions()
    report['seconds'] = time.perf_counter() - start
    report['status'] = 'all passed'
    report['case_digest'] = hashlib.sha256(json.dumps(report['case_hashes']).encode()).hexdigest()
    Path(args.out).write_text(json.dumps(report, indent=1) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'case_hashes'}, indent=1))


if __name__ == '__main__':
    main()
