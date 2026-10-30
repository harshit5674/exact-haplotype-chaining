"""Input-domain experiments (PROOF.md section 13).

D1  negative weights: the specification and fast paths with the weight check relaxed (a
    subclass of Prepared defined HERE; the package and working copy are unchanged), checked
    against the brute force and both pointer checkers on RS-style and hostile instances.
    The native engine must reject negative weights with an error, never answer.
D2  native int64 domain: inputs whose explicit bound is just below 2^60 must be solved
    exactly (compared with the arbitrary-precision specification); at or above 2^60 they
    must raise NativeRangeError. Includes coordinates near 2^60 and gamma near 2^60.
D3  S >= 2^29: the native class engine's int64 range-scan branch (never used on MHC) with
    length-variant alleles of size ~2^30, against the specification and brute force.
D4  per-haplotype lane bound: keys >= 2^38 must be refused by the per-haplotype engine
    while the class engine still answers exactly.
D5  wide panels with mid-graph haplotypes (H in 200..1024, native W = 4..16): native compact
    records identical to the specification, native expanded/packed/per-haplotype C, B, R equal
    to direct_oracle, verify_pointers on every native pointer (too large for the brute force).
Writes results/domain.json. Single process. `--only=D4` reruns one part and keeps the others.
"""
from __future__ import annotations
import collections
import json
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import harness as H  # noqa: E402
import families as fam  # noqa: E402
import hostile as ho  # noqa: E402
gc, gn = H.gc, H.gn


class RelaxedPrepared(gc.Prepared):
    """Identical index; check_anchors allows negative integer weights (D1 only)."""

    def check_anchors(self, anchors, gamma):
        neg = [gc.Anchor(a.v, a.x, a.y, a.c, a.d, abs(a.weight)) for a in anchors]
        return super().check_anchors(neg, gamma)


def d1_negative(count):
    stats = collections.Counter()
    failures = []
    native_rejections = collections.Counter()
    for seed in range(count):
        if seed % 2:
            name, (L, E, P), A, gamma = fam.random_stress_instance(50_000 + seed)
        else:
            gen = list(ho.GENERATORS)[seed // 2 % len(ho.GENERATORS)]
            name, (L, E, P), A, gamma = ho.GENERATORS[gen](10_000 + seed)
        rng = random.Random(seed)
        A = [(v, x, y, c, d, w - rng.choice([0, 1, 3, 7, 20, 1000])) for (v, x, y, c, d, w) in A]
        stats['negative_weight_anchors'] += sum(a[5] < 0 for a in A)
        g = H.GraphCache(L, E, P, prepared_cls=RelaxedPrepared)
        f = H.check_instance(g, A, gamma, stats, native=False)
        if f:
            failures.append(H.failure_record(L, E, P, A, gamma, f, list(map(str, name))))
        if g.native is not None and any(a[5] < 0 for a in A):
            try:
                g.native.solve_compact([gc.Anchor(*a) for a in A], gamma)
                native_rejections['ANSWERED (unexpected)'] += 1
            except ValueError as e:
                native_rejections[str(e)[:60]] += 1
    return dict(instances=count, stats=dict(stats), failing_instances=len(failures),
                failures=failures[:30], native_on_negative=dict(native_rejections))


def native_vs_spec(G, p, anchors, gamma):
    """Native compact/expanded/packed outputs against the arbitrary-precision specification."""
    spec = gc.solve(p, anchors, gamma)
    C, ptr = spec.expand()
    comp = G.solve_compact(anchors, gamma)
    ok = comp.records() == spec.records and list(comp.B) == spec.B and comp.R_values() == spec.R
    for res in (comp.expand(), G.solve_packed(anchors, gamma, mode='class')):
        nC, nptr, nB, nR = res.as_dicts(anchors)
        ok &= (nC, nB, nR) == (C, spec.B, spec.R)
        gc.verify_pointers(p, anchors, gamma, nC, nB, nptr)
    return ok


def d2_int64():
    out = []
    LIMIT = 1 << 60
    g = gc.Graph([3, 5, 2, 4], [(0, 1), (0, 2), (1, 3), (2, 3)], [[0, 1, 3], [0, 2, 3], [0, 1, 3]])
    p = g.prepare()
    G = gn.GapGraph(p)
    S = p.S
    for slack in (-2, -1, 0, 1, 2, 1000):
        # bound = sum(w) + 2S + max x + max y + max c + max d + max w + gamma + 16
        base = [gc.Anchor(0, 1, 2, 1, 2, 0), gc.Anchor(1, 1, 3, 4, 6, 0), gc.Anchor(3, 2, 4, 8, 9, 0)]
        mx, my, mc, md = 2, 4, 8, 9
        gamma = 7
        target = LIMIT - 1 - slack  # desired bound value (bound < LIMIT accepted)
        rest = target - (2 * S + mx + my + mc + md + gamma + 16)
        # weights: w0 = w1 = small, w2 large; sum(w) + max(w) = rest
        w_small = 5
        w2 = (rest - 2 * w_small) // 2
        anchors = [gc.Anchor(a.v, a.x, a.y, a.c, a.d, w) for a, w in zip(base, (w_small, w_small, w2))]
        bound = sum(a.weight for a in anchors) + 2 * S + mx + my + mc + md + max(a.weight for a in anchors) + gamma + 16
        try:
            ok = native_vs_spec(G, p, anchors, gamma)
            out.append(dict(case='weights', bound_minus_limit=bound - LIMIT, result='answered', exact=ok))
        except gn.NativeRangeError as e:
            out.append(dict(case='weights', bound_minus_limit=bound - LIMIT, result='NativeRangeError'))
    for big in (LIMIT // 8, LIMIT // 4 - 100, LIMIT // 2):
        anchors = [gc.Anchor(0, 1, 1, big, big + 1, 9), gc.Anchor(3, 1, 1, big + 3, big + 3, 9),
                   gc.Anchor(1, 2, 2, big + 2, big + 2, 4)]
        for gamma in (0, big // 3):
            bound = sum(a.weight for a in anchors) + 2 * S + 2 + 2 + (big + 3) + (big + 3) + 9 + gamma + 16
            try:
                ok = native_vs_spec(G, p, anchors, gamma)
                out.append(dict(case='coordinates', big=big, gamma=gamma, bound_minus_limit=bound - LIMIT,
                                result='answered', exact=ok))
            except gn.NativeRangeError:
                out.append(dict(case='coordinates', big=big, gamma=gamma, bound_minus_limit=bound - LIMIT,
                                result='NativeRangeError'))
    for gamma in (LIMIT - 300, LIMIT - 207, LIMIT - 206, LIMIT - 205, LIMIT):
        anchors = [gc.Anchor(0, 1, 1, 1, 1, 50), gc.Anchor(3, 1, 1, 5, 5, 50)]
        try:
            ok = native_vs_spec(G, p, anchors, gamma)
            out.append(dict(case='gamma', gamma_minus_limit=gamma - LIMIT, result='answered', exact=ok))
        except gn.NativeRangeError:
            out.append(dict(case='gamma', gamma_minus_limit=gamma - LIMIT, result='NativeRangeError'))
        except ValueError as e:
            out.append(dict(case='gamma', gamma_minus_limit=gamma - LIMIT, result='ValueError: ' + str(e)[:80]))
    return out


def d3_large_S(count):
    stats = collections.Counter()
    failures = []
    rng = random.Random(4242)
    branch = None
    for k in range(count):
        big = (1 << 30) + rng.randint(0, 1 << 20)
        lengths = [rng.randint(1, 3), big, big + rng.choice([0, 1, 7, 1 << 25]), rng.randint(1, 3), 2]
        edges = [(0, 1), (0, 2), (1, 3), (2, 3), (3, 4), (0, 3)]
        paths = rng.sample([[0, 1, 3, 4], [0, 2, 3, 4], [0, 3, 4], [1, 3], [0, 2, 3], [3, 4]], rng.randint(2, 5))
        verts = sorted({v for q in paths for v in q})
        A = []
        for _ in range(rng.randint(3, 6)):
            v = rng.choice(verts)
            L = lengths[v]
            x = rng.choice([1, L, rng.randint(1, L)])
            y = rng.choice([x, L])
            c = rng.randint(1, 6)
            A.append((v, x, y, c, c + rng.choice([0, 1]), rng.choice([1 << 31, 1 << 32, (1 << 31) + 5])))
        gamma = rng.choice([0, 3, 1 << 30])
        g = H.GraphCache(lengths, edges, paths)
        if branch is None and g.native is not None:
            branch = 'int64 range scans (S >= 2^29)' if g.p.S >= (1 << 29) else 'int32'
        f = H.check_instance(g, A, gamma, stats)
        if f:
            failures.append(H.failure_record(lengths, edges, paths, A, gamma, f, ['D3', k]))
    return dict(instances=count, native_branch=branch, stats=dict(stats), failing_instances=len(failures),
                failures=failures[:20])


def d4_lane_bound():
    g = gc.Graph([4, 2, 6, 3], [(0, 1), (0, 2), (1, 3), (2, 3)], [[0, 1, 3], [0, 2, 3]])
    p = g.prepare()
    G = gn.GapGraph(p)
    out = []
    # max_key = 2w + 3 + S(15) + max y(2) + max d(6) = 2w + 26: 2^37 - 14 gives 2^38 - 2 (accepted), 2^37 - 13 gives 2^38.
    for w in (1 << 36, (1 << 37) - 14, (1 << 37) - 13, 1 << 37, 1 << 40):
        anchors = [gc.Anchor(0, 1, 1, 1, 1, w), gc.Anchor(1, 1, 2, 3, 4, 3), gc.Anchor(3, 1, 1, 6, 6, w)]
        b = G.batch(anchors, 5)
        rec = dict(weight=w, max_key=b.max_key)
        rec['class_exact'] = native_vs_spec(G, p, anchors, 5)
        try:
            res = G.solve_packed(b, mode='per_haplotype')
            nC, nptr, nB, nR = res.as_dicts(anchors)
            eC, eB, eR, _ = gc.direct_oracle(p, anchors, 5)
            rec['per_haplotype'] = 'answered, exact=%s' % ((nC, nB, nR) == (eC, eB, eR))
            rec['consistent'] = b.max_key < (1 << 38) and (nC, nB, nR) == (eC, eB, eR)
        except gn.NativeRangeError:
            rec['per_haplotype'] = 'NativeRangeError'
            rec['consistent'] = b.max_key >= (1 << 38)
        out.append(rec)
    return out


def d5_wide_mid(count):
    rng = random.Random(5150)
    stats = collections.Counter()
    failures = []
    for k in range(count):
        V = rng.randint(3, 8)
        lengths = [rng.randint(1, 6) for _ in range(V)]
        fams = []
        for _ in range(rng.randint(2, 6)):
            s0 = rng.choice([0, rng.randrange(V)])
            e0 = rng.choice([V - 1, rng.randrange(s0, V)])
            fams.append([s0] + [v for v in range(s0 + 1, e0) if rng.random() < .6] + ([e0] if e0 > s0 else []))
        H = rng.choice([200, 256, 257, 511, 513, 700, 1000, 1024])
        paths = [list(rng.choice(fams)) for _ in range(H)]
        edges = sorted({(u, v) for p in paths for u, v in zip(p, p[1:])} |
                       {(u, v) for u in range(V) for v in range(u + 1, V) if rng.random() < .1})
        verts = sorted({v for p in paths for v in p})
        A = []
        for _ in range(rng.randint(3, 12)):
            v = rng.choice(verts)
            x = rng.randint(1, lengths[v]); y = rng.randint(x, lengths[v])
            c = rng.randint(1, 20)
            A.append(gc.Anchor(v, x, y, c, c + rng.choice([0, 1]), rng.choice([0, 1, 3, 9, 40])))
        gamma = rng.choice([0, 1, 4, 30])
        p = gc.Graph(lengths, edges, paths).prepare()
        G = gn.GapGraph(p)
        spec = gc.solve(p, A, gamma)
        sC, _ = spec.expand()
        eC, eB, eR, _ = gc.direct_oracle(p, A, gamma)
        probs = []
        if (sC, spec.B, spec.R) != (eC, eB, eR):
            probs.append('spec vs oracle')
        comp = G.solve_compact(A, gamma)
        if comp.records() != spec.records:
            probs.append('native records')
        outs = [comp.expand(), G.solve_packed(A, gamma, mode='class'), G.solve_packed(A, gamma, mode='class', warm=True)]
        for cut in (gn.DENSE_CUT, 0, 257):
            outs.append(G.solve_packed(A, gamma, mode='per_haplotype', dense_cut=cut))
        for res in outs:
            nC, nptr, nB, nR = res.as_dicts(A)
            if (nC, nB, nR) != (eC, eB, eR):
                probs.append('values ' + res.stats['mode'])
            try:
                gc.verify_pointers(p, A, gamma, nC, nB, nptr)
            except AssertionError as e:
                probs.append('pointers ' + res.stats['mode'] + ' ' + repr(e)[:200])
        stats['instances'] += 1
        stats['states'] += sum(map(len, eC))
        stats['W%d' % G.W] += 1
        stats['missing_without_last2reach'] += sum(1 for z in A for h in range(p.H) if not p.on[z.v][h] and p.last[z.v][h] < 0)
        if probs:
            failures.append(dict(k=k, lengths=lengths, edges=edges, paths=paths, H=H,
                                 anchors=[list(vars(z).values()) for z in A], gamma=gamma, problems=probs))
    return dict(stats=dict(stats), failing_instances=len(failures), failures=failures[:10])


def main():
    assert sys.flags.optimize == 0
    runs = dict(D1_negative_weights=lambda: d1_negative(4000), D2_int64=d2_int64,
                D3_large_S=lambda: d3_large_S(400), D4_lane_bound=d4_lane_bound, D5_wide_mid=lambda: d5_wide_mid(120))
    only = [a.split('=', 1)[1] for a in sys.argv[1:] if a.startswith('--only=')]
    path = HERE / 'results/domain.json'
    report = json.loads(path.read_text()) if only and path.exists() else {}
    for key, fn in runs.items():
        if not only or key.split('_')[0] in only:
            report[key] = fn()
    (HERE / 'results').mkdir(exist_ok=True)
    (HERE / 'results/domain.json').write_text(json.dumps(report, indent=1, default=str) + '\n')
    print(json.dumps({k: (v if not isinstance(v, dict) else {kk: vv for kk, vv in v.items() if kk != 'failures'})
                      for k, v in report.items()}, indent=1, default=str))


if __name__ == '__main__':
    main()
