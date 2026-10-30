"""Independent checks of the native emulation (giraffe_emul.cpp), written separately from it.

1. Shortest distances: native D(u,v) against Python Dijkstra on random vertex pairs.
2. Full instances: pure-Python transitions (own rule), pure-Python Giraffe rule and pure-Python
   exact per-haplotype DP, compared record by record with the native outputs. Uses real instances
   on the 60-haplotype panel and on random small haplotype subsets (which force switches).
3. Brute force: small sub-instances; every chain enumerated, every labeling enumerated, giving
   the optimum of chains ending at each anchor; compared with the native exact B_i, and the
   heuristic's stored scores are checked to be values of their traced chains and <= the optimum.

  python3.13 experiments/exactness_payoff/check_independent.py --out experiments/exactness_payoff/checks_pruned.json
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
import random
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import emul  # noqa: E402

ROOT = HERE.parents[1]
DATA = {'hprc60': ROOT / 'experiments/hprc_mhc_panel/data/hprc60', 'mhc60': ROOT / 'experiments/mhc_pilot/data'}


# ------------------------------------------------------------ pure-Python reference implementations

def py_transitions(panel, sa, qlen):
    """Transitions into each destination of the sorted anchors, own implementation of SPEC section 3."""
    gl = max(20000, int(0.10501002120802233 * qlen))
    il = max(5000, int(2.45 * qlen))
    out = [[] for _ in sa]
    for i, a in enumerate(sa):
        for j in range(i):
            j0 = emul.py_transition(panel, sa[j], a, gl, il)
            if j0 is not None:
                out[i].append((j, j0))
    return out


def py_giraffe(items, masks, trans, rho, cb):
    """Literal re-implementation of chain_items_dp (vg v1.77.0 lines 417-582), Python integers."""
    n = len(items)
    INF_SRC = float('inf')  # nowhere = SIZE_MAX
    score = list(items)
    src = [INF_SRC] * n
    paths = list(masks)
    rec = [0] * n
    evb = [cb] * n
    for i in range(n):
        for j, j0 in trans[i]:
            # from nowhere
            e_now, e_cur = items[i] + cb, score[i] + evb[i]
            if e_now > e_cur or (e_now == e_cur and (items[i], INF_SRC) > (score[i], src[i])):
                score[i], src[i], paths[i], rec[i], evb[i] = items[i], INF_SRC, masks[i], 0, cb
            common = paths[j] & masks[i]
            jump = j0 - (rho if common == 0 else 0)
            cs = score[j] + jump + items[i]
            cp = masks[i] if common == 0 else common
            cr = rec[j] + (1 if common == 0 else 0)
            bonus = 0
            if cb > 0:
                pre = bin(paths[j]).count('1')
                if pre > 0 and common != 0:
                    bonus = (cb * bin(cp).count('1')) // pre
            e_from, e_best = cs + bonus, score[i] + evb[i]
            if e_from > e_best or (e_from == e_best and (cs, j) > (score[i], src[i])):
                score[i], src[i], paths[i], rec[i], evb[i] = cs, j, cp, cr, bonus
    return score, [-1 if s == INF_SRC else s for s in src], paths, rec


def py_exact(items, masks, trans, rho, H):
    """Per-haplotype DP: C[i][h] = item + max(0, max_j max(C[j][h], B_j - rho) + jump0)."""
    n = len(items)
    C, B = [], []
    for i in range(n):
        labs = [h for h in range(H) if masks[i] >> h & 1]
        row = {}
        for h in labs:
            best = 0
            for j, j0 in trans[i]:
                v = max(C[j].get(h, -10 ** 18), B[j] - rho) + j0
                if v > best:
                    best = v
            row[h] = items[i] + best
        C.append(row)
        B.append(max(row.values()))
    return C, B


# ------------------------------------------------------------ brute force

def brute_force(items, masks, trans, rho, H):
    """Optimum of chains ending at each anchor: all chains, all labelings (H small)."""
    n = len(items)
    succ = [[] for _ in range(n)]
    for i in range(n):
        for j, j0 in trans[i]:
            succ[j].append((i, j0))
    best = [None] * n
    chains = 0
    stats = dict(exhaustive=0, label_dp=0)

    def labelings_min(chain):
        lab = [[h for h in range(H) if masks[k] >> h & 1] for k in chain]
        size = 1
        for x in lab:
            size *= len(x)
        if size <= 4096:  # enumerate every labeling
            stats['exhaustive'] += 1
            m = None
            for combo in itertools.product(*lab):
                c = sum(1 for a, b in zip(combo, combo[1:]) if a != b)
                m = c if m is None or c < m else m
            return m
        stats['label_dp'] += 1  # too many labelings: forward label DP (independent of the native code)
        return emul.min_switch_labels([masks[k] for k in chain], H)[0]

    def dfs(chain, base):
        nonlocal chains
        chains += 1
        end = chain[-1]
        v = base - rho * labelings_min(chain)
        if best[end] is None or v > best[end]:
            best[end] = v
        for nxt, j0 in succ[end]:
            dfs(chain + [nxt], base + j0 + items[nxt])

    for s in range(n):
        dfs([s], items[s])
    return best, chains, stats


# ------------------------------------------------------------ checks

def check_distances(panel, rng, pairs=2000):
    bad = 0
    for _ in range(pairs):
        u = rng.randrange(panel.V)
        v = rng.randrange(u, panel.V)
        if panel.native_dist(u, v) != panel.py_dist(u, v):
            bad += 1
    return dict(pairs=pairs, mismatches=bad)


def compare_full(panel, case, rho, bonuses=(0, 12)):
    r = emul.solve(panel, case['anchors'], case['qlen'], rho=rho, bonuses=bonuses)
    sa = r.anchors
    items = [emul.item_points(a) for a in sa]
    masks = [panel.A[a[0]] for a in sa]
    trans = py_transitions(panel, sa, case['qlen'])
    T = sum(len(t) for t in trans)
    assert T == r.stats['transitions'], (T, r.stats['transitions'])
    n = r.n
    mism = 0
    for b, cb in enumerate(bonuses):
        s, src, paths, rec = py_giraffe(items, masks, trans, rho, cb)
        for i in range(n):
            k = b * n + i
            if (s[i], src[i], paths[i], rec[i]) != (r.hs[k], r.hsrc[k], r.hmask[k], r.hrec[k]):
                mism += 1
    C, B = py_exact(items, masks, trans, rho, panel.H)
    for i in range(n):
        if B[i] != r.B[i]:
            mism += 1
        for h in range(panel.H):
            want = C[i].get(h)
            got = r.C[i * panel.H + h]
            if (want is None and got != emul.NEG) or (want is not None and want != got):
                mism += 1
    opt = max(B)
    heur = [max(r.hs[b * n:(b + 1) * n]) for b in range(len(bonuses))]
    return dict(id=case['id'], H=panel.H, rho=rho, N=n, transitions=T, mismatches=mism, opt=opt,
                heuristic_best=heur, exact_switches=None)


def sub_instance(rng, sa, chain_idx, m):
    """m anchors: a strided sample of an exact chain plus a few random others, re-sorted."""
    stride = rng.randint(3, 60)
    start = rng.randrange(max(1, len(chain_idx) - stride * m))
    picked = set(chain_idx[start::stride][:m - 3])
    while len(picked) < m:
        picked.add(rng.randrange(len(sa)))
    return [sa[k] for k in sorted(picked)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--seed', type=int, default=20261006)
    ap.add_argument('--full', type=int, default=6, help='full real instances per dataset')
    ap.add_argument('--subset_cases', type=int, default=6, help='full instances on small haplotype subsets')
    ap.add_argument('--brute', type=int, default=150, help='brute-force sub-instances')
    args = ap.parse_args()
    out = Path(args.out)
    if out.exists():
        raise SystemExit(f'{out} exists')
    rng = random.Random(args.seed)
    t0 = time.perf_counter()
    report = dict(seed=args.seed, library=emul.build().identity, distances={}, full=[], subset_full=[], brute=[])
    for name, d in DATA.items():
        panel = emul.load_panel(d / 'graph.json', name)
        report['distances'][name] = check_distances(panel, rng)
        cases = [json.loads(s) for s in (d / 'instances.jsonl').read_text().splitlines()]
        # Smallest instances whose transition count keeps pure Python fast, plus a mid-size one.
        sized = sorted(cases, key=lambda c: len(c['anchors']))
        chosen = sized[:args.full - 1] + [sized[len(sized) // 3]]
        for c in chosen:
            for rho in (2, 10):
                report['full'].append(dict(dataset=name, **compare_full(panel, c, rho)))
        # Small haplotype subsets force recombination; full instances, native vs pure Python.
        for k in range(args.subset_cases):
            haps = sorted(rng.sample(range(panel.H), rng.choice((3, 5, 8))))
            sub, vmap = emul.subset_panel(panel, haps, f'{name}-sub{k}')
            c = rng.choice(sized[:len(sized) // 2])
            anchors = [[vmap[a[0]], *a[1:]] for a in c['anchors'] if vmap[a[0]] >= 0]
            if len(anchors) < 2:
                continue
            case = dict(c, anchors=anchors)
            for rho in (2, 28):
                row = compare_full(sub, case, rho)
                row.update(dataset=name, haplotypes=haps)
                report['subset_full'].append(row)
        # Brute force on sub-instances drawn from mid-size true-strand instances.
        pool = [c for c in cases if 1500 <= len(c['anchors']) <= 6000]
        done = 0
        while done < args.brute // len(DATA):
            haps = sorted(rng.sample(range(panel.H), rng.choice((2, 3, 4))))
            sub, vmap = emul.subset_panel(panel, haps)
            c = rng.choice(pool)
            anchors = [[vmap[a[0]], *a[1:]] for a in c['anchors'] if vmap[a[0]] >= 0]
            if len(anchors) < 12:
                continue
            full = emul.solve(sub, anchors, c['qlen'], rho=2, bonuses=(0,))
            chain_idx = [i for i, _ in emul.exact_chain(full)]
            if len(chain_idx) < 6:
                continue
            m = rng.randint(6, 11)
            sa = sub_instance(rng, full.anchors, chain_idx, m)
            rho = rng.choice((1, 2, 5, 10, 28))
            r = emul.solve(sub, sa, c['qlen'], rho=rho, bonuses=(0, 12))
            items = list(r.item)
            masks = [sub.A[a[0]] for a in r.anchors]
            trans = py_transitions(sub, r.anchors, c['qlen'])
            if sum(len(t) for t in trans) == 0:
                continue
            best, nchains, lstats = brute_force(items, masks, trans, rho, sub.H)
            ok_exact = all(best[i] == r.B[i] for i in range(r.n))
            ok_heur = True
            for b in range(2):
                for i in range(r.n):
                    s = r.hs[b * r.n + i]
                    if s > best[i]:
                        ok_heur = False
                    ch = emul.heuristic_chain(r, b, end=i)
                    val = sum(items[k] for k in ch)
                    for x, y in zip(ch, ch[1:]):
                        val += dict(trans[y])[x]
                    val -= rho * emul.greedy_resets([masks[k] for k in ch])[0]
                    if val != s:
                        ok_heur = False
            gaps = [max(best) - max(r.hs[b * r.n:(b + 1) * r.n]) for b in range(2)]
            report['brute'].append(dict(dataset=name, id=c['id'], haplotypes=haps, rho=rho, n=r.n,
                                        transitions=sum(len(t) for t in trans), chains=nchains,
                                        labelings_exhaustive=lstats['exhaustive'], labelings_dp=lstats['label_dp'],
                                        exact_ok=ok_exact, heuristic_ok=ok_heur, opt=max(best), gaps=gaps))
            done += 1
    report['seconds'] = round(time.perf_counter() - t0, 1)
    report['summary'] = dict(
        distance_mismatches=sum(v['mismatches'] for v in report['distances'].values()),
        full_instances=len(report['full']), full_mismatches=sum(x['mismatches'] for x in report['full']),
        subset_instances=len(report['subset_full']),
        subset_mismatches=sum(x['mismatches'] for x in report['subset_full']),
        subset_with_heuristic_gap=sum(1 for x in report['subset_full'] if max(x['heuristic_best']) < x['opt']
                                      or min(x['heuristic_best']) < x['opt']),
        brute_cases=len(report['brute']), brute_exact_failures=sum(1 for x in report['brute'] if not x['exact_ok']),
        brute_heuristic_failures=sum(1 for x in report['brute'] if not x['heuristic_ok']),
        brute_chains=sum(x['chains'] for x in report['brute']),
        brute_chain_labelings_exhaustive=sum(x['labelings_exhaustive'] for x in report['brute']),
        brute_chain_labelings_dp=sum(x['labelings_dp'] for x in report['brute']),
        brute_with_gap=sum(1 for x in report['brute'] if max(x['gaps']) > 0))
    out.write_text(json.dumps(report, indent=1) + '\n')
    print(json.dumps(report['summary'], indent=1))


if __name__ == '__main__':
    main()
