"""Hostile synthetic check: random small instances where Giraffe's single record often loses.

Random DAG transitions, random masks over few haplotypes, random items and gap terms; the native DP
core (ge_solve_explicit, the same code as the geometric solver) against exhaustive enumeration of
every chain and every labeling, plus the pure-Python Giraffe rule from check_independent.py.

  python3.13 experiments/exactness_payoff/check_synthetic.py --out experiments/exactness_payoff/checks_pruned_synthetic.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import random
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import emul  # noqa: E402
import check_independent as ci  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--cases', type=int, default=3000)
    ap.add_argument('--seed', type=int, default=6406)
    args = ap.parse_args()
    out = Path(args.out)
    if out.exists():
        raise SystemExit(f'{out} exists')
    rng = random.Random(args.seed)
    stats = dict(cases=0, exact_failures=0, heuristic_rule_mismatches=0, heuristic_above_opt=0,
                 heuristic_value_mismatch=0, with_gap={'0': 0, '12': 0}, chains=0)
    examples = []
    for case in range(args.cases):
        n = rng.randint(3, 9)
        H = rng.randint(2, 4)
        items = [rng.randint(5, 20) for _ in range(n)]
        masks = []
        for _ in range(n):
            m = 0
            while m == 0:
                m = rng.getrandbits(H)
            masks.append(m)
        p = rng.choice((0.3, 0.5, 0.8))
        trans = [[(j, -rng.randint(0, 4)) for j in range(i) if rng.random() < p] for i in range(n)]
        rho = rng.choice((1, 2, 5, 10))
        flat = [(j, i, j0) for i in range(n) for j, j0 in trans[i]]
        r = emul.solve_explicit(n, H, items, masks, flat, rho=rho, bonuses=(0, 12))
        best, nchains, _ = ci.brute_force(items, masks, trans, rho, H)
        stats['chains'] += nchains
        stats['cases'] += 1
        if any(best[i] != r.B[i] for i in range(n)):
            stats['exact_failures'] += 1
        for b, cb in enumerate((0, 12)):
            s, src, paths, rec = ci.py_giraffe(items, masks, trans, rho, cb)
            if any((s[i], src[i], paths[i], rec[i]) != (r.hs[b * n + i], r.hsrc[b * n + i], r.hmask[b * n + i],
                                                        r.hrec[b * n + i]) for i in range(n)):
                stats['heuristic_rule_mismatches'] += 1
            for i in range(n):
                if r.hs[b * n + i] > best[i]:
                    stats['heuristic_above_opt'] += 1
                ch = emul.heuristic_chain(r, b, end=i)
                val = sum(items[k] for k in ch) + sum(dict(trans[y])[x] for x, y in zip(ch, ch[1:]))
                val -= rho * emul.greedy_resets([masks[k] for k in ch])[0]
                if val != r.hs[b * n + i]:
                    stats['heuristic_value_mismatch'] += 1
            gap = max(best) - max(r.hs[b * n:(b + 1) * n])
            if gap > 0:
                stats['with_gap'][str(cb)] += 1
                if len(examples) < 5:
                    examples.append(dict(n=n, H=H, items=items, masks=masks, transitions=flat, rho=rho, bonus=cb,
                                         exact=max(best), heuristic=max(r.hs[b * n:(b + 1) * n])))
    out.write_text(json.dumps(dict(seed=args.seed, library=emul.build().identity, stats=stats,
                                   examples=examples), indent=1) + '\n')
    print(json.dumps(stats))


if __name__ == '__main__':
    main()
