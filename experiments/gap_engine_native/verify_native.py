"""Stage B exactness on real MHC reads for the native engines (gap_native.cpp).

For each selected read (Minichain scaling: weights x 200) and panel it requires:
  * native class compact records identical to the specification gap_chaining.solve
    (value, origin, mask and witness of every record), with identical B and R;
  * identical C, B, R arrays from every native output: compact expanded, class packed
    (cold and warm pair cache) and per-haplotype;
  * the same C, B, R as real_check_mosaic.py's independent graph-major per-haplotype
    implementation `mine` (its own offsets, last2reach and distances);
  * gap_chaining.verify_pointers on the class pointers and on the per-haplotype pointers.
Defaults: all 40 primary reads at H=60, six reads at H=960 (spread by size), and six
H=60 reads with gamma 20000.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import common as cm  # noqa: E402
import gap_chaining as gc  # noqa: E402
import gap_native as gn  # noqa: E402


def check_level(H, count, gamma):
    cases = cm.primary_cases()
    order = sorted(range(len(cases)), key=lambda k: len(cases[k]['anchors']))
    picks = order if count >= len(order) else order[::max(1, len(order) // count)][:count]
    p = cm.prepared(H)
    G = gn.GapGraph(p)
    ns = cm.independent_mine(H)
    ns['GAMMA'] = gamma
    rows = []
    for k in picks:
        case = cases[k]
        A = cm.anchors(case)
        b = G.batch(A, gamma)
        spec = gc.solve(p, A, gamma)
        comp = G.solve_compact(b)
        row = dict(H=H, gamma=gamma, id=case['id'], N=b.n, K=b.K,
                   record_identity=comp.records() == spec.records and list(comp.B) == spec.B
                   and comp.R_values() == spec.R)
        outs = [comp.expand(), G.solve_packed(b, mode='class'), G.solve_packed(b, mode='class', warm=True),
                G.solve_packed(b, mode='per_haplotype')]
        row['native_arrays_identical'] = len({o.arrays_digest() for o in outs}) == 1
        row['digest'] = outs[0].arrays_digest()
        mC, mB, mR = ns['mine'](A)
        row['pointer_counts'] = {}
        agree = True
        for res in (outs[0], outs[3]):
            C, ptr, B, R = res.as_dicts(A)
            agree &= (C == mC and B == mB and R == mR)
            row['pointer_counts'][res.stats['mode']] = gc.verify_pointers(p, A, gamma, C, B, ptr)
            del C, ptr
        row['independent_agree'] = agree
        row['passed'] = row['record_identity'] and row['native_arrays_identical'] and agree
        rows.append(row)
        print(json.dumps(row), flush=True)
    G.warm_cache_reset()
    return dict(H=H, gamma=gamma, reads=len(rows), passed=sum(r['passed'] for r in rows),
                states=sum(r['K'] for r in rows), rows=rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--real60', type=int, default=40)
    ap.add_argument('--real960', type=int, default=6)
    ap.add_argument('--gamma20000', type=int, default=6)
    args = ap.parse_args()
    if args.out.exists():
        raise SystemExit('refusing to overwrite ' + str(args.out))
    start = time.perf_counter()
    gn.build()
    levels = []
    if args.real60:
        levels.append(check_level(60, args.real60, cm.GAMMA))
    if args.gamma20000:
        levels.append(check_level(60, args.gamma20000, 20000))
    if args.real960:
        levels.append(check_level(960, args.real960, cm.GAMMA))
    failed = sum(lv['reads'] - lv['passed'] for lv in levels)
    out = dict(script='experiments/gap_engine_native/verify_native.py', argv=sys.argv[1:], python=sys.version,
               levels=levels, status='all passed' if not failed else 'FAILED', failures=failed,
               seconds=round(time.perf_counter() - start, 1))
    args.out.write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(dict(status=out['status'], seconds=out['seconds'],
                          states={f"H{lv['H']}_g{lv['gamma']}": lv['states'] for lv in levels})))
    return 0 if not failed else 1


if __name__ == '__main__':
    sys.exit(main())
