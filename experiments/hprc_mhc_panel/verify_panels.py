"""Exactness of the native engines on the real HPRC MHC panels.

For each selected primary instance (Minichain scaling, weights x 200) it requires:
  * native class compact records identical to gap_chaining.solve (value, origin, mask and
    witness of every record), with identical B and R;
  * identical C, B, R from compact-expanded, class packed (cold and warm) and per-haplotype;
  * the same C, B, R as real_check_mosaic.py's independent graph-major per-haplotype `mine`,
    run on the panel's graph with its own index;
  * gap_chaining.verify_pointers on the class pointers and the per-haplotype pointers.
Defaults: all 40 primary instances on the 463-haplotype panel, 10 per subset panel, 10 on
the 463-haplotype panel with gamma 20000, and 10 on the leave-CHM13-out panel.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panel_common as pc  # noqa: E402
import gap_chaining as gc  # noqa: E402
import gap_native as gn  # noqa: E402


def check(name, count, gamma):
    p = pc.prepared(name)
    G = gn.GapGraph(p)
    ns = pc.independent_mine(name, gamma)
    cases = sorted(pc.primary_cases(name), key=lambda c: len(c['anchors']))
    picks = cases if count >= len(cases) else cases[::max(1, len(cases) // count)][:count]
    rows = []
    for case in picks:
        A = pc.anchors(case)
        b = G.batch(A, gamma)
        spec = gc.solve(p, A, gamma)
        comp = G.solve_compact(b)
        row = dict(panel=name, H=p.H, gamma=gamma, id=case['id'], N=b.n, K=b.K,
                   record_identity=comp.records() == spec.records and list(comp.B) == spec.B
                   and comp.R_values() == spec.R)
        outs = [comp.expand(), G.solve_packed(b, mode='class'), G.solve_packed(b, mode='class', warm=True),
                G.solve_packed(b, mode='per_haplotype')]
        row['native_arrays_identical'] = len({o.arrays_digest() for o in outs}) == 1
        row['digest'] = outs[0].arrays_digest()
        mC, mB, mR = ns['mine'](A)
        agree, row['pointer_counts'] = True, {}
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
    return dict(panel=name, H=p.H, gamma=gamma, reads=len(rows), passed=sum(r['passed'] for r in rows),
                states=sum(r['K'] for r in rows), rows=rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--full', type=int, default=40, help='instances on the 463-haplotype panel')
    ap.add_argument('--subset', type=int, default=10, help='instances per subset panel')
    ap.add_argument('--gamma20000', type=int, default=10)
    ap.add_argument('--leaveout', type=int, default=10, help='instances on hprc462_noCHM13')
    args = ap.parse_args()
    if args.out.exists():
        raise SystemExit('refusing to overwrite ' + str(args.out))
    start = time.perf_counter()
    gn.build()
    levels = [check('hprc463', args.full, pc.GAMMA)]
    for name in ('hprc60', 'hprc120', 'hprc240'):
        levels.append(check(name, args.subset, pc.GAMMA))
    if args.gamma20000:
        levels.append(check('hprc463', args.gamma20000, 20000))
    if args.leaveout:
        levels.append(check('hprc462_noCHM13', args.leaveout, pc.GAMMA))
    failed = sum(lv['reads'] - lv['passed'] for lv in levels)
    out = dict(script='experiments/hprc_mhc_panel/verify_panels.py', argv=sys.argv[1:], python=sys.version,
               levels=levels, status='all passed' if not failed else 'FAILED', failures=failed,
               seconds=round(time.perf_counter() - start, 1))
    args.out.write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(dict(status=out['status'], seconds=out['seconds'],
                          states={f"{lv['panel']}_g{lv['gamma']}": lv['states'] for lv in levels})))
    return 0 if not failed else 1


if __name__ == '__main__':
    sys.exit(main())
