"""Exactness of the native engines on a new locus (adapted from hprc_mhc_panel/verify_panels.py).

For each selected primary instance (Minichain scaling, weights x 200) it requires:
  * native class compact records identical to gap_chaining.solve (value, origin, mask and
    witness of every record), with identical B and R;
  * identical C, B, R from compact-expanded, class packed (cold and warm) and per-haplotype;
  * the same C, B, R as real_check_mosaic.py's independent graph-major per-haplotype `mine`,
    run on the panel's graph with its own index;
  * gap_chaining.verify_pointers on the class pointers and the per-haplotype pointers.
Defaults (as verify_panels.py): all 40 primary instances on the all-haplotype panel, 10 per
subset panel, 10 on the all-haplotype panel with gamma 20000, and 10 on the leave-out panel.
Instances are spread by size (sorted by anchor count, every k-th). Every failure is recorded
with its instance id; the run continues after a failure and the status says FAILED.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import locus_panels as lp  # noqa: E402
import gap_chaining as gc  # noqa: E402

lp.native_library_ready()
import gap_native as gn  # noqa: E402


def load_done(logs):
    """Rows already verified and printed in full by earlier (stopped) runs: (panel, gamma, id) -> (row, log)."""
    done = {}
    for log in logs:
        for line in log.read_text().splitlines():
            if line.startswith('{"locus"'):
                r = json.loads(line)
                key = (r['panel'], r['gamma'], r['id'])
                if key in done:
                    assert done[key][0]['passed'] == r['passed'] and done[key][0].get('digest') == r.get('digest'), key
                    continue
                done[key] = (r, log.name)
    return done


def check(locus, name, count, gamma, groups, done=None):
    p = lp.prepared(locus, name)
    G = gn.GapGraph(p)
    ns = lp.independent_mine(locus, name, gamma)
    cases = sorted(lp.primary_cases(locus, name), key=lambda c: len(c['anchors']))
    picks = cases if count >= len(cases) else cases[::max(1, len(cases) // count)][:count]
    rows = []
    for case in picks:
        if done and (name, gamma, case['id']) in done:  # verified by an earlier run; reused as logged
            prior, log = done[(name, gamma, case['id'])]
            row = dict(prior, resumed_from=log)
            rows.append(row)
            print(json.dumps(row), flush=True)
            continue
        row = dict(locus=locus, panel=name, H=p.H, gamma=gamma, id=case['id'], group=groups[case['read']]['group'],
                   source=groups[case['read']]['source'])
        try:
            A = lp.anchors(case)
            b = G.batch(A, gamma)
            row.update(N=b.n, K=b.K)
            spec = gc.solve(p, A, gamma)
            comp = G.solve_compact(b)
            row['record_identity'] = (comp.records() == spec.records and list(comp.B) == spec.B
                                      and comp.R_values() == spec.R)
            del spec
            outs = [comp.expand(), G.solve_packed(b, mode='class'), G.solve_packed(b, mode='class', warm=True),
                    G.solve_packed(b, mode='per_haplotype')]
            row['native_arrays_identical'] = len({o.arrays_digest() for o in outs}) == 1
            row['digest'] = outs[0].arrays_digest()
            keep = (outs[0], outs[3])
            del outs, comp
            mC, mB, mR = ns['mine'](A)
            agree, row['pointer_counts'] = True, {}
            for res in keep:
                C, ptr, B, R = res.as_dicts(A)
                agree &= (C == mC and B == mB and R == mR)
                row['pointer_counts'][res.stats['mode']] = gc.verify_pointers(p, A, gamma, C, B, ptr)
                del C, ptr
            del mC, mB, mR, keep
            row['independent_agree'] = agree
            row['passed'] = row['record_identity'] and row['native_arrays_identical'] and agree
        except Exception:  # record the failure with its input and continue
            row['passed'] = False
            row['error'] = traceback.format_exc()[-3000:]
        rows.append(row)
        print(json.dumps(row), flush=True)
    G.warm_cache_reset()
    return dict(locus=locus, panel=name, H=p.H, gamma=gamma, reads=len(rows), passed=sum(r['passed'] for r in rows),
                states=sum(r.get('K', 0) for r in rows), anchors=sum(r.get('N', 0) for r in rows),
                failures=[dict(id=r['id'], panel=name, gamma=gamma, error=r.get('error'),
                               checks={k: r.get(k) for k in ('record_identity', 'native_arrays_identical', 'independent_agree')})
                          for r in rows if not r['passed']], rows=rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('locus')
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--full', type=int, default=40, help='instances on the all-haplotype panel')
    ap.add_argument('--subset', type=int, default=10, help='instances per subset panel')
    ap.add_argument('--gamma20000', type=int, default=10)
    ap.add_argument('--leaveout', type=int, default=10)
    ap.add_argument('--levels', default='full,subsets,gamma20000,leaveout', help='which levels to run')
    ap.add_argument('--resume', type=Path, nargs='*', default=[],
                    help='logs of earlier stopped runs: instances printed there in full are reused, not recomputed')
    args = ap.parse_args()
    if args.out.exists():
        raise SystemExit('refusing to overwrite ' + str(args.out))
    run = set(args.levels.split(','))
    start = time.perf_counter()
    lib = gn.build()._name
    info = lp.panels(args.locus)
    groups = lp.read_groups(args.locus)
    done = load_done(args.resume)
    levels = []
    if 'full' in run and args.full:
        levels.append(check(args.locus, info['all'], args.full, lp.GAMMA, groups, done))
    if 'subsets' in run and args.subset:
        for name in info['subsets']:
            levels.append(check(args.locus, name, args.subset, lp.GAMMA, groups, done))
    if 'gamma20000' in run and args.gamma20000:
        levels.append(check(args.locus, info['all'], args.gamma20000, 20000, groups, done))
    if 'leaveout' in run and args.leaveout:
        levels.append(check(args.locus, info['leaveout'], args.leaveout, lp.GAMMA, groups, done))
    failed = sum(lv['reads'] - lv['passed'] for lv in levels)
    out = dict(script='experiments/second_locus/verify_locus.py', argv=sys.argv[1:], python=sys.version, library=lib,
               levels=levels, status='all passed' if not failed else 'FAILED', failures=failed,
               instances=sum(lv['reads'] for lv in levels), states=sum(lv['states'] for lv in levels),
               seconds=round(time.perf_counter() - start, 1),
               resumed_rows=sum('resumed_from' in r for lv in levels for r in lv['rows']),
               resumed_from=sorted(str(x.name) for x in args.resume))
    args.out.write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(dict(status=out['status'], seconds=out['seconds'], instances=out['instances'], states=out['states'],
                          per_level={f"{lv['panel']}_g{lv['gamma']}": [lv['reads'], lv['passed'], lv['states']] for lv in levels})))
    return 0 if not failed else 1


if __name__ == '__main__':
    sys.exit(main())
