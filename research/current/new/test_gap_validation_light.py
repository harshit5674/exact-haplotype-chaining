"""Cross-check the memory-light validator (gap_validation_light) against gap_validation (2026-10-20).

On the families of test_gap_validation.py (fixed cases and 350 wide random cases), both validators must
accept every engine output and the direct oracle must agree; on single-field corruptions of C, B, R and
pointers, both must reject. With --panel, the same agreement is also checked on real panel instances.

  python3.13 -B research/current/new/test_gap_validation_light.py --out <fresh-result.json> \
      [--panel experiments/hprc_mhc_panel/data/hprc60 --instances 40]
"""
import argparse
import json
import sys
import time
from pathlib import Path

import gap_chaining as gc
import gap_native as gn
import gap_validation as gv
import gap_validation_light as gvl
from test_gap_chaining import random_case, fixed_cases


def corrupt_and_check(result, batches):
    rejected = 0
    for field in ('C', 'B', 'R', 'back'):
        arr = getattr(result, field)
        if not len(arr):
            continue
        old = arr[0]
        arr[0] = len(result.C) if field == 'back' else old + 1
        try:
            for vb in batches:
                try:
                    vb.check(result)
                except AssertionError:
                    rejected += 1
                else:
                    raise AssertionError(('corruption accepted', field, type(vb).__module__))
        finally:
            arr[0] = old
    return rejected


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--panel', type=Path)
    ap.add_argument('--instances', type=int, default=40)
    args = ap.parse_args()
    if args.out.exists():
        raise SystemExit('refusing overwrite')
    t = time.perf_counter()
    passed = states = skipped = rejected = 0
    cases = [case for _, case in fixed_cases()] + [random_case(10000 + s, wide=True) for s in range(350)]
    for graph, A, gamma in cases:
        p = graph.prepare()
        try:
            G = gn.GapGraph(p)
            b = G.batch(A, gamma)
        except (gn.NativeRangeError, OverflowError):
            skipped += 1
            continue
        if not b.max_key < 1 << 38:
            skipped += 1
            continue
        gdict = dict(lengths=graph.lengths, edges=graph.edges, paths=graph.paths)
        v0, v1 = gv.Validator(gdict), gvl.Validator(gdict)
        b0, b1 = v0.batch(A, gamma), v1.batch(A, gamma)
        expected = gc.direct_oracle(p, A, gamma)
        for mode in ('per_haplotype', 'class', 'auto'):
            result = G.solve_packed(b, mode=mode)
            b0.check(result)
            b1.check(result)
            C, ptr, B, R = result.as_dicts(A)
            assert (C, B, R) == expected[:3]
        if b.K:
            rejected += corrupt_and_check(result, (b0, b1))
        for x in (b0, b1, v0, v1):
            x.close()
        passed += 1
        states += b.K
    panel = None
    if args.panel:
        sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'experiments/gap_fallback'))
        import auto_real as data
        graph, instances = data.load(args.panel)
        p = gc.Graph(graph['lengths'], [tuple(e) for e in graph['edges']], graph['paths']).prepare()
        G = gn.GapGraph(p)
        v0, v1 = gv.Validator(graph), gvl.Validator(graph)
        n = pstates = 0
        for c in instances[:args.instances]:
            A = [gc.Anchor(*x[:5], x[5] * data.SCALE) for x in c['anchors']]
            b0, b1 = v0.batch(A, data.GAMMA), v1.batch(A, data.GAMMA)
            for mode in ('per_haplotype', 'class', 'auto'):
                result = G.solve_packed(A, data.GAMMA, mode=mode)
                b0.check(result)
                b1.check(result)
            rejected += corrupt_and_check(result, (b0, b1))
            b0.close(); b1.close()
            n += 1
            pstates += len(result.C)
        v0.close(); v1.close()
        panel = dict(path=str(args.panel), instances=n, states=pstates)
    out = dict(passed=passed, states=states, skipped_outside_native_domain=skipped,
               corruptions_rejected_by_both=rejected, panel=panel, seconds=time.perf_counter() - t, argv=sys.argv)
    args.out.write_text(json.dumps(out, indent=2) + '\n')
    print(args.out.read_text())


if __name__ == '__main__':
    main()
