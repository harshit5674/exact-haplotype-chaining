"""Regression cases for charged cache work, atomic warm-cache growth, and budget domains."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess

import gap_chaining as gc
import gap_native as gn


def collision_case(m):
    """C28: H=1, descending origins, a shared target, and colliding pair-cache keys.

    Source r restarts at 8m+4r. Target r uniquely chooses that source, at 8m+3r+1.
    The old auto counter admitted >10M uncharged probes at m=4096 without falling back.
    """
    stride = max(1, m // 2048)
    target = m * stride + 1
    path = [i * stride for i in range(1, m + 1)] + [target]
    g = gc.Graph([1] * (target + 1), list(zip(path, path[1:])), [path])
    a = []
    for r in range(m):
        a.extend([gc.Anchor((m-r)*stride, 1, 1, 2*r+1, 2*r+1, 8*m+4*r),
                  gc.Anchor(target, 1, 1, 2*r+2, 2*r+2, 1)])
    return g, a, 0


def run_budget_regressions():
    here = Path(__file__).resolve().parent
    compiler = os.environ.get('CXX', 'clang++')
    flags = ['-std=c++20', '-O2']
    identity = (here/'test_gap_budget.cpp').read_bytes() + (here/'gap_native.cpp').read_bytes() + repr((compiler, flags)).encode()
    exe = gn.ROOT/'results/native_build'/('gap-budget-test-'+hashlib.sha256(identity).hexdigest()[:16])
    exe.parent.mkdir(parents=True, exist_ok=True)
    if not exe.exists():
        temp = exe.with_name(exe.name+'.'+str(os.getpid())+'.tmp')
        subprocess.run([compiler, *flags, str(here/'test_gap_budget.cpp'), '-o', str(temp)], check=True)
        temp.replace(exe)
    unit = json.loads(subprocess.check_output([str(exe)], text=True))
    rows = []
    for m in (16, 128, 1024, 4096):
        g, a, gamma = collision_case(m)
        p = g.prepare(); G = gn.GapGraph(p); b = G.batch(a, gamma)
        ref = G.solve_packed(b, mode='per_haplotype')
        assert all(ref.C[2*r] == 8*m+4*r and ref.C[2*r+1] == 8*m+3*r+1 for r in range(m))
        if m == 16:
            C, ptr, B, R = ref.as_dicts(a)
            eC, eB, eR, _ = gc.direct_oracle(p, a, gamma)
            assert (C, B, R) == (eC, eB, eR)
        for warm in (False, True):
            if warm:
                assert G.solve_packed(b, mode='class', warm=True).arrays_digest() == ref.arrays_digest()
            for factor in (0.0, 0.25, 1.0, 4.0):
                z = G.solve_packed(b, mode='auto', warm=warm, auto_factor=factor)
                assert z.arrays_digest() == ref.arrays_digest(), (m, warm, factor)
                C, ptr, B, R = z.as_dicts(a)
                gc.verify_pointers(p, a, gamma, C, B, ptr)
                st = z.stats
                assert st['auto_work'] <= st['auto_budget']
                assert st['auto_cache_probes'] <= st['auto_work']
                if m == 4096 and not warm and factor == 1:
                    assert st['fallback'], 'collision family must trigger the guard'
                rows.append(dict(m=m, warm=warm, factor=factor, **{k:v for k,v in st.items()
                    if k.startswith('auto_') or k == 'fallback'}))
            # The same persistent cache must remain valid after partial runs.
            assert G.solve_packed(b, mode='class', warm=warm).arrays_digest() == ref.arrays_digest()
    p = gc.Graph([2], [], [[0], [0]]).prepare(); G = gn.GapGraph(p)
    a = [gc.Anchor(0, 1, 1, 1, 1, 1 << 40), gc.Anchor(0, 2, 2, 2, 2, 1)]
    z = G.solve_packed(a, 0, mode='auto', auto_factor=0)
    assert not z.stats['auto_fallback_available'] and math.isinf(z.stats['auto_budget'])
    assert not z.stats['fallback']
    C, ptr, B, R = z.as_dicts(a); eC, eB, eR, _ = gc.direct_oracle(p, a, 0)
    assert (C, B, R) == (eC, eB, eR)
    gc.verify_pointers(p, a, 0, C, B, ptr)
    small = [gc.Anchor(0, 1, 1, 1, 1, 1)]
    for factor in (-1.0, float('inf'), float('nan'), 1e308):
        try: G.solve_packed(small, 0, mode='auto', auto_factor=factor)
        except ValueError: pass
        else: raise AssertionError(('invalid factor accepted', factor))
    empty = G.solve_packed([], 0, mode='auto')
    assert not empty.C and empty.stats['auto_work'] == empty.stats['auto_budget'] == 0
    return dict(cache_unit=unit, collisions=rows, outside_fallback_domain='correct and explicitly unguarded',
                invalid_factors_rejected=4, empty='passed', status='all passed')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--out', default='verification_budget.json')
    args = ap.parse_args()
    report = run_budget_regressions()
    Path(args.out).write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k != 'collisions'}, indent=2))
