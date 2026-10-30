"""Greedy minimizer for a failing instance recorded by run_families.py / hostile.py.

Usage: python3.13 -B minimize.py results/<file>.json [index] [--native] [--relaxed]
Takes failure record `index` (default 0) and shrinks it while harness.check_instance still
reports a failure of the same (solver, kind): drop anchors, drop haplotypes, drop edges not on
haplotypes, zero/halve weights, shift query coordinates down, lower gamma. Writes
results/minimized_<source>_<index>.json with the original and minimized inputs.
"""
from __future__ import annotations
import collections
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import harness as H  # noqa: E402


def failing_kinds(inst, native, relaxed):
    L, E, P, A, gamma = inst
    try:
        prepared_cls = None
        if relaxed:
            import domain_tests
            prepared_cls = domain_tests.RelaxedPrepared
        g = H.GraphCache(L, E, P, prepared_cls=prepared_cls)
        f = H.check_instance(g, A, gamma, collections.Counter(), native=native)
    except Exception as e:
        return {('exception', type(e).__name__)}
    return {(x['solver'], x['kind']) for x in f}


def candidates(inst):
    L, E, P, A, gamma = inst
    for k in range(len(A)):
        yield (L, E, P, A[:k] + A[k + 1:], gamma)
    for h in range(len(P)):
        P2 = P[:h] + P[h + 1:]
        covered = {v for p in P2 for v in p}
        if P2 and all(a[0] in covered for a in A):
            yield (L, E, P2, A, gamma)
    on_path = {(u, v) for p in P for u, v in zip(p, p[1:])}
    for e in E:
        if tuple(e) not in on_path:
            yield (L, [x for x in E if x != e], P, A, gamma)
    for k, a in enumerate(A):
        v, x, y, c, d, w = a
        for w2 in (0, w // 2):
            if w2 != w:
                yield (L, E, P, A[:k] + [(v, x, y, c, d, w2)] + A[k + 1:], gamma)
        if c > 1:
            yield (L, E, P, A[:k] + [(v, x, y, c - 1, d - 1, w)] + A[k + 1:], gamma)
    for g2 in (0, 1, gamma // 2):
        if g2 < gamma:
            yield (L, E, P, A, g2)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    native = '--native' in sys.argv
    relaxed = '--relaxed' in sys.argv
    src = Path(args[0])
    idx = int(args[1]) if len(args) > 1 else 0
    data = json.loads(src.read_text())
    recs = data.get('failures') or [f for g in data.get('generators', {}).values() for f in g['failures']]
    rec = recs[idx]
    inst = (rec['lengths'], [tuple(e) for e in rec['edges']], rec['paths'],
            [tuple(a) for a in rec['anchors']], int(rec['gamma']))
    target = failing_kinds(inst, native, relaxed)
    assert target, 'record does not fail any more'
    changed = True
    while changed:
        changed = False
        for cand in candidates(inst):
            if failing_kinds(cand, native, relaxed) & target:
                inst, changed = cand, True
                break
    out = dict(source=str(src), index=idx, failure_kinds=sorted(map(list, target)),
               original=rec, minimized=dict(lengths=inst[0], edges=inst[1], paths=inst[2],
                                            anchors=inst[3], gamma=inst[4]))
    path = HERE / 'results' / ('minimized_%s_%d.json' % (src.stem, idx))
    path.write_text(json.dumps(out, indent=1, default=str) + '\n')
    print(json.dumps(out['minimized'], default=str))


if __name__ == '__main__':
    main()
