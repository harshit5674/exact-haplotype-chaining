"""Mutation test: does the adversarial harness detect deliberately broken solvers?

For each mutant in mutants.py, run it on a deterministic subset of each family and count
instances where its C/B/R differ from the brute force or either pointer check fails.
Non-equivalent mutants must be detected; equivalent ones must produce no report (no false
alarms). Writes results/mutation.json. Usage:
  python3.13 -B experiments/proof_review/mutation_test.py [--per-family 3000]
"""
from __future__ import annotations
import argparse
import collections
import itertools
import json
import multiprocessing as mp
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def subset(family, k):
    import families as fam
    if family == 'EXA':
        it = itertools.chain.from_iterable(fam.exa_instances(V, 3, fam.ANCHOR_SETS_A, 0, 7) for V in (2, 3))
    elif family == 'EXB2':
        it = itertools.chain.from_iterable(fam.exb2_instances(n, 0, 37) for n in fam.HOSTILE_GRAPHS)
    elif family == 'EXB3':
        it = itertools.chain.from_iterable(fam.exb3_instances(n, 0, 23) for n in fam.HOSTILE_GRAPHS_3)
    else:
        it = (fam.random_stress_instance(s) for s in range(k))
    return list(itertools.islice(it, k))


def run(job):
    import brute_chains as bc
    import mutants as mu
    import gap_chaining as gc
    mutant, family, k = job
    det = collections.Counter()
    example = None
    facts_cache = {}
    for name, (L, E, P), A, gamma in subset(family, k):
        key = (tuple(L), tuple(map(tuple, E)), tuple(map(tuple, P)))
        if key not in facts_cache:
            facts_cache.clear()
            facts_cache[key] = bc.Facts(L, E, P)
        facts = facts_cache[key]
        C, Cnot, B, R = bc.chain_dag_optimum(facts, A, gamma)
        det['instances'] += 1
        try:
            p, anchors, mC, mB, mR, mptr = mu.mutant_outputs(L, E, P, A, gamma, mutant)
        except Exception as e:
            det['crash'] += 1
            example = example or dict(name=name, why='crash ' + repr(e)[:200])
            continue
        fC = {(i, h): v for i, row in enumerate(mC) for h, v in row.items()}
        why = []
        if fC != C:
            why.append('C')
        if list(mB) != B:
            why.append('B')
        if list(mR) != R:
            why.append('R')
        try:
            gc.verify_pointers(p, anchors, gamma, mC, mB, mptr)
        except AssertionError:
            why.append('verify_pointers')
        fptr = {}
        for i, row in enumerate(mptr):
            for h, w in row.items():
                fptr[(i, h)] = (w.kind, w.j, w.g if w.kind in ('remote', 'same') else -1,
                                w.w if w.kind == 'remote' else -1)
        if bc.trace_and_score(facts, A, gamma, fC, fptr):
            why.append('traceback')
        if why:
            det['detected'] += 1
            for wh in why:
                det['by_' + wh] += 1
            if example is None:
                example = dict(name=list(map(str, name)), lengths=L, edges=E, paths=P, anchors=A,
                               gamma=gamma, why=why)
    return mutant, family, dict(det), example


def main():
    import mutants as mu
    ap = argparse.ArgumentParser()
    ap.add_argument('--per-family', type=int, default=3000)
    ap.add_argument('--workers', type=int, default=3)
    args = ap.parse_args()
    assert args.workers <= 3 and sys.flags.optimize == 0
    fams = ('EXA', 'EXB2', 'EXB3', 'RS')
    jobs = [(m, f, args.per_family) for m in mu.NON_EQUIVALENT + mu.EQUIVALENT for f in fams]
    with mp.get_context('spawn').Pool(args.workers) as pool:
        res = pool.map(run, jobs, chunksize=1)
    out = collections.defaultdict(dict)
    for m, f, det, ex in res:
        out[m][f] = dict(det, example=ex)
    summary = {}
    for m in mu.NON_EQUIVALENT + mu.EQUIVALENT:
        tot = sum(out[m][f].get('detected', 0) for f in fams)
        summary[m] = dict(equivalent=m in mu.EQUIVALENT, detected_instances=tot,
                          verdict=('OK' if (tot > 0) != (m in mu.EQUIVALENT) else 'PROBLEM'))
    report = dict(per_family=args.per_family, summary=summary, detail=out)
    (HERE / 'results').mkdir(exist_ok=True)
    (HERE / 'results/mutation.json').write_text(json.dumps(report, indent=1, default=str) + '\n')
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
