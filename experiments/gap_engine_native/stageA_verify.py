"""Stage A exactness gates for gap_fastpath.FastPathSolver.

Every gate runs two plans of the fast-path solver - 'prompt' (certified tests i-iii,
range floors/switches/minima, per-label remote queries) and 'full' (plus the grouped
remote step and the range floor cut) - and requires:
  * the same per-anchor state records, B, R and Fenwick node contents as the
    specification gap_chaining.OriginSolver (record-level identity, stronger than C/B/R)
    from both plans;
  * identical expanded C, B, R to an independent reference;
  * gap_chaining.verify_pointers on the fast-path pointers (legality, exact score,
    acyclic traceback) for every state.
Gates:
  random       test_gap_chaining.random_case (wide label counts every 4th seed) plus its
               7 fixed cases, against direct_oracle and scalar_fenwick;
  wide         127-1024 labels around word boundaries, against direct_oracle;
  independent  experiments/gap_solution_review/independent_check.py's generator
               (seed 777) and oracle;
  real         the 40 primary MHC reads at H=60, and mosaic panels at H=960, against
               real_check_mosaic.py's independent graph-major per-haplotype `mine`.
A subset of the random gate also runs the ablations tests='none'|'i'|'ii'.
Failures are recorded with their inputs before the script exits nonzero.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from multiprocessing import get_context
from pathlib import Path
import random
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import common as cm  # noqa: E402
import gap_chaining as gc  # noqa: E402
import gap_fastpath as fp  # noqa: E402
from test_gap_chaining import random_case, fixed_cases  # noqa: E402


def node_sets(solver) -> list:
    return [sorted((r.coef, r.origin, r.mask, r.witness) for r in node) for node in solver.nodes]


PLANS = dict(prompt=dict(tests='iii', remote='labels', cut='none'),
             full=dict(tests='iii', remote='grouped', cut='range'))


def spec_and_fast(p, a, gamma, plans=('prompt', 'full'), **cfg):
    """Spec result, the last fast result, and whether every fast run is record-identical."""
    spec = gc.OriginSolver(p, a, gamma)
    rs = spec.run()
    target = node_sets(spec)
    same, rf = True, None
    for name in plans:
        fast = fp.FastPathSolver(p, a, gamma, **(PLANS[name] if name else cfg))
        rf = fast.run()
        same &= (rf.records == rs.records and rf.B == rs.B and rf.R == rs.R
                 and node_sets(fast) == target)
    return rs, rf, same


def digest(C, B, R) -> str:
    return hashlib.sha256(json.dumps([[sorted(r.items()) for r in C], B, R]).encode()).hexdigest()


def failure(kind, name, case, detail):
    g, a, gamma = case
    return dict(kind=kind, name=name, detail=str(detail)[:2000], lengths=g.lengths, edges=g.edges,
                paths=g.paths, anchors=[list(vars(z).values()) for z in a], gamma=str(gamma))


def check_small(name, case, ablate, report, oracle='direct'):
    g, a, gamma = case
    p = g.prepare()
    rs, rf, same = spec_and_fast(p, a, gamma)
    if not same:
        report['failures'].append(failure('records', name, case, 'record mismatch with spec'))
        return
    if ablate:
        for tests in ('none', 'i', 'ii'):
            if not spec_and_fast(p, a, gamma, plans=(None,), tests=tests)[2]:
                report['failures'].append(failure('records_' + tests, name, case, 'ablation mismatch'))
                return
        report['ablation_cases'] += 1
    C, ptr = rf.expand()
    eC, eB, eR, eptr = gc.direct_oracle(p, a, gamma)
    if (C, rf.B, rf.R) != (eC, eB, eR):
        report['failures'].append(failure('direct_oracle', name, case, 'C/B/R mismatch'))
        return
    if oracle == 'direct+scalar':
        sC, sB, sR, _ = gc.scalar_fenwick(p, a, gamma)
        if (sC, sB, sR) != (eC, eB, eR):
            report['failures'].append(failure('scalar_fenwick', name, case, 'reference disagreement'))
            return
    counts = gc.verify_pointers(p, a, gamma, C, rf.B, ptr)
    for k, v in counts.items():
        report['pointer_counts'][k] += v
    report['cases'] += 1
    report['states'] += sum(map(len, C))
    report['max_H'] = max(report['max_H'], p.H)
    for k in ('sliced_pairs', 'sliced_thresholds', 'sliced_minima', 'pairs_level3', 'pairs_level4',
              'remote_single', 'remote_sliced', 'remote_pruned', 'cut_records'):
        report['fast_counters'][k] = report['fast_counters'].get(k, 0) + rf.cost[k]
    report['hashes'].append([name, digest(C, rf.B, rf.R)])


def new_report():
    return dict(cases=0, states=0, max_H=0, ablation_cases=0, failures=[], fast_counters={},
                pointer_counts=dict(restart=0, local=0, remote=0, same=0), hashes=[])


def gate_random(count, ablate_every):
    report = new_report()
    for name, case in fixed_cases():
        check_small(name, case, True, report, 'direct+scalar')
    for seed in range(count):
        check_small(f'random_{seed}', random_case(seed, wide=seed % 4 == 0),
                    seed % ablate_every == 0, report, 'direct+scalar')
    return report


def wide_case(seed):
    """Few anchors, many labels: 127..1024 haplotypes around 64-bit word boundaries."""
    rng = random.Random(10_000 + seed)
    H = rng.choice([96, 127, 128, 129, 191, 192, 255, 256, 257, 383, 384, 511, 512, 513, 767, 960, 1000, 1023, 1024])
    V = rng.randint(2, 9)
    lengths = [rng.randint(1, 9) for _ in range(V)]
    # A few haplotype "families" with private indels, so origins and route lengths vary.
    families = [[0] + [v for v in range(1, V - 1) if rng.random() < rng.uniform(.2, .9)] + [V - 1]
                for _ in range(rng.randint(1, 6))]
    paths = []
    for h in range(H):
        base = rng.choice(families)
        if rng.random() < .15:
            base = [0] + [v for v in range(1, V - 1) if rng.random() < .5] + [V - 1]
        paths.append(list(base))
    edges = {(u, v) for path in paths for u, v in zip(path, path[1:])}
    for u in range(V):
        for v in range(u + 1, V):
            if rng.random() < .1:
                edges.add((u, v))
    graph = gc.Graph(lengths, sorted(edges), paths)
    vertices = sorted({v for path in paths for v in path})
    anchors = []
    for _ in range(rng.randint(0, 22)):
        v = rng.choice(vertices)
        x = rng.randint(1, lengths[v]); y = rng.randint(x, lengths[v])
        c = rng.randint(1, 40); d = c + rng.randint(0, 5)
        anchors.append(gc.Anchor(v, x, y, c, d, rng.choice([0, 1, 7, 200, 3400]) * rng.randint(1, 4)))
    return graph, anchors, rng.choice([0, 1, 5, 100, 10000])


def gate_wide(count):
    report = new_report()
    for seed in range(count):
        check_small(f'wide_{seed}', wide_case(seed), seed % 5 == 0, report)
    return report


def gate_independent(count):
    ns = cm.independent_oracle()
    rng = random.Random(777)
    report = new_report()
    for k in range(count):
        lengths, edges, paths, A, gamma = ns['gen'](rng)
        case = (gc.Graph(lengths, edges, paths), A, gamma)
        p = case[0].prepare()
        rs, rf, same = spec_and_fast(p, A, gamma)
        if not same:
            report['failures'].append(failure('records', f'independent_{k}', case, 'record mismatch'))
            continue
        C, ptr = rf.expand()
        eC, eB, eR = ns['oracle'](lengths, edges, paths, A, gamma)
        if (C, rf.B, rf.R) != (eC, eB, eR):
            report['failures'].append(failure('independent_oracle', f'independent_{k}', case, 'C/B/R mismatch'))
            continue
        counts = gc.verify_pointers(p, A, gamma, C, rf.B, ptr)
        for key, v in counts.items():
            report['pointer_counts'][key] += v
        report['cases'] += 1
        report['states'] += sum(map(len, C))
        report['max_H'] = max(report['max_H'], p.H)
    report.pop('hashes')
    return report


_NS, _PREP = {}, {}


def real_one(job):
    H, idx, gamma = job
    if H not in _NS:
        _NS[H] = cm.independent_mine(H)
        _PREP[H] = cm.prepared(H)
    ns, p = _NS[H], _PREP[H]
    case = cm.primary_cases()[idx]
    A = cm.anchors(case)
    t = time.perf_counter()
    rs, rf, same = spec_and_fast(p, A, gamma)
    solve_seconds = time.perf_counter() - t
    row = dict(H=H, id=case['id'], N=len(A), gamma=gamma, record_identity=same)
    if same:
        C, ptr = rf.expand()
        ns['GAMMA'] = gamma
        mC, mB, mR = ns['mine'](A)
        row['independent_agree'] = (C == mC and rf.B == mB and rf.R == mR)
        row['pointer_counts'] = gc.verify_pointers(p, A, gamma, C, rf.B, ptr)
        row['states'] = sum(map(len, C))
        row['digest'] = digest(C, rf.B, rf.R)
        del C, ptr, mC
    row['solve_seconds'] = round(solve_seconds, 2)
    row['passed'] = bool(same and row.get('independent_agree'))
    return row


def gate_real(H, count, gamma, workers):
    n = len(cm.primary_cases())
    order = sorted(range(n), key=lambda k: len(cm.primary_cases()[k]['anchors']))
    picks = order if count >= n else order[::max(1, n // count)][:count]
    jobs = [(H, k, gamma) for k in picks]
    if workers > 1:
        with get_context('spawn').Pool(workers) as pool:
            rows = pool.map(real_one, jobs, chunksize=1)
    else:
        rows = [real_one(j) for j in jobs]
    for r in rows:
        print(json.dumps(r), flush=True)
    return dict(H=H, gamma=gamma, reads=len(rows), passed=sum(r['passed'] for r in rows),
                states=sum(r.get('states', 0) for r in rows),
                pointer_counts={k: sum(r.get('pointer_counts', {}).get(k, 0) for r in rows)
                                for k in ('restart', 'local', 'remote', 'same')},
                rows=rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--random', type=int, default=3000)
    ap.add_argument('--wide', type=int, default=400)
    ap.add_argument('--independent', type=int, default=3000)
    ap.add_argument('--real60', type=int, default=40)
    ap.add_argument('--real960', type=int, default=6)
    ap.add_argument('--gamma20000', type=int, default=6, help='H=60 reads rerun with gamma 20000')
    ap.add_argument('--workers', type=int, default=4)
    args = ap.parse_args()
    if args.out.exists():
        raise SystemExit('refusing to overwrite ' + str(args.out))
    start = time.perf_counter()
    out = dict(script='experiments/gap_engine_native/stageA_verify.py', argv=sys.argv[1:],
               python=sys.version, gates={})
    out['gates']['random'] = gate_random(args.random, 10)
    print('random', {k: v for k, v in out['gates']['random'].items() if k != 'hashes'}, flush=True)
    out['gates']['wide'] = gate_wide(args.wide)
    print('wide', {k: v for k, v in out['gates']['wide'].items() if k != 'hashes'}, flush=True)
    out['gates']['independent'] = gate_independent(args.independent)
    print('independent', out['gates']['independent'], flush=True)
    if args.real60:
        out['gates']['real_H60'] = gate_real(60, args.real60, cm.GAMMA, args.workers)
    if args.gamma20000:
        out['gates']['real_H60_gamma20000'] = gate_real(60, args.gamma20000, 20000, args.workers)
    if args.real960:
        out['gates']['real_H960'] = gate_real(960, args.real960, cm.GAMMA, min(2, args.workers))
    fails = sum(len(out['gates'][k]['failures']) for k in ('random', 'wide', 'independent'))
    fails += sum(g['reads'] - g['passed'] for k, g in out['gates'].items() if k.startswith('real'))
    out['failures_total'] = fails
    out['status'] = 'all passed' if not fails else 'FAILED'
    out['seconds'] = round(time.perf_counter() - start, 1)
    args.out.write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(dict(status=out['status'], seconds=out['seconds'])), flush=True)
    return 0 if not fails else 1


if __name__ == '__main__':
    sys.exit(main())
