"""Run the adversarial families through harness.check_instance with at most 3 workers.

Usage (from the repository root; never with -O):
  python3.13 -B experiments/proof_review/run_families.py --family EXA --workers 3
  families: EXA (V<=3, H<=3), EXA4 (V=4, H<=2), EXB2, EXB3, RS (--rs N seeds), ALL
Writes experiments/proof_review/results/<family>.json with counts, solver-check counters and
every failure (full input, solver, kind, detail). Exits nonzero if any failure occurred.
"""
from __future__ import annotations
import argparse
import collections
import json
import multiprocessing as mp
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
OUT = HERE / 'results'

MAX_FAILURES_KEPT = 300


def make_iter(job):
    import families as fam
    kind = job['family']
    s, ns = job['shard'], job['nshards']
    if kind == 'EXA':
        for V in (1, 2, 3):
            yield from fam.exa_instances(V, 3, fam.ANCHOR_SETS_A, s, ns)
    elif kind == 'EXA4':
        yield from fam.exa_instances(4, 2, fam.ANCHOR_SETS_A4, s, ns)
    elif kind == 'EXB2':
        for name in fam.HOSTILE_GRAPHS:
            yield from fam.exb2_instances(name, s, ns)
    elif kind == 'EXB3':
        for name in fam.HOSTILE_GRAPHS_3:
            yield from fam.exb3_instances(name, s, ns)
    elif kind == 'RS':
        for seed in range(job['rs_start'] + s, job['rs_start'] + job['rs'], ns):
            yield fam.random_stress_instance(seed)
    else:
        raise ValueError(kind)


def worker(job):
    import harness as H
    stats = collections.Counter()
    failures = []
    nfail = 0
    cache_key, gcache = None, None
    t0 = time.perf_counter()
    for name, graph, anchors, gamma in make_iter(job):
        lengths, edges, paths = graph
        key = (tuple(lengths), tuple(map(tuple, edges)), tuple(map(tuple, paths)))
        try:
            if key != cache_key:
                gcache = H.GraphCache(lengths, edges, paths)
                cache_key = key
            f = H.check_instance(gcache, anchors, gamma, stats)
        except Exception as e:  # a crash is a failure with its input
            f = [dict(solver='harness', kind='exception', detail=traceback.format_exc()[-2000:])]
            cache_key = None
        if f:
            nfail += 1
            if len(failures) < MAX_FAILURES_KEPT:
                failures.append(H.failure_record(lengths, edges, paths, anchors, gamma, f, list(map(str, name))))
    stats['seconds'] = time.perf_counter() - t0
    return dict(job=job, stats=dict(stats), failing_instances=nfail, failures=failures)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--family', required=True)
    ap.add_argument('--workers', type=int, default=3)
    ap.add_argument('--shards', type=int, default=0, help='default: 4 x workers')
    ap.add_argument('--rs', type=int, default=20000)
    ap.add_argument('--rs-start', type=int, default=0)
    args = ap.parse_args()
    assert args.workers <= 3, 'shared machine: at most 3 workers'
    assert sys.flags.optimize == 0, 'assertions must stay enabled'
    fams = ['EXA', 'EXA4', 'EXB2', 'EXB3', 'RS'] if args.family == 'ALL' else [args.family]
    OUT.mkdir(exist_ok=True)
    any_fail = False
    for family in fams:
        nshards = args.shards or 4 * args.workers
        jobs = [dict(family=family, shard=s, nshards=nshards, rs=args.rs, rs_start=args.rs_start)
                for s in range(nshards)]
        t = time.perf_counter()
        ctx = mp.get_context('spawn')
        with ctx.Pool(args.workers) as pool:
            parts = pool.map(worker, jobs, chunksize=1)
        stats = collections.Counter()
        failures, nfail = [], 0
        for part in parts:
            stats.update({k: v for k, v in part['stats'].items() if k != 'seconds'})
            nfail += part['failing_instances']
            failures.extend(part['failures'])
        report = dict(family=family, wall_seconds=round(time.perf_counter() - t, 1),
                      worker_seconds=round(sum(p['stats'].get('seconds', 0) for p in parts), 1),
                      workers=args.workers, shards=nshards,
                      rs_seeds=[args.rs_start, args.rs_start + args.rs] if family == 'RS' else None,
                      failing_instances=nfail, stats=dict(sorted(stats.items())),
                      failures=failures[:MAX_FAILURES_KEPT])
        name = family if family != 'RS' else 'RS_%d_%d' % (args.rs_start, args.rs_start + args.rs)
        (OUT / (name + '.json')).write_text(json.dumps(report, indent=1, default=str) + '\n')
        print(json.dumps({k: v for k, v in report.items() if k != 'failures'}, default=str))
        any_fail |= nfail > 0
    sys.exit(1 if any_fail else 0)


if __name__ == '__main__':
    main()
