#!/usr/bin/env python3
"""Deterministic differential, arithmetic, indexing, and traceback tests."""
import argparse, json, random, time, hashlib
from dataclasses import asdict
from pathlib import Path
from gap_chaining import *


def random_case(seed, wide=False):
    rng = random.Random(seed)
    V = rng.randint(2, 15)
    H = rng.choice([1, 2, 3, 5, 8, 16, 31, 63, 64, 65, 70]) if wide else rng.randint(1, 10)
    lengths = [rng.randint(1, 12) for _ in range(V)]
    paths = []
    for h in range(H):
        prob = rng.uniform(.15, .9)
        paths.append([0]+[v for v in range(1,V-1) if rng.random() < prob]+[V-1])
    edges = {(u,v) for path in paths for u,v in zip(path,path[1:])}
    for u in range(V):
        for v in range(u+1,V):
            if rng.random()<.07:
                edges.add((u,v))
    graph = Graph(lengths, sorted(edges), paths)
    vertices = sorted(set(v for path in paths for v in path))
    N = rng.randint(0, 35)
    anchors = []
    for i in range(N):
        v = rng.choice(vertices)
        x = rng.randint(1,lengths[v]); y = rng.randint(x,lengths[v])
        c = rng.randint(1,60); d = c+rng.randint(0,7)
        weight = rng.choice([0,1,2,10,200,1000]) * rng.randint(1,5)
        if seed % 29 == 0:
            weight += (1 << 250) * rng.randint(1,7)
        anchors.append(Anchor(v,x,y,c,d,weight))
    gamma = rng.choice([0,1,2,10,100,10000,1<<270])
    return graph, anchors, gamma


def check_index(p):
    # Independent all-pairs DAG shortest paths.
    n = p.V
    dd = [[None]*n for _ in range(n)]
    for s in range(n):
        dd[s][s] = 0
        for v in range(s+1,n):
            values = [dd[s][u]+p.graph.lengths[u] for u in p.pred[v] if dd[s][u] is not None]
            if values:
                dd[s][v] = min(values)
    for v in range(n):
        for h,path in enumerate(p.graph.paths):
            ws = [w for w in path if dd[w][v] is not None]
            w = max(ws, default=-1)
            assert p.last[v][h] == w
            assert p.distance[v][h] == (None if w < 0 else dd[w][v])


def remote_only_case():
    # All graph edges lie on a source-to-sink haplotype. Target labels have no
    # published anchors, but a third (remote) label supplies a better chain.
    g = Graph([1]*5, [(0,1),(1,2),(2,4),(0,3),(3,4),(0,2),(2,3)],
              [[0,1,2,4],[0,3,4],[0,2,3,4]])
    return g, [Anchor(1,1,1,1,1,200),Anchor(3,1,1,2,2,200)], 10


def fixed_cases():
    yield 'remote_only', remote_only_case()
    # Zero penalty, length-variant routes, exact coordinate ties, very large scores.
    g = Graph([2,1,7,2], [(0,1),(0,2),(1,3),(2,3)], [[0,1,3],[0,2,3]])
    for gamma in [0,1,5,10000]:
        a = [Anchor(0,1,1,1,1,200),Anchor(3,1,1,10,10,200),
             Anchor(3,1,2,10,11,400),Anchor(3,2,2,11,11,200)]
        yield f'diamond_{gamma}', (g,a,gamma)
    # All 130 labels share one path, so an anchor class may have >64 labels.
    g = Graph([12], [], [[0] for _ in range(130)])
    a = [Anchor(0,i,i,i,i,200) for i in range(1,13)]
    yield '130_identical', (g,a,0)
    g = Graph([1,1000000000000000000000000000000,1,1],
              [(0,1),(0,2),(1,3),(2,3)], [[0,1,3],[0,2,3]])
    a = [Anchor(0,1,1,1,1,10**100),Anchor(3,1,1,3,3,10**100)]
    yield 'large_coordinates_and_scores', (g,a,10**101)


def arithmetic_tests(count=1000):
    rng = random.Random(12345)
    checks = 0
    for _ in range(count):
        H = rng.choice([1,2,7,63,64,65,129])
        maxoff = rng.randint(1,10000)
        # Arbitrary nonnegative coordinate arrays are tested directly. This
        # checks the numeric primitive on a superset of realizable DAG data.
        fake = Graph([maxoff+1], [], [[0] for h in range(H)]).prepare()
        anchors = [Anchor(0,1,1,1,1,1) for i in range(6)]
        solver = OriginSolver(fake, anchors, 0)
        vals = [[rng.randint(0,maxoff) for h in range(H)] for v in range(6)]
        fake.phibits = [[sum(1<<h for h in range(H) if (row[h]>>b)&1)
                         for b in range(fake.b)] for row in vals]
        fake.phi = vals
        # Map each origin to a separate synthetic coordinate row.
        solver.a = [Anchor(v,1,1,1,1,1) for v in range(6)]
        left, right = [], []
        for dest in [left,right]:
            buckets = {}
            for h in range(H):
                if rng.random()<.15:
                    continue
                idx = rng.randrange(6)
                buckets[idx] = buckets.get(idx,0)|(1<<h)
            for idx,mask in buckets.items():
                coef = rng.choice([rng.randint(-20000,20000),10**100+rng.randint(-10,10),-10**100])
                dest.append(Record(coef,idx,mask,Witness('local',idx)))
        out = solver.merge(left,right)
        for h in range(H):
            candidates = [r.coef+vals[r.origin][h] for r in left+right if r.mask>>h&1]
            got = [r.coef+vals[r.origin][h] for r in out if r.mask>>h&1]
            assert got == ([] if not candidates else [max(candidates)])
            checks += 1
    return checks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cases',type=int,default=1200)
    ap.add_argument('--arithmetic',type=int,default=1500)
    ap.add_argument('--out',default='verification.json')
    args = ap.parse_args()
    start=time.perf_counter()
    report = dict(random_cases=0,fixed_cases=0,states=0,pointer_counts=dict(restart=0,local=0,remote=0,same=0),
                  max_H=0,total_lambda=0,total_coded_records=0,total_route_calls=0,
                  arithmetic_lane_checks=arithmetic_tests(args.arithmetic),case_hashes=[])
    def run_one(name, case):
        g,a,gamma = case
        p = g.prepare(); check_index(p)
        expected,B,R,op = direct_oracle(p,a,gamma)
        dense,db,dr,dp = scalar_fenwick(p,a,gamma)
        assert (dense,db,dr)==(expected,B,R), name
        compact = solve(p,a,gamma)
        C,ptr = compact.expand()
        assert (C,compact.B,compact.R)==(expected,B,R), (name,C,expected,compact.R,R)
        counts = verify_pointers(p,a,gamma,C,compact.B,ptr)
        verify_pointers(p,a,gamma,dense,db,dp)
        verify_pointers(p,a,gamma,expected,B,op)
        for k,v in counts.items(): report['pointer_counts'][k]+=v
        report['states'] += sum(map(len,C))
        report['max_H'] = max(report['max_H'],p.H)
        report['total_lambda'] += compact.stats['lambda_records']
        report['total_coded_records'] += compact.stats['coded_records']
        report['total_route_calls'] += compact.stats['route_calls']
        digest = hashlib.sha256(json.dumps([C,B,R],sort_keys=True).encode()).hexdigest()
        report['case_hashes'].append([name,digest])
    for name,case in fixed_cases():
        run_one(name,case); report['fixed_cases']+=1
    for seed in range(args.cases):
        run_one(f'random_{seed}',random_case(seed,wide=seed%4==0));report['random_cases']+=1
    report['seconds'] = time.perf_counter()-start
    report['status'] = 'all passed'
    Path(args.out).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='case_hashes'},indent=2))

if __name__ == '__main__': main()
