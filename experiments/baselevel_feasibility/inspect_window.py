"""Extract a base-level MHC window from HPRC release 2 and test whether it is a faithful DAG (2026-10-20).

Feasibility spike. For a CHM13 chr6 interval, gbz-base (pinned, vendor/) extracts
the subgraph of the release-2 Minigraph-Cactus graph v2.1 (hprc-v2.1-mc-chm13.gbz.db) with all haplotype
walks, extended to every top-level snarl overlapping the interval (--extend-snarls; with --snarls or no
extension, alternative alleles outside the reference path are missing and every haplotype breaks into
hundreds of fragments, as the first 250 kb test showed). The script then orients every walk as CHM13 runs (reversing a
walk that runs against it, as the MICB pilot did), and checks, without removing anything:
  - which oriented segments (handles) the walks use, and whether any segment is used in both orientations;
  - whether a walk visits a handle twice (a repeat, which no DAG path can do);
  - whether the union of walk steps and GFA links between used handles is acyclic;
  - node-length distribution, walk counts and spans.
Nothing is pruned to obtain a DAG; a failure is reported as a failure.

  python3.13 -B experiments/baselevel_feasibility/inspect_window.py --start 31200000 --end 31450000 \
      --out experiments/baselevel_feasibility/windows/w31200k
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import heapq
import json
from pathlib import Path
import re
import statistics
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DB = HERE / 'data/hprc-v2.1-mc-chm13.gbz.db'
GBZ_BASE = HERE / 'vendor/gbz-base/target/release/gbz-base'
RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def flip(h):
    return h[0], not h[1]


def parse(gfa: Path):
    seq, links, walks = {}, [], []
    for line in gfa.open():
        a = line.rstrip('\n').split('\t')
        if a[0] == 'S':
            seq[int(a[1])] = a[2]
        elif a[0] == 'L':
            links.append(((int(a[1]), a[2] == '-'), (int(a[3]), a[4] == '-')))
        elif a[0] == 'W':
            steps = [(int(n), o == '<') for o, n in re.findall(r'([><])(\d+)', a[6])]
            walks.append(dict(sample=a[1], hap=a[2], contig=a[3], start=a[4], end=a[5], steps=steps))
    return seq, links, walks


def toposort(handles, adj):
    indeg = dict.fromkeys(handles, 0)
    for u in handles:
        for v in adj[u]:
            indeg[v] += 1
    queue = [u for u in handles if not indeg[u]]
    heapq.heapify(queue)
    order = []
    while queue:
        u = heapq.heappop(queue)
        order.append(u)
        for v in sorted(adj[u]):
            indeg[v] -= 1
            if not indeg[v]:
                heapq.heappush(queue, v)
    return order, [u for u in handles if indeg[u] > 0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--start', type=int, required=True, help='CHM13 chr6 start, 0-based')
    ap.add_argument('--end', type=int, required=True, help='CHM13 chr6 end, half-open')
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--limit', default='5M')
    args = ap.parse_args()
    args.out.mkdir(parents=True)
    gfa = args.out / 'window.gfa'
    cmd = [str(GBZ_BASE), 'query', '--sample', 'CHM13', '--contig', 'chr6', '--interval', f'{args.start}..{args.end}',
           '--context', '0', '--extend-snarls', '--limit', args.limit, str(DB)]
    with gfa.open('w') as fh:
        run = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True)
    (args.out / 'query.log').write_text(run.stderr)
    assert run.returncode == 0, run.stderr[-2000:]
    seq, links, walks = parse(gfa)
    adj = defaultdict(set)
    for u, v in links:
        adj[u].add(v)
        adj[flip(v)].add(flip(u))
    ref = [w for w in walks if w['sample'] == 'CHM13']
    assert len(ref) == 1, [w['sample'] for w in ref]
    ref_steps = ref[0]['steps']
    ref_pos = {h: i for i, h in enumerate(ref_steps)}
    oriented, report = [], dict(walks=len(walks), reversed=0, undecided=0, repeated_handle_walks=0)
    for w in walks:
        steps = w['steps']
        fwd = [ref_pos[h] for h in steps if h in ref_pos]
        rev = [ref_pos[flip(h)] for h in steps if flip(h) in ref_pos]
        if len(rev) > len(fwd) or (len(rev) == len(fwd) and rev and rev[0] > rev[-1]):
            steps = [flip(h) for h in reversed(steps)]
            report['reversed'] += 1
        elif not fwd and not rev:
            report['undecided'] += 1
        if len(set(steps)) != len(steps):
            report['repeated_handle_walks'] += 1
        oriented.append(steps)
    handles = {h for s in oriented for h in s}
    segs_both = {n for n, r in handles if (n, not r) in handles}
    edges = set()
    for s in oriented:
        for u, v in zip(s, s[1:]):
            assert v in adj[u], ('walk step without a GFA link', u, v)
            edges.add((u, v))
    for u in handles:
        for v in adj[u]:
            if v in handles:
                edges.add((u, v))
    dag_adj = defaultdict(set)
    for u, v in edges:
        dag_adj[u].add(v)
    order, stuck = toposort(handles, dag_adj)
    lengths = [len(seq[n]) for n, _ in handles]
    span_full = sum(1 for s in oriented if s and s[0] == ref_steps[0] and s[-1] == ref_steps[-1])
    samples = Counter((w['sample'], w['hap']) for w in walks)
    report.update(
        command=cmd, gbz_base=GBZ_BASE.name, db=DB.name, start=args.start, end=args.end,
        gfa_sha256=hashlib.sha256(gfa.read_bytes()).hexdigest(), gfa_bytes=gfa.stat().st_size,
        segments=len(seq), links=len(links), handles=len(handles), segments_in_both_orientations=len(segs_both),
        edges_used=len(edges), acyclic=not stuck, handles_in_cycles_or_after=len(stuck),
        walks_spanning_window=span_full, distinct_sample_hap=len(samples),
        sample_hap_with_several_walks=sum(1 for c in samples.values() if c > 1),
        unknown_named_walks=sum(1 for w in walks if 'unknown' in w['sample'].lower()),
        reference_steps=len(ref_steps), reference_bases=sum(len(seq[n]) for n, _ in ref_steps),
        node_length=dict(median=statistics.median(lengths), mean=statistics.fmean(lengths),
                         one_base=sum(1 for x in lengths if x == 1), lt17=sum(1 for x in lengths if x < 17),
                         ge17_bases=sum(x for x in lengths if x >= 17), total_bases=sum(lengths)),
        walk_steps=sum(len(s) for s in oriented))
    (args.out / 'inspect.json').write_text(json.dumps(report, indent=1) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'command'}, indent=1))


if __name__ == '__main__':
    main()
