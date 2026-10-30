"""Anchors for each HPRC MHC panel, exported by the pinned Minichain and converted exactly.

For data/<panel>/graph.gfa.gz this runs the pinned Minichain v1.3 build with the research
seed export (experiments/mhc_pilot/minichain-export.patch: it dumps every seed hit after
collection and skips chaining) on the same 100 sampled CHM13 PacBio reads as the pilot
(experiments/mhc_pilot/data/sample100.fq). It then converts seeds to anchors with the
pilot's rules (experiments/mhc_pilot/prepare.py `convert`): 1-based closed coordinates,
reverse-strand hits mapped to the forward graph by reverse-complementing both sides,
weight = seed span, and an exact sequence check of every anchor against the graph and read.
Outputs per panel: anchors.tsv, graph.json (vertices in topological order) and
instances.jsonl (one instance per read and strand), plus conversion.json.
"""
from __future__ import annotations
import argparse
import collections
import gzip
import hashlib
import heapq
import json
import os
from pathlib import Path
import re
import subprocess

HERE = Path(__file__).resolve().parent
PILOT = HERE.parents[0] / 'mhc_pilot'
MINICHAIN = PILOT / 'vendor/minichain/minichain'
READS = PILOT / 'data/sample100.fq'
RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def load_graph(path: Path):
    seqs, edges, paths, names = {}, [], [], []
    with gzip.open(path, 'rt') as f:
        for line in f:
            a = line.rstrip('\n').split('\t')
            if a[0] == 'S':
                seqs[a[1]] = a[2]
            elif a[0] == 'L':
                assert a[2] == a[4] == '+' and a[5] == '0M'
                edges.append((a[1], a[3]))
            elif a[0] == 'W':
                assert '<' not in a[6]
                paths.append(re.findall(r'>([^><]+)', a[6]))
                names.append('|'.join(a[1:4]))
    assert set(seqs) == {v for p in paths for v in p}, 'every segment must lie on a haplotype'
    indeg = dict.fromkeys(seqs, 0)
    adj = collections.defaultdict(list)
    for u, v in edges:
        adj[u].append(v)
        indeg[v] += 1
    ready = [(int(v[1:]), v) for v in seqs if not indeg[v]]
    heapq.heapify(ready)
    order = []
    while ready:
        _, u = heapq.heappop(ready)
        order.append(u)
        for v in adj[u]:
            indeg[v] -= 1
            if not indeg[v]:
                heapq.heappush(ready, (int(v[1:]), v))
    assert len(order) == len(seqs), 'graph has a cycle'
    ids = {v: i for i, v in enumerate(order)}
    graph = dict(lengths=[len(seqs[v]) for v in order], edges=[[ids[u], ids[v]] for u, v in edges],
                 paths=[[ids[v] for v in p] for p in paths], names=order, haplotypes=names)
    for p in graph['paths']:
        assert all(u < v for u, v in zip(p, p[1:])), 'paths must follow the topological order'
    return graph, seqs, ids


def export(panel: Path):
    out = panel / 'anchors.tsv'
    if out.exists():
        raise SystemExit(f'{out} exists: the exporter appends, use a fresh path')
    env = dict(os.environ, RECOMB_DUMP_ANCHORS=str(out), OMP_NUM_THREADS='1')
    cmd = [str(MINICHAIN), '-x', 'lr', '-t', '1', str(panel / 'graph.gfa.gz'), str(READS)]
    proc = subprocess.run(cmd, env=env, capture_output=True, text=True)
    (panel / 'export.log').write_text(' '.join(cmd) + '\n' + proc.stderr)
    assert proc.returncode == 0, proc.stderr[-2000:]
    return cmd


def convert(panel: Path, cmd):
    graph, seqs, ids = load_graph(panel / 'graph.gfa.gz')
    reads = {}
    with READS.open() as f:
        while h := f.readline():
            reads[h.split()[0][1:]] = f.readline().strip()
            f.readline()
            f.readline()
    cases, seen, seeds, current, reverse = [], set(), 0, None, {}
    with (panel / 'anchors.tsv').open() as f:
        for line in f:
            a = line.rstrip('\n').split('\t')
            if a[0] == 'R':
                if current:
                    assert sum(len(x['anchors']) for x in current) == expected
                    cases.extend(current)
                name, qlen, expected = a[1], int(a[2]), int(a[3])
                assert name not in seen
                seen.add(name)
                assert len(reads[name]) == qlen
                current = [dict(id=name + ':' + s, read=name, strand=s, qlen=qlen, anchors=[]) for s in ('+', '-')]
            else:
                assert a[0] == 'A'
                seg, strand = a[1], int(a[2])
                y, d, span = map(int, a[3:])
                x, c = y - span + 1, d - span + 1
                if strand and seg not in reverse:
                    reverse[seg] = seqs[seg].translate(RC)[::-1]
                ref = seqs[seg] if strand == 0 else reverse[seg]
                assert 0 <= x <= y < len(ref) and 0 <= c <= d < qlen
                assert ref[x:y + 1].upper() == reads[name][c:d + 1].upper(), (name, seg, x, y, c, d)
                if strand:
                    x, y = len(ref) - 1 - y, len(ref) - 1 - x
                    c, d = qlen - 1 - d, qlen - 1 - c
                current[strand]['anchors'].append([ids[seg], x + 1, y + 1, c + 1, d + 1, span])
                seeds += 1
        if current:
            assert sum(len(x['anchors']) for x in current) == expected
            cases.extend(current)
    assert seen == set(reads)
    (panel / 'graph.json').write_text(json.dumps(graph))
    with (panel / 'instances.jsonl').open('w') as f:
        for c in cases:
            f.write(json.dumps(c) + '\n')
    meta = dict(panel=panel.name, H=len(graph['paths']), vertices=len(graph['lengths']), edges=len(graph['edges']),
                reads=len(seen), orientation_instances=len(cases), nonempty_instances=sum(bool(c['anchors']) for c in cases),
                anchors=seeds, exact_sequence_checks=seeds, command=cmd,
                minichain_sha256=hashlib.sha256(MINICHAIN.read_bytes()).hexdigest(),
                reads_sha256=hashlib.sha256(READS.read_bytes()).hexdigest(),
                anchor_weight='minimizer span', coordinates='1-based closed',
                reverse_policy='reverse-complement both graph segment and query to the forward graph',
                filter='unchanged Minichain lr seed collection; no anchor filtering')
    (panel / 'conversion.json').write_text(json.dumps(meta, indent=1) + '\n')
    return meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('panels', nargs='+')
    args = ap.parse_args()
    for name in args.panels:
        panel = HERE / 'data' / name
        cmd = export(panel)
        meta = convert(panel, cmd)
        print(json.dumps({k: v for k, v in meta.items() if k not in ('command',)}))


if __name__ == '__main__':
    main()
