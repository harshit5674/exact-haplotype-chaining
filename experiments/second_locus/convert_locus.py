"""Anchors for each new-locus panel, exported by the pinned Minichain and converted exactly.

Adapted from experiments/hprc_mhc_panel/convert.py (copied, with the read set and paths as
parameters; the conversion rules and checks are unchanged): the pinned Minichain v1.3 build
with the research seed export (experiments/mhc_pilot/minichain-export.patch: it dumps every
seed hit after collection and skips chaining) is run with `-x lr -t 1` on the panel graph and
the simulated reads; seeds become anchors with 1-based closed coordinates, reverse-strand hits
are mapped to the forward graph by reverse-complementing both sides, weight = seed span, and
every anchor is checked against the graph segment and read sequence.
Outputs per panel directory: anchors.tsv (the raw dump), export.log, graph.json (vertices in
topological order), instances.jsonl (one instance per read and strand) and conversion.json.

Usage: convert_locus.py kir hprc459 hprc60 ...   (graph and reads from data/kir/)
       convert_locus.py mhc_sim hprc463 hprc60   (graphs read from experiments/hprc_mhc_panel/data,
                                                  outputs in data/mhc_sim/<panel>/)
The exporter appends to its destination, so an existing anchors.tsv is refused.
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
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import locus_common as lc  # noqa: E402

RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def load_graph(path: Path):
    """convert.py's load_graph, unchanged."""
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


def export(out_dir: Path, gfa: Path, reads: Path):
    out = out_dir / 'anchors.tsv'
    if out.exists():
        raise SystemExit(f'{out} exists: the exporter appends, use a fresh path')
    out_dir.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, RECOMB_DUMP_ANCHORS=str(out), OMP_NUM_THREADS='1')
    cmd = [str(lc.MINICHAIN), '-x', 'lr', '-t', '1', str(gfa), str(reads)]
    proc = subprocess.run(cmd, env=env, capture_output=True, text=True)
    (out_dir / 'export.log').write_text(' '.join(cmd) + '\n' + proc.stderr)
    assert proc.returncode == 0, proc.stderr[-2000:]
    return cmd


def convert(out_dir: Path, gfa: Path, reads_path: Path, cmd):
    """convert.py's convert, with the graph and read paths as parameters; checks unchanged."""
    graph, seqs, ids = load_graph(gfa)
    reads = {}
    with reads_path.open() as f:
        while h := f.readline():
            reads[h.split()[0][1:]] = f.readline().strip()
            f.readline()
            f.readline()
    cases, seen, seeds, current, reverse = [], set(), 0, None, {}
    with (out_dir / 'anchors.tsv').open() as f:
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
    (out_dir / 'graph.json').write_text(json.dumps(graph))
    with (out_dir / 'instances.jsonl').open('w') as f:
        for c in cases:
            f.write(json.dumps(c) + '\n')
    meta = dict(panel=out_dir.name, graph=str(gfa.relative_to(lc.ROOT)), graph_sha256=lc.digest(gfa),
                H=len(graph['paths']), vertices=len(graph['lengths']), edges=len(graph['edges']),
                reads=len(seen), orientation_instances=len(cases), nonempty_instances=sum(bool(c['anchors']) for c in cases),
                anchors=seeds, exact_sequence_checks=seeds, command=[str(Path(x).relative_to(lc.ROOT)) if x.startswith('/') else x for x in cmd],
                minichain_sha256=hashlib.sha256(lc.MINICHAIN.read_bytes()).hexdigest(),
                reads_file=str(reads_path.relative_to(lc.ROOT)), reads_sha256=hashlib.sha256(reads_path.read_bytes()).hexdigest(),
                anchor_weight='minimizer span', coordinates='1-based closed',
                reverse_policy='reverse-complement both graph segment and query to the forward graph',
                filter='unchanged Minichain lr seed collection; no anchor filtering')
    (out_dir / 'conversion.json').write_text(json.dumps(meta, indent=1) + '\n')
    return meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('locus', help="kir, chr20 or mhc_sim")
    ap.add_argument('panels', nargs='+')
    args = ap.parse_args()
    reads = lc.DATA / args.locus / 'reads.fq'
    for name in args.panels:
        out_dir = lc.DATA / args.locus / name
        gfa = (lc.MHC / 'data' / name if args.locus == 'mhc_sim' else out_dir) / 'graph.gfa.gz'
        cmd = export(out_dir, gfa, reads)
        meta = convert(out_dir, gfa, reads, cmd)
        print(json.dumps({k: v for k, v in meta.items() if k not in ('command',)}), flush=True)


if __name__ == '__main__':
    main()
