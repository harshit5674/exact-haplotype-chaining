"""Turn an extracted base-level window into a project panel with real-read anchors (2026-10-20, feasibility spike).

Input: windows/<name>/window.gfa from inspect_window.py (HPRC release 2 Minigraph-Cactus v2.1, CHM13 chr6
interval, every top-level snarl overlapping it, all haplotype walks). Nothing is pruned: every segment, link
and walk of the extraction is kept. Walks are oriented as CHM13 runs (inspect_window.py found none to reverse);
oriented segments are renumbered s1..sV in a deterministic topological order. Haplotypes are named H001.. in
the order gbz-base lists them (it reports them as unknown haplotypes; only CHM13 keeps its name).

Reads: the real CHM13 MHC HiFi reads that lie entirely inside the window, from read_placement.tsv (k-mer
placement on CHM13; support >= 0.5 of indexed k-mers). Seeds: the pinned Minichain seed exporter
(experiments/mhc_pilot/vendor/minichain), converted and sequence-checked by the unchanged
experiments/hprc_mhc_panel/convert.py, exactly as the SV-level panels and the MICB pilot were.

  python3.13 -B experiments/baselevel_feasibility/prepare_window.py --window windows/w31250k_250kb
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import gzip
import hashlib
import heapq
import importlib.util
import json
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FASTQ = ROOT / 'experiments/mhc_pilot/data/MHC_CHM13_PacBio_filt.fq.gz'
RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def flip(h):
    return h[0], not h[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--window', type=Path, required=True, help='directory with window.gfa and inspect.json')
    ap.add_argument('--min-support', type=float, default=0.5)
    args = ap.parse_args()
    win = args.window if args.window.is_absolute() else HERE / args.window
    info = json.loads((win / 'inspect.json').read_text())
    assert info['acyclic'] and not info['repeated_handle_walks'] and not info['segments_in_both_orientations']
    out = win / 'panel'
    out.mkdir()                                     # refuse to overwrite
    seq, links, walks = {}, [], []
    for line in (win / 'window.gfa').open():
        a = line.rstrip('\n').split('\t')
        if a[0] == 'S':
            seq[int(a[1])] = a[2]
        elif a[0] == 'L':
            assert a[5] in ('0M', '*'), a[5]
            links.append(((int(a[1]), a[2] == '-'), (int(a[3]), a[4] == '-')))
        elif a[0] == 'W':
            walks.append((a[1:6], [(int(n), o == '<') for o, n in re.findall(r'([><])(\d+)', a[6])]))
    ref = next(s for m, s in walks if m[0] == 'CHM13')
    ref_pos = {h: i for i, h in enumerate(ref)}
    oriented = []
    for m, s in walks:
        fwd = sum(h in ref_pos for h in s)
        rev = sum(flip(h) in ref_pos for h in s)
        oriented.append((m, [flip(h) for h in reversed(s)] if rev > fwd else s))
    handles = {h for _, s in oriented for h in s}
    assert {n for n, _ in handles} == set(seq), 'every segment lies on a haplotype in this extraction'
    adj = defaultdict(set)
    for u, v in links:
        adj[u].add(v)
        adj[flip(v)].add(flip(u))
    edges = {(u, v) for u in handles for v in adj[u] if v in handles}
    for _, s in oriented:
        assert all((u, v) in edges for u, v in zip(s, s[1:])), 'walk step without a link'
    succ = defaultdict(list)
    indeg = dict.fromkeys(handles, 0)
    for u, v in edges:
        succ[u].append(v)
        indeg[v] += 1
    ready = [u for u in handles if not indeg[u]]
    heapq.heapify(ready)
    order = []
    while ready:
        u = heapq.heappop(ready)
        order.append(u)
        for v in sorted(succ[u]):
            indeg[v] -= 1
            if not indeg[v]:
                heapq.heappush(ready, v)
    assert len(order) == len(handles), 'cycle'
    ids = {h: i for i, h in enumerate(order)}
    sequences = [seq[n].translate(RC)[::-1] if r else seq[n] for n, r in order]
    names = ['CHM13' if m[0] == 'CHM13' else f'H{k:03d}' for k, (m, _) in enumerate(oriented)]
    assert len(set(names)) == len(names)
    with gzip.open(out / 'graph.gfa.gz', 'wt') as fh:
        fh.write('H\tVN:Z:1.1\n')
        for i, s in enumerate(sequences):
            fh.write(f'S\ts{i + 1}\t{s}\n')
        for u, v in sorted((ids[u], ids[v]) for u, v in edges):
            fh.write(f'L\ts{u + 1}\t+\ts{v + 1}\t+\t0M\n')
        for name, (_, s) in zip(names, oriented):
            p = [ids[h] for h in s]
            assert all(x < y for x, y in zip(p, p[1:]))
            fh.write(f"W\t{name}\t1\twindow\t0\t{sum(len(sequences[v]) for v in p)}\t" + ''.join(f'>s{v + 1}' for v in p) + '\n')
    # Real reads lying inside the window (CHM13 chr6 coordinates of the extraction).
    lo, hi = info['start'], info['end']
    keep = []
    for r in csv.DictReader((HERE / 'read_placement.tsv').open(), delimiter='\t'):
        if r['strand'] == '.':
            continue
        support = int(r['votes']) / (int(r['length']) / 4)
        if support >= args.min_support and lo <= int(r['chm13_chr6_start']) and int(r['chm13_chr6_end']) <= hi:
            keep.append(r['read'])
    keep_set = set(keep)
    with gzip.open(FASTQ, 'rt') as src, (out / 'reads.fq').open('w') as dst:
        while True:
            block = [src.readline() for _ in range(4)]
            if not block[0]:
                break
            if block[0][1:].split()[0] in keep_set:
                dst.writelines(block)
    spec = importlib.util.spec_from_file_location('panel_convert', ROOT / 'experiments/hprc_mhc_panel/convert.py')
    cv = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cv)
    cv.READS = out / 'reads.fq'
    cmd = cv.export(out)
    meta = cv.convert(out, cmd)
    graph = json.loads((out / 'graph.json').read_text())
    lengths = graph['lengths']
    anchored = defaultdict(int)
    for line in (out / 'instances.jsonl').open():
        for a in json.loads(line)['anchors']:
            anchored[a[0]] += 1
    prov = dict(window=info['start'], end=info['end'], source_gfa_sha256=info['gfa_sha256'], reads=len(keep),
                min_support=args.min_support, segments=len(order), edges=len(edges), haplotypes=len(oriented),
                reversed_walks=sum(1 for (m, s0), (_, s) in zip(walks, oriented) if s is not s0),
                anchored_vertices=len(anchored), anchored_one_base_vertices=sum(1 for v in anchored if lengths[v] == 1),
                vertices_ge17=sum(1 for x in lengths if x >= 17), conversion=meta,
                fastq_sha256=hashlib.sha256((out / 'reads.fq').read_bytes()).hexdigest())
    (out / 'provenance.json').write_text(json.dumps(prov, indent=1) + '\n')
    print(json.dumps({k: v for k, v in prov.items() if k != 'conversion'}, indent=1))
    print(json.dumps({k: v for k, v in meta.items() if k not in ('command',)}, indent=1))


if __name__ == '__main__':
    main()
