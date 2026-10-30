"""Export anchors for simulated reads with the pinned Minichain seed exporter, converted exactly as
experiments/hprc_mhc_panel/convert.py does (same binary, same flags, same conversion and checks).

The pinned binary (experiments/mhc_pilot/vendor/minichain/minichain, research export patch) dumps every
seed hit and skips chaining when RECOMB_DUMP_ANCHORS is set; it appends, so the dump path must be new.
Conversion: 1-based closed coordinates, reverse-strand hits mapped to the forward graph by
reverse-complementing both segment and read, weight = seed span, exact sequence check of every anchor.

  python3.13 experiments/exactness_payoff/export_anchors.py --panel hprc60 --prefix experiments/exactness_payoff/sim/hprc60
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'experiments/hprc_mhc_panel'))
import convert  # noqa: E402

MINICHAIN = ROOT / 'experiments/mhc_pilot/vendor/minichain/minichain'
PANELS = {name: ROOT / 'experiments/hprc_mhc_panel/data' / name for name in ('hprc60', 'hprc463')}
RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def read_fastq(path: Path) -> dict:
    reads = {}
    with gzip.open(path, 'rt') as f:
        while h := f.readline():
            reads[h.split()[0][1:]] = f.readline().strip()
            f.readline()
            f.readline()
    return reads


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--panel', choices=sorted(PANELS), required=True)
    ap.add_argument('--prefix', required=True, help='simulation prefix (<prefix>.fq.gz exists)')
    args = ap.parse_args()
    prefix = Path(args.prefix)
    fq = Path(str(prefix) + '.fq.gz')
    dump = Path(str(prefix) + '.anchors.tsv')
    inst_path = Path(str(prefix) + '.instances.jsonl')
    if dump.exists() or inst_path.exists():
        raise SystemExit('outputs exist: the exporter appends, use a fresh prefix')
    panel = PANELS[args.panel]
    env = dict(os.environ, RECOMB_DUMP_ANCHORS=str(dump), OMP_NUM_THREADS='1')
    cmd = [str(MINICHAIN), '-x', 'lr', '-t', '1', str(panel / 'graph.gfa.gz'), str(fq)]
    proc = subprocess.run(cmd, env=env, capture_output=True, text=True)
    Path(str(prefix) + '.export.log').write_text(' '.join(cmd) + '\n' + proc.stderr)
    assert proc.returncode == 0, proc.stderr[-2000:]

    graph, seqs, ids = convert.load_graph(panel / 'graph.gfa.gz')
    stored = json.loads((panel / 'graph.json').read_text())
    assert graph['lengths'] == stored['lengths'] and graph['paths'] == stored['paths'], 'vertex numbering differs'
    reads = read_fastq(fq)
    cases, seen, seeds, current, reverse, expected = [], set(), 0, None, {}, 0
    with dump.open() as f:
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
    missing = sorted(set(reads) - seen)
    with inst_path.open('x') as f:
        for c in cases:
            f.write(json.dumps(c) + '\n')
    meta = dict(panel=args.panel, reads=len(reads), exported_reads=len(seen), reads_without_record=missing,
                orientation_instances=len(cases), anchors=seeds, exact_sequence_checks=seeds, command=cmd,
                minichain_sha256=hashlib.sha256(MINICHAIN.read_bytes()).hexdigest(),
                fastq_sha256=hashlib.sha256(fq.read_bytes()).hexdigest())
    Path(str(prefix) + '.conversion.json').write_text(json.dumps(meta, indent=1) + '\n')
    # The raw dump is reproducible from the FASTQ and the pinned binary; keep it compressed.
    with dump.open('rb') as src, gzip.open(str(dump) + '.gz', 'wb') as dst:
        dst.write(src.read())
    dump.unlink()
    print(json.dumps({k: v for k, v in meta.items() if k != 'command'}))


if __name__ == '__main__':
    main()
