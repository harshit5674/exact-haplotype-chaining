"""Place the 2,104 real CHM13 MHC HiFi reads on CHM13's MHC sequence by exact k-mer voting (2026-10-20).

Feasibility spike for a base-level MHC window. The reads are CHM13 HiFi reads,
so exact k-mer matches against CHM13's own MHC sequence locate them. The reference is the CHM13 walk of
Minichain's 60-haplotype MHC graph (experiments/mhc_pilot/data/MHC-CHM13.0.gfa.gz), which is CHM13v2.0
chr6:28,381,554-33,301,857 (experiments/hprc_mhc_panel/REPORT.md). Every 4th reference k-mer that occurs
once in the reference is indexed; every read k-mer is looked up in both orientations, and the read is placed
on the most supported diagonal (reference offset minus read offset). Output: one TSV row per read with the
strand, the implied CHM13 chr6 interval (0-based, half-open) and the vote support.

  python3.13 -B experiments/baselevel_feasibility/place_reads.py --out experiments/baselevel_feasibility/read_placement.tsv
"""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
GRAPH = ROOT / 'experiments/mhc_pilot/data/MHC-CHM13.0.gfa.gz'
READS = ROOT / 'experiments/mhc_pilot/data/MHC_CHM13_PacBio_filt.fq.gz'
CHR6_START = 28_381_554          # CHM13v2.0 chr6 offset of the walk's first base
K, STEP = 24, 4
RC = str.maketrans('ACGTacgt', 'TGCAtgca')


def chm13_sequence() -> str:
    seq, walk = {}, None
    with gzip.open(GRAPH, 'rt') as fh:
        for line in fh:
            a = line.rstrip('\n').split('\t')
            if a[0] == 'S':
                seq[a[1]] = a[2]
            elif a[0] == 'W' and a[1] == 'CHM13#0':
                walk = re.findall(r'([><])([^><]+)', a[6])
    assert walk
    return ''.join(seq[n] if o == '>' else seq[n].translate(RC)[::-1] for o, n in walk).upper()


def reads():
    with gzip.open(READS, 'rt') as fh:
        while True:
            head = fh.readline()
            if not head:
                return
            s = fh.readline().strip().upper()
            fh.readline(); fh.readline()
            yield head[1:].split()[0], s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    ref = chm13_sequence()
    assert len(ref) == 33_301_857 - CHR6_START, len(ref)
    counts = Counter(ref[i:i + K] for i in range(len(ref) - K + 1))
    index = {ref[i:i + K]: i for i in range(0, len(ref) - K + 1, STEP) if counts[ref[i:i + K]] == 1}
    del counts
    rows = []
    for name, s in reads():
        best = (0, None, None)
        for strand, q in (('+', s), ('-', s.translate(RC)[::-1])):
            votes = Counter()
            for j in range(len(q) - K + 1):
                i = index.get(q[j:j + K])
                if i is not None:
                    votes[i - j] += 1
            if votes:
                diag, n = votes.most_common(1)[0]
                # Indels shift the diagonal slightly; pool votes within 200 bp of the mode.
                n = sum(c for d, c in votes.items() if abs(d - diag) <= 200)
                if n > best[0]:
                    best = (n, strand, diag)
        n, strand, diag = best
        if strand is None:
            rows.append((name, len(s), '.', -1, -1, 0))
        else:
            lo = max(0, diag)
            hi = min(len(ref), diag + len(s))
            rows.append((name, len(s), strand, CHR6_START + lo, CHR6_START + hi, n))
    with args.out.open('x') as fh:
        fh.write('read\tlength\tstrand\tchm13_chr6_start\tchm13_chr6_end\tvotes\n')
        for r in rows:
            fh.write('\t'.join(map(str, r)) + '\n')
    placed = sum(1 for r in rows if r[2] != '.')
    print(f'{placed} of {len(rows)} reads placed; reference {len(ref):,} bp, {len(index):,} indexed {K}-mers')


if __name__ == '__main__':
    main()
