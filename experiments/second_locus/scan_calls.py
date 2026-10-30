"""One streaming pass over the 464-haplotype bubble calls, for locus selection.

Input (read only, shared): experiments/hprc_mhc_panel/data/CHM13-464.call.bed.gz, the
`minigraph call` + `mgutils.js merge` output of H. Li and G. Hickey (Zenodo 16728828,
v4.1); its published MD5 is checked first.

Outputs (data/scan/):
  bubbles.tsv.gz        one row per bubble, genome-wide: CHM13 chromosome, START, END,
                        VS, VE, number of alleles, NS, LBUBBLE, LSEGDUP, ANNO, the number
                        of allele walks containing a reverse step, and the raw CHM13 and
                        GRCh38 calls (GRCh38's call gives its GRCh38 coordinates, which
                        are used to lift well-known GRCh38 locus coordinates to CHM13).
  calls_<chrom>.bed.gz  the raw call lines of the candidate chromosomes, for threading.
  scan.json             counts, the input checksum and the script's arguments.
Nothing is filtered or changed; the raw lines are copied byte for byte.
"""
from __future__ import annotations
import argparse
import collections
import gzip
import hashlib
import json
from pathlib import Path
import time

HERE = Path(__file__).resolve().parent
SRC = HERE.parents[0] / 'hprc_mhc_panel/data/CHM13-464.call.bed.gz'
MD5 = 'c478805f722684ea7c0190bc8fc50205'
DEFAULT_CHROMS = 'chr1,chr2,chr6,chr14,chr19,chr20,chr22'


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1 << 22), b''):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--chroms', default=DEFAULT_CHROMS)
    args = ap.parse_args()
    out = HERE / 'data/scan'
    if (out / 'scan.json').exists():
        raise SystemExit('refusing to overwrite data/scan')
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    assert md5(SRC) == MD5, 'call BED does not match the published MD5'
    keep = {f'CHM13#0#{c}' for c in args.chroms.split(',')}
    writers = {c: gzip.open(out / f'calls_{c.split("#")[-1]}.bed.gz', 'wt', compresslevel=3) for c in keep}
    per_chrom = collections.Counter()
    header = None
    with gzip.open(SRC, 'rt') as f, gzip.open(out / 'bubbles.tsv.gz', 'wt', compresslevel=6) as summ:
        summ.write('chrom\tbeg\tend\tvs\tve\tNA\tNS\tLBUBBLE\tLSEGDUP\tANNO\trev_alleles\tmax_alen\tchm13_call\tgrch38_call\n')
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                header = line.rstrip('\n').split('\t')
                assert header[5:7] == ['CHM13.pri', 'GRCh38.pri']
                for w in writers.values():
                    w.write(line)
                continue
            chrom, beg, end, info, fmt, rest = line.split('\t', 5)
            assert fmt == 'GT:CSTRAND:CTG:CS:CE'
            kv = dict(x.split('=', 1) for x in info.split(';') if '=' in x)
            walks = [w for w in kv['AWALK'].split(',') if w]  # empty for NA=0 bubbles (e.g. alpha satellite)
            alen = [int(x) for x in kv['ALEN'].split(',') if x]
            chm13, grch38, _ = rest.split('\t', 2)
            summ.write('\t'.join([chrom, beg, end, kv['VS'], kv['VE'], kv['NA'], kv['NS'], kv.get('LBUBBLE', ''),
                                  kv.get('LSEGDUP', ''), kv.get('ANNO', ''), str(sum('<' in w for w in walks)),
                                  str(max(alen, default=0)), chm13, grch38]) + '\n')
            per_chrom[chrom] += 1
            if chrom in writers:
                writers[chrom].write(line)
    for w in writers.values():
        w.close()
    meta = dict(script='experiments/second_locus/scan_calls.py', argv=vars(args), source=str(SRC.relative_to(HERE.parents[1])),
                source_md5=MD5, samples=len(header) - 5, bubbles=sum(per_chrom.values()),
                bubbles_per_chrom=dict(sorted(per_chrom.items())), seconds=round(time.time() - t0, 1))
    (out / 'scan.json').write_text(json.dumps(meta, indent=1) + '\n')
    print(json.dumps({k: v for k, v in meta.items() if k != 'bubbles_per_chrom'}))


if __name__ == '__main__':
    main()
