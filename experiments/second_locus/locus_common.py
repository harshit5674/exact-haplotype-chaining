"""Shared helpers for the second-locus experiment (paths, bubble calls, threading, lifting).

The threading rule is build_panel.py's (experiments/hprc_mhc_panel), unchanged: bubbles fully
inside the window must form a chain (each sink is the next source); a haplotype's walk is
the first source, then for every bubble its called allele walk and the bubble's sink. A
haplotype is kept iff it has a call at every bubble, only forward steps, and all calls on
one contig and one strand with monotone contig coordinates. Here the exclusion reasons are
also counted separately (build_panel.py records only the first failing rule).
"""
from __future__ import annotations
import bisect
import collections
import gzip
import hashlib
import json
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATA = HERE / 'data'
MHC = ROOT / 'experiments/hprc_mhc_panel'
GFA = MHC / 'data/CHM13-464.gfa.gz'
CALLS = MHC / 'data/CHM13-464.call.bed.gz'
SOURCES = {
    'CHM13-464.gfa.gz': dict(path=GFA, url='https://zenodo.org/records/16728828/files/CHM13-464.gfa.gz',
                             md5='a6f058f96b14f07d1cf88f71a5dcd518'),
    'CHM13-464.call.bed.gz': dict(path=CALLS, url='https://zenodo.org/records/16728828/files/CHM13-464.call.bed.gz',
                                  md5='c478805f722684ea7c0190bc8fc50205'),
}
MINICHAIN = ROOT / 'experiments/mhc_pilot/vendor/minichain/minichain'
NEW = ROOT / 'research/current/new'
SCALE, GAMMA = 200, 10000
STEP = re.compile(r'[><][^><]+')

# Chosen loci (survey_loci.py, survey.json). Coordinates: CHM13v2.0, 0-based half-open; the
# bubbles fully inside [start, end) are threaded. KIR's window is the GRCh38 LRC/KIR span
# chr19:54,025,634-55,084,318 lifted through GRCh38's own bubble calls (build_locus.py
# recomputes and asserts the lift).
LOCI = {
    'kir': dict(chrom='chr19', start=57104214, end=58179183, walk_name='LRC_KIR', grch38=[54025634, 55084318],
                label='LRC/KIR (leukocyte receptor complex with the KIR gene cluster), chr19',
                note='SV-rich and medically relevant; GRCh38 LRC/KIR span lifted to CHM13 via GRCh38 calls'),
    'chr20': dict(chrom='chr20', start=50000000, end=55000000, walk_name='CHR20_50_55', grch38=None,
                  label='chr20:50-55 Mb (typical euchromatic region)',
                  note='typical: 5 Mb window with genome-typical bubble density and little segmental duplication'),
}


def digest(path: Path, algo='sha256') -> str:
    h = hashlib.new(algo)
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1 << 22), b''):
            h.update(block)
    return h.hexdigest()


def haplotype_name(sample: str) -> str:
    """build_panel.py's naming: HG002.pat -> HG002#1, CHM13.pri -> CHM13#0."""
    base, tag = sample.rsplit('.', 1)
    return f"{base}#{dict(pri=0, pat=1, mat=2, hap1=1, hap2=2)[tag]}"


def load_bubble_table(chrom: str) -> list[list[str]]:
    """Rows of data/scan/bubbles.tsv.gz for one CHM13 chromosome (e.g. 'chr19'), in file order."""
    out = []
    key = f'CHM13#0#{chrom}\t'
    with gzip.open(DATA / 'scan/bubbles.tsv.gz', 'rt') as f:
        f.readline()
        for line in f:
            if line.startswith(key):
                out.append(line.rstrip('\n').split('\t'))
    return out


def grch38_points(rows, chrom: str):
    """(GRCh38 contig start, CHM13 bubble start) for bubbles GRCh38 calls on the same chromosome, + strand."""
    pts = []
    for a in rows:
        if a[13] == '.':
            continue
        gt, strand, ctg, cs, ce = a[13].split(':')
        if ctg == f'GRCh38#0#{chrom}' and strand == '+':
            pts.append((int(cs), int(a[1])))
    pts.sort()
    return pts


def lift(pts, g: int) -> int:
    """CHM13 coordinate for GRCh38 coordinate g: offset of the nearest GRCh38-called bubble at or before g."""
    k = bisect.bisect_right(pts, (g, float('inf'))) - 1
    k = max(k, 0)
    gs, cs = pts[k]
    return cs + (g - gs)


def read_calls(chrom: str, start: int, end: int):
    """Bubbles with start <= START and END <= end on a CHM13 chromosome, as build_panel.read_calls."""
    bubbles, samples = [], None
    key = f'CHM13#0#{chrom}\t'
    with gzip.open(DATA / f'scan/calls_{chrom}.bed.gz', 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                samples = line.rstrip('\n').split('\t')[5:]
                continue
            if not line.startswith(key):
                continue
            a = line.rstrip('\n').split('\t')
            beg, stop = int(a[1]), int(a[2])
            if beg < start or stop > end:
                continue
            info = dict(kv.split('=', 1) for kv in a[3].split(';') if '=' in kv)
            assert a[4] == 'GT:CSTRAND:CTG:CS:CE'
            calls = []
            for x in a[5:]:
                if x == '.':
                    calls.append(None)
                else:
                    gt, strand, ctg, cs, ce = x.split(':')
                    calls.append((int(gt), strand, ctg, int(cs), int(ce)))
            bubbles.append(dict(beg=beg, end=stop, vs=info['VS'], ve=info['VE'],
                                awalk=info['AWALK'].split(','), calls=calls, anno=info.get('ANNO', ''),
                                lbubble=int(info.get('LBUBBLE') or 0),
                                alen=[int(x) for x in info.get('ALEN', '').split(',') if x]))
    bubbles.sort(key=lambda b: (b['beg'], b['end']))
    assert samples and all(len(b['calls']) == len(samples) for b in bubbles)
    return samples, bubbles


def chain_breaks(bubbles) -> list[int]:
    """Indices k with bubbles[k].ve != bubbles[k+1].vs."""
    return [k for k in range(len(bubbles) - 1) if bubbles[k]['ve'] != bubbles[k + 1]['vs']]


def thread(samples, bubbles):
    """build_panel.thread, plus a full per-haplotype diagnosis (all rules, not only the first)."""
    assert not chain_breaks(bubbles), 'bubbles must form a chain'
    walks, status = {}, {}
    for h, sample in enumerate(samples):
        steps, contigs, strands, coords = [bubbles[0]['vs']], [], set(), []
        missing = reverse_alleles = 0
        for b in bubbles:
            call = b['calls'][h]
            if call is None:
                missing += 1
                continue
            gt, strand, ctg, cs, ce = call
            allele = b['awalk'][gt]
            if allele != '*':
                st = STEP.findall(allele)
                reverse_alleles += any(s[0] != '>' for s in st)
                steps.extend(st)
            steps.append(b['ve'])
            if not contigs or contigs[-1] != ctg:
                contigs.append(ctg)
            strands.add(strand)
            coords.append(cs)
        distinct = set(contigs)
        monotone = (all(x <= y for x, y in zip(coords, coords[1:])) if strands == {'+'} else
                    all(x >= y for x, y in zip(coords, coords[1:])) if strands == {'-'} else False)
        # build_panel.py order: missing call, then reverse step, then contiguity
        why = None
        if missing:
            why = 'missing call'
        elif any(s[0] != '>' for s in steps):
            why = 'reverse step'
        elif not (len(distinct) == 1 and len(strands) == 1 and monotone):
            why = 'not one contiguous contig'
        status[sample] = dict(kept=why is None, reason=why, missing_calls=missing, reverse_alleles=reverse_alleles,
                              contigs=len(distinct), contig_switches=len(contigs) - 1, strands=''.join(sorted(strands)),
                              monotone=monotone)
        if why is None:
            walks[haplotype_name(sample)] = [s[1:] for s in steps]
    return walks, status


def summarize_status(status) -> dict:
    s = list(status.values())
    return dict(samples=len(s), kept=sum(v['kept'] for v in s),
                excluded=dict(collections.Counter(v['reason'] for v in s if not v['kept'])),
                with_missing_calls=sum(v['missing_calls'] > 0 for v in s),
                with_reverse_alleles=sum(v['reverse_alleles'] > 0 for v in s),
                multi_contig=sum(v['contigs'] > 1 for v in s),
                mixed_strand=sum(len(v['strands']) > 1 for v in s),
                non_monotone_single_contig=sum(v['contigs'] == 1 and len(v['strands']) == 1 and not v['monotone'] for v in s))


def dump(path: Path, obj):
    path.write_text(json.dumps(obj, indent=1) + '\n')
