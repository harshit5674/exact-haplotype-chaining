"""Build real MHC haplotype panels from HPRC release 2 (464-haplotype minigraph set).

Inputs (data/, public downloads; sizes and checksums are recorded in panels.json):
  CHM13-464.gfa.gz       minigraph graph of 464 haplotypes (462 HPRC release 2 assemblies,
                         GRCh38, CHM13) on CHM13v2.0; Zenodo record 16728828 (v4.1), H. Li
                         and G. Hickey
  CHM13-464.call.bed.gz  `minigraph call` + `mgutils.js merge`: for every bubble, each
                         haplotype's allele walk, contig, strand and contig coordinates
  ../mhc_pilot/data/MHC-CHM13.0.gfa.gz  Minichain's 60-haplotype MHC graph, used only to
                         locate the same region

Region: the CHM13 MHC sequence of the 60-haplotype graph, located by a unique exact match
in the CHM13 chr6 backbone of the 464-haplotype graph. The bubbles fully inside it form a
chain (each sink is the next source). A haplotype's walk is source, called allele walk,
sink, along the chain. This mirrors Minichain's own recipe (haplotypes threaded through a
minigraph SV graph), with the threading taken from the published calls.

A haplotype is kept if it has a call at every bubble and only forward steps. Contiguity
(all calls from one contig and strand, contig coordinates monotone) is checked and
required for the panels. Every walk step must be a forward link of the source graph.
Output per panel: data/<panel>/graph.gfa.gz in the style of Minichain's MHC graph (S lines,
L lines +/+ 0M, W lines `W sample#hap 1 MHC 1 1 walk`), restricted to segments that the
panel's haplotypes use. Subsets are seeded random samples that always contain CHM13.
"""
from __future__ import annotations
import argparse
import collections
import gzip
import hashlib
import json
from pathlib import Path
import random
import re

HERE = Path(__file__).resolve().parent
DATA = HERE / 'data'
PILOT_GFA = HERE.parents[0] / 'mhc_pilot/data/MHC-CHM13.0.gfa.gz'
SOURCES = {
    'CHM13-464.gfa.gz': dict(url='https://zenodo.org/records/16728828/files/CHM13-464.gfa.gz',
                             md5='a6f058f96b14f07d1cf88f71a5dcd518'),
    'CHM13-464.call.bed.gz': dict(url='https://zenodo.org/records/16728828/files/CHM13-464.call.bed.gz',
                                  md5='c478805f722684ea7c0190bc8fc50205'),
}
CHROM = 'CHM13#0#chr6'


def digest(path: Path, algo='sha256') -> str:
    h = hashlib.new(algo)
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1 << 22), b''):
            h.update(block)
    return h.hexdigest()


def locate_region():
    """[start, end) on CHM13 chr6 of the 60-haplotype graph's CHM13 MHC walk."""
    seqs, walk = {}, None
    with gzip.open(PILOT_GFA, 'rt') as f:
        for line in f:
            a = line.rstrip('\n').split('\t')
            if a[0] == 'S':
                seqs[a[1]] = a[2]
            elif a[0] == 'W' and a[1] == 'CHM13#0':
                walk = re.findall(r'>([^><]+)', a[6])
    mhc = ''.join(seqs[s] for s in walk).upper()
    backbone = []
    with gzip.open(DATA / 'CHM13-464.gfa.gz', 'rt') as f:
        for line in f:
            if line[0] != 'S' or f'\tSN:Z:{CHROM}\t' not in line:
                continue
            a = line.rstrip('\n').split('\t')
            tags = dict(t.split(':', 2)[0::2] for t in a[3:])
            if tags.get('SR') == '0':
                backbone.append((int(tags['SO']), a[1], a[2]))
    backbone.sort()
    assert all(backbone[k][0] + len(backbone[k][2]) == backbone[k + 1][0] for k in range(len(backbone) - 1))
    chrom = ''.join(s for _, _, s in backbone).upper()
    start = chrom.find(mhc)
    assert start >= 0 and chrom.find(mhc, start + 1) < 0, 'MHC sequence must occur exactly once'
    return dict(chrom=CHROM, start=start, end=start + len(mhc), length=len(mhc),
                note='0-based half-open CHM13v2.0 coordinates; exact unique match of the 60-haplotype graph CHM13 walk')


def read_calls(region):
    bubbles, samples = [], None
    with gzip.open(DATA / 'CHM13-464.call.bed.gz', 'rt') as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                samples = line.rstrip('\n').split('\t')[5:]
                continue
            if not line.startswith(CHROM + '\t'):
                continue
            a = line.rstrip('\n').split('\t')
            beg, end = int(a[1]), int(a[2])
            if beg < region['start'] or end > region['end']:
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
            bubbles.append(dict(beg=beg, end=end, vs=info['VS'], ve=info['VE'],
                                awalk=info['AWALK'].split(','), calls=calls))
    bubbles.sort(key=lambda b: (b['beg'], b['end']))
    assert samples and all(len(b['calls']) == len(samples) for b in bubbles)
    return samples, bubbles


def haplotype_name(sample: str) -> str:
    base, tag = sample.rsplit('.', 1)
    return f"{base}#{dict(pri=0, pat=1, mat=2, hap1=1, hap2=2)[tag]}"


def thread(samples, bubbles):
    """Walks along the bubble chain, and per-haplotype status."""
    for a, b in zip(bubbles, bubbles[1:]):
        assert a['ve'] == b['vs'], 'bubbles must form a chain'
    walks, status = {}, {}
    for h, sample in enumerate(samples):
        steps, why, contigs, strands, coords = [bubbles[0]['vs']], None, set(), set(), []
        for b in bubbles:
            call = b['calls'][h]
            if call is None:
                why = 'missing call'
                break
            gt, strand, ctg, cs, ce = call
            allele = b['awalk'][gt]
            if allele != '*':
                steps.extend(re.findall(r'[><][^><]+', allele))
            steps.append(b['ve'])
            contigs.add(ctg)
            strands.add(strand)
            coords.append(cs)
        if why is None and any(s[0] != '>' for s in steps):
            why = 'reverse step'
        contiguous = None
        if why is None:
            monotone = (all(x <= y for x, y in zip(coords, coords[1:])) if strands == {'+'} else
                        all(x >= y for x, y in zip(coords, coords[1:])) if strands == {'-'} else False)
            contiguous = len(contigs) == 1 and len(strands) == 1 and monotone
            if not contiguous:
                why = 'not one contiguous contig'
        status[sample] = dict(kept=why is None, reason=why, contigs=len(contigs), contiguous=contiguous)
        if why is None:
            walks[haplotype_name(sample)] = [s[1:] for s in steps]
    return walks, status


def load_segments(names):
    seqs, edges = {}, set()
    with gzip.open(DATA / 'CHM13-464.gfa.gz', 'rt') as f:
        for line in f:
            if line[0] == 'S':
                a = line.split('\t', 3)
                if a[1] in names:
                    seqs[a[1]] = a[2].rstrip('\n')
            elif line[0] == 'L':
                a = line.rstrip('\n').split('\t')
                if a[1] in names and a[3] in names:
                    edges.add((a[1], a[2], a[3], a[4], a[5]))
    return seqs, edges


def write_panel(name, walks, seqs, links):
    out = DATA / name
    out.mkdir(parents=True, exist_ok=True)
    used = {s for w in walks.values() for s in w}
    steps = {(u, v) for w in walks.values() for u, v in zip(w, w[1:])}
    for u, v in steps:
        assert (u, '+', v, '+', '0M') in links, ('walk step is not a forward link', u, v)
    path = out / 'graph.gfa.gz'
    with gzip.open(path, 'wt') as f:
        for s in sorted(used, key=lambda x: int(x[1:])):
            f.write(f'S\t{s}\t{seqs[s]}\n')
        for u, v in sorted(steps, key=lambda e: (int(e[0][1:]), int(e[1][1:]))):
            f.write(f'L\t{u}\t+\t{v}\t+\t0M\n')
        for hap, w in walks.items():
            f.write(f'W\t{hap}\t1\tMHC\t1\t1\t' + ''.join('>' + s for s in w) + '\n')
    lens = [sum(len(seqs[s]) for s in w) for w in walks.values()]
    return dict(panel=name, H=len(walks), segments=len(used), edges=len(steps),
                total_segment_length=sum(len(seqs[s]) for s in used),
                haplotype_length_min=min(lens), haplotype_length_max=max(lens),
                haplotypes=list(walks), gfa=str(path.relative_to(HERE)), gfa_sha256=digest(path))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--subsets', type=lambda s: [int(x) for x in s.split(',')], default=[60, 120, 240])
    ap.add_argument('--seed', type=int, default=4640)
    args = ap.parse_args()
    summary_path = HERE / 'panels.json'
    if summary_path.exists():
        raise SystemExit('refusing to overwrite panels.json')
    for n, src in SOURCES.items():
        assert digest(DATA / n, 'md5') == src['md5'], f'{n} does not match the published MD5'
    region = locate_region()
    samples, bubbles = read_calls(region)
    walks, status = thread(samples, bubbles)
    names = {s for w in walks.values() for s in w}
    seqs, links = load_segments(names)
    assert names <= set(seqs), 'missing segment sequences'
    reasons = collections.Counter(v['reason'] for v in status.values() if not v['kept'])
    panels = [write_panel(f'hprc{len(walks)}', walks, seqs, links)]
    rng = random.Random(args.seed)
    others = sorted(h for h in walks if h != 'CHM13#0')
    assert 'CHM13#0' in walks
    for size in args.subsets:
        chosen = ['CHM13#0'] + rng.sample(others, size - 1)
        panels.append(write_panel(f'hprc{size}', {h: walks[h] for h in chosen}, seqs, links))
    first, last = bubbles[0], bubbles[-1]
    summary = dict(
        script='experiments/hprc_mhc_panel/build_panel.py', argv=vars(args), region=region,
        bubbles=len(bubbles), chain_source=first['vs'], chain_sink=last['ve'],
        chain_span=[first['beg'], last['end']], samples=len(samples),
        kept=len(walks), excluded=dict(reasons),
        excluded_samples={s: v['reason'] for s, v in status.items() if not v['kept']},
        multi_contig_samples=sum(1 for v in status.values() if v['contigs'] > 1),
        inputs={n: dict(bytes=(DATA / n).stat().st_size, md5=digest(DATA / n, 'md5'), sha256=digest(DATA / n),
                        url=SOURCES[n]['url'], published_md5=SOURCES[n]['md5']) for n in SOURCES},
        pilot_gfa_md5=digest(PILOT_GFA, 'md5'), panels=panels)
    summary_path.write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps({k: v for k, v in summary.items() if k not in ('excluded_samples', 'panels')}, indent=1))
    for pnl in panels:
        print(pnl['panel'], {k: v for k, v in pnl.items() if k not in ('haplotypes', 'panel')})


if __name__ == '__main__':
    main()
