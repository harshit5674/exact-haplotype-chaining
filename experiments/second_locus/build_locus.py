"""Build HPRC release 2 haplotype panels for new loci, in the style of build_panel.py.

Inputs (read only): the 464-haplotype minigraph graph CHM13-464.gfa.gz and the bubble calls
CHM13-464.call.bed.gz in experiments/hprc_mhc_panel/data (Zenodo 16728828, v4.1; published
MD5s checked), via data/scan/calls_<chrom>.bed.gz written by scan_calls.py (raw lines).

Per locus (locus_common.LOCI) the bubbles fully inside the window are threaded exactly as in
build_panel.py: a chain of bubbles; a haplotype is kept iff it has a call at every bubble,
only forward steps, and one contig, one strand, monotone coordinates; every walk step must be
a forward link of the source graph. Panels (data/<locus>/<panel>/graph.gfa.gz: S lines,
L lines +/+ 0M, W lines `W sample#hap 1 <locus> 1 1 walk`), each pruned to its own segments:
  hprc<N>            all usable haplotypes (the hprc463 equivalent)
  hprc60/120/240     seeded random subsets that always contain CHM13 (build_panel.py's
                     procedure and seed 4640: one Random(seed), `others` sorted, sequential
                     rng.sample per size, so the subsets are independent, not nested)
  hprc<M>_leaveout   all usable haplotypes minus the held-out samples (both haplotypes)
Held-out samples (read sources absent from the panel): samples with both haplotypes usable
and in none of the 60/120/240 subsets, drawn with --holdout-seed. Present sources: CHM13 and
haplotypes in all three subsets (so present in every panel, including the leave-out one).
Output: data/<locus>/panels.json with sizes, checksums, exclusions, held-out and present sets.
"""
from __future__ import annotations
import argparse
import collections
import gzip
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import locus_common as lc  # noqa: E402


def load_segments(names):
    """build_panel.load_segments: sequences of the named segments and the links among them."""
    seqs, edges = {}, set()
    with gzip.open(lc.GFA, 'rt') as f:
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


def write_panel(locus, name, walks, seqs, links):
    """build_panel.write_panel, with the locus name in the W lines and a per-locus folder."""
    out = lc.DATA / locus / name
    out.mkdir(parents=True, exist_ok=False)
    used = {s for w in walks.values() for s in w}
    steps = {(u, v) for w in walks.values() for u, v in zip(w, w[1:])}
    for u, v in steps:
        assert (u, '+', v, '+', '0M') in links, ('walk step is not a forward link', u, v)
    path = out / 'graph.gfa.gz'
    walk_name = lc.LOCI[locus]['walk_name']
    with gzip.open(path, 'wt') as f:
        for s in sorted(used, key=lambda x: int(x[1:])):
            f.write(f'S\t{s}\t{seqs[s]}\n')
        for u, v in sorted(steps, key=lambda e: (int(e[0][1:]), int(e[1][1:]))):
            f.write(f'L\t{u}\t+\t{v}\t+\t0M\n')
        for hap, w in walks.items():
            f.write(f'W\t{hap}\t1\t{walk_name}\t1\t1\t' + ''.join('>' + s for s in w) + '\n')
    lens = [sum(len(seqs[s]) for s in w) for w in walks.values()]
    return dict(panel=name, H=len(walks), segments=len(used), edges=len(steps),
                total_segment_length=sum(len(seqs[s]) for s in used),
                haplotype_length_min=min(lens), haplotype_length_max=max(lens),
                haplotype_length_median=sorted(lens)[len(lens) // 2],
                haplotypes=list(walks), gfa=str(path.relative_to(HERE)), gfa_sha256=lc.digest(path))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--loci', default='kir,chr20')
    ap.add_argument('--subsets', type=lambda s: [int(x) for x in s.split(',')], default=[60, 120, 240])
    ap.add_argument('--seed', type=int, default=4640, help="build_panel.py's subset seed")
    ap.add_argument('--holdout', type=int, default=3, help='held-out samples (both haplotypes) per locus')
    ap.add_argument('--holdout-seed', type=int, default=2027)
    ap.add_argument('--present', type=int, default=2, help='present sources besides CHM13 (in all three subsets)')
    args = ap.parse_args()
    loci = args.loci.split(',')
    for locus in loci:
        if (lc.DATA / locus / 'panels.json').exists():
            raise SystemExit(f'refusing to overwrite data/{locus}/panels.json')
    inputs = {}
    for n, src in lc.SOURCES.items():
        assert lc.digest(src['path'], 'md5') == src['md5'], f'{n} does not match the published MD5'
        inputs[n] = dict(bytes=src['path'].stat().st_size, md5=src['md5'], sha256=lc.digest(src['path']), url=src['url'])
    plans, names = {}, set()
    for locus in loci:
        spec = lc.LOCI[locus]
        if spec['grch38']:
            pts = lc.grch38_points(lc.load_bubble_table(spec['chrom']), spec['chrom'])
            assert [lc.lift(pts, g) for g in spec['grch38']] == [spec['start'], spec['end']], 'lift changed'
        samples, bubbles = lc.read_calls(spec['chrom'], spec['start'], spec['end'])
        walks, status = lc.thread(samples, bubbles)
        assert 'CHM13#0' in walks
        plans[locus] = (samples, bubbles, walks, status)
        names |= {s for w in walks.values() for s in w}
        print(locus, lc.summarize_status(status), flush=True)
    seqs, links = load_segments(names)
    assert names <= set(seqs), 'missing segment sequences'
    for locus in loci:
        samples, bubbles, walks, status = plans[locus]
        spec = lc.LOCI[locus]
        panels = [write_panel(locus, f'hprc{len(walks)}', walks, seqs, links)]
        rng = random.Random(args.seed)
        others = sorted(h for h in walks if h != 'CHM13#0')
        subsets = {}
        for size in args.subsets:
            chosen = ['CHM13#0'] + rng.sample(others, size - 1)
            subsets[size] = set(chosen)
            panels.append(write_panel(locus, f'hprc{size}', {h: walks[h] for h in chosen}, seqs, links))
        in_any = set().union(*subsets.values())
        in_all = set.intersection(*subsets.values())
        by_sample = collections.defaultdict(list)
        for h in walks:
            by_sample[h.split('#')[0]].append(h)
        eligible = sorted(s for s, hs in by_sample.items() if len(hs) == 2 and not set(hs) & in_any)
        hrng = random.Random(args.holdout_seed)
        held = sorted(hrng.sample(eligible, args.holdout))
        held_haps = sorted(h for s in held for h in by_sample[s])
        present = ['CHM13#0'] + sorted(hrng.sample(sorted(in_all - {'CHM13#0'}), args.present))
        leave = {h: w for h, w in walks.items() if h not in held_haps}
        panels.append(write_panel(locus, f'hprc{len(leave)}_leaveout', leave, seqs, links))
        first, last = bubbles[0], bubbles[-1]
        summary = dict(
            script='experiments/second_locus/build_locus.py', argv=vars(args), locus=locus,
            region=dict(chrom=f'CHM13#0#{spec["chrom"]}', start=spec['start'], end=spec['end'],
                        length=spec['end'] - spec['start'], grch38=spec['grch38'], label=spec['label'], note=spec['note'],
                        coordinates='0-based half-open CHM13v2.0; bubbles fully inside are threaded'),
            bubbles=len(bubbles), chain_source=first['vs'], chain_sink=last['ve'], chain_span=[first['beg'], last['end']],
            samples=len(samples), kept=len(walks), status_summary=lc.summarize_status(status),
            excluded=dict(collections.Counter(v['reason'] for v in status.values() if not v['kept'])),
            excluded_samples={s: v for s, v in status.items() if not v['kept']},
            multi_contig_samples=sum(1 for v in status.values() if v['contigs'] > 1),
            subsets_nested=False, subset_intersection=sorted(in_all), eligible_holdout_samples=len(eligible),
            holdout_samples=held, holdout_haplotypes=held_haps, present_sources=present,
            inputs=inputs, panels=panels)
        lc.dump(lc.DATA / locus / 'panels.json', summary)
        print(locus, {k: v for k, v in summary.items() if k in ('bubbles', 'kept', 'excluded', 'holdout_samples', 'present_sources')})
        for pnl in panels:
            print(' ', pnl['panel'], {k: v for k, v in pnl.items() if k not in ('haplotypes', 'panel')}, flush=True)


if __name__ == '__main__':
    main()
