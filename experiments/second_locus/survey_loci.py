"""Survey candidate loci: how many of the 464 haplotypes thread cleanly, and why others fail.

Candidate SV-rich, medically relevant loci are given in GRCh38 coordinates (well-known
gene-cluster spans, stated below) and lifted to CHM13v2.0 through GRCh38's own bubble calls
(the GRCh38 haplotype is one of the 464; its call at every bubble gives GRCh38 contig
coordinates). Typical regions are 5 Mb CHM13 windows on chr20 and chr22 with genome-typical
bubble density and little segmental duplication. The MHC window of the earlier experiment is
included as a check that this code reproduces build_panel.py (463 kept, NA20129.pat missing).

For every window the bubbles fully inside it are threaded with build_panel.py's rule
(locus_common.thread). A signature check prints the largest allele-length differences, to
confirm from the data that the lifted window contains the expected structural variation
(for example the common ~60 kb RHD deletion or KIR gene-content haplotypes).
Output: survey.json (no files outside this folder are written).
"""
from __future__ import annotations
import argparse
import collections
import json
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import locus_common as lc  # noqa: E402

# name: (chrom, GRCh38 start, GRCh38 end, note); spans are approximate gene-cluster coordinates
GRCH38_LOCI = {
    'MHC_check': ('chr6', 28510120, 33480577, 'GRCh38 MHC alt-locus span; the earlier experiment used CHM13 chr6:28,381,554-33,301,857'),
    'KIR_LRC': ('chr19', 54025634, 55084318, 'GRCh38 LRC/KIR alt-locus span (LILR and KIR clusters)'),
    'RHD_RHCE': ('chr1', 25250000, 25450000, 'RHD, TMEM50A, RHCE'),
    'AMY': ('chr1', 103500000, 103800000, 'AMY2B, AMY2A, AMY1A-C'),
    'CYP2D6': ('chr22', 42100000, 42160000, 'CYP2D6, CYP2D7, CYP2D8P'),
    'IGH': ('chr14', 105586437, 106879844, 'IGH constant, D, J and V genes'),
    'IGK': ('chr2', 88857361, 90235368, 'IGK proximal and distal V clusters'),
    'IGL': ('chr22', 22026076, 22922913, 'IGL V, J and C genes'),
}
# extra flanked windows around small loci, as (name, chrom, CHM13 start, CHM13 end) after lifting
FLANKS = {'RHD_RHCE': 1_000_000, 'AMY': 1_000_000, 'CYP2D6': 1_000_000, 'KIR_LRC': 500_000}
TYPICAL = [('chr20', 33_000_000, 38_000_000), ('chr20', 44_000_000, 49_000_000), ('chr20', 50_000_000, 55_000_000),
           ('chr20', 38_000_000, 43_000_000), ('chr22', 26_000_000, 31_000_000), ('chr22', 30_000_000, 35_000_000),
           ('chr22', 33_000_000, 38_000_000)]


def window_stats(samples, bubbles, walks, status):
    out = lc.summarize_status(status)
    kept_idx = [h for h, s in enumerate(samples) if status[s]['kept']]
    called = [len({b['calls'][h][0] for h in kept_idx}) for b in bubbles] if kept_idx else []
    spread = sorted(((max(b['alen']) - min(b['alen']) if b['alen'] else 0), b['beg'], b['end'], len(b['alen']))
                    for b in bubbles)[::-1]
    out.update(bubbles=len(bubbles), span=[bubbles[0]['beg'], bubbles[-1]['end']] if bubbles else None,
               bubbles_per_mb=round(len(bubbles) / max(1e-9, (bubbles[-1]['end'] - bubbles[0]['beg']) / 1e6), 1) if bubbles else 0,
               alleles_called_per_bubble_kept=round(statistics.mean(called), 2) if called else None,
               bubbles_with_allele_spread_ge_10kb=sum(1 for s in spread if s[0] >= 10_000),
               sv_bases_max_allele=sum(max(b['alen']) for b in bubbles if b['alen']),
               top_allele_spreads=[dict(spread=s[0], beg=s[1], end=s[2], alleles=s[3]) for s in spread[:5]],
               anno=dict(collections.Counter(b['anno'] or '-' for b in bubbles).most_common(6)))
    if walks:
        lens = sorted(len(w) for w in walks.values())
        out['walk_steps_min_max'] = [lens[0], lens[-1]]
    out['excluded_samples'] = {s: dict(reason=v['reason'], missing_calls=v['missing_calls'],
                                       reverse_alleles=v['reverse_alleles'], contigs=v['contigs'],
                                       strands=v['strands'], monotone=v['monotone'])
                               for s, v in status.items() if not v['kept']}
    return out


def survey(chrom, start, end):
    samples, bubbles = lc.read_calls(chrom, start, end)
    if not bubbles:
        return dict(bubbles=0)
    breaks = lc.chain_breaks(bubbles)
    if breaks:
        return dict(bubbles=len(bubbles), chain_breaks=len(breaks))
    walks, status = lc.thread(samples, bubbles)
    return window_stats(samples, bubbles, walks, status)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, default=HERE / 'survey.json')
    args = ap.parse_args()
    if args.out.exists():
        raise SystemExit('refusing to overwrite ' + str(args.out))
    results, tables = [], {}
    for name, (chrom, g0, g1, note) in GRCH38_LOCI.items():
        if chrom not in tables:
            rows = lc.load_bubble_table(chrom)
            tables[chrom] = lc.grch38_points(rows, chrom)
        pts = tables[chrom]
        c0, c1 = lc.lift(pts, g0), lc.lift(pts, g1)
        windows = [('core', c0, c1)]
        if name in FLANKS:
            f = FLANKS[name]
            windows.append((f'core+-{f // 1000}kb', c0 - f, c1 + f))
        for label, s, e in windows:
            res = dict(name=name, kind='medical' if name != 'MHC_check' else 'check', window=label, chrom=chrom,
                       grch38=[g0, g1], chm13=[s, e], note=note, **survey(chrom, s, e))
            results.append(res)
            print(json.dumps({k: v for k, v in res.items() if k not in ('excluded_samples', 'top_allele_spreads', 'anno')}), flush=True)
    for chrom, s, e in TYPICAL:
        res = dict(name=f'{chrom}_{s // 1_000_000}_{e // 1_000_000}Mb', kind='typical', window='5Mb', chrom=chrom,
                   chm13=[s, e], **survey(chrom, s, e))
        results.append(res)
        print(json.dumps({k: v for k, v in res.items() if k not in ('excluded_samples', 'top_allele_spreads', 'anno')}), flush=True)
    lc.dump(args.out, dict(script='experiments/second_locus/survey_loci.py', results=results))


if __name__ == '__main__':
    main()
