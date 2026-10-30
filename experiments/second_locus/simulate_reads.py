"""Simulate HiFi-like reads from panel haplotypes (standard library only, seeded).

Haplotype sequence: the concatenated segment sequences of the haplotype's walk in the
locus's all-haplotype panel graph (graph-implied sequence: SV-level alleles on the CHM13
backbone; small variants are not represented, as in the graph itself).

Model, per read (one random.Random(seed) stream; draws in this order):
  source     reads alternate between the 'absent' group (held-out samples, both haplotypes:
             absent from the 60/120/240 subsets and the leave-out panel, present only in the
             all-haplotype panel) and the 'present' group (CHM13 and haplotypes in all three
             subsets: present in every panel); the haplotype is uniform within its group
  length     lognormal: exp(N(ln 18000, 0.15^2)), redrawn until within [10 kb, 30 kb]
             (real CHM13 MHC HiFi reads: median 20.2 kb, log-sd 0.134)
  start      uniform over the haplotype sequence (reads lie fully inside the region)
  strand     '+' or '-' with probability 1/2; '-' reads are reverse-complemented
  errors     per-read rate e ~ U(0.001, 0.003) per template base (real reads: QV-implied 0.22%);
             each event is a substitution (1/2, uniform other base), an insertion (1/4) or a
             deletion (1/4); indel length 1 (p 0.7) or 2 (p 0.3); an inserted base repeats the
             previous base with probability 1/2 (homopolymer-like), else is uniform
Primary instances: the first 20 reads in both orientations (40 instances: 10 absent-source
and 10 present-source reads), as the MHC pilot used the first 20 sampled reads.
Outputs: data/<name>/reads.fq (constant quality), reads.tsv (truth), simulation.json.
"""
from __future__ import annotations
import argparse
import gzip
import json
import math
import random
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import locus_common as lc  # noqa: E402

RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
BASES = 'ACGT'
MODEL = dict(length_median=18000, length_log_sd=0.15, length_min=10000, length_max=30000,
             error_min=0.001, error_max=0.003, p_sub=0.5, p_ins=0.25, p_del=0.25,
             indel_len_probs={1: 0.7, 2: 0.3}, p_homopolymer_insert=0.5, primary_reads=20)


def walk_sequences(gfa: Path, wanted):
    seqs, walks = {}, {}
    with gzip.open(gfa, 'rt') as f:
        for line in f:
            a = line.rstrip('\n').split('\t')
            if a[0] == 'S':
                seqs[a[1]] = a[2]
            elif a[0] == 'W' and a[1] in wanted:
                walks[a[1]] = re.findall(r'>([^><]+)', a[6])
    assert set(walks) == set(wanted), set(wanted) - set(walks)
    return {h: ''.join(seqs[s] for s in w).upper() for h, w in walks.items()}


def mutate(seq: str, e: float, rng: random.Random):
    """Apply substitution / insertion / deletion events at rate e per template base."""
    out, i, n, counts = [], 0, len(seq), dict(sub=0, ins=0, dele=0)
    log1me = math.log1p(-e)
    while True:
        gap = int(math.log(1.0 - rng.random()) / log1me)  # geometric number of error-free bases
        if i + gap >= n:
            out.append(seq[i:])
            break
        out.append(seq[i:i + gap])
        i += gap
        u = rng.random()
        k = 1 if rng.random() < MODEL['indel_len_probs'][1] else 2
        if u < MODEL['p_sub']:
            out.append(rng.choice([b for b in BASES if b != seq[i]]))
            i += 1
            counts['sub'] += 1
        elif u < MODEL['p_sub'] + MODEL['p_ins']:
            out.append(seq[i])
            for _ in range(k):
                out.append(seq[i] if rng.random() < MODEL['p_homopolymer_insert'] else rng.choice(BASES))
            i += 1
            counts['ins'] += k
        else:
            i += k
            counts['dele'] += k
    return ''.join(out), counts


def simulate(name, gfa, groups, count, seed, prefix):
    out = lc.DATA / name
    out.mkdir(parents=True, exist_ok=True)
    if (out / 'reads.fq').exists():
        raise SystemExit(f'refusing to overwrite data/{name}/reads.fq')
    wanted = sorted({h for hs in groups.values() for h in hs})
    seqs = walk_sequences(gfa, wanted)
    rng = random.Random(seed)
    order = list(groups)
    rows, mu, sd = [], math.log(MODEL['length_median']), MODEL['length_log_sd']
    with (out / 'reads.fq').open('w') as fq:
        for r in range(count):
            group = order[r % len(order)]
            src = rng.choice(groups[group])
            seq = seqs[src]
            while True:
                length = round(math.exp(rng.gauss(mu, sd)))
                if MODEL['length_min'] <= length <= min(MODEL['length_max'], len(seq)):
                    break
            start = rng.randrange(0, len(seq) - length + 1)
            frag = seq[start:start + length]
            strand = rng.choice('+-')
            if strand == '-':
                frag = frag.translate(RC)[::-1]
            e = rng.uniform(MODEL['error_min'], MODEL['error_max'])
            read, counts = mutate(frag, e, rng)
            rid = f'{prefix}{r:04d}'
            fq.write(f'@{rid} src={src} group={group} start={start} end={start + length} strand={strand} err={e:.5f}\n'
                     f'{read}\n+\n{"I" * len(read)}\n')
            rows.append(dict(read=rid, source=src, group=group, start=start, end=start + length, strand=strand,
                             template_length=length, read_length=len(read), error_rate=round(e, 6), **counts,
                             primary=r < MODEL['primary_reads']))
    with (out / 'reads.tsv').open('w') as f:
        keys = list(rows[0])
        f.write('\t'.join(keys) + '\n')
        for row in rows:
            f.write('\t'.join(str(row[k]) for k in keys) + '\n')
    meta = dict(script='experiments/second_locus/simulate_reads.py', name=name, graph=str(gfa.relative_to(lc.ROOT)),
                graph_sha256=lc.digest(gfa), seed=seed, reads=count, model=MODEL, groups=groups,
                source_lengths={h: len(s) for h, s in seqs.items()},
                group_counts={g: sum(r['group'] == g for r in rows) for g in groups},
                primary_ids=[f"{r['read']}:{s}" for r in rows if r['primary'] for s in '+-'],
                observed=dict(template_length_median=sorted(r['template_length'] for r in rows)[count // 2],
                              mean_error_events_per_kb=round(1000 * sum(r['sub'] + r['ins'] + r['dele'] for r in rows)
                                                             / sum(r['template_length'] for r in rows), 3)),
                reads_sha256=lc.digest(out / 'reads.fq'))
    lc.dump(out / 'simulation.json', meta)
    print(json.dumps({k: v for k, v in meta.items() if k not in ('primary_ids', 'model')}))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--loci', default='kir,chr20')
    ap.add_argument('--reads', type=int, default=100)
    ap.add_argument('--seed', type=int, default=7000, help='locus k (0-based in --loci) uses seed + k; the MHC check uses seed + 100')
    ap.add_argument('--mhc-check', type=int, default=100, help='reads from CHM13 MHC for the realism check (0 to skip)')
    args = ap.parse_args()
    for k, locus in enumerate(args.loci.split(',')):
        meta = json.loads((lc.DATA / locus / 'panels.json').read_text())
        allp = meta['panels'][0]
        groups = dict(absent=meta['holdout_haplotypes'], present=meta['present_sources'])
        simulate(locus, HERE / allp['gfa'], groups, args.reads, args.seed + k, f'{locus}_r')
    if args.mhc_check:
        simulate('mhc_sim', lc.MHC / 'data/hprc463/graph.gfa.gz', dict(present=['CHM13#0']), args.mhc_check,
                 args.seed + 100, 'mhcsim_r')


if __name__ == '__main__':
    main()
