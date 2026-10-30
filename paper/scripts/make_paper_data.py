#!/usr/bin/env python3
"""Generate every number, table row and plot table that the paper reads from saved results.

Inputs are saved result files under experiments/; nothing is re-run. Outputs go to
paper/generated/ and are overwritten on each run:

  numbers.tex           \\dataset{key}{value} definitions; the text uses \\data{key}
  rows_identity.tex     rows of Table 1 (identity with Minichain's DP tables)
  rows_loci.tex         rows of the original 40-instance class table (no longer printed)
  rows_exactness.tex    rows of Table 3 (single-record heuristic vs exact optimum)
  rows_supp_speed.tex   full rows of the supplementary per-panel speed table
  rows_supp_times.tex, rows_supp_ratios.tex  its two readable table panels
  rows_supp_counters.tex rows of the supplementary work-counter table
  rows_supp_guard.tex   paired guarded-engine panel measurements (Table S4, upper panel)
  rows_supp_adversarial.tex  adversarial instances under the guard (Table S4, lower panel)
  rows_b100_main.tex    rows of Table 2 (guarded engine vs per-haplotype, 100 reads, three loci)
  rows_b100_times.tex, rows_b100_ratios.tex  the two panels of Table S1 (every panel, 100 reads)
  rows_b100_counters.tex rows of Table S2 (work counters, 100 reads)
  rows_b100_stages.tex  rows of Table S5 (stage breakdown behind the loss at small H)
  rows_memory.tex       native heap per engine, index and output sizes (experiments/memory_ablation)
  rows_groupings.tex    states, keys, score classes and origin records per anchor (experiments/memory_ablation)
  b100_speedup_<locus>.dat  guarded and compact per-read speedups per panel size (Figure 1b)
  time_vs_h.dat         matched timing against Minichain, medians per read (Figure 1a)
  speedup_<locus>.dat   original 40-instance class speedups per panel size (no longer plotted)
  scatter_<locus>.dat   original 40-instance speedups and records per anchor (no longer plotted)
  b100_scatter_<locus>.dat, b100_scatter_clipped.dat  per-read guarded speedup against records per
                        anchor at the largest panels, and the reads beyond the x axis (Figure S1)
  provenance.json       SHA-256 of each input and of this script

The 100-read benchmark's speedups come from its idle-host rerun; the per-read analysis behind
Figure S1 keeps the first run, on which it was specified (BENCH_RUNS).

Summary values are cross-checked against the per-instance rows they came from, and the
identity summaries must report zero mismatches. Values that exist only in report text,
not in a saved data file, are listed in TEXT_SOURCED with their location.

Run from the repository root:  python3.13 paper/scripts/make_paper_data.py
Python 3.10+, standard library only, assertions enabled.
"""
from __future__ import annotations

import hashlib
import json
import math
import platform
import random
import re
import statistics
import sys
from datetime import datetime, timezone
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

if not __debug__:
    sys.exit('make_paper_data.py: run with assertions enabled (no -O or PYTHONOPTIMIZE)')

ROOT = Path(__file__).resolve().parents[2]
EXP = ROOT / 'experiments'
OUT = ROOT / 'paper' / 'generated'
GAMMA = 10000           # primary switch penalty (Minichain -R10000 with weights x200)
GAMMA_ALT = 20000       # sensitivity setting

# locus -> (summary JSON, per-instance JSONL, panel sizes in order, leave-out panel size)
LOCI = {
    'mhc': (EXP / 'hprc_mhc_panel/benchmark_v2_summary.json',
            EXP / 'hprc_mhc_panel/benchmark_v2.jsonl', (60, 120, 240, 463), 462),
    'kir': (EXP / 'second_locus/benchmark_kir_summary.json',
            EXP / 'second_locus/benchmark_kir.jsonl', (60, 120, 240, 459), 453),
    'chr20': (EXP / 'second_locus/benchmark_chr20_summary.json',
              EXP / 'second_locus/benchmark_chr20.jsonl', (60, 120, 240, 463), 457),
}
LOCUS_NAME = {'mhc': 'MHC', 'kir': 'LRC/KIR', 'chr20': 'chr20'}
LOCUS_READS = {'mhc': 'real', 'kir': 'simulated', 'chr20': 'simulated'}
MATCHED = EXP / 'minichain_matched/results'
HPRC_PANELS = (60, 120, 240, 463)
GUARD = EXP / 'gap_safeguard_repair'
FALLBACK = GUARD / 'real_validation.json'
EXACT = EXP / 'exactness_payoff'

# Numbers that exist only in report text (computed in an earlier session, not saved as data).
TEXT_SOURCED = {
    'exact.sharedseq.lo': ('90', 'experiments/exactness_payoff/RESULTS.md, "Answer in brief"'),
    'exact.sharedseq.hi': ('95', 'experiments/exactness_payoff/RESULTS.md, "Answer in brief"'),
    'baseline.grouping.gain': ('23', 'experiments/hprc_mhc_panel/REPORT.md, "Correction to the 2026-10-04 numbers"'),
    'proofreview.instances': ('1.17', 'experiments/proof_review/REPORT.md, "Test counts" (1,174,300 instances)'),
    'kir.region.mb': ('1.04', 'experiments/second_locus/RESULTS.md section 1 (bubble chain 57,130,780-58,173,439)'),
    'chr20.region.mb': ('4.84', 'experiments/second_locus/RESULTS.md section 1 (bubble chain 50,154,364-54,992,990)'),
    'kir.bubbles.permb': ('157', 'experiments/second_locus/RESULTS.md section 1 (MHC 58, genome median 38)'),
    'mhc.bubbles.permb': ('58', 'experiments/second_locus/RESULTS.md section 1'),
    'genome.bubbles.permb': ('38', 'experiments/second_locus/RESULTS.md section 1 (median over 1-Mb bins)'),
}
TIMING_META = [EXP / 'hprc_mhc_panel/benchmark_v2.meta.json',
               EXP / 'second_locus/benchmark_kir.meta.json',
               EXP / 'second_locus/benchmark_chr20.meta.json',
               *[MATCHED / f'bench_{label}_R10000.meta.json'
                 for label in ('mhc2104', 'hprc60', 'hprc120', 'hprc240', 'hprc463')]]

inputs_used: set[Path] = set()
numbers: dict[str, tuple[str, str]] = {}
KEY_RE = re.compile(r'^[A-Za-z0-9.\-]+$')


def load(path: Path):
    inputs_used.add(path)
    with open(path) as fh:
        return json.load(fh)


def load_jsonl(path: Path):
    inputs_used.add(path)
    with open(path) as fh:
        return [json.loads(line) for line in fh if line.strip()]


def put(key: str, value, source: str) -> None:
    assert KEY_RE.match(key), key
    assert key not in numbers, f'duplicate key {key}'
    numbers[key] = (str(value), source)


def _half_up(x: float, places: str) -> str:
    # Round half up from the decimal repr, as the reports did by hand (1.315 -> 1.32).
    return str(Decimal(repr(x)).quantize(Decimal(places), rounding=ROUND_HALF_UP))


def f2(x: float) -> str:
    return _half_up(x, '0.01')


def f1(x: float) -> str:
    return _half_up(x, '0.1')


def fint(n: int) -> str:
    return f'{n:,}'


def fmega(n: int) -> str:
    """Millions with one decimal, e.g. 986545746 -> 986.5."""
    return f'{n / 1e6:.1f}'


def pct(a: float, b: float) -> str:
    return f'{round(100 * a / b)}'


def rng3(xs):
    return [round(min(xs), 3), round(statistics.median(xs), 3), round(max(xs), 3)]


def avg_ranks(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    ranks = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        for k in range(i, j + 1):
            ranks[order[k]] = (i + j) / 2
        i = j + 1
    return ranks


def spearman(a, b) -> float:
    ra, rb = avg_ranks(a), avg_ranks(b)
    ma, mb = statistics.fmean(ra), statistics.fmean(rb)
    cov = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    va = sum((x - ma) ** 2 for x in ra)
    vb = sum((y - mb) ** 2 for y in rb)
    return cov / (va * vb) ** 0.5


# ---------------------------------------------------------------------------------------------
# Identity with Minichain's DP tables (Table 1)
# ---------------------------------------------------------------------------------------------
def identity():
    rows = []

    def check(summary, label):
        assert summary['mismatch_count'] == 0 and summary['mismatches'] == [], label
        pairs = summary['totals']['pairs']
        assert set(summary['engines']) == {'per_haplotype', 'class'}, label
        for engine in ('per_haplotype', 'class'):
            for field, count in summary['match'][engine].items():
                assert count == pairs, (label, engine, field)
        return summary['totals']

    mhc = load(MATCHED / 'identity_mhc2104_R10000_minichain.summary.json')
    assert mhc['distance'] == 'minichain'
    t = check(mhc, 'mhc2104')
    put('id.mhc.reads', fint(t['reads']), 'identity_mhc2104 totals.reads')
    put('id.mhc.pairs', fint(t['pairs']), 'identity_mhc2104 totals.pairs')
    put('id.mhc.anchors', fmega(t['anchors']), 'identity_mhc2104 totals.anchors (millions)')
    put('id.mhc.states', fmega(t['states']), 'identity_mhc2104 totals.states (millions)')
    rows.append(('Minichain benchmark (MHC graph)', '60', fint(t['reads']), fint(t['pairs']),
                 fmega(t['anchors']), fmega(t['states'])))
    total_states, total_pairs = t['states'], t['pairs']
    hprc_states = hprc_pairs = 0
    for H in HPRC_PANELS:
        s = load(MATCHED / f'identity_hprc{H}_R10000_minichain.summary.json')
        assert s['distance'] == 'minichain'
        t = check(s, f'hprc{H}')
        rows.append(('HPRC release 2, MHC', str(H), fint(t['reads']), fint(t['pairs']),
                     fmega(t['anchors']), fmega(t['states'])))
        hprc_states += t['states']
        hprc_pairs += t['pairs']
        put(f'id.hprc{H}.states', fmega(t['states']), f'identity_hprc{H} totals.states (millions)')
    total_states += hprc_states
    total_pairs += hprc_pairs
    put('id.hprc.states', f'{hprc_states / 1e6:.0f}', 'sum of identity_hprc* totals.states (millions)')
    put('id.hprc.pairs', fint(hprc_pairs), 'sum of identity_hprc* totals.pairs')
    put('id.total.states', f'{total_states / 1e9:.2f}', 'identity: MHC 2,104 reads + HPRC panels, states (billions)')
    put('id.total.pairs', fint(total_pairs), 'identity: MHC 2,104 reads + HPRC panels, pairs')
    rows.append(('Total', '', '', fint(total_pairs), '', f'{total_states / 1e6:,.1f}'))

    alg2 = load(MATCHED / 'identity_mhc2104_R10000_truedist_vs_alg2.summary.json')
    assert alg2['distance'] == 'true' and alg2['index_header']['truedist'] == '1'
    t = check(alg2, 'truedist_vs_alg2')
    put('id.alg2.pairs', fint(t['pairs']), 'identity_mhc2104_truedist_vs_alg2 totals.pairs')
    put('id.alg2.states', fmega(t['states']), 'identity_mhc2104_truedist_vs_alg2 totals.states (millions)')

    dev = load(MATCHED / 'variant_truedist_R10000.json')
    c = dev['counts']
    assert c['best_changed'] == 0 and dev['reads_best_changed'] == 0
    put('dev.dist.tables', fint(c['hc_changed']), 'variant_truedist counts.hc_changed')
    put('dev.dist.tables.pct', pct(c['hc_changed'], c['pairs']), 'variant_truedist hc_changed / pairs (%)')
    put('dev.dist.best', fint(c['best_changed']), 'variant_truedist counts.best_changed')
    gate = load(MATCHED / 'variant_gate_R10000.json')
    assert gate['counts']['best_changed'] == 0
    put('dev.gate.tables', fint(gate['counts']['hc_changed']), 'variant_gate counts.hc_changed')

    with open(OUT / 'rows_identity.tex', 'w') as fh:
        fh.write('% Generated by paper/scripts/make_paper_data.py; do not edit.\n')
        fh.write('% Data & H & Reads & (read, component) pairs & Anchors (M) & States compared (M)\n')
        for i, r in enumerate(rows):
            if r[0] == 'Total':
                fh.write('\\midrule\n')
                fh.write(f'Total & & & {r[3]} & & {r[5]} \\\\\n')
            else:
                fh.write(' & '.join(r) + ' \\\\\n')


# ---------------------------------------------------------------------------------------------
# Matched timing against Minichain (Figure 1a, text)
# ---------------------------------------------------------------------------------------------
TIMING_MODES = ('minichain_chaining', 'per_haplotype', 'class_expanded', 'class_compact')


def matched_timing():
    mhc = load(MATCHED / 'timing_mhc2104_R10000.json')
    assert mhc['reads'] == 2104 and mhc['unmatched_reads'] == 0
    for m in TIMING_MODES:
        assert mhc['per_read_ms'][m]['n'] == 2104
    tot = mhc['totals_seconds']
    put('mc.mhc.total.minichain', f'{tot["minichain_chaining"]:.0f}', 'timing_mhc2104 totals_seconds.minichain_chaining (s)')
    put('mc.mhc.total.perhap', f1(tot['per_haplotype']), 'timing_mhc2104 totals_seconds.per_haplotype (s)')
    put('mc.mhc.total.compact', f1(tot['class_compact']), 'timing_mhc2104 totals_seconds.class_compact (s)')
    put('mc.mhc.total.expanded', f1(tot['class_expanded']), 'timing_mhc2104 totals_seconds.class_expanded (s)')
    pr = mhc['per_read_ms']
    put('mc.mhc.med.minichain', f1(pr['minichain_chaining']['median']), 'timing_mhc2104 per_read_ms.minichain_chaining.median')
    put('mc.mhc.med.perhap', f2(pr['per_haplotype']['median']), 'timing_mhc2104 per_read_ms.per_haplotype.median')
    put('mc.mhc.med.compact', f2(pr['class_compact']['median']), 'timing_mhc2104 per_read_ms.class_compact.median')
    rr = mhc['paired_ratio']
    for key, name in (('minichain_chaining/per_haplotype', 'perhap'), ('minichain_chaining/class_compact', 'compact'),
                      ('minichain_chaining/class_expanded', 'expanded')):
        put(f'mc.mhc.x.{name}', f'{rr[key]["median"]:.0f}', f'timing_mhc2104 paired_ratio.{key}.median')
        put(f'mc.mhc.x.{name}.p10', f'{rr[key]["p10"]:.0f}', f'timing_mhc2104 paired_ratio.{key}.p10')
        put(f'mc.mhc.x.{name}.p90', f'{rr[key]["p90"]:.0f}', f'timing_mhc2104 paired_ratio.{key}.p90')
    put('mc.mhc.cls.compact', f2(rr['per_haplotype/class_compact']['median']), 'timing_mhc2104 paired_ratio.per_haplotype/class_compact.median')
    put('mc.mhc.cls.expanded', f2(rr['per_haplotype/class_expanded']['median']), 'timing_mhc2104 paired_ratio.per_haplotype/class_expanded.median')

    lines = ['H minichain minichain_p10 minichain_p90 perhap perhap_p10 perhap_p90 '
             'expanded expanded_p10 expanded_p90 compact compact_p10 compact_p90']
    medians_by_h = {}
    for H in HPRC_PANELS:
        t = load(MATCHED / f'timing_hprc{H}_R10000.json')
        assert t['reads'] == 100 and t['unmatched_reads'] == 0
        pr, rr = t['per_read_ms'], t['paired_ratio']
        for m in TIMING_MODES:
            assert pr[m]['n'] == 100
        cols = [str(H)]
        for m in TIMING_MODES:
            cols += [f'{pr[m]["median"]:.4f}', f'{pr[m]["p10"]:.4f}', f'{pr[m]["p90"]:.4f}']
        lines.append(' '.join(cols))
        put(f'mc.h{H}.med.minichain', fint(round(pr['minichain_chaining']['median'])), f'timing_hprc{H} per_read_ms.minichain_chaining.median')
        put(f'mc.h{H}.med.perhap', f2(pr['per_haplotype']['median']), f'timing_hprc{H} per_read_ms.per_haplotype.median')
        put(f'mc.h{H}.med.compact', f2(pr['class_compact']['median']), f'timing_hprc{H} per_read_ms.class_compact.median')
        put(f'mc.h{H}.med.expanded', f2(pr['class_expanded']['median']), f'timing_hprc{H} per_read_ms.class_expanded.median')
        for key, name in (('minichain_chaining/per_haplotype', 'perhap'), ('minichain_chaining/class_compact', 'compact'),
                          ('minichain_chaining/class_expanded', 'expanded')):
            put(f'mc.h{H}.x.{name}', f'{rr[key]["median"]:.0f}', f'timing_hprc{H} paired_ratio.{key}.median')
        put(f'mc.h{H}.cls.compact', f2(rr['per_haplotype/class_compact']['median']), f'timing_hprc{H} paired_ratio.per_haplotype/class_compact.median')
        put(f'mc.h{H}.cls.expanded', f2(rr['per_haplotype/class_expanded']['median']), f'timing_hprc{H} paired_ratio.per_haplotype/class_expanded.median')
        medians_by_h[H] = {m: pr[m]['median'] for m in TIMING_MODES}
    with open(OUT / 'time_vs_h.dat', 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    lo, hi = HPRC_PANELS[0], HPRC_PANELS[-1]
    put('mc.growth.h', f1(hi / lo), f'HPRC panel sizes {hi}/{lo}')
    for m, name in (('minichain_chaining', 'minichain'), ('per_haplotype', 'perhap'), ('class_compact', 'compact')):
        put(f'mc.growth.{name}', f1(medians_by_h[hi][m] / medians_by_h[lo][m]),
            f'timing_hprc{hi} / timing_hprc{lo} per_read_ms.{m}.median')


# ---------------------------------------------------------------------------------------------
# Class vs per-haplotype on three loci (Figure 1b, 1c, Table 2, supplement)
# ---------------------------------------------------------------------------------------------
MODES = ('class_compact', 'class_expanded')


def recompute(rows, H, gamma):
    rs = [r for r in rows if r['H'] == H and r['gamma'] == gamma]
    assert len(rs) == 40, (H, gamma, len(rs))
    ph = [r['kernel_median']['per_haplotype'] for r in rs]
    out = {'rows': rs}
    for m in MODES:
        ratios = [a / r['kernel_median'][m] for a, r in zip(ph, rs)]
        out[m] = (rng3(ratios), sum(x > 1 for x in ratios), ratios)
    return out


def loci():
    loci_rows, supp_rows, counter_rows = [], [], []
    pooled_recs, pooled_cmp = [], []
    for locus, (summary_path, rows_path, sizes, leaveout) in LOCI.items():
        summary = load(summary_path)['levels']
        rows = load_jsonl(rows_path)
        dat = ['H compact compact_min compact_max expanded expanded_min expanded_max']
        for H in sizes + (leaveout,):
            for gamma in (GAMMA, GAMMA_ALT):
                lv = summary[f'H{H}_gamma{gamma}']
                rc = recompute(rows, H, gamma)
                assert lv['reads'] == 40
                assert lv['states'] == sum(r['K'] for r in rc['rows'])
                assert lv['anchors'] == sum(r['N'] for r in rc['rows'])
                for m in MODES:
                    assert rc[m][0] == lv['kernel_speedup'][m], (locus, H, gamma, m, rc[m][0], lv['kernel_speedup'][m])
                    assert rc[m][1] == lv['kernel_wins'][m], (locus, H, gamma, m)
                c, e = lv['kernel_speedup']['class_compact'], lv['kernel_speedup']['class_expanded']
                wc, we = lv['kernel_wins']['class_compact'], lv['kernel_wins']['class_expanded']
                km = lv['kernel_ms_median']
                spa = lv['states'] / lv['anchors']
                cpa = lv['counters_per_anchor']
                label = f'{H}' + (r'$^{\ast}$' if H == leaveout else '')
                if gamma == GAMMA:
                    supp_rows.append((LOCUS_NAME[locus], label, f'{spa:.0f}', f2(km['per_haplotype']),
                                      f2(km['class_compact']), f2(km['class_expanded']),
                                      f'{f2(c[1])} [{f2(c[0])}, {f2(c[2])}]', f'{f2(e[1])} [{f2(e[0])}, {f2(e[2])}]',
                                      f'{wc}/{we}'))
                    counter_rows.append((LOCUS_NAME[locus], label,
                                         f'{lv["pair_decision_shares"]["pairs_same_vertex"][1]:.3f}',
                                         f'{lv["pair_decision_shares"]["pairs_loose"][1]:.3f}',
                                         f'{lv["pair_decision_shares"]["pairs_range"][1]:.3f}',
                                         f'{lv["pair_decision_shares"]["pairs_sliced"][1]:.3f}',
                                         f'{cpa["pairs_sliced"][1]:.3f}', f'{cpa["output_records"][1]:.3f}'))
                    if H != leaveout:
                        dat.append(f'{H} {c[1]} {c[0]} {c[2]} {e[1]} {e[0]} {e[2]}')
                    key = f'{locus}.h{H}'
                    put(f'{key}.x.compact', f2(c[1]), f'{summary_path.name} H{H}_gamma{gamma}.kernel_speedup.class_compact[median]')
                    put(f'{key}.x.expanded', f2(e[1]), f'{summary_path.name} H{H}_gamma{gamma}.kernel_speedup.class_expanded[median]')
                    put(f'{key}.x.compact.min', f2(c[0]), f'{summary_path.name} H{H}_gamma{gamma}.kernel_speedup.class_compact[min]')
                    put(f'{key}.x.compact.max', f2(c[2]), f'{summary_path.name} H{H}_gamma{gamma}.kernel_speedup.class_compact[max]')
                    put(f'{key}.x.expanded.min', f2(e[0]), f'{summary_path.name} H{H}_gamma{gamma}.kernel_speedup.class_expanded[min]')
                    put(f'{key}.x.expanded.max', f2(e[2]), f'{summary_path.name} H{H}_gamma{gamma}.kernel_speedup.class_expanded[max]')
                    put(f'{key}.wins.compact', str(wc), f'{summary_path.name} H{H}_gamma{gamma}.kernel_wins.class_compact')
                    put(f'{key}.wins.expanded', str(we), f'{summary_path.name} H{H}_gamma{gamma}.kernel_wins.class_expanded')
                    put(f'{key}.ms.perhap', f2(km['per_haplotype']), f'{summary_path.name} H{H}_gamma{gamma}.kernel_ms_median.per_haplotype')
                    put(f'{key}.ms.compact', f2(km['class_compact']), f'{summary_path.name} H{H}_gamma{gamma}.kernel_ms_median.class_compact')
                    put(f'{key}.ms.expanded', f2(km['class_expanded']), f'{summary_path.name} H{H}_gamma{gamma}.kernel_ms_median.class_expanded')
                    put(f'{key}.sliced', f'{cpa["pairs_sliced"][1]:.3f}', f'{summary_path.name} H{H}_gamma{gamma}.counters_per_anchor.pairs_sliced[median]')
                    put(f'{key}.records', f'{cpa["output_records"][1]:.2f}', f'{summary_path.name} H{H}_gamma{gamma}.counters_per_anchor.output_records[median]')
                    put(f'{key}.statesperanchor', f'{spa:.0f}', f'{summary_path.name} H{H}_gamma{gamma} states/anchors')
                    if H in (sizes[0], sizes[-1]):
                        loci_rows.append((LOCUS_NAME[locus] if H == sizes[0] else '',
                                          LOCUS_READS[locus] if H == sizes[0] else '', str(H), f'{spa:.0f}',
                                          f'{cpa["output_records"][1]:.2f}', f2(km['per_haplotype']),
                                          f2(km['class_compact']), f2(km['class_expanded']),
                                          f2(c[1]), f2(e[1]), f'{wc}/{we}'))
                    if H == sizes[-1]:
                        # Per-instance scatter at the largest all-haplotype panel.
                        recs = [r['stats']['class_compact']['output_records'] / r['N'] for r in rc['rows']]
                        sl = [r['stats']['class_compact']['pairs_sliced'] / r['N'] for r in rc['rows']]
                        cmp_, exp_ = rc['class_compact'][2], rc['class_expanded'][2]
                        with open(OUT / f'scatter_{locus}.dat', 'w') as fh:
                            fh.write('records sliced compact expanded\n')
                            for row in sorted(zip(recs, sl, cmp_, exp_)):
                                fh.write(' '.join(f'{v:.4f}' for v in row) + '\n')
                        pooled_recs += recs
                        pooled_cmp += cmp_
                        put(f'{locus}.hmax', str(H), f'largest all-haplotype panel of {locus}')
                        put(f'{locus}.hmax.x.compact', f2(c[1]), f'{summary_path.name} H{H}_gamma{gamma} compact median')
                        put(f'{locus}.hmax.x.expanded', f2(e[1]), f'{summary_path.name} H{H}_gamma{gamma} expanded median')
                        put(f'{locus}.hmax.wins.compact', str(wc), f'{summary_path.name} H{H}_gamma{gamma} compact wins')
                        put(f'{locus}.hmax.wins.expanded', str(we), f'{summary_path.name} H{H}_gamma{gamma} expanded wins')
                        put(f'{locus}.hmax.sliced', f'{cpa["pairs_sliced"][1]:.3f}', f'{summary_path.name} H{H}_gamma{gamma} pairs_sliced median')
                        put(f'{locus}.hmax.records', f'{cpa["output_records"][1]:.2f}', f'{summary_path.name} H{H}_gamma{gamma} output_records median')
                        if locus == 'mhc':
                            sh = lv['pair_decision_shares']
                            for name, field in (('same', 'pairs_same_vertex'), ('loose', 'pairs_loose'),
                                                ('range', 'pairs_range'), ('sliced', 'pairs_sliced')):
                                put(f'mhc.hmax.share.{name}', f'{100 * sh[field][1]:.1f}', f'{summary_path.name} H{H} pair_decision_shares.{field}[median] (%)')
                            st = lv['stages_ms']['class_expanded']
                            put('mhc.hmax.expand.ms', f2(st['export_seconds']), f'{summary_path.name} H{H} stages_ms.class_expanded.export_seconds')
                            put('mhc.hmax.expand.pct', pct(st['export_seconds'], km['class_expanded']),
                                f'{summary_path.name} H{H} export / kernel_ms_median.class_expanded (%)')
                else:
                    # gamma sensitivity: largest change of the median ratio against GAMMA
                    base = summary[f'H{H}_gamma{GAMMA}']['kernel_speedup']
                    d = max(abs(lv['kernel_speedup'][m][1] - base[m][1]) for m in MODES)
                    numbers.setdefault(f'{locus}.gammadelta', ('0.00', ''))
                    prev = float(numbers[f'{locus}.gammadelta'][0])
                    numbers[f'{locus}.gammadelta'] = (f2(max(prev, d)), f'max |median(gamma={GAMMA_ALT}) - median(gamma={GAMMA})| over panels, {summary_path.name}')
        with open(OUT / f'speedup_{locus}.dat', 'w') as fh:
            fh.write('\n'.join(dat) + '\n')
    gd = max(float(numbers[f'{loc}.gammadelta'][0]) for loc in LOCI)
    put('all.gammadelta', f2(gd), 'max over loci of the gamma sensitivity')
    put('scatter.spearman', f2(spearman(pooled_recs, pooled_cmp)),
        'Spearman correlation, records per anchor vs compact speedup, 120 instances at the largest panels (computed here)')
    put('scatter.n', str(len(pooled_recs)), 'instances in the pooled scatter')

    with open(OUT / 'rows_loci.tex', 'w') as fh:
        fh.write('% Generated by paper/scripts/make_paper_data.py; do not edit.\n')
        fh.write('% Locus & H & records/anchor & per-hap ms & compact ms & expanded ms & x compact & x expanded\n')
        for i, r in enumerate(loci_rows):
            if i and r[0]:
                fh.write('\\addlinespace\n')
            fh.write(' & '.join(r[j] for j in (0, 2, 4, 5, 6, 7, 8, 9)) + ' \\\\\n')
    with open(OUT / 'rows_supp_speed.tex', 'w') as fh:
        fh.write('% Generated by paper/scripts/make_paper_data.py; do not edit.\n')
        prev = None
        for r in supp_rows:
            if prev is not None and r[0] != prev:
                fh.write('\\midrule\n')
            fh.write(' & '.join(r if r[0] != prev else ('',) + r[1:]) + ' \\\\\n')
            prev = r[0]
    # Two readable panels of Table S1 retain the original nine-column data above.
    for filename, columns in [('rows_supp_times.tex', (0, 1, 2, 3, 4, 5)),
                              ('rows_supp_ratios.tex', (0, 1, 6, 7, 8))]:
        with open(OUT / filename, 'w') as fh:
            fh.write('% Generated by paper/scripts/make_paper_data.py; do not edit.\n')
            prev = None
            for r in supp_rows:
                if prev is not None and r[0] != prev:
                    fh.write('\\midrule\n')
                row = r if r[0] != prev else ('',) + r[1:]
                fh.write(' & '.join(row[j] for j in columns) + ' \\\\\n')
                prev = r[0]
    with open(OUT / 'rows_supp_counters.tex', 'w') as fh:
        fh.write('% Generated by paper/scripts/make_paper_data.py; do not edit.\n')
        prev = None
        for r in counter_rows:
            if prev is not None and r[0] != prev:
                fh.write('\\midrule\n')
            fh.write(' & '.join(r if r[0] != prev else ('',) + r[1:]) + ' \\\\\n')
            prev = r[0]


# ---------------------------------------------------------------------------------------------
# Work-budget fallback (auto mode): fallbacks on real panels
# ---------------------------------------------------------------------------------------------
def fallback():
    fb = load(FALLBACK)
    assert fb['status'] == 'all passed' and fb['phase'] == 'validate'
    panels = fb['panels']
    total = sum(p['instances'] for p in panels)
    by_locus = {}
    for p in panels:
        assert p['instances'] == len(p['rows'])
        assert p['states'] == sum(r['K'] for r in p['rows'])
        assert p['fallbacks'] == sum(r['guard']['fallback'] for r in p['rows'])
        for r in p['rows']:
            assert r['guard']['auto_fallback_available']
            assert r['guard']['auto_work'] <= r['guard']['auto_budget']
            assert len(r['digest']) == len(r['pointer_digest']) == 64
        by_locus.setdefault(p['locus'], 0)
        by_locus[p['locus']] += p['fallbacks']
    src = 'gap_safeguard_repair/real_validation.json'
    put('auto.panels', str(len(panels)), src + ' number of panels')
    put('auto.instances', fint(total), src + ' instances over all panels')
    put('auto.fallbacks', str(sum(by_locus.values())), src + ' summed fallbacks')
    for locus in ('mhc', 'kir', 'chr20'):
        put(f'auto.fallbacks.{locus}', str(by_locus[locus]), f'{src} fallbacks summed over {locus} panels')
    put('guard.real.states.million', fmega(sum(p['states'] for p in panels)), src + ' total valid states in millions')


def safeguard_timing():
    """Repaired auto timings, kept separate from the historical class-only experiment."""
    bench = load(GUARD / 'benchmarks.json')
    assert bench['status'] == 'all passed' and bench['phase'] == 'benchmark'
    assert bench['limit'] == 0
    repetitions = bench['repetitions']
    panels = bench['panels']
    put('guard.repetitions', repetitions, 'guard benchmark repetitions per mode/instance')
    put('guard.timed.panels', len(panels), 'guard benchmark panels')
    put('guard.timed.instances', sum(p['instances'] for p in panels), 'guard benchmark panel-instance count (overlapping corpus)')
    put('guard.native.library', Path(bench['libraries']['new']).stem, 'guard benchmark native library, distinct from original timing build')

    def checked_medians(row):
        med = {}
        for mode, ts in row['kernel_seconds'].items():
            assert len(ts) == repetitions and all(t > 0 for t in ts)
            med[mode] = statistics.median(ts)
            assert math.isclose(med[mode], row['medians'][mode], rel_tol=1e-12)
        assert row['guards']['new_auto']['auto_fallback_available']
        assert row['guards']['new_auto']['auto_work'] <= row['guards']['new_auto']['auto_budget']
        assert row['guards']['new_auto']['auto_factor'] == 1.0
        return med

    rows = []
    all_loads = []
    for p in panels:
        locus, H = p['locus'], p['H']
        assert len(p['rows']) == p['instances']
        medians = [checked_medians(r) for r in p['rows']]
        speedups = [m['per_haplotype'] / m['new_auto'] for m in medians]
        overheads = [m['new_auto'] / m['new_class'] for m in medians]
        for key, values in [('ph_over_new_auto', speedups), ('new_auto_over_class', overheads)]:
            assert math.isclose(statistics.median(values), p['ratios'][key]['median'], rel_tol=1e-12)
            for r, value in zip(p['rows'], values):
                assert math.isclose(value, r['ratios'][key], rel_tol=1e-12)
        fallbacks = sum(r['guards']['new_auto']['fallback'] for r in p['rows'])
        assert fallbacks == p['fallbacks']['new_auto']
        src = f'gap_safeguard_repair/benchmarks.json {locus} H={H}'
        prefix = f'guard.{locus}.h{H}'
        put(prefix + '.speedup', f2(statistics.median(speedups)), src + ' median paired per-haplotype / repaired auto')
        put(prefix + '.overhead.pct', f1(100 * (statistics.median(overheads) - 1)), src + ' repaired auto / current class minus one, percent')
        put(prefix + '.fallbacks', fallbacks, src + ' fallback count')
        rows.append((LOCUS_NAME[locus], str(H), str(p['instances']),
                     *[f2(1000 * statistics.median(m[mode] for m in medians))
                       for mode in ('per_haplotype', 'new_class', 'new_auto')],
                     f2(statistics.median(speedups)), str(fallbacks)))
        all_loads.extend([p['load_before'][0], p['load_after'][0]])
    put('guard.load.lo', f1(min(all_loads)), 'minimum one-minute load before/after guard panel timings')
    put('guard.load.hi', f1(max(all_loads)), 'maximum one-minute load before/after guard panel timings')
    with open(OUT / 'rows_supp_guard.tex', 'w') as fh:
        fh.write('% Generated from the repaired safeguard benchmark; expanded output throughout.\n')
        for row in rows:
            fh.write(' & '.join(row) + ' \\\\\n')

    # Adversarial rows are checked against raw repeated timings too, never mixed with panels.
    hostile_ratios, class_ratios, adversarial_rows = [], [], []
    for row in bench['adversarial']:
        med = checked_medians(row)
        assert row['guards']['new_auto']['fallback']
        hostile_ratios.append(med['new_auto'] / med['per_haplotype'])
        class_ratios.append((med['new_class'] / med['per_haplotype'], row['family']))
        adversarial_rows.append((row['family'].capitalize(), fint(row['size']), fint(row['N']), fint(row['H']),
                                 *[f'{1000 * med[mode]:.3f}' for mode in ('per_haplotype', 'new_class', 'new_auto')],
                                 f2(med['new_auto'] / med['per_haplotype'])))
    put('guard.hostile.instances', len(hostile_ratios), 'guard benchmark adversarial rows')
    put('guard.hostile.maxslowdown', f2(max(hostile_ratios)), 'maximum new auto / per-haplotype among saved hostile rows')
    worst, family = max(class_ratios)
    assert family == 'star'  # the text says the class engine's worst case is on the star family
    put('guard.hostile.class.maxslowdown', f'{worst:.0f}', 'maximum unguarded class / per-haplotype among saved hostile rows (star family)')
    with open(OUT / 'rows_supp_adversarial.tex', 'w') as fh:
        fh.write('% Generated from the guard benchmark adversarial rows; every row falls back.\n')
        prev = None
        for r in adversarial_rows:
            if prev is not None and r[0] != prev:
                fh.write('\\midrule\n')
            fh.write(' & '.join(r if r[0] != prev else ('',) + r[1:]) + ' \\\\\n')
            prev = r[0]

    # Source version and counters from the completed native suite, not a fresh test claim.
    manifest = load(GUARD / 'validation_manifest.json')
    run = ROOT / manifest['gap_run']
    run_meta = load(run)
    assert run_meta['status'] == 'passed'
    test = load(run.parent / 'gap_native' / 'verification_native.json')
    put('guard.test.instances', fint(test['cases']), str(run.parent.relative_to(ROOT)) + ' native accepted instances')
    put('guard.test.states', fint(test['states']), str(run.parent.relative_to(ROOT)) + ' native valid states')


# ---------------------------------------------------------------------------------------------
# 100-read guarded benchmark (Table 2, Figure 1b, Sections 6.4-6.6, Tables S1, S2, S5)
# ---------------------------------------------------------------------------------------------
BENCH = EXP / 'benchmark_100reads'
BENCH_PANELS = {  # locus -> all-haplotype panels in order of H, plus the leave-out panel if any
    'mhc': ((('hprc60', 60), ('hprc120', 120), ('hprc240', 240), ('hprc463', 463)), None),
    'kir': ((('hprc60', 60), ('hprc120', 120), ('hprc240', 240), ('hprc459', 459)), ('hprc453_leaveout', 453)),
    'chr20': ((('hprc60', 60), ('hprc120', 120), ('hprc240', 240), ('hprc463', 463)), ('hprc457_leaveout', 457)),
}
BENCH_MODES = ('per_haplotype', 'class', 'auto', 'compact')
BENCH_RESAMPLES = 2000
# Timing runs of the 100-read benchmark, each replayed against the one shared validation. The paper's speedups
# come from the second rerun on an idle host (experiments/benchmark_100reads_rerun, 2026-10-21: same inputs,
# engines, validation and protocol). The per-read analysis of Section 6.5 was specified on the first run and
# keeps it.
BENCH_RUNS = {'first': (BENCH / 'timings', BENCH / 'summary.json'),
              'idle': (EXP / 'benchmark_100reads_rerun/timings2', EXP / 'benchmark_100reads_rerun/timings2-summary.json')}
BENCH_HEADLINE = 'idle'


def _bench_quantile(values, fraction):
    # Nearest rank on the sorted list, as experiments/benchmark_100reads/summarize.py.
    ordered = sorted(values)
    return ordered[round(fraction * (len(ordered) - 1))]


def _bench_bootstrap(pairs, seed):
    """Replays summarize.py's paired whole-read bootstrap for one mode (fresh RNG per mode)."""
    rng = random.Random(seed)
    medians, totals = [], []
    for _ in range(BENCH_RESAMPLES):
        selected = [pairs[rng.randrange(len(pairs))] for _ in pairs]
        medians.append(statistics.median(a / b for a, b in selected))
        totals.append(sum(a for a, b in selected) / sum(b for a, b in selected))
    return ([_bench_quantile(medians, .025), _bench_quantile(medians, .975)],
            [_bench_quantile(totals, .025), _bench_quantile(totals, .975)])


def _close(a, b):
    return math.isclose(a, b, rel_tol=1e-12, abs_tol=0.0)


def bench_panel(locus, name, H, summary_panel, repetitions, timings):
    """One panel: replay the read-level summary from raw rows and derive counters and stages."""
    key = f'{locus}-{name}'
    trows = load_jsonl(timings / f'{key}.jsonl')
    vrows = {r['id']: r for r in load_jsonl(BENCH / 'validation' / f'{key}.jsonl')}
    tmeta, vmeta = load(timings / f'{key}.json'), load(BENCH / 'validation' / f'{key}.json')
    assert tmeta['status'] == vmeta['status'] == 'passed'
    assert tmeta['H'] == vmeta['H'] == summary_panel['H'] == H
    assert tmeta['input_hashes'] == vmeta['input_hashes']
    assert len(trows) == len(vrows) == 200 and {r['id'] for r in trows} == set(vrows)
    groups = {}
    for r in trows:
        v = vrows[r['id']]
        assert (r['N'], r['K']) == (v['N'], v['K'])
        # Timing rows carry the validation run's guard and class statistics.
        assert r['guard'] == v['stats']['auto'] and r['class_stats'] == v['stats']['class']
        assert r['compact_records'] == v['stats']['compact']['compact_records']
        assert r['guard']['auto_fallback_available'] and r['guard']['auto_work'] <= r['guard']['auto_budget']
        for m in BENCH_MODES:
            assert len(r['kernel_seconds'][m]) == repetitions and min(r['kernel_seconds'][m]) > 0
        groups.setdefault(r['read'], []).append(r)
    assert len(groups) == 100 and all(len(c) == 2 for c in groups.values())
    srows = summary_panel['read_rows']
    assert [s['read'] for s in srows] == list(groups)  # same order, so the bootstrap replays exactly
    assert sum(r['K'] for r in trows) == summary_panel['states'] == tmeta['states'] == vmeta['states']
    fallbacks = sum(r['guard']['fallback'] for r in trows)
    assert fallbacks == summary_panel['fallback_components']

    reads = []
    for s in srows:
        comps = groups[s['read']]
        t = {m: statistics.median(sum(c['kernel_seconds'][m][j] for c in comps) for j in range(repetitions))
             for m in BENCH_MODES}
        for m in BENCH_MODES:
            assert _close(t[m], s['kernel_seconds'][m]), (key, s['read'], m)
        N, K = sum(c['N'] for c in comps), sum(c['K'] for c in comps)
        assert (N, K) == (s['N'], s['K']) and N > 0

        def total(field):
            return sum(c['class_stats'][field] for c in comps)

        def stage(mode, field):  # single validation run of each mode, both components
            return sum(vrows[c['id']]['stats'][mode][field] for c in comps)
        reads.append(dict(read=s['read'], prior=s['prior_primary'], t=t, N=N, K=K, comps=comps,
                          fallback=any(c['guard']['fallback'] for c in comps),
                          records=sum(c['compact_records'] for c in comps),
                          sliced=total('pairs_sliced'), pairs=total('pairs'), same=total('pairs_same_vertex'),
                          loose=total('pairs_loose'), range=total('pairs_range'),
                          ph_kernel=stage('per_haplotype', 'kernel_seconds'),
                          class_kernel=stage('class', 'kernel_seconds'),
                          class_export=stage('class', 'export_seconds')))

    out = dict(H=H, states=summary_panel['states'], fallbacks=fallbacks)
    for m in ('class', 'auto', 'compact'):
        pairs = [(r['t']['per_haplotype'], r['t'][m]) for r in reads]
        ratios = [a / b for a, b in pairs]
        ref = summary_panel['kernel'][m]
        med, tot = statistics.median(ratios), sum(a for a, _ in pairs) / sum(b for _, b in pairs)
        assert _close(med, ref['median_speedup']) and _close(tot, ref['total_speedup']), (key, m)
        assert sum(x > 1 for x in ratios) == ref['wins'] and _close(1 / min(ratios), ref['max_slowdown'])
        out[m] = dict(median=med, total=tot, wins=ref['wins'], slowest=1 / min(ratios),
                      ms=1000 * statistics.median(b for _, b in pairs))
    out['perhap_ms'] = 1000 * statistics.median(r['t']['per_haplotype'] for r in reads)
    ci, total_ci = _bench_bootstrap([(r['t']['per_haplotype'], r['t']['auto']) for r in reads],
                                    20261017 + sum(key.encode()))
    ref = summary_panel['kernel']['auto']
    assert all(_close(a, b) for a, b in zip(ci + total_ci, ref['median_ci95'] + ref['total_ci95'])), key
    out['auto'].update(ci=ci, total_ci=total_ci)
    for label, prior in (('previous_20', True), ('additional_80', False)):
        sub = [r['t']['per_haplotype'] / r['t']['auto'] for r in reads if r['prior'] == prior]
        assert len(sub) == (20 if prior else 80)
        assert _close(statistics.median(sub), summary_panel['subgroups'][label]['auto']['median_speedup'])
        out[label] = statistics.median(sub)

    out['records'] = statistics.median(r['records'] / r['N'] for r in reads)
    out['sliced'] = statistics.median(r['sliced'] / r['N'] for r in reads)
    out['statesperanchor'] = statistics.median(r['K'] / r['N'] for r in reads)
    with_pairs = [r for r in reads if r['pairs'] > 0]
    assert len(with_pairs) >= 95, (key, len(with_pairs))
    out['pair_reads'] = len(with_pairs)
    for field in ('same', 'loose', 'range', 'sliced'):
        out['share.' + field] = statistics.median(r[field] / r['pairs'] for r in with_pairs)
    # Stage breakdown. Ratios of a read's class time to its per-haplotype time use one validation
    # run per mode; the guard's premium and every speedup use the five timed repetitions.
    out['nsperstate'] = statistics.median(1e9 * r['t']['per_haplotype'] / r['K'] for r in reads)
    out['recwork'] = statistics.median((r['class_kernel'] - r['class_export']) / r['ph_kernel'] for r in reads)
    out['expand'] = statistics.median(r['class_export'] / r['ph_kernel'] for r in reads)
    out['expand_share'] = statistics.median(r['class_export'] / r['class_kernel'] for r in reads)
    out['premium'] = statistics.median(r['t']['auto'] / r['t']['class'] for r in reads)
    single = statistics.median(r['ph_kernel'] / r['class_kernel'] for r in reads)
    out['stage_agreement'] = abs(single - out['class']['median']) / out['class']['median']
    out['load'] = [tmeta['load_before'][0], tmeta['load_after'][0]]
    out['reads'] = reads  # per-read rows for mechanism()
    return out


def bench_results(run):
    """Replay one timing run of the 100-read benchmark (BENCH_RUNS) against the shared validation."""
    timings, summary_path = BENCH_RUNS[run]
    summary = load(summary_path)
    assert summary['status'] == 'passed' and summary['source_unchanged']
    tman, vman = load(timings / 'manifest.json'), load(BENCH / 'validation' / 'manifest.json')
    assert tman['status'] == vman['status'] == 'passed'
    assert (tman['phase'], vman['phase']) == ('benchmark', 'validate')
    assert tman['source_hashes'] == vman['source_hashes'] and tman['library'] == vman['library']
    assert tman['gamma'] == GAMMA and tman['weight_scale'] == 200 and tman['auto_factor'] == 1.0
    assert tman['pair_cache'] == 'cold' and tman['threads'] == 1 and tman['assertions'] and vman['assertions']
    vsum = load(BENCH / 'validation_summary.json')
    assert vsum['status'] == 'passed' and vsum['panels'] == len(summary['panels']) == 14
    assert vsum['components'] == summary['components'] == 2800 and vsum['states'] == summary['states']
    by_panel = {(p['locus'], p['panel']): p for p in summary['panels']}
    assert len(by_panel) == 14

    results, loads = {}, []
    for locus, (panels, leaveout) in BENCH_PANELS.items():
        for name, H in panels + ((leaveout,) if leaveout else ()):
            res = bench_panel(locus, name, H, by_panel[(locus, name)], tman['repetitions'], timings)
            results[(locus, H)] = res
            loads += res['load']
    assert len(results) == 14
    assert sum(r['states'] for r in results.values()) == vsum['states']
    return results, loads, tman, vsum


def bench100():
    results, loads, tman, vsum = bench_results(BENCH_HEADLINE)
    repetitions = tman['repetitions']
    vsrc, src = 'benchmark_100reads', str(BENCH_RUNS[BENCH_HEADLINE][0].relative_to(EXP))
    put('b100.reads', '100', vsrc + ' reads per panel')
    put('b100.panels', str(len(results)), vsrc + ' panels')
    put('b100.components', fint(vsum['components']), vsrc + '/validation_summary.json components (read orientations)')
    put('b100.states.billion', f'{vsum["states"] / 1e9:.2f}', vsrc + '/validation_summary.json states, counted once per panel (billions)')
    put('b100.direct.checks', str(vsum['direct_checks']), vsrc + '/validation_summary.json direct-recurrence cross-checks')
    put('b100.historical', str(vsum['historical_components_matched']), vsrc + '/validation_summary.json earlier components reproduced')
    put('b100.reps', str(repetitions), src + '/manifest.json repetitions')
    put('b100.resamples', fint(BENCH_RESAMPLES), vsrc + '/summarize.py bootstrap resamples (replayed here)')
    put('b100.load.lo', f1(min(loads)), src + ' minimum one-minute load before/after panel timings')
    put('b100.load.hi', f1(max(loads)), src + ' maximum one-minute load before/after panel timings')
    put('b100.native.library', Path(tman['library']).stem, src + '/manifest.json native library')
    # The reproducibility appendix says this is the repaired-guard benchmark's library.
    assert numbers['b100.native.library'][0] == numbers['guard.native.library'][0]
    by_locus = {}
    for (locus, _), r in results.items():
        by_locus[locus] = by_locus.get(locus, 0) + r['fallbacks']
    put('b100.fallbacks', str(sum(by_locus.values())), src + ' fallback components over all panels')
    assert by_locus['mhc'] == by_locus['chr20'] == 0  # the text says all fallbacks are on LRC/KIR
    put('b100.fallbacks.kir', str(by_locus['kir']), src + ' fallback components on LRC/KIR panels')

    for (locus, H), r in results.items():
        p = f'b100.{locus}.h{H}'
        s = f'{src} {locus}-H{H}'
        a = r['auto']
        put(p + '.x.guard', f2(a['median']), s + ' median per-read per-haplotype/auto kernel ratio')
        put(p + '.x.guard.lo', f2(a['ci'][0]), s + ' 95% bootstrap interval of that median, low')
        put(p + '.x.guard.hi', f2(a['ci'][1]), s + ' 95% bootstrap interval of that median, high')
        put(p + '.total.guard', f2(a['total']), s + ' summed per-haplotype / summed auto time')
        put(p + '.total.guard.lo', f2(a['total_ci'][0]), s + ' 95% bootstrap interval of the total ratio, low')
        put(p + '.total.guard.hi', f2(a['total_ci'][1]), s + ' 95% bootstrap interval of the total ratio, high')
        put(p + '.wins.guard', str(a['wins']), s + ' reads where auto is faster')
        put(p + '.slowest', f2(a['slowest']), s + ' slowest read, auto time / per-haplotype time')
        put(p + '.x.class', f2(r['class']['median']), s + ' median per-read per-haplotype/class (expanded, unguarded)')
        put(p + '.x.compact', f2(r['compact']['median']), s + ' median per-read per-haplotype/compact (unguarded)')
        put(p + '.ms.perhap', f2(r['perhap_ms']), s + ' median per-read per-haplotype kernel ms')
        put(p + '.ms.guard', f2(a['ms']), s + ' median per-read auto kernel ms')
        put(p + '.records', f2(r['records']), s + ' median per-read compact records / anchors')
        put(p + '.sliced', f'{r["sliced"]:.3f}', s + ' median per-read bit-sliced pair decisions / anchors')
        put(p + '.statesperanchor', f'{r["statesperanchor"]:.0f}', s + ' median per-read states / anchors')
        put(p + '.fallbacks', str(r['fallbacks']), s + ' fallback components')
        put(p + '.prev20', f2(r['previous_20']), s + ' median per-read auto speedup, 20 reads of the earlier benchmark')
        put(p + '.add80', f2(r['additional_80']), s + ' median per-read auto speedup, 80 added reads')
        put(p + '.nsperstate', f1(r['nsperstate']), s + ' median per-read per-haplotype ns per state (timed)')
        put(p + '.recwork', f2(r['recwork']), s + ' median per-read class time without expansion / per-haplotype time (validation run)')
        put(p + '.expand', f2(r['expand']), s + ' median per-read class expansion time / per-haplotype time (validation run)')
        put(p + '.premium.pct', f'{100 * (r["premium"] - 1):.0f}', s + ' median per-read auto / class time minus one, percent (timed)')
    mhc = results[('mhc', 463)]
    for field in ('same', 'loose', 'range', 'sliced'):
        put(f'b100.mhc.h463.share.{field}', f'{100 * mhc["share." + field]:.1f}',
            f'{src} mhc-H463 median per-read share of pair decisions: {field} (%)')
    put('b100.mhc.h463.expand.share', f'{100 * mhc["expand_share"]:.0f}',
        f'{src} mhc-H463 median per-read class expansion time / class kernel time (%, validation run)')

    largest = [results[(locus, panels[-1][1])] for locus, (panels, _) in BENCH_PANELS.items()]
    smallest = [results[(locus, panels[0][1])] for locus, (panels, _) in BENCH_PANELS.items()]
    allh = [r for r in results.values()]
    for r in smallest:
        assert r['auto']['wins'] == 0  # the text says no read is faster at 60 haplotypes
    kir = results[('kir', 459)]['auto']['total_ci']
    assert kir[0] < 1 < kir[1]  # the text says this interval includes parity
    for loc in ('mhc', 'chr20'):
        assert results[(loc, 463)]['auto']['total_ci'][0] > 1
    put('b100.hmax.guard.min', f2(min(r['auto']['median'] for r in largest)), src + ' smallest guarded median at the largest panels')
    put('b100.hmax.guard.max', f2(max(r['auto']['median'] for r in largest)), src + ' largest guarded median at the largest panels')
    put('b100.h60.guard.min', f2(min(r['auto']['median'] for r in smallest)), src + ' smallest guarded median at H=60')
    put('b100.h60.guard.max', f2(max(r['auto']['median'] for r in smallest)), src + ' largest guarded median at H=60')
    put('b100.premium.h60.lo', f'{100 * (min(r["premium"] for r in smallest) - 1):.0f}', src + ' guard premium at H=60, min over loci (%)')
    put('b100.premium.h60.hi', f'{100 * (max(r["premium"] for r in smallest) - 1):.0f}', src + ' guard premium at H=60, max over loci (%)')
    put('b100.premium.hmax.lo', f'{100 * (min(r["premium"] for r in largest) - 1):.0f}', src + ' guard premium at the largest panels, min (%)')
    put('b100.premium.hmax.hi', f'{100 * (max(r["premium"] for r in largest) - 1):.0f}', src + ' guard premium at the largest panels, max (%)')
    put('b100.expand.lo', f'{100 * min(r["expand"] for r in allh):.0f}', src + ' class expansion / per-haplotype time, min over panels (%)')
    put('b100.expand.hi', f'{100 * max(r["expand"] for r in allh):.0f}', src + ' class expansion / per-haplotype time, max over panels (%)')
    put('b100.pairreads.min', str(min(r['pair_reads'] for r in allh)),
        src + ' fewest reads with at least one record-pair decision on any panel (share medians use these)')
    agreement = max(r['stage_agreement'] for r in allh)
    assert agreement < 0.10, agreement
    put('b100.stage.agree.pct', f'{100 * agreement:.0f}',
        src + ' max relative difference, single validation run vs timed median class speedup (%)')

    def label(locus, H):
        return f'{H}' + (r'$^{\ast}$' if BENCH_PANELS[locus][1] and H == BENCH_PANELS[locus][1][1] else '')

    def ci(x, lohi):
        return f'{f2(x)} [{f2(lohi[0])}, {f2(lohi[1])}]'
    main_rows, times, ratios, counters, stages = [], [], [], [], []
    for locus, (panels, leaveout) in BENCH_PANELS.items():
        hs = [H for _, H in panels] + ([leaveout[1]] if leaveout else [])
        lines = ['H guard guard_lo guard_hi class compact']
        for H in hs:
            r = results[(locus, H)]
            a = r['auto']
            name = LOCUS_NAME[locus]
            times.append((name, label(locus, H), f'{r["statesperanchor"]:.0f}', f2(r['perhap_ms']),
                          f2(r['class']['ms']), f2(a['ms']), f2(r['compact']['ms'])))
            ratios.append((name, label(locus, H), ci(a['median'], a['ci']), ci(a['total'], a['total_ci']),
                           f2(r['class']['median']), f2(r['compact']['median']), str(a['wins']), str(r['fallbacks'])))
            counters.append((name, label(locus, H), *(f'{r["share." + f]:.3f}' for f in ('same', 'loose', 'range', 'sliced')),
                             f'{r["sliced"]:.3f}', f'{r["records"]:.3f}'))
            stages.append((name, label(locus, H), f1(r['nsperstate']), f2(r['recwork']), f2(r['expand']),
                           f2(r['premium']), f2(a['median'])))
            if H in (panels[0][1], panels[-1][1]):
                main_rows.append((name if H == panels[0][1] else '', str(H), f2(r['records']), f2(r['perhap_ms']),
                                  f2(a['ms']), ci(a['median'], a['ci']), f2(a['total']), f2(r['compact']['median'])))
            if not (leaveout and H == leaveout[1]):
                lines.append(f'{H} {a["median"]:.4f} {a["ci"][0]:.4f} {a["ci"][1]:.4f} '
                             f'{r["class"]["median"]:.4f} {r["compact"]["median"]:.4f}')
        with open(OUT / f'b100_speedup_{locus}.dat', 'w') as fh:
            fh.write('\n'.join(lines) + '\n')

    def write_rows(filename, rows):
        with open(OUT / filename, 'w') as fh:
            fh.write(f'% Generated by paper/scripts/make_paper_data.py from experiments/{src}; do not edit.\n')
            prev = None
            for r in rows:
                if prev is not None and r[0] and r[0] != prev:
                    fh.write('\\midrule\n' if filename != 'rows_b100_main.tex' else '\\addlinespace\n')
                fh.write(' & '.join(r if r[0] != prev else ('',) + r[1:]) + ' \\\\\n')
                prev = r[0] or prev
    write_rows('rows_b100_main.tex', main_rows)
    write_rows('rows_b100_times.tex', times)
    write_rows('rows_b100_ratios.tex', ratios)
    write_rows('rows_b100_counters.tex', counters)
    write_rows('rows_b100_stages.tex', stages)
    return results


def bench_first_run(results, first, first_loads):
    """The first run of the 100-read benchmark, on a loaded host, against the reported run: its load range and
    the largest relative difference between the two runs' per-panel guarded medians (rounded up)."""
    assert results.keys() == first.keys()
    diff = max(abs(r['auto']['median'] / first[k]['auto']['median'] - 1) for k, r in results.items())
    src = str(BENCH_RUNS['first'][0].relative_to(EXP)) + ' (first run)'
    put('b100.first.load.lo', f1(min(first_loads)), src + ' minimum one-minute load before/after panel timings')
    put('b100.first.load.hi', f1(max(first_loads)), src + ' maximum one-minute load before/after panel timings')
    put('b100.first.maxdiff.pct', str(math.ceil(100 * diff)),
        src + ' largest |reported / first - 1| over per-panel guarded medians (%, rounded up)')


# ---------------------------------------------------------------------------------------------
# Native memory per engine and record counts under alternative groupings (2026-10-20,
# experiments/memory_ablation)
# ---------------------------------------------------------------------------------------------
MEMABL = EXP / 'memory_ablation/results'


def memory_ablation(b100):
    """Supplementary tables S8 (memory) and S9 (groupings); main-text keys mem.* and abl.*.

    Records per anchor are recomputed from the compact output of a separate pass with the same engine and
    inputs, so they must equal the 100-read benchmark's counters exactly (asserted)."""
    s = load(MEMABL / 'summary.json')
    man = load(MEMABL / 'run1/manifest.json')
    assert man['status'] == 'passed' and man['gamma'] == GAMMA and man['weight_scale'] == 200
    assert man['pair_cache'] == 'cold' and man['auto_factor'] == 1.0 and man['assertions']
    by = {(p['locus'], p['H']): p for p in s['panels']}
    assert len(by) == 14
    mem_rows, abl_rows = [], []
    for locus, (panels, leaveout) in BENCH_PANELS.items():
        for name, H in panels + ((leaveout,) if leaveout else ()):
            p = by[(locus, H)]
            assert p['panel'] == name
            assert abs(p['per_anchor.records'] - b100[(locus, H)]['records']) < 1e-9, (locus, H)
            lab = f'{H}' + (r'$^{\ast}$' if leaveout and H == leaveout[1] else '')
            mem_rows.append((LOCUS_NAME[locus], lab, f1(p['index_native_mb']),
                             *(f2(p[f'{m}.native_peak.median_mb']) for m in ('per_haplotype', 'class', 'auto', 'compact')),
                             f2(p['per_haplotype.output.median_mb']), f2(p['compact.output.median_mb']),
                             f2(p['auto.total.max_mb'])))
            abl_rows.append((LOCUS_NAME[locus], lab, f'{p["per_anchor.states"]:.0f}', f1(p['per_anchor.keys']),
                             f2(p['per_anchor.scores']), f2(p['per_anchor.records']),
                             f1(p['ratio.keys_over_records']), f2(p['ratio.scores_over_records'])))
            if H in (panels[0][1], panels[-1][1]):
                k, src = f'{locus}.h{H}', f'memory_ablation summary {locus}-{name}'
                put(f'mem.{k}.index', f1(p['index_native_mb']), src + ' native index (MiB)')
                for m, short in (('per_haplotype', 'ph'), ('class', 'class'), ('auto', 'guard'), ('compact', 'compact')):
                    put(f'mem.{k}.{short}', f2(p[f'{m}.native_peak.median_mb']), src + f' median per-read native peak, {m} (MiB)')
                    put(f'mem.{k}.{short}.max', f2(p[f'{m}.native_peak.max_mb']), src + f' max per-read native peak, {m} (MiB)')
                put(f'mem.{k}.out.expanded', f2(p['per_haplotype.output.median_mb']), src + ' median expanded output (MiB)')
                put(f'mem.{k}.out.compact', f2(p['compact.output.median_mb']), src + ' median compact output (MiB)')
                put(f'mem.{k}.out.ratio', f'{p["ratio.expanded_over_compact_output"]:.0f}', src + ' median expanded/compact output')
                put(f'abl.{k}.states', f'{p["per_anchor.states"]:.0f}', src + ' median per-read states per anchor')
                put(f'abl.{k}.keys', f1(p['per_anchor.keys']), src + ' median per-read distinct keys per anchor')
                put(f'abl.{k}.scores', f2(p['per_anchor.scores']), src + ' median per-read (score, witness) classes per anchor')
                put(f'abl.{k}.records', f2(p['per_anchor.records']), src + ' median per-read origin records per anchor')
                put(f'abl.{k}.keys.ratio', f'{p["ratio.keys_over_records"]:.0f}', src + ' median per-read keys/records')
                put(f'abl.{k}.scores.ratio', f1(p['ratio.scores_over_records']), src + ' median per-read scores/records')
    worst = max(1 - p['per_anchor.score_values'] / p['per_anchor.scores'] for p in s['panels'])
    put('abl.witness.maxpct', f1(100 * worst), 'memory_ablation summary: largest relative drop of the per-panel median '
        'score classes when the witness rule is dropped (%)')
    rss = s['minichain_peak_rss']
    for h in ('60', '463'):
        put(f'mem.minichain.h{h}', f'{rss[h]["median_gb"]:.2f}', f'minichain_matched time_hprc{h} logs, median peak RSS (GB, 100 reads)')

    def write(filename, rows):
        with open(OUT / filename, 'w') as fh:
            fh.write('% Generated by paper/scripts/make_paper_data.py from experiments/memory_ablation; do not edit.\n')
            prev = None
            for r in rows:
                if prev is not None and r[0] != prev:
                    fh.write('\\midrule\n')
                fh.write(' & '.join(r if r[0] != prev else ('',) + r[1:]) + ' \\\\\n')
                prev = r[0]
    write('rows_memory.tex', mem_rows)
    write('rows_groupings.tex', abl_rows)


# ---------------------------------------------------------------------------------------------
# Minichain with the engines as its chainer (2026-10-20, experiments/minichain_integrated)
# ---------------------------------------------------------------------------------------------
INTEG = EXP / 'minichain_integrated/results'
INTEG_DATASETS = ('mhc2104', 'hprc60', 'hprc120', 'hprc240', 'hprc463')


def integration():
    """Correctness of the drop-in chainer: DP-table identity for every component, GAF agreement per read.

    The engine-off binary must reproduce the pinned binary's saved GAF byte for byte, and every engine must
    reproduce every DP hash of the engine-off run; otherwise nothing is reported."""
    m = load(INTEG / 'correctness/manifest.json')
    assert m['status'] == 'passed' and set(m['comparisons']) == set(INTEG_DATASETS)
    comps = m['comparisons']
    assert all(c['off_gaf_equals_pinned'] for c in comps.values())

    def tot(v, k):
        return sum(c[v][k] for c in comps.values())
    for v in ('ph', 'class', 'auto'):
        assert tot(v, 'dp_hash_mismatches') == 0, v
        assert tot(v, 'components') == tot('ph', 'components') and tot(v, 'reads') == tot('ph', 'reads')
    src = 'minichain_integrated/results/correctness/manifest.json'
    reads, components = tot('ph', 'reads'), tot('ph', 'components')
    put('int.datasets', str(len(comps)), src + ' datasets')
    put('int.reads', fint(reads), src + ' reads over all datasets')
    put('int.components', fint(components), src + ' (read, component) DP tables compared per engine')
    for v, name in (('ph', 'ph'), ('class', 'class'), ('auto', 'guard')):
        same, diff = tot(v, 'gaf_identical_reads'), tot(v, 'gaf_different_reads')
        assert same + diff == reads
        put(f'int.{name}.identical', fint(same), src + f' reads with byte-identical GAF, {v}')
        put(f'int.{name}.identical.pct', pct(same, reads), src + f' reads with byte-identical GAF, {v} (%)')
        put(f'int.{name}.different', fint(diff), src + f' reads whose GAF differs, {v}')
        put(f'int.{name}.sameinterval', fint(tot(v, 'different_reads_with_same_path_and_target_interval')),
            src + f' differing reads with the same path name and target interval, {v}')
        if diff:
            put(f'int.{name}.sameinterval.pct', pct(tot(v, 'different_reads_with_same_path_and_target_interval'), diff),
                src + f' differing reads with the same path name and target interval, {v} (% of differing)')
    s = INTEG / 'summary.json'
    if s.exists() and 'timing' in load(s):
        tm = load(INTEG / 'timing/manifest.json')
        assert tm['status'] == 'passed' and len(tm['runs']) == 3 * 4 * len(INTEG_DATASETS)
        loads = [x for r in tm['runs'] for x in (r['host_before']['load'][0], r['host_after']['load'][0])]
        put('int.load.hi', f'{max(loads):.1f}', 'minichain_integrated timing manifest: largest 1-minute load around any run')
        put('int.reps', '3', 'minichain_integrated timing repetitions per variant')
        t = load(s)['timing']
        x = t['hprc463']
        put('int.e2e.class.over.ph.hprc463', f2(x['ph']['real_seconds'] / x['class']['real_seconds']),
            'minichain_integrated summary hprc463: per-haplotype over class end-to-end real time')
        put('int.dp.hprc463.class', f1(x['class']['dp_seconds']), 'minichain_integrated summary hprc463 class: summed DP-section seconds')
        put('int.rest.hprc463.class', f1(x['class']['real_seconds'] - x['class']['dp_seconds']),
            'minichain_integrated summary hprc463 class: real time outside the engine DP section (s)')
        put('int.postdp.hprc463.off', f1(x['off']['chaining_seconds'] - x['off']['dp_seconds']),
            "minichain_integrated summary hprc463 off: Minichain's chaining time after its DP (backtracking) (s)")
        for d in ('mhc2104', 'hprc463'):
            off = t[d]['off']
            put(f'int.time.{d}.off', f'{off["real_seconds"]:.0f}', f'minichain_integrated summary {d} off real seconds (median)')
            put(f'int.rss.{d}.off', f'{off["peak_rss_gb"]:.2f}', f'minichain_integrated summary {d} off peak RSS (GB)')
            for v, name in (('ph', 'ph'), ('class', 'class'), ('auto', 'guard')):
                x = t[d][v]
                put(f'int.time.{d}.{name}', f1(x['real_seconds']), f'minichain_integrated summary {d} {v} real seconds (median)')
                put(f'int.x.{d}.{name}', f1(x['end_to_end_speedup']), f'minichain_integrated summary {d} {v} end-to-end speedup')
                put(f'int.rss.{d}.{name}', f'{x["peak_rss_gb"]:.2f}', f'minichain_integrated summary {d} {v} peak RSS (GB)')


# ---------------------------------------------------------------------------------------------
# Base-level windows of HPRC release 2 (2026-10-20, experiments/baselevel_feasibility; analysis
# fixed in its REPORT.md before timing or exactness ran)
# ---------------------------------------------------------------------------------------------
BASE = EXP / 'baselevel_feasibility'
BASE_WINDOWS = (('w31250k_250kb', 'class I'), ('w31700k_250kb', 'class III'), ('w32375k_100kb', 'class II'))


def baselevel():
    """Correctness and counters per window; timing and exactness once their summaries exist."""
    timing = load(BASE / 'timing/summary.json') if (BASE / 'timing/summary.json').exists() else None
    exact = load(BASE / 'exactness/summary.json') if (BASE / 'exactness/summary.json').exists() else None
    tw = {x['window']: x for x in timing['windows']} if timing else {}
    rows, counters = [], []
    tot = dict(reads=0, anchors=0, states=0)
    loads, guarded, compact = [], [], []
    spans, best_losses = [], []  # extracted reference spans; per-read losses in the better orientation (2026-10-30)
    for w, region in BASE_WINDOWS:
        d = BASE / 'windows' / w
        ins, prov = load(d / 'inspect.json'), load(d / 'panel/provenance.json')
        eng, val = load(d / 'engine_check.json')['summary'], load(d / 'validation.json')
        assert ins['acyclic'] and not ins['repeated_handle_walks'] and not ins['segments_in_both_orientations']
        assert ins['walks'] == ins['walks_spanning_window'] == eng['H']
        assert val['status'] == 'passed' and val['states'] == eng['states'] and val['instances'] == eng['instances']
        assert prov['conversion']['anchors'] == eng['anchors'] and prov['conversion']['exact_sequence_checks'] == eng['anchors']
        spans.append(ins['reference_bases'])
        k, src = f'base.{w[:7]}', f'baselevel_feasibility windows/{w}'
        put(f'{k}.start', f'{ins["start"] / 1e6:.3f}', src + '/inspect.json start (Mb)')
        put(f'{k}.end', f'{ins["end"] / 1e6:.3f}', src + '/inspect.json end (Mb)')
        put(f'{k}.vertices', fint(ins['segments']), src + '/inspect.json segments')
        put(f'{k}.reads', str(prov['reads']), src + '/panel/provenance.json reads inside the window')
        put(f'{k}.anchors', fint(eng['anchors']), src + '/engine_check.json anchors')
        put(f'{k}.states.million', f'{eng["states"] / 1e6:.1f}', src + '/engine_check.json states (millions)')
        put(f'{k}.records', f2(eng['records_per_anchor_median']), src + '/engine_check.json median per-read records per anchor')
        put(f'{k}.sliced', f'{eng["sliced_per_anchor_median"]:.3f}', src + '/engine_check.json median per-read sliced pairs per anchor')
        put(f'{k}.fallbacks', str(eng['fallback_instances']), src + '/engine_check.json guard fallbacks (instances)')
        put(f'{k}.instances', str(eng['instances']), src + '/engine_check.json instances')
        for t in ('reads', 'anchors', 'states'):
            tot[t] += {'reads': prov['reads'], 'anchors': eng['anchors'], 'states': eng['states']}[t]
        counters.append((region, f'{ins["start"] / 1e6:.2f}--{ins["end"] / 1e6:.2f}', fint(ins['segments']), str(prov['reads']),
                         fint(eng['anchors']), f'{eng["states_per_anchor_median"]:.0f}', f2(eng['records_per_anchor_median']),
                         f'{eng["sliced_per_anchor_median"]:.3f}', f'{eng["fallback_instances"]}/{eng["instances"]}'))
        row = [region]
        if w in tw:
            x = tw[w]
            assert x['instances'] == eng['instances'] and x['states'] == eng['states']
            a = x['kernel']['auto']
            assert a['reads'] == prov['reads'] and a['median_ci95'][0] > 1  # the text: faster in every window
            loads += [x['host_before']['load'][0], x['host_after']['load'][0]]
            guarded.append(a['median_speedup'])
            compact.append(x['kernel']['compact']['median_speedup'])
            put(f'{k}.x.guard', f2(a['median_speedup']), f'baselevel timing summary {w} median per-read guarded speedup')
            put(f'{k}.x.guard.lo', f2(a['median_ci95'][0]), f'baselevel timing summary {w} 95% CI low')
            put(f'{k}.x.guard.hi', f2(a['median_ci95'][1]), f'baselevel timing summary {w} 95% CI high')
            put(f'{k}.total.guard', f2(a['total_speedup']), f'baselevel timing summary {w} total guarded speedup')
            put(f'{k}.wins.guard', str(a['wins']), f'baselevel timing summary {w} reads faster guarded')
            put(f'{k}.x.compact', f2(x['kernel']['compact']['median_speedup']), f'baselevel timing summary {w} compact')
            row += [f'{f2(a["median_speedup"])} [{f2(a["median_ci95"][0])}, {f2(a["median_ci95"][1])}]',
                    f2(a['total_speedup']), f2(x['kernel']['compact']['median_speedup'])]
        else:
            row += ['--', '--', '--']
        if exact:
            reps = [v for name, v in exact.items() if name.startswith(w)]
            below = {rho: [max(r['rhos'][rho]['h0']['instances_below'], r['rhos'][rho]['h12']['instances_below']) for r in reps]
                     for rho in ('2', '56')}
            inst = sum(r['rhos']['2']['instances'] for r in reps)
            put(f'{k}.giraffe.instances', fint(inst), f'baselevel exactness summary {w}: instances over replicates')
            for rho in ('2', '56'):
                put(f'{k}.giraffe.below.rho{rho}', str(sum(below[rho])), f'baselevel exactness {w} rho {rho}: instances below exact (max over bonus), summed over replicates')
                loss = max(r['rhos'][rho][b]['max_rel_gap'] for r in reps for b in ('h0', 'h12'))
                put(f'{k}.giraffe.maxloss.rho{rho}', f'{100 * loss:.1f}', f'baselevel exactness {w} rho {rho}: largest relative loss (%)')
            put(f'{k}.giraffe.below.rho2.pct', pct(sum(below['2']), inst), f'baselevel exactness {w} rho 2: instances below exact over instances (%)')
            sw = [r['rhos']['2']['exact_with_switch'] for r in reps]
            put(f'{k}.giraffe.switch.rho2', str(sum(sw)), f'baselevel exactness {w} rho 2: instances whose exact optimum switches, summed over samples')
            # Distinct reads behind the rho-2 instances below the optimum (2026-10-25), from the per-instance rows,
            # which must reproduce the summary's counts. The text attributes the counted instances to these reads,
            # so the instances below under either bonus must be exactly those counted (the larger bonus set).
            seen, below_reads = set(), set()
            for name in sorted(n for n in exact if n.startswith(w)):
                rows2 = [x for x in load_jsonl(BASE / 'exactness' / f'{name}.jsonl') if str(x['rho']) == '2']
                assert len(rows2) == exact[name]['rhos']['2']['instances'], name
                for b in ('h0', 'h12'):
                    assert sum(x[b]['gap'] > 0 for x in rows2) == exact[name]['rhos']['2'][b]['instances_below'], (name, b)
                either = [x for x in rows2 if x['h0']['gap'] > 0 or x['h12']['gap'] > 0]
                assert len(either) == max(exact[name]['rhos']['2'][b]['instances_below'] for b in ('h0', 'h12')), name
                seen |= {x['read'] for x in rows2}
                below_reads |= {x['read'] for x in either}
            assert len(seen) == prov['reads'], (w, len(seen))
            put(f'{k}.giraffe.below.rho2.reads', str(len(below_reads)), f'baselevel exactness {w} rho 2: distinct reads with an instance below exact (either bonus), from the per-instance rows')
            # Each read's better orientation (2026-10-30). A read's two orientations are separate instances; the text
            # also counts read-sample pairs whose better orientation stays below the optimum, from the same rows.
            pairs = {}
            for name in sorted(n for n in exact if n.startswith(w)):
                for x in load_jsonl(BASE / 'exactness' / f'{name}.jsonl'):
                    if str(x['rho']) == '2':
                        pairs.setdefault((name, x['read']), {})[x['strand']] = x
            assert 2 * len(pairs) == inst and all(sorted(s) == ['+', '-'] for s in pairs.values()), w
            weaker, worst, short = 0, (0.0, False), {'h0': set(), 'h12': set()}
            for key, s in pairs.items():
                for x, y in ((s['+'], s['-']), (s['-'], s['+'])):
                    rel = max(x[b]['rel_gap'] for b in ('h0', 'h12'))
                    weaker += rel > 0 and x['opt'] < y['opt']
                    if rel > worst[0]:
                        worst = (rel, x['opt'] < y['opt'])
                top = max(s['+']['opt'], s['-']['opt'])
                for b in short:
                    loss = top - max(s['+'][b]['best'], s['-'][b]['best'])
                    assert loss >= 0, (w, key, b)
                    if loss > 0:
                        short[b].add(key)
                        best_losses.append((loss, loss / top))
            assert worst[0] == max(r['rhos']['2'][b]['max_rel_gap'] for r in reps for b in ('h0', 'h12')), w
            assert worst[1] or not w.startswith('w32375k'), 'the text: the largest class II loss is in a weaker orientation'
            either = short['h0'] | short['h12']
            assert len(either) == max(len(v) for v in short.values()), w  # one bonus's set holds the other's
            put(f'{k}.giraffe.best.pairs', str(len(pairs)), f'baselevel exactness {w} rho 2: read-sample pairs')
            put(f'{k}.giraffe.best.below.rho2', str(len(either)), f'baselevel exactness {w} rho 2: read-sample pairs whose better orientation is below exact (either bonus)')
            put(f'{k}.giraffe.best.reads.rho2', str(len({r for _, r in either})), f'baselevel exactness {w} rho 2: distinct reads of those pairs')
            put(f'{k}.giraffe.weaker.rho2', str(weaker), f'baselevel exactness {w} rho 2: instances below exact in the read\'s weaker orientation')
            for b in ('h0', 'h12'):
                put(f'{k}.giraffe.below.rho2.{b}', str(sum(r['rhos']['2'][b]['instances_below'] for r in reps)),
                    f'baselevel exactness {w} rho 2: instances below exact with consistency bonus {b[1:]}')
            row += [f'{sum(below["2"])}/{inst}', f'{sum(below["56"])}/{inst}']
        else:
            row += ['--', '--']
        rows.append(row)
    if tw:
        put('base.load.hi', f1(max(loads)), 'baselevel timing summary: maximum one-minute load before/after window timings')
        put('base.x.guard.lo', f2(min(guarded)), 'baselevel timing summary: smallest median per-read guarded speedup over windows')
        put('base.x.guard.hi', f2(max(guarded)), 'baselevel timing summary: largest such speedup')
        put('base.x.compact.lo', f2(min(compact)), 'baselevel timing summary: smallest median per-read compact speedup over windows')
        put('base.x.compact.hi', f2(max(compact)), 'baselevel timing summary: largest such speedup')
    if exact:
        allrows = list(exact.values())
        put('base.giraffe.instances', fint(sum(r['rhos']['2']['instances'] for r in allrows)),
            'baselevel exactness summary: instances over windows and 64-haplotype samples')
        for rho in ('2', '56'):
            put(f'base.giraffe.below.rho{rho}', str(sum(max(r['rhos'][rho]['h0']['instances_below'], r['rhos'][rho]['h12']['instances_below'])
                                                       for r in allrows)),
                f'baselevel exactness summary rho {rho}: instances below exact, max over bonus, summed')
        worst = max(r['rhos'][rho][b]['max_rel_gap'] for r in allrows for rho in r['rhos'] for b in ('h0', 'h12'))
        put('base.giraffe.maxloss', f'{100 * worst:.2f}', 'baselevel exactness summary: largest relative loss over all rho and bonuses (%)')
        put('base.giraffe.samples', str(len(allrows)), 'baselevel exactness summary: window samples')
        points = max(p for p, _ in best_losses)
        assert points == 1, 'the text says "at most 1 point"'
        put('base.giraffe.best.maxloss.points', str(points), 'baselevel exactness rho 2: largest loss in a read\'s better orientation (score points)')
        put('base.giraffe.best.maxloss.pct', f'{100 * max(r for _, r in best_losses):.2f}',
            'baselevel exactness rho 2: largest relative loss in a read\'s better orientation (%)')
    put('base.span.kb.lo', f'{min(spans) / 1000:.0f}', 'baselevel windows: smallest extracted CHM13 reference span (kb, inspect.json reference_bases)')
    put('base.span.kb.hi', f'{max(spans) / 1000:.0f}', 'baselevel windows: largest extracted CHM13 reference span (kb)')
    put('base.reads', str(tot['reads']), 'baselevel windows reads (sum)')
    put('base.anchors', fint(tot['anchors']), 'baselevel windows anchors (sum)')
    put('base.states.million', f'{tot["states"] / 1e6:.0f}', 'baselevel windows states (sum, millions)')
    for name, table in (('rows_baselevel_counters.tex', counters), ('rows_baselevel_results.tex', rows)):
        with open(OUT / name, 'w') as fh:
            fh.write('% Generated by paper/scripts/make_paper_data.py from experiments/baselevel_feasibility; do not edit.\n')
            for r in table:
                fh.write(' & '.join(r) + ' \\\\\n')


# ---------------------------------------------------------------------------------------------
# Pre-specified per-read association analysis (Section 6.5, Supplementary Figure S1)
# ---------------------------------------------------------------------------------------------
MECH = EXP / 'mechanism_100reads'
MECH_LARGEST = {'mhc': 463, 'kir': 459, 'chr20': 463}
MECH_CLIP = 4.2  # x-axis limit of figures/fig-scatter.tex; reads beyond it are drawn at the edge


def _signed(x: float) -> str:
    s = f2(x)
    return '$-$' + s[1:] if s.startswith('-') else s


def _auroc(a, b):
    # Mann-Whitney U / (n_a n_b), ties counted one half.
    return sum((x > y) + 0.5 * (x == y) for x in a for y in b) / (len(a) * len(b))


def mechanism(results):
    """Recompute the analysis's point estimates from the raw rows of the first timing run, on which it was
    specified, and check them, with the Holm-tested verdicts the text relies on, against its saved results."""
    plan = MECH / 'PLAN.md'
    inputs_used.update((plan, MECH / 'PLAN.sha256'))
    res = load(MECH / 'results.json')
    frozen = (MECH / 'PLAN.sha256').read_text().split()[0]
    assert res['status'] == 'passed' and res['plan_sha256'] == frozen == hashlib.sha256(plan.read_bytes()).hexdigest()
    assert res['decisions']['section_6_5_sentence'] == 'supported at all three largest panels'
    tests, src = res['primary']['tests'], 'mechanism_100reads'
    names = {(locus, H): name for locus, (panels, leaveout) in BENCH_PANELS.items()
             for name, H in panels + ((leaveout,) if leaveout else ())}

    def columns(reads):
        y = [math.log(r['t']['per_haplotype'] / r['t']['auto']) for r in reads]
        return y, [r['records'] / r['N'] for r in reads], [r['sliced'] / r['N'] for r in reads]
    mhckir, chr20big = [], []
    for (locus, H), r in results.items():
        key = f'{locus}-{names[(locus, H)]}'
        y, rec, sl = columns(r['reads'])
        for x, field in ((rec, 'x_rec'), (sl, 'x_sl')):
            assert _close(spearman(x, y), tests[f'P1/guard/{key}/{field}']['rho']), (key, field)
        verdict = {f: tests[f'P1/guard/{key}/{f}']['verdict'] for f in ('x_rec', 'x_sl')}
        if locus in ('mhc', 'kir'):  # the text: both counters track at every panel size
            assert verdict == {'x_rec': 'tracks', 'x_sl': 'tracks'}, key
            mhckir.append(spearman(rec, y))
        else:  # the text: records track from 240 haplotypes, bit-sliced decisions never
            assert verdict['x_sl'] != 'tracks' and (verdict['x_rec'] == 'tracks') == (H >= 240), key
            if H >= 240:
                chr20big.append(spearman(rec, y))
            else:  # "most reads have about one record per anchor"
                assert statistics.median(rec) < 1.02 and _bench_quantile(rec, .75) < 1.1, key
    put('mech.mhckir.rec.weak', _signed(max(mhckir)), src + ' weakest per-read Spearman rho, guarded speedup vs records/anchor, MHC and KIR panels')
    put('mech.mhckir.rec.strong', _signed(min(mhckir)), src + ' strongest such rho, MHC and KIR panels')
    put('mech.chr20big.rec.weak', _signed(max(chr20big)), src + ' weakest such rho, chr20 panels with H >= 240')
    put('mech.chr20big.rec.strong', _signed(min(chr20big)), src + ' strongest such rho, chr20 panels with H >= 240')

    pooled_y, pooled_rec = [], []
    for locus, H in MECH_LARGEST.items():
        reads = results[(locus, H)]['reads']
        y, rec, _ = columns(reads)
        n = [r['N'] for r in reads]
        rho = spearman(rec, y)
        r_yn, r_xn = spearman(y, n), spearman(rec, n)
        partial = (rho - r_yn * r_xn) / ((1 - r_yn ** 2) * (1 - r_xn ** 2)) ** 0.5
        assert _close(partial, res['exploratory']['E2'][f'E2/{locus}/y_guard~x_rec|N']['rho']), locus
        put(f'mech.{locus}.rec', _signed(rho), f'{src} {locus}-H{H} per-read Spearman rho, guarded speedup vs records/anchor (Holm-tested)')
        put(f'mech.{locus}.rec.partialN', _signed(partial), f'{src} {locus}-H{H} that rho with anchors per read held fixed (partial Spearman)')
        pooled_y += y
        pooled_rec += rec
        with open(OUT / f'b100_scatter_{locus}.dat', 'w') as fh:
            fh.write('records sliced guard compact fallback\n')
            for r in sorted(reads, key=lambda r: (r['records'] / r['N'], r['read'])):
                fh.write(f"{r['records'] / r['N']:.4f} {r['sliced'] / r['N']:.4f} "
                         f"{r['t']['per_haplotype'] / r['t']['auto']:.4f} "
                         f"{r['t']['per_haplotype'] / r['t']['compact']:.4f} {int(r['fallback'])}\n")
    assert _close(spearman(pooled_rec, pooled_y), tests['P2/guard/x_rec']['rho'])
    assert tests['P2/guard/x_rec']['verdict'] == tests['P3/guard/x_rec']['verdict'] == 'tracks'
    put('mech.pooled.rec', _signed(spearman(pooled_rec, pooled_y)), src + ' pooled rho over the three largest panels')
    put('mech.n', str(len(pooled_y)), src + ' reads at the three largest panels')
    clipped = sorted((r for locus, H in MECH_LARGEST.items() for r in results[(locus, H)]['reads']
                      if r['records'] / r['N'] > MECH_CLIP), key=lambda r: r['records'] / r['N'])
    assert len(clipped) == 2 and all(r['fallback'] and r in results[('kir', 459)]['reads'] for r in clipped)
    assert f'xmax={MECH_CLIP}' in (ROOT / 'paper/figures/fig-scatter.tex').read_text()
    with open(OUT / 'b100_scatter_clipped.dat', 'w') as fh:
        fh.write('records guard\n')
        for r in clipped:
            fh.write(f"{MECH_CLIP - 0.06:.2f} {r['t']['per_haplotype'] / r['t']['auto']:.4f}\n")
    put('mech.clip.a', f1(clipped[0]['records'] / clipped[0]['N']), src + ' records/anchor of the first LRC/KIR read beyond the Figure S1 axis')
    put('mech.clip.b', f1(clipped[1]['records'] / clipped[1]['N']), src + ' records/anchor of the second such read')

    # The guard's fallbacks: all on LRC/KIR (asserted in bench100).
    ids, at_output, aurocs = {}, 0, {}
    for (locus, H), r in results.items():
        if locus != 'kir':
            continue
        comps = [c for rd in r['reads'] for c in rd['comps'] if c['N'] > 0]
        fb = [c for c in comps if c['guard']['fallback']]
        aurocs[f'{locus}-{names[(locus, H)]}'] = _auroc([c['compact_records'] / c['N'] for c in fb],
                                                        [c['compact_records'] / c['N'] for c in comps if not c['guard']['fallback']])
        for c in fb:
            ids.setdefault(c['id'], set()).add(H)
            g = c['guard']
            at_output += g['auto_denied_work'] == c['N'] + g['W'] * c['compact_records'] + c['K']
    f = res['fallback']
    assert len(ids) == f['F4']['distinct_components'] and at_output == f['F1_summary']['stage']['output']
    assert aurocs.keys() == f['F2'].keys() and all(_close(a, f['F2'][k]['records_per_anchor']['auroc']) for k, a in aurocs.items())
    assert min(aurocs.values()) > 0.95  # the text: fallbacks are the instances with the most records
    late = sum(row['excess'] >= row['predicted_excess_validation'] for row in f['F1'])
    assert late >= 0.75 * len(f['F1'])  # the text: usually only after most of the class work
    put('mech.fb.distinct', str(len(ids)), src + ' distinct LRC/KIR instances among the fallbacks')
    put('mech.fb.allpanels', str(sum(len(v) == 5 for v in ids.values())), src + ' of these, instances falling back on all five LRC/KIR panels')
    put('mech.fb.output', str(at_output), src + ' fallbacks refused only at the output charge (refused units = N + W records + K)')
    put('mech.fb.auroc.lo', f2(min(aurocs.values())), src + ' AUROC of records/anchor separating fallback instances, minimum over LRC/KIR panels')
    put('mech.fb.auroc.hi', f2(max(aurocs.values())), src + ' that AUROC, maximum over LRC/KIR panels')
    # The fallback costs on this run must match the saved analysis; the paper reports them for the reported run
    # (fallback_cost()).
    ratios, total_ph = _fallback_cost(results)
    rows459 = [row['guarded_over_ph'] for row in f['F1'] if row['panel'] == 'kir-hprc459']
    assert _close(min(ratios), min(rows459)) and _close(max(ratios), max(rows459))
    assert _close(total_ph, f['F3']['kir-hprc459']['total_fallbacks_at_per_haplotype_time'])


def _fallback_cost(results):
    """LRC/KIR fallback instances at 459 haplotypes: guarded / per-haplotype time (timed medians), and the total
    speedup with each fallback instance charged its own per-haplotype time in the same repetition. The total is
    a reconstruction from measured times, not a measurement of another guard."""
    reads = results[('kir', 459)]['reads']
    reps = len(reads[0]['comps'][0]['kernel_seconds']['auto'])
    ratios = [statistics.median(c['kernel_seconds']['auto']) / statistics.median(c['kernel_seconds']['per_haplotype'])
              for rd in reads for c in rd['comps'] if c['guard']['fallback']]
    alt = sum(statistics.median(sum(c['kernel_seconds']['per_haplotype' if c['guard']['fallback'] else 'auto'][j]
                                    for c in rd['comps']) for j in range(reps)) for rd in reads)
    total_ph = sum(rd['t']['per_haplotype'] for rd in reads) / alt
    assert total_ph > results[('kir', 459)]['auto']['total']
    return ratios, total_ph


def fallback_cost(results):
    """Section 6.5's fallback costs, on the reported run of the 100-read benchmark."""
    ratios, total_ph = _fallback_cost(results)
    src = str(BENCH_RUNS[BENCH_HEADLINE][0].relative_to(EXP))
    put('b100.fb.h459.lo', f2(min(ratios)), src + ' kir-H459 fallback instances: guarded / per-haplotype time, minimum (timed medians)')
    put('b100.fb.h459.hi', f2(max(ratios)), src + ' that ratio, maximum')
    put('b100.kir.h459.total.fbph', f2(total_ph), src + ' kir-H459 total speedup with fallbacks charged their per-haplotype time (reconstruction)')


PILOT = EXP / 'baselevel_mhc_pilot'


def micb_pilot():
    """Base-level feasibility pilot (Supplement C): input conversion and correctness, no timing."""
    v = load(PILOT / 'validation.json')
    assert v['status'] == 'passed' and v['modes'] == ['per_haplotype', 'class', 'auto']
    conv = load(PILOT / 'conversion_provenance.json')
    seeds = load(PILOT / 'seed_diagnostics.json')
    data = load(PILOT / 'data' / 'conversion.json')
    rows = load_jsonl(PILOT / 'validation_rows.jsonl')
    assert conv['nodes'] == seeds['nodes'] == data['vertices'] and conv['edges'] == data['edges']
    assert conv['haplotypes'] == data['H'] and conv['excluded_micb_edges'] == 0
    assert conv['excluded_biological_labels'] == 1 and len(conv['excluded_ambiguous_walks']) == 2  # 'both of its walks'
    assert conv['sequence_lengths']['1'] == seeds['one_base_nodes'] and seeds['anchored_one_base_nodes'] == 0
    assert len(rows) == v['instances'] == data['orientation_instances']
    assert sum(r['N'] for r in rows) == v['anchors'] == data['anchors'] == data['exact_sequence_checks']
    assert sum(r['K'] for r in rows) == v['states']
    for r in rows:
        assert r['guard']['auto_fallback_available'] and r['guard']['auto_work'] <= r['guard']['auto_budget']
    nonempty = [r for r in rows if r['N'] > 0]
    assert len(nonempty) == v['nonempty_instances'] == data['nonempty_instances']
    fallbacks = [r for r in rows if r['guard']['fallback']]
    assert len(fallbacks) == v['fallbacks']
    tiny = [r for r in nonempty if r['N'] < 10]
    small = [r for r in nonempty if 10 <= r['N'] < 30]
    assert all(r['guard']['fallback'] for r in tiny)  # the text says all of these fall back
    assert not any(r['guard']['fallback'] for r in nonempty if r['N'] >= 30)  # and none from 30 anchors
    small_fb = sum(r['guard']['fallback'] for r in small)
    assert len(fallbacks) == len(tiny) + small_fb
    charge = 4096
    fixed = [r for r in tiny if r['guard']['auto_denied_work'] == charge
             and r['guard']['auto_budget'] - r['guard']['auto_work'] < charge]
    m = re.search(r'chr6:(\d+)-(\d+)', conv['locus'])
    assert m
    src = 'baselevel_mhc_pilot'
    put('pilot.kb', f'{(int(m.group(2)) - int(m.group(1))) / 1000:.0f}', src + '/conversion_provenance.json locus interval (kb)')
    put('pilot.nodes', fint(conv['nodes']), src + '/conversion_provenance.json nodes')
    put('pilot.edges', fint(conv['edges']), src + '/conversion_provenance.json edges')
    put('pilot.onebase', fint(seeds['one_base_nodes']), src + '/seed_diagnostics.json one_base_nodes')
    put('pilot.labels', str(conv['haplotypes']), src + '/conversion_provenance.json haplotypes (biological labels kept)')
    put('pilot.excluded.walks', str(len(conv['excluded_ambiguous_walks'])), src + '/conversion_provenance.json excluded_ambiguous_walks')
    put('pilot.reads', fint(data['reads']), src + '/data/conversion.json reads')
    put('pilot.nonempty', fint(v['nonempty_instances']), src + '/validation.json nonempty_instances')
    put('pilot.anchors', fint(v['anchors']), src + '/validation.json anchors (all sequence-checked)')
    put('pilot.states', fint(v['states']), src + '/validation.json states')
    put('pilot.anchored.nodes', str(seeds['anchored_nodes']), src + '/seed_diagnostics.json anchored_nodes')
    put('pilot.fallbacks', fint(len(fallbacks)), src + '/validation_rows.jsonl fallback instances')
    put('pilot.fb.tiny', fint(len(tiny)), src + '/validation_rows.jsonl nonempty instances with fewer than 10 anchors (all fall back)')
    put('pilot.fb.tiny.fixed', fint(len(fixed)), src + f'/validation_rows.jsonl of these, refused one {charge}-unit charge exceeding the remaining budget')
    put('pilot.fixed.charge', fint(charge), src + ' the refused charge (auto_denied_work)')
    put('pilot.small', fint(len(small)), src + '/validation_rows.jsonl nonempty instances with 10-29 anchors')
    put('pilot.fb.small', fint(small_fb), src + '/validation_rows.jsonl of these, fallbacks')


def origin_example():
    """Calculated illustration in Figure 1, the introduction and Section 3; not a measurement.

    Point anchors of equal weight sit on unit-length vertices u and v, on either side of alleles
    of lengths 2 (haplotype h1) and 5 (h2). The competing record of Figure 1 is a restart at a
    further point anchor on v at read position 1. results/paper-framing-draft-20261019/
    check_framing_example.py runs this instance through the specification solver and the direct
    recurrence.
    """
    source_weight = target_weight = 20
    source_query, target_query = 1, 8
    gamma = 100
    competitor_weight, competitor_query = 38, 1
    alpha = target_weight - 1 - target_query + 2
    normalized = source_weight + alpha + 1 + source_query
    key = normalized + 1 + target_query               # published coefficient at the second anchor
    competitor_key = competitor_weight + 1 + competitor_query  # restart at v, so origin v
    routes = {}
    for label, allele_length in [('short', 2), ('long', 5)]:
        route = 1 + allele_length
        gap = route + target_query - source_query - 2
        score = source_weight + target_weight - gap
        assert score == normalized - route and score > target_weight
        routes[label] = route
        put('example.allele.' + label, allele_length, 'Figure 1 / Section 3 illustrative allele length')
        put('example.route.' + label, route, 'Section 3 illustrative graph, calculated from vertex lengths')
        put('example.score.' + label, score, 'Section 3 illustrative chain, weights minus exact gap')
    assert gamma > 2 * target_weight  # the chosen switch penalty cannot improve the example
    # Keys offered to a later anchor on label h: the shared record gives key + phi_h(u) = key,
    # the competitor competitor_key + phi_h(v) = competitor_key + route.
    lead = {label: key - (competitor_key + route) for label, route in routes.items()}
    assert lead['short'] > 0 > lead['long']  # the winning record depends on the haplotype
    put('example.weight', source_weight, 'Figure 1 / Section 3 illustrative anchor weight')
    put('example.read.first', source_query, 'Figure 1 / Section 3 illustrative read position, first anchor')
    put('example.read.second', target_query, 'Figure 1 / Section 3 illustrative read position, second anchor')
    put('example.gamma', gamma, 'Figure 1 / Section 3 illustrative switch penalty')
    put('example.normalized', normalized, 'Section 3 illustrative normalized coefficient')
    put('example.key', key, 'Section 3 illustrative publication coefficient')
    put('example.compete.margin', key - competitor_key,
        'Figure 1 / introduction: shared record minus competing restart at v, published coefficients')
    put('example.compete.short', lead['short'], 'Figure 1 / introduction: shared record lead on h1')
    put('example.compete.long', -lead['long'], 'Figure 1 / introduction: shared record deficit on h2')


# ---------------------------------------------------------------------------------------------
# Exactness against the single-record heuristic (Table 3)
# ---------------------------------------------------------------------------------------------
def giraffe_preset() -> dict:
    """The vg v1.77.0 hifi parameters as the emulator defines them (the HIFI literal in emul.py)."""
    import ast
    path = EXACT / 'emul.py'
    inputs_used.add(path)
    for node in ast.parse(path.read_text()).body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'HIFI' for t in node.targets):
            assert isinstance(node.value, ast.Call) and node.value.func.id == 'dict'
            return {kw.arg: ast.literal_eval(kw.value) for kw in node.value.keywords}
    raise AssertionError('HIFI preset not found in emul.py')


def exactness():
    p2 = load(EXACT / 'phase2_summary.json')
    datasets = list(p2)
    assert len(datasets) == 10
    # Both settings of Giraffe's consistency bonus, and the hifi preset they belong to (Table S3).
    bonuses = {tuple(load(EXACT / f'phase2_{d}.jsonl.meta.json')['bonuses']) for d in datasets}
    assert bonuses == {(0, 12)}, bonuses
    hifi = giraffe_preset()
    assert hifi['consistency_bonus'] == 12 and hifi['match'] == 1
    put('exact.bonus.off', '0', 'phase2_<dataset>.jsonl.meta.json bonuses[0]')
    put('exact.bonus.hifi', str(hifi['consistency_bonus']), 'emul.py HIFI consistency_bonus (= phase2 meta bonuses[1])')
    put('giraffe.hifi.rho', str(hifi['rho']), 'emul.py HIFI rho')
    put('giraffe.hifi.itembonus', str(hifi['item_bonus']), 'emul.py HIFI item_bonus')
    put('giraffe.hifi.gapscale', str(hifi['gap_scale']), 'emul.py HIFI gap_scale')
    rows = []
    for rho in ('2', '28'):
        inst = sum(p2[d][rho]['instances'] for d in datasets)
        reads = p2[datasets[0]][rho]['reads']
        per_bonus = [sum(p2[d][rho][b]['instances_with_gap'] for d in datasets) for b in ('h0', 'h12')]
        gaps = max(per_bonus)
        rel = [p2[d][rho][b]['gap_rel'] for d in datasets for b in ('h0', 'h12') if p2[d][rho][b]['gap_rel'] is not None]
        maxrel = max((x['max'] if isinstance(x, dict) else x) for x in rel) if rel else 0.0
        put(f'exact.real.rho{rho}.instances', fint(inst), f'phase2_summary sum of instances, rho={rho}')
        put(f'exact.real.rho{rho}.below', str(gaps), f'phase2_summary max over bonus of summed instances_with_gap, rho={rho}')
        put(f'exact.real.rho{rho}.maxloss', f'{100 * maxrel:.3f}', f'phase2_summary max gap_rel, rho={rho} (%)')
        rows.append(('Real reads' if rho == '2' else '', rho + (' (default)' if rho == '2' else ''), fint(inst),
                     ' / '.join(map(str, per_bonus)), f'{100 * maxrel:.3f}\\%' if gaps else '0', '--', '--'))
    put('exact.real.reads', str(reads), 'phase2_summary reads per dataset')
    put('exact.real.datasets', str(len(datasets)), 'phase2_summary datasets')
    sw = {p2[d]['2']['exact_with_switch'] for d in ('hprc60', 'mhc60')}
    assert len(sw) == 1
    put('exact.real.switch', str(sw.pop()), 'phase2_summary exact_with_switch (60-haplotype panels, rho=2)')
    put('exact.real.perpanel', str(p2['hprc60']['2']['instances']), 'phase2_summary instances per dataset')

    p3 = load(EXACT / 'phase3_summary.json')['phase3_hprc60_inf']
    allm = p3['all']['methods']
    put('exact.sim.reads', str(p3['all']['reads']), 'phase3 informative set reads')
    put('exact.sim.needswitch', str(p3['all']['truth_needs_switch']), 'phase3 informative set truth_needs_switch')
    for rho in ('2', '56'):
        subs = [p3['paired'][f'rho{rho}_b{b}']['suboptimal']['reads'] for b in (0, 12)]
        sub, sub_lo = max(subs), min(subs)
        we = allm[f'exact_rho{rho}']['walk_accuracy_pooled']
        wg = [allm[f'giraffe_b{b}_rho{rho}']['walk_accuracy_pooled'] for b in (0, 12)]
        put(f'exact.sim.rho{rho}.below', str(sub), f'phase3 paired rho{rho} suboptimal reads (max over bonus)')
        put(f'exact.sim.rho{rho}.below.lo', str(sub_lo), f'phase3 paired rho{rho} suboptimal reads (min over bonus)')
        put(f'exact.sim.rho{rho}.walk.exact', f'{we:.3f}', f'phase3 exact_rho{rho} walk_accuracy_pooled')
        put(f'exact.sim.rho{rho}.walk.heur', f'{min(wg):.3f}', f'phase3 giraffe rho{rho} walk_accuracy_pooled (min over bonus)')
        rows.append(('Simulated, switch needed' if rho == '2' else '', rho + (' (default)' if rho == '2' else ''),
                     str(p3['all']['reads']), ' / '.join(map(str, subs)), '--', f'{we:.3f}',
                     ' / '.join(f'{w:.3f}' for w in wg)))
    put('exact.sim.rho56.below.pct', pct(int(numbers['exact.sim.rho56.below'][0]), p3['all']['reads']), 'phase3 rho56 suboptimal / reads (%)')
    # Coverage of the chains in Table S3: walk accuracy is pooled over truly placed chain anchors, so show
    # that exact and heuristic chains select comparable anchor sets.
    shown = [f'exact_rho{r}' for r in ('2', '56')] + [f'giraffe_b{b}_rho{r}' for r in ('2', '56') for b in (0, 12)]
    lens = [allm[m]['chain_len']['mean'] for m in shown]
    spans = [allm[m]['span_cov']['mean'] for m in shown]
    placed = [allm[m]['placed_frac']['mean'] for m in shown]
    put('exact.sim.chainlen.lo', fint(round(min(lens))), 'phase3 all methods chain_len mean, min over Table S3 methods')
    put('exact.sim.chainlen.hi', fint(round(max(lens))), 'phase3 all methods chain_len mean, max over Table S3 methods')
    put('exact.sim.spancov', f1(100 * min(spans)), 'phase3 all methods span_cov mean, min over Table S3 methods (%)')
    put('exact.sim.placed.lo', f1(100 * min(placed)), 'phase3 all methods placed_frac mean, min over Table S3 methods (%)')
    put('exact.sim.placed.hi', f1(100 * max(placed)), 'phase3 all methods placed_frac mean, max over Table S3 methods (%)')

    p4 = load(EXACT / 'phase4_summary.json')
    drops = [p4[k]['paired_vs_full']['exact_rho2']['reads_with_drop'] for k in p4 if k != 'full']
    walks = [p4[k]['paired_vs_full']['exact_rho2']['walk_acc_subset_mean'] for k in p4 if k != 'full']
    full_walk = {p4[k]['paired_vs_full']['exact_rho2']['walk_acc_full_mean'] for k in p4 if k != 'full'}
    assert full_walk == {1.0}
    put('exact.sample.reads', str(p4['s16_r0']['reads']), 'phase4 reads')
    put('exact.sample.drop.lo', str(min(drops)), 'phase4 exact_rho2 reads_with_drop, min over subsets')
    put('exact.sample.drop.hi', str(max(drops)), 'phase4 exact_rho2 reads_with_drop, max over subsets')
    put('exact.sample.walk.lo', f'{min(walks):.3f}', 'phase4 exact_rho2 walk_acc_subset_mean, min')
    put('exact.sample.walk.hi', f'{max(walks):.3f}', 'phase4 exact_rho2 walk_acc_subset_mean, max')

    with open(OUT / 'rows_exactness.tex', 'w') as fh:
        fh.write('% Generated by paper/scripts/make_paper_data.py; do not edit.\n')
        fh.write('% Data (panel sizes) & rho & instances & below exact & max loss & walk acc exact & walk acc heuristic\n')
        for i, r in enumerate(rows):
            if i and r[0]:
                fh.write('\\addlinespace\n')
            fh.write(' & '.join(r) + ' \\\\\n')


def semantics():
    """Small instances on which other readings of the problem change the optimum (Supplement A.9)."""
    # Since 2026-10-20 the prose example is the covered instance: a third haplotype carries the edge that
    # lay on no haplotype, as CGJ's graph model requires. The original instance's result is unchanged and
    # must give the same two values.
    prose = load(MATCHED / 'problem5_prose_example_covered.json')
    assert prose['uncovered_edges'] == [] and prose['uncovered_vertices'] == []
    assert prose['contract_B'][1] == prose['contract_R'][1] == prose['direct_oracle_B'][1]
    original = load(MATCHED / 'problem5_prose_example.json')
    assert original['contract_B'][1] == prose['contract_B'][1]
    assert original['prose_problem5_value_at_i_h1'] == prose['prose_problem5_value_at_i_h1']
    put('prose.alg2', str(prose['contract_B'][1]), 'problem5_prose_example_covered.json contract_B[1] (Algorithm 2)')
    put('prose.literal', str(prose['prose_problem5_value_at_i_h1']), 'problem5_prose_example_covered.json prose_problem5_value_at_i_h1')
    examples = {e['example'][:3]: e for e in load(EXP / 'proof_review/results/semantics_examples.json')}
    exit_, between = examples['(d)'], examples['(c)']
    assert exit_['gamma'] == between['gamma']
    put('sem.gamma', str(exit_['gamma']), 'semantics_examples.json gamma')
    put('sem.exit.contract', str(exit_['contract_C_of_anchor1']), 'semantics_examples.json (d) contract value')
    put('sem.exit.shortest', str(exit_['shortest_exit_value']), 'semantics_examples.json (d) shortest-exit reading')
    put('sem.between.contract', str(between['contract_C']['(1, 1)']), 'semantics_examples.json (c) contract C(i,h)')
    put('sem.between.anywhere', str(between['switch_anywhere_value']), 'semantics_examples.json (c) switch-anywhere reading')


def realism():
    """Simulated against real CHM13 reads on the MHC (Section 6.1)."""
    r = load(EXP / 'second_locus/realism_check.json')
    assert set(r['panels']) == {'hprc60', 'hprc463'}
    distances = []
    for p in r['panels'].values():
        assert p['reads_real'] == p['reads_simulated'] == 100
        for section in ('per_read', 'class_counters_major_orientation'):
            distances += [v['ks_distance'] for k, v in p[section].items() if k != 'qlen']
    put('realism.ks.max', f'{max(distances):.2f}', 'second_locus/realism_check.json max KS distance over both panels, read length excluded')
    q = r['panels']['hprc463']['per_read']['qlen']
    assert q == r['panels']['hprc60']['per_read']['qlen']
    put('realism.readlen.real', f1(q['real'][1] / 1000), 'second_locus/realism_check.json median real read length (kb)')
    put('realism.readlen.sim', f1(q['simulated'][1] / 1000), 'second_locus/realism_check.json median simulated read length (kb)')
    put('realism.readlen.ks', f'{q["ks_distance"]:.2f}', 'second_locus/realism_check.json read-length KS distance')


def context():
    """Panel source and Minichain's end-to-end time, for the introduction."""
    panels = load(EXP / 'hprc_mhc_panel/panels.json')
    assert panels['kept'] == 463 and panels['samples'] - panels['kept'] == sum(panels['excluded'].values())
    put('graph.haplotypes', str(panels['samples']), 'hprc_mhc_panel/panels.json samples (haplotypes in the source graph)')
    put('mhc.region.mb', f'{panels["region"]["length"] / 1e6:.2f}', 'hprc_mhc_panel/panels.json region.length (Mb)')
    log = EXP / 'minichain_matched/data/pinned_end_to_end_R10000.log'
    inputs_used.add(log)
    m = re.search(r'\[M::main\] Real time: ([0-9.]+) sec', log.read_text())
    assert m, 'no Real time line in the Minichain log'
    total = float(m.group(1))
    chaining = load(MATCHED / 'timing_mhc2104_R10000.json')['totals_seconds']['minichain_chaining']
    assert chaining < total
    put('mc.endtoend.total', f1(total), 'pinned Minichain end-to-end on 2,104 reads, single thread (log "Real time")')
    put('mc.chain.share', f'{100 * chaining / total:.0f}',
        'Minichain Chaining() total (timing runs) / end-to-end Real time (separate run), percent')


def machine():
    """Host and load during every timing run used in the paper (from the runs' meta files)."""
    loads, libs, compilers = [], set(), set()
    for path in TIMING_META:
        meta = load(path)
        loads += [meta['load_average_before'][0], meta['load_average_after'][0]]
        libs.add(Path(meta['library']).name)
        compilers.add(meta['compiler'])
        assert meta['machine'] == 'arm64', path
    assert len(libs) == 1, f'timing runs used different native builds: {libs}'
    assert len(compilers) == 1, compilers
    put('load.lo', f1(min(loads)), 'min 1-minute load average before/after the timing runs (meta.json)')
    put('load.hi', f1(max(loads)), 'max 1-minute load average before/after the timing runs (meta.json)')
    put('native.library', libs.pop().replace('.dylib', ''), 'native library used by the original class timing runs (meta.json)')


def write_numbers():
    for key, (value, source) in TEXT_SOURCED.items():
        put(key, value, 'TEXT-SOURCED: ' + source)
    with open(OUT / 'numbers.tex', 'w') as fh:
        fh.write('% Generated by paper/scripts/make_paper_data.py from saved results; do not edit.\n')
        fh.write('% Each line: \\dataset{key}{value} % source. Use \\data{key} in the text.\n')
        for key in sorted(numbers):
            value, source = numbers[key]
            assert '%' not in value and '\\' not in value.replace('\\,', ''), (key, value)
            fh.write(f'\\dataset{{{key}}}{{{value}}}% {source}\n')


def write_provenance():
    def sha(p: Path) -> str:
        return hashlib.sha256(p.read_bytes()).hexdigest()
    prov = {
        'script': str(Path(__file__).resolve().relative_to(ROOT)),
        'script_sha256': sha(Path(__file__).resolve()),
        'python': sys.version.split()[0],
        'platform': platform.platform(),
        'generated_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'inputs': {str(p.relative_to(ROOT)): sha(p) for p in sorted(inputs_used)},
        'text_sourced': {k: v[1] for k, v in TEXT_SOURCED.items()},
        'outputs': sorted(p.name for p in OUT.iterdir() if p.name != 'provenance.json'),
    }
    with open(OUT / 'provenance.json', 'w') as fh:
        json.dump(prov, fh, indent=1)
        fh.write('\n')


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    identity()
    matched_timing()
    loci()
    fallback()
    safeguard_timing()
    b100 = bench100()
    fallback_cost(b100)
    first, first_loads, _, _ = bench_results('first')
    bench_first_run(b100, first, first_loads)
    mechanism(first)
    memory_ablation(b100)
    integration()
    baselevel()
    micb_pilot()
    origin_example()
    exactness()
    semantics()
    realism()
    context()
    machine()
    write_numbers()
    write_provenance()
    print(f'wrote {len(numbers)} numbers and {len(list(OUT.iterdir()))} files to {OUT.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
