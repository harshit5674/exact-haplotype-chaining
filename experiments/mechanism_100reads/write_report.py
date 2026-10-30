"""Write REPORT.md and the Figure S1 preview TeX from results.json and posthoc.json.

Every number below comes from those files; each qualitative claim in the prose is asserted.
  python3.13 -B experiments/mechanism_100reads/write_report.py
"""
from decimal import ROUND_HALF_UP, Decimal
import hashlib
import json
from pathlib import Path
import re
import statistics as st

HERE = Path(__file__).resolve().parent
R = json.loads((HERE/'results.json').read_text())
PH = json.loads((HERE/'posthoc.json').read_text())
AUDIT = json.loads((HERE/'archive_audit.json').read_text())
EDIT = json.loads((HERE/'paper_edit'/'verification.json').read_text())
assert EDIT['status'] == 'passed' and EDIT['body_end_page'] <= 10 and EDIT['archive_audit_after_edit'] == 'passed'
assert AUDIT['status'] == 'passed' and not AUDIT['working_code_changed_since_import']
SUMMARY = json.loads((HERE.parents[0]/'benchmark_100reads'/'summary.json').read_text())
mine = {k.split('timings/')[-1]: v for k, v in R['inputs_sha256'].items() if '/timings/' in k}
theirs = {k.split('timings/')[-1]: v for k, v in SUMMARY['input_sha256'].items()}
assert set(mine) == set(theirs) and all(mine[k] == theirs[k] for k in mine)
PAPER = dict(re.findall(r'\\dataset\{([^}]*)\}\{([^}]*)\}', (HERE.parents[1]/'paper/generated/numbers.tex').read_text()))
LOCI = ('mhc', 'kir', 'chr20')
NAME = {'mhc': 'MHC', 'kir': 'KIR', 'chr20': 'chr20'}
LARGEST = {'mhc': 'mhc-hprc463', 'kir': 'kir-hprc459', 'chr20': 'chr20-hprc463'}
PANELS = [c['panel'] for c in R['checks']]
H = {k: v['H'] for k, v in R['secondary']['S5_reliability'].items()}
MINUS = '\u2212'
CLIP = 4.2  # x-axis limit of the proposed figure; reads beyond it are drawn at the edge


def num(x, d=2):
    return f'{x:.{d}f}'.replace('-', MINUS)


def half_up(x, d=2):
    # The paper generator's rounding (half up from the decimal repr).
    return str(Decimal(repr(x)).quantize(Decimal(1).scaleb(-d), rounding=ROUND_HALF_UP))


def tex(x, d=2):
    s = half_up(x, d)
    return '$-$' + s[1:] if s.startswith('-') else s


def disp(x, d=2):
    # The generator's rounding, shown with a Unicode minus for Markdown.
    return half_up(x, d).replace('-', MINUS)


def ci(e, d=2):
    return f"[{num(e['ci95'][0], d)}, {num(e['ci95'][1], d)}]"


def pv(p):
    return '<0.001' if p < 0.001 else f'{p:.3f}'


def label(key):
    locus, name = key.split('-', 1)
    return NAME[locus], f"{H[key]}{' LO' if 'leaveout' in name else ''}"


def verdict(e):
    return {'tracks': 'tracks', 'no detectable association': 'none detected',
            'opposite to prediction': 'opposite', 'undefined': 'undefined'}[e['verdict']]


P = R['primary']['tests']
S1, S2 = R['secondary']['S1_compact']['tests'], R['secondary']['S2_class']['tests']
S3, S4, S5 = (R['secondary'][k] for k in ('S3_kir_without_fallback_reads', 'S4_bridge', 'S5_reliability'))
E = R['exploratory']
F = R['fallback']


def t(fam, outcome, key, x):
    return fam[f'P1/{outcome}/{key}/{x}']


def big(locus, x='x_rec', fam=P, outcome='guard'):
    return t(fam, outcome, LARGEST[locus], x)


# ---------------------------------------------------------------- claims made in the prose
mhckir = [k for k in PANELS if not k.startswith('chr20')]
chr20 = [k for k in PANELS if k.startswith('chr20')]
chr20_small = ['chr20-hprc60', 'chr20-hprc120']
chr20_big = [k for k in chr20 if k not in chr20_small]
assert all(t(P, 'guard', k, x)['verdict'] == 'tracks' for k in mhckir for x in ('x_rec', 'x_sl'))
assert all(t(P, 'guard', k, 'x_rec')['verdict'] == 'tracks' for k in chr20_big)
assert all(t(P, 'guard', k, 'x_rec')['verdict'] != 'tracks' for k in chr20_small)
assert all(t(P, 'guard', k, 'x_sl')['verdict'] == 'no detectable association' for k in chr20)
assert all(v['verdict'] == 'tracks' for k, v in P.items() if k.startswith(('P2', 'P3')))
assert all(t(P, 'guard', k, 'x_rec')['magnitude'] == 'moderate' for k in chr20_big)
assert R['decisions']['section_6_5_sentence'] == 'supported at all three largest panels'
assert all(v['verdict'] == 'tracks' for v in S1.values())
assert [k for k, v in S2.items() if v['verdict'] != 'tracks'] == ['P1/class/chr20-hprc120/x_sl']
assert all(v['p'] < 0.001 for v in S3.values())
assert all(v['discarded'] == 0 for fam in (P, S1, S2, S3) for v in fam.values())
assert all(v['rho'] < 0 for fam in (P, S1, S2, S3) for v in fam.values() if v['p'] is not None and v['p'] < .05)
min_p = 1/(1 + R['permutations'])
assert all(abs(v['p_holm'] - 32*min_p) < 1e-12 for k, v in P.items() if v['verdict'] == 'tracks')
rel_g = [S5[k]['split_half_guard'] for k in PANELS]
rel_c = [S5[k]['split_half_compact'] for k in PANELS]
partial_rec_N = {s: E['E2'][f'E2/{s}/y_guard~x_rec|N'] for s in LOCI + ('pooled',)}
assert all(abs(v['rho']) >= .75 for v in partial_rec_N.values())
assert all(abs(E['E2'][f'E2/{s}/y_guard~x_rec|{z}']['rho']) >= .5 for s in LOCI + ('pooled',) for z in ('N', 's'))
assert E['E1']['E1/chr20/y_guard~N']['rho'] > 0 and E['E1']['E1/chr20/x_rec~N']['rho'] > 0
assert abs(partial_rec_N['chr20']['rho']) > abs(t(P, 'guard', 'chr20-hprc463', 'x_rec')['rho'])
assert all(E['E1'][f'E1/{s}/x_rec~N']['rho'] > 0 and E['E1'][f'E1/{s}/x_rec~s']['rho'] < 0 for s in LOCI)
assert all(abs(S3[f'S3/guard/{k}/x_rec']['rho']) >= .69 and
           0 <= abs(t(P, 'guard', k, 'x_rec')['rho']) - abs(S3[f'S3/guard/{k}/x_rec']['rho']) <= .06
           for k in (c['panel'] for c in R['checks'] if c['panel'].startswith('kir')))
w_guard = {s: E['E3'][f'E3/{s}/y_guard~w'] for s in LOCI + ('pooled',)}
assert all(v['rho'] <= -.9 for v in w_guard.values())
fs = F['F1_summary']
rows = F['F1']
assert fs['components'] == len(rows) == 59 and fs['stage'] == {'output': 18, 'sweep': 41}
by_panel = {}
for r in rows:
    by_panel.setdefault(r['panel'], []).append(r)
kir_panels = list(by_panel)
h459 = by_panel['kir-hprc459']
full_sweep = [r for r in h459 if r['excess'] >= r['predicted_excess_validation']]
early = sorted((r for r in h459 if r['excess'] < r['predicted_excess_validation']), key=lambda r: r['id'])
assert len(full_sweep) == 10 and [r['id'] for r in early] == ['kir_r0024:+', 'kir_r0095:-']
assert all(r['records_per_anchor'] > 7 for r in early)
heavy = {'kir_r0024:+', 'kir_r0095:-'}
early_ratio = [r['excess']/r['predicted_excess_validation'] for r in early]
assert all(.3 <= x <= .4 for x in early_ratio)
ratio = lambda r: r['excess']/r['predicted_excess_validation']
low = [r for r in rows if ratio(r) < .8]
low_other = [r for r in low if r['id'] not in heavy]
assert sum(r['id'] in heavy for r in low) == 10 == sum(r['id'] in heavy for r in rows)
assert [(r['id'], r['panel']) for r in low_other] == [('kir_r0016:-', 'kir-hprc453_leaveout')]
assert all(r['records_per_anchor'] > 4.5 for r in rows if r['id'] in heavy)
others_min = min(ratio(r) for r in rows if r not in low)
nonheavy_min = min(ratio(r) for r in rows if r['id'] not in heavy)
assert others_min >= .995
aur = {k: F['F2'][k]['records_per_anchor']['auroc'] for k in kir_panels}
aur_sl = {k: F['F2'][k]['sliced_per_anchor']['auroc'] for k in kir_panels}
aur_s = {k: F['F2'][k]['states_per_anchor']['auroc'] for k in kir_panels}
aur_c = {k: F['F2'][k]['class_over_ph']['auroc'] for k in kir_panels}
assert all(v < .5 for v in aur_s.values())
rec_ratio = [F['F2'][k]['records_per_anchor']['fallback']['median']/F['F2'][k]['records_per_anchor']['other']['median'] for k in kir_panels]
small_slow = {k: F['F2'][k]['two_by_two'] for k in ('kir-hprc60', 'kir-hprc120')}
assert all(v['other_class_slower'] > v['other_class_faster'] for v in small_slow.values())
unneeded = [r['class_over_ph'] for r in rows if r['class_over_ph'] <= 1]
assert len(unneeded) == 59 - fs['class_slower_than_ph']
imp = F['F3']
rec = F['F4']
assert rec['h453_h459_shared'] == 12 == len(h459) == len(by_panel['kir-hprc453_leaveout'])
assert all(imp[k]['fallback_reads'] == imp[k]['fallback_components'] for k in kir_panels)
assert all(r['admitted_over_budget'] > .9 for r in rows)
bridge = S4
assert bridge['a_old_benchmark_120_instances']['record_count_mismatches'] == 0
clipped = [r for r in R['figure']['kir'] if r['records'] > CLIP]
assert len(clipped) == 2 and all(r['fallback'] for r in clipped)
assert all(r['records'] <= CLIP for locus in ('mhc', 'chr20') for r in R['figure'][locus])
lo_hi = lambda xs: (min(xs), max(xs))
drift = [v['spread_new_over_old_compact_speedup'] for v in PH['loci'].values()]
assert all(abs(v['median_new_over_old_per_haplotype_time'] - 1) > abs(v['median_new_over_old_compact_time'] - 1) for v in PH['loci'].values())
ph_drift = [v['median_new_over_old_per_haplotype_time'] for v in PH['loci'].values()]
totals_other = {l: float(PAPER[f'b100.{l}.h463.total.guard']) for l in ('mhc', 'chr20')}
assert F['F3']['kir-hprc459']['total_without_fallback_reads'] < min(totals_other.values())


# ---------------------------------------------------------------- tables
def primary_table():
    lines = ['| Locus | H | Records/anchor: ρ [95% CI] | Holm p | | Bit-sliced/anchor: ρ [95% CI] | Holm p | |',
             '| --- | ---: | ---: | ---: | --- | ---: | ---: | --- |']
    for k in PANELS:
        a, b = t(P, 'guard', k, 'x_rec'), t(P, 'guard', k, 'x_sl')
        lines.append(f"| {' | '.join(label(k))} | {num(a['rho'])} {ci(a)} | {pv(a['p_holm'])} | {verdict(a)} | "
                     f"{num(b['rho'])} {ci(b)} | {pv(b['p_holm'])} | {verdict(b)} |")
    for kind, name in (('P2', 'Pooled, 300 reads'), ('P3', 'Mean within locus')):
        a, b = P[f'{kind}/guard/x_rec'], P[f'{kind}/guard/x_sl']
        lines.append(f"| {name} | 463/459/463 | {num(a['rho'])} {ci(a)} | {pv(a['p_holm'])} | {verdict(a)} | "
                     f"{num(b['rho'])} {ci(b)} | {pv(b['p_holm'])} | {verdict(b)} |")
    return '\n'.join(lines)


def secondary_table():
    def cell(fam, outcome, k, x):
        e = t(fam, outcome, k, x)
        return num(e['rho']) + ('' if e['verdict'] == 'tracks' else ' (n.d.)')
    lines = ['| Locus | H | Compact: records | Compact: bit-sliced | Class: records | Class: bit-sliced |',
             '| --- | ---: | ---: | ---: | ---: | ---: |']
    for k in PANELS:
        lines.append(f"| {' | '.join(label(k))} | {cell(S1, 'compact', k, 'x_rec')} | {cell(S1, 'compact', k, 'x_sl')} | "
                     f"{cell(S2, 'class', k, 'x_rec')} | {cell(S2, 'class', k, 'x_sl')} |")
    for kind, name in (('P2', 'Pooled'), ('P3', 'Mean within locus')):
        lines.append(f"| {name} | | {num(S1[f'{kind}/compact/x_rec']['rho'])} | {num(S1[f'{kind}/compact/x_sl']['rho'])} | "
                     f"{num(S2[f'{kind}/class/x_rec']['rho'])} | {num(S2[f'{kind}/class/x_sl']['rho'])} |")
    return '\n'.join(lines)


def s3_table():
    lines = ['| H | Reads | Records/anchor: ρ [95% CI] | Bit-sliced/anchor: ρ [95% CI] | All-read records ρ |',
             '| ---: | ---: | ---: | ---: | ---: |']
    for k in kir_panels:
        a, b = S3[f'S3/guard/{k}/x_rec'], S3[f'S3/guard/{k}/x_sl']
        lines.append(f"| {label(k)[1]} | {a['n']} | {num(a['rho'])} {ci(a)} | {num(b['rho'])} {ci(b)} | {num(t(P, 'guard', k, 'x_rec')['rho'])} |")
    return '\n'.join(lines)


def bridge_table():
    names = [('a_old_benchmark_120_instances', 'Old benchmark: 120 instances, per instance, compact'),
             ('b_new_timings_earlier_20_reads', 'New timings: the same 120 instances'),
             ('c_new_timings_all_components', 'New timings: all 600 instances'),
             ('d_read_level_compact', 'New timings: 300 reads, compact'),
             ('e_read_level_guarded', 'New timings: 300 reads, guarded (primary)')]
    lines = ['| Unit and outcome | Pooled ρ | MHC | KIR | chr20 |', '| --- | ---: | ---: | ---: | ---: |']
    for key, name in names:
        b = bridge[key]
        lines.append(f"| {name} | {num(b['pooled'])} | {' | '.join(num(b['within'][l]) for l in LOCI)} |")
    return '\n'.join(lines)


def posthoc_table():
    lines = ['| Locus | Instances | Spearman, old vs new compact speedup | New/old speedup: median [min, max] | New/old per-hap. time | New/old compact time |',
             '| --- | ---: | ---: | ---: | ---: | ---: |']
    for l in LOCI:
        v = PH['loci'][l]
        sp = v['spread_new_over_old_compact_speedup']
        lines.append(f"| {NAME[l]} {v['H']} | {v['components']} | {num(v['spearman_old_vs_new_compact_speedup'])} | "
                     f"{num(sp['median'])} [{num(sp['min'])}, {num(sp['max'])}] | {num(v['median_new_over_old_per_haplotype_time'])} | "
                     f"{num(v['median_new_over_old_compact_time'])} |")
    return '\n'.join(lines)


def reliability_table():
    lines = ['| Locus | H | Split-half ρ, guarded | Split-half ρ, compact | Guarded speedup: median [IQR] | Records/anchor: median [IQR] | Fallback reads |',
             '| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for k in PANELS:
        v = S5[k]
        g, x = v['speedup_guard'], v['x_rec']
        lines.append(f"| {' | '.join(label(k))} | {num(v['split_half_guard'], 3)} | {num(v['split_half_compact'], 3)} | "
                     f"{num(g['median'])} [{num(g['q25'])}, {num(g['q75'])}] | {num(x['median'])} [{num(x['q25'])}, {num(x['q75'])}] | {v['fallback_reads']} |")
    return '\n'.join(lines)


def exploratory_table():
    def row(name, section, pattern):
        cells = []
        for s in LOCI + ('pooled',):
            e = E[section][pattern.format(s=s)]
            cells.append(f"{num(e['rho'])} {ci(e)}")
        return f"| {name} | {' | '.join(cells)} |"
    lines = ['| Statistic | MHC 463 | KIR 459 | chr20 463 | Pooled |', '| --- | ---: | ---: | ---: | ---: |']
    lines += [row('ρ(speedup, anchors N)', 'E1', 'E1/{s}/y_guard~N'),
              row('ρ(speedup, states/anchor)', 'E1', 'E1/{s}/y_guard~s'),
              row('ρ(records, bit-sliced)', 'E1', 'E1/{s}/x_rec~x_sl'),
              row('ρ(records, states/anchor)', 'E1', 'E1/{s}/x_rec~s'),
              row('ρ(bit-sliced, states/anchor)', 'E1', 'E1/{s}/x_sl~s'),
              row('ρ(records, N)', 'E1', 'E1/{s}/x_rec~N'),
              row('ρ(bit-sliced, N)', 'E1', 'E1/{s}/x_sl~N'),
              row('Partial ρ(speedup, records ∣ N)', 'E2', 'E2/{s}/y_guard~x_rec|N'),
              row('Partial ρ(speedup, records ∣ states/anchor)', 'E2', 'E2/{s}/y_guard~x_rec|s'),
              row('Partial ρ(speedup, bit-sliced ∣ N)', 'E2', 'E2/{s}/y_guard~x_sl|N'),
              row('Partial ρ(speedup, bit-sliced ∣ states/anchor)', 'E2', 'E2/{s}/y_guard~x_sl|s'),
              row('ρ(guarded speedup, work ratio w)', 'E3', 'E3/{s}/y_guard~w'),
              row('ρ(class speedup, work ratio w)', 'E3', 'E3/{s}/y_class~w')]
    return '\n'.join(lines)


def fallback_panel_table():
    lines = ['| H | Fallbacks | Refused at output | Guarded / per-hap. time: median [min, max] | Class / per-hap. | Excess, measured | Excess, class sweep (validation run) |',
             '| ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    for k in kir_panels:
        rs = by_panel[k]
        g = [r['guarded_over_ph'] for r in rs]
        c = [r['class_over_ph'] for r in rs]
        lines.append(f"| {label(k)[1]} | {len(rs)} | {sum(r['stage'] == 'output' for r in rs)} | {num(st.median(g))} [{num(min(g))}, {num(max(g))}] | "
                     f"{num(st.median(c))} [{num(min(c))}, {num(max(c))}] | {num(st.median(r['excess'] for r in rs))} | "
                     f"{num(st.median(r['predicted_excess_validation'] for r in rs))} |")
    return '\n'.join(lines)


def separation_table():
    lines = ['| H | Records/anchor: fallback vs other (median) | AUROC records | AUROC bit-sliced | AUROC states/anchor | AUROC N | AUROC class/per-hap. | Class slower than per-hap.: fallback / other |',
             '| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    for k in kir_panels:
        v = F['F2'][k]
        tb = v['two_by_two']
        nf, no = tb['fallback_class_slower'] + tb['fallback_class_faster'], tb['other_class_slower'] + tb['other_class_faster']
        lines.append(f"| {label(k)[1]} | {num(v['records_per_anchor']['fallback']['median'])} vs {num(v['records_per_anchor']['other']['median'])} | "
                     f"{num(v['records_per_anchor']['auroc'], 3)} | {num(v['sliced_per_anchor']['auroc'], 3)} | {num(v['states_per_anchor']['auroc'], 3)} | "
                     f"{num(v['N']['auroc'], 3)} | {num(v['class_over_ph']['auroc'], 3)} | {tb['fallback_class_slower']}/{nf} vs {tb['other_class_slower']}/{no} |")
    return '\n'.join(lines)


def impact_table():
    lines = ['| H | Fallback reads | Median speedup: fallback / other reads | Share of guarded time | Total: measured | Total: other reads only | Total: fallbacks at per-hap. time | Total: fallbacks at class time |',
             '| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    for k in kir_panels:
        v = imp[k]
        lines.append(f"| {label(k)[1]} | {v['fallback_reads']} | {num(v['median_speedup_fallback_reads'])} / {num(v['median_speedup_other_reads'])} | "
                     f"{100*v['fallback_share_of_guarded_time']:.0f}% | {num(v['total_measured'])}× | {num(v['total_without_fallback_reads'])}× | "
                     f"{num(v['total_fallbacks_at_per_haplotype_time'])}× | {num(v['total_fallbacks_at_class_time'])}× |")
    return '\n'.join(lines)


def fallback_rows_table():
    lines = ['| H | Instance | N | States/anchor | Records/anchor | Bit-sliced/anchor | Refused at | Admitted / budget | Per-hap. ms | Guarded / per-hap. | Class / per-hap. | Excess | Class sweep / per-hap. |',
             '| ---: | --- | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for r in rows:
        lines.append(f"| {label(r['panel'])[1]} | `{r['id']}` | {r['N']:,} | {r['states_per_anchor']:.0f} | {r['records_per_anchor']:.2f} | "
                     f"{r['sliced_per_anchor']:.2f} | {r['stage']} | {r['admitted_over_budget']:.3f} | {r['per_haplotype_ms']:.1f} | "
                     f"{r['guarded_over_ph']:.2f} | {r['class_over_ph']:.2f} | {r['excess']:.2f} | {r['predicted_excess_validation']:.2f} |")
    return '\n'.join(lines)


# ---------------------------------------------------------------- figure preview
def figure_body(xmax, clip):
    cap = (f"Per-read guarded speedup over the per-haplotype engine, both with expanded output, at the largest panels, "
           f"against compact records per anchor; each read combines its two instances. There are {P['P2/guard/x_rec']['n']} reads. "
           f"Open squares are LRC/KIR reads on which the guard fell back")
    if clip:
        cap += (f"; the two drawn as triangles at the right edge have {tex(clipped[0]['records'], 1)} and "
                f"{tex(clipped[1]['records'], 1)} records per anchor")
    cap += (f". Spearman correlations are {tex(big('kir')['rho'])} on LRC/KIR, {tex(big('mhc')['rho'])} on the MHC, "
            f"{tex(big('chr20')['rho'])} on chr20 and {tex(P['P2/guard/x_rec']['rho'])} pooled. "
            f"The analysis was specified before it was run.")
    restrict = f', restrict x to domain=0:{xmax}' if clip else ''
    edge = ''
    if clip:
        pts = ' '.join(f"({xmax - 0.06},{r['guard']:.4f})" for r in clipped)
        edge = (f"\\addplot[kir, only marks, mark=triangle, mark size=2.2pt, mark options={{rotate=-90, fill=white}}, forget plot] "
                f"coordinates {{{pts}}};\n")
    ticks = ','.join(str(x) for x in ([1, 1.5, 2, 2.5, 3, 3.5, 4] if clip else [1, 2, 3, 4, 5, 6, 7]))
    return f"""% Generated by write_report.py from results.json; preview only.
\\begin{{figure}}[h]
\\centering
\\begin{{tikzpicture}}
\\begin{{axis}}[
  paperaxis, width=10cm, height=6cm,
  title={{Guarded speedup per read}},
  xmin=0.9, xmax={xmax}, ymin=0.5, ymax=3.2, ymode=log,
  xtick={{{ticks}}}, ytick={{0.5,1,2,3}}, yticklabels={{0.5,1,2,3}},
  xlabel={{compact records per anchor}}, ylabel={{speedup, expanded output}},
  legend pos=north east,
]
\\addplot[color=axisline, densely dotted, line width=0.8pt, mark=none, forget plot] coordinates {{(0.9,1) ({xmax},1)}};
\\addplot[chr, only marks, mark size=1.6pt, fill opacity=0.75, draw opacity=0.9] table[x=records, y=guard] {{b100_scatter_chr20.dat}};
\\addlegendentry{{chr20}}
\\addplot[mhc, only marks, mark size=1.3pt, fill opacity=0.75, draw opacity=0.9] table[x=records, y=guard] {{b100_scatter_mhc.dat}};
\\addlegendentry{{MHC}}
\\addplot[kir, only marks, mark size=1.2pt, fill opacity=0.75, draw opacity=0.9, restrict expr to domain={{\\thisrow{{fallback}}}}{{0:0}}] table[x=records, y=guard] {{b100_scatter_kir.dat}};
\\addlegendentry{{LRC/KIR}}
\\addplot[kir, only marks, mark size=1.6pt, mark options={{fill=white}}, restrict expr to domain={{\\thisrow{{fallback}}}}{{1:1}}{restrict}] table[x=records, y=guard] {{b100_scatter_kir.dat}};
\\addlegendentry{{LRC/KIR, fallback}}
{edge}\\end{{axis}}
\\end{{tikzpicture}}
\\caption{{{cap}}}
\\end{{figure}}
"""


(HERE/'preview'/'fig-s1-body.tex').write_text(figure_body(CLIP, True))
(HERE/'preview'/'fig-s1-body-fullrange.tex').write_text(figure_body(7.6, False))


# ---------------------------------------------------------------- report
def recurrence_text():
    counts = rec['components_by_number_of_panels']
    assert list(counts) == ['1', '2', '3', '4', '5']
    words = ['one KIR panel', 'two', 'three', 'four', 'all five']
    parts = [f"{counts[k]} {'instance falls' if counts[k] == 1 else 'instances fall'} back on {words[0]}" if k == '1'
             else f"{counts[k]} on {words[int(k) - 1]}" for k in counts]
    return ', '.join(parts[:-1]) + ' and ' + parts[-1]


def rng_text(keys, x='x_rec', fam=P, outcome='guard'):
    vals = [t(fam, outcome, k, x)['rho'] for k in keys]
    lo, hi = min(vals), max(vals)
    return f'{num(hi)} to {num(lo)}'


fb_ratio = {k: lo_hi([r['guarded_over_ph'] for r in by_panel[k]]) for k in kir_panels}
p2, p3 = P['P2/guard/x_rec'], P['P3/guard/x_rec']
p2s, p3s = P['P2/guard/x_sl'], P['P3/guard/x_sl']
i459 = imp['kir-hprc459']
plan_hash = R['plan_sha256']
chr20_q75 = {k: S5[k]['x_rec']['q75'] for k in chr20_small}
chr20_iqr60 = S5['chr20-hprc60']['speedup_guard']
mean_excess_ratio = fs['excess_over_predicted']

report = f"""# Mechanism analysis of the 100-read benchmark — 2026-10-18

**Status: completed.** This is a pre-specified association analysis of saved data. [PLAN.md](PLAN.md) was frozen before any association between a speedup and a counter was computed ([PLAN.sha256](PLAN.sha256): `{plan_hash[:12]}…`). [analysis.py](analysis.py) checks that hash, and was run once; [results.json](results.json) holds every value. This report and the figure preview are generated from it by [write_report.py](write_report.py). No timing was run. No solver, benchmark file or archive changed. One post hoc diagnostic, written after the results were seen, is marked as such. After the analysis, the user approved the proposed paper edits; they were applied and checked (Section 6.6).

## Answer in brief

- **Speedup tracks records per anchor read by read, within every locus at about 460 haplotypes.** The guarded per-read speedup has Spearman ρ = {num(big('mhc')['rho'])} on MHC, {num(big('kir')['rho'])} on KIR and {num(big('chr20')['rho'])} on chr20 at the largest panels (Holm-adjusted p = {pv(big('mhc')['p_holm'])} each, the floor of the test). By the pre-specified rule, the Section 6.5 sentence "within loci the speedup fell as records per anchor grew" is **supported at all three largest panels**.
- **Across panel sizes:** strong on MHC and KIR at every H ({rng_text([k for k in PANELS if k.startswith('mhc')])} and {rng_text([k for k in PANELS if k.startswith('kir')])}). Moderate on chr20 from H=240 ({rng_text(chr20_big)}). Not detectable on chr20 at H=60 and 120, where the upper quartile of records per anchor is {num(chr20_q75['chr20-hprc60'])} and {num(chr20_q75['chr20-hprc120'])}.
- **Bit-sliced decisions per anchor** track on MHC and KIR at every H ({rng_text(mhckir, 'x_sl')}), but not on chr20 at any H.
- **Across loci,** the pooled correlation at the largest panels is {num(p2['rho'])} {ci(p2)} for records and {num(p2s['rho'])} for bit-sliced decisions. The mean within-locus correlations, {num(p3['rho'])} and {num(p3s['rho'])}, are also significant. So the pooled association is not only a difference between the three loci.
- **Read size does not explain it (exploratory).** Controlling for the number of anchors, the records correlation is {num(partial_rec_N['mhc']['rho'])}, {num(partial_rec_N['kir']['rho'])} and {num(partial_rec_N['chr20']['rho'])}. On chr20 it gets stronger, because longer reads there have both more records and higher speedups. The guard's own work ratio (admitted class work over the per-haplotype estimate) predicts the guarded speedup closely: ρ = {num(w_guard['mhc']['rho'])}, {num(w_guard['kir']['rho'])}, {num(w_guard['chr20']['rho'])} and {num(w_guard['pooled']['rho'])} pooled.
- **The 59 KIR fallbacks are {rec['distinct_components']} read orientations,** one per read, and {rec['components_by_number_of_panels']['5']} of them fall back on all five KIR panels.
  - **What they are.** They are the instances with the most records. Records per anchor separate them from the other instances with AUROC {num(min(aur.values()), 2)}–{num(max(aur.values()), 2)} on every KIR panel.
  - **What they cost.** The guard usually refuses late: apart from the two most record-heavy instances, the time lost (estimate pass plus abandoned attempt) is at least {num(nonheavy_min)} of the class engine's full sweep time, and usually more. So a fallback costs {num(fs['guarded_over_ph']['min'])}–{num(fs['guarded_over_ph']['max'])}× the per-haplotype time: {num(fb_ratio['kir-hprc459'][0])}–{num(fb_ratio['kir-hprc459'][1])}× at H=459 and {num(fb_ratio['kir-hprc60'][0])}–{num(fb_ratio['kir-hprc60'][1])}× at H=60. On {fs['stage']['output']} of the 59 it refused only the final output charge, after the whole sweep.
  - **Effect on KIR 459.** The measured total speedup is {num(i459['total_measured'])}×. Charging each fallback only its own per-haplotype time would give {num(i459['total_fallbacks_at_per_haplotype_time'])}×; charging its unguarded class time would give {num(i459['total_fallbacks_at_class_time'])}×. These are reconstructions from measured times.
- **Bridge to the old figure.** The old −0.78 reproduces exactly. With the new timings, the same 120 instances give {num(bridge['b_new_timings_earlier_20_reads']['pooled'])}, all 600 give {num(bridge['c_new_timings_all_components']['pooled'])}, and the 300 reads give {num(bridge['e_read_level_guarded']['pooled'])}. Post hoc: old and new compact speedups of the same instances have rank correlation only {num(min(v['spearman_old_vs_new_compact_speedup'] for v in PH['loci'].values()))}–{num(max(v['spearman_old_vs_new_compact_speedup'] for v in PH['loci'].values()))} between the two runs, although the new run's split-half reliability for guarded speedup is {num(min(rel_g), 2)}–{num(max(rel_g), 3)}.

## 1. Inputs, checks and methods

- **Data.** The 14 panels of [benchmark_100reads](../benchmark_100reads/REPORT.md), 100 reads each, both orientations summed, five paired repetitions, at the native-kernel boundary. The SHA-256 of every input is in `results.json`.
- **Checks.** All [plan checks](PLAN.md#integrity-checks-the-run-stops-on-any-failure) passed:
  - the manifests, source hashes and native library agree;
  - the timing rows' guard, class and compact statistics equal the validation rows;
  - per-read times, medians, totals and wins replay `summary.json`;
  - the per-panel predictor medians reproduce every paper key `b100.*.records` and `b100.*.sliced`;
  - the fallback count is 59, all on KIR;
  - the old pooled −0.78 reproduces from the original files, with all 120 old record counts equal to the new ones;
  - all {len(mine)} timing files are byte-identical to those `summary.json` was built from;
  - the archive audit passed afterwards ([archive_audit.json](archive_audit.json): {AUDIT['archive_files_checked']} archive files, {AUDIT['published_checksums_checked']} published checksums, no working-code change).
- **Methods.** Spearman ρ uses average ranks. Intervals are percentile bootstraps over reads ({R['resamples']:,} resamples), stratified by locus when pooled. p-values come from {R['permutations']:,} permutations, Holm-adjusted within each family of {R['primary']['holm_family_size']} tests. The smallest attainable p is 1/{R['permutations'] + 1:,}, so the smallest Holm-adjusted p is {pv(32*min_p)}; every primary association that tracks sits at that floor. No bootstrap resample was degenerate.

## 2. Primary result: guarded speedup per read

Spearman ρ of ln(per-haplotype time / guarded time) against each counter, 100 reads per panel. "Tracks" means ρ < 0 with Holm p < 0.05; "none detected" is not evidence of no association. LO is a leave-out panel.

{primary_table()}

**Pre-specified decisions.**
- **Section 6.5 sentence: supported at all three largest panels.** It needs to say "per read" and to name the new figure.
- **Across loci: the pooled and within-locus associations are both significant.** The pooled ρ ({num(p2['rho'])}) is stronger than the within-locus mean ({num(p3['rho'])}), because the loci also differ in level, which the counters order (Section 6.5).

**Reading.**
- **MHC and KIR.** Records per anchor account for {100*min(t(P, 'guard', k, 'x_rec')['rho']**2 for k in mhckir):.0f}–{100*max(t(P, 'guard', k, 'x_rec')['rho']**2 for k in mhckir):.0f}% of the rank variance (ρ²) at every panel size, from H=60, where every read is slower than per-haplotype, to H≈460.
- **chr20.** chr20 is the locus with the fewest records. Its reads cluster at about one record per anchor, and the association is weaker.
  - At H=60 and 120 the guarded speedups also vary little (IQR {num(chr20_iqr60['q25'])}–{num(chr20_iqr60['q75'])} at H=60), and no association is detectable.
  - There the unguarded class and compact speedups do track records (Section 3.1). So the remaining variation in guarded time is mostly something these counters do not measure.
  - One candidate is the guard's own estimate and metering cost, which adds {PAPER['b100.premium.h60.lo']}–{PAPER['b100.premium.h60.hi']}% to class time at H=60 (paper Table S5). This reading is untested.

## 3. Secondary analyses

### 3.1 Compact and class outcomes

Each column is its own Holm family. "(n.d.)" marks no detectable association; every other entry tracks.

{secondary_table()}

- **Compact speedup** tracks both counters on every panel, chr20 at H=60 and 120 included.
- **Unguarded class speedup** tracks everywhere except bit-sliced decisions on chr20 at H=120.
- **Guard dilution.** The guarded outcome is thus the one that loses the association on chr20 at small H. Pooled correlations hardly depend on the outcome ({num(S1['P2/compact/x_rec']['rho'])} compact, {num(S2['P2/class/x_rec']['rho'])} class, {num(p2['rho'])} guarded).

### 3.2 KIR without fallback reads

Removing the reads with a fallback weakens KIR's association only a little. The fallbacks are not what makes KIR's correlation strong.

{s3_table()}

### 3.3 Bridge to the old Figure S1

{bridge_table()}

- **Unit and outcome.** Going from instances to reads and from compact to guarded speedup changes the pooled correlation little: {num(bridge['c_new_timings_all_components']['pooled'])} → {num(bridge['e_read_level_guarded']['pooled'])}.
- **Timing run.** Most of the change from the old −0.78 comes from the timing run. The predictor values are identical, so only the speedups differ.
- **Old within-locus correlations.** The old within-locus ρ was {num(bridge['a_old_benchmark_120_instances']['within']['mhc'])} on MHC, {num(bridge['a_old_benchmark_120_instances']['within']['kir'])} on KIR and {num(bridge['a_old_benchmark_120_instances']['within']['chr20'])} on chr20. The old figure's "within loci" reading therefore held on MHC and KIR, but not clearly on chr20.

**Post hoc diagnostic (not in the plan).** It compares the same 40 instances per locus between the old benchmark and the new one ([posthoc.py](posthoc.py), [posthoc.json](posthoc.json)):

{posthoc_table()}

**Between-run drift.** Per-instance speedups changed by {num(100*(min(d['min'] for d in drift) - 1), 0)}% to +{100*(max(d['max'] for d in drift) - 1):.0f}% between the two runs. On every locus per-haplotype times changed more than compact times (medians {num(100*(min(ph_drift) - 1), 0)}% to +{100*(max(ph_drift) - 1):.0f}%). Within-run split-half reliability (Section 3.4) is near 1, so it overstates how reproducible a single instance's speedup is from one day to the next. The predictors are exact counts. If this timing noise is unrelated to the counters, it can only weaken the correlations, not create them.

### 3.4 Reliability and spread

The split-half correlation is between per-read log speedups from repetitions 1, 3 and 5 and from repetitions 2 and 4 of the same run.

{reliability_table()}

## 4. Exploratory analyses (largest panels; point estimates and 95% bootstrap intervals, no tests)

{exploratory_table()}

- **Size and states per anchor.** Records per anchor rise with read length and fall with states per anchor. Neither explains the speedup association: the partial correlations stay strong.
- **States per anchor across loci.** The pooled ρ of speedup with states per anchor, {num(E['E1']['E1/pooled/y_guard~s']['rho'])}, mostly reflects differences between loci. Within loci it ranges from {num(E['E1']['E1/chr20/y_guard~s']['rho'])} on chr20 to {num(E['E1']['E1/kir/y_guard~s']['rho'])} on KIR.
- **Work ratio.** The guard's work ratio, which counts all charged class work, is a better per-read predictor than either counter. It is close to a model of the time ratio itself, so this mainly shows that the accounting units are a faithful cost model. Two caveats apply: it is computed from the guarded call of the same validation run, and it is undefined for fallback reads.

## 5. The 59 KIR fallbacks

### 5.1 Which instances fall back

{separation_table()}

- **Records and bit-sliced decisions.** Fallback instances have {num(min(rec_ratio), 1)}–{num(max(rec_ratio), 1)} times the median records per anchor of the other instances, and many more bit-sliced decisions.
- **States per anchor.** They also have fewer states per anchor (AUROC below 0.5), so the per-haplotype engine is relatively cheap on them.
- **Class time.** The guard catches the instances on which class is slowest relative to per-haplotype (AUROC {num(min(aur_c.values()), 2)}–{num(max(aur_c.values()), 3)}).
- **Small H.** At H=60 and 120, most non-fallback instances are also slower with class than per-haplotype ({small_slow['kir-hprc60']['other_class_slower']} of {small_slow['kir-hprc60']['other_class_slower'] + small_slow['kir-hprc60']['other_class_faster']} and {small_slow['kir-hprc120']['other_class_slower']} of {small_slow['kir-hprc120']['other_class_slower'] + small_slow['kir-hprc120']['other_class_faster']}). The guard's budget is in work units and admits them; this is the small-H loss of paper Section 6.5.
- **Unneeded fallbacks.** On 453 LO and 459, {F['F2']['kir-hprc453_leaveout']['two_by_two']['fallback_class_faster']} and {F['F2']['kir-hprc459']['two_by_two']['fallback_class_faster']} fallbacks were instances on which unguarded class would have been faster, by at most {100*(1 - min(unneeded)):.0f}%.

### 5.2 What a fallback costs

The excess is (guarded − per-haplotype) / per-haplotype time, from the five timed repetitions. The class-sweep column is (class kernel − class export) / per-haplotype kernel from the single validation run. It is what a fallback after the whole sweep would add, apart from the estimate pass.

{fallback_panel_table()}

- **Admitted work.** Every refused attempt had admitted {num(fs['admitted_over_budget']['min'], 3)}–{num(fs['admitted_over_budget']['max'], 3)} of its budget: {fs['stage']['output']} were refused at the output charge, {fs['stage']['sweep']} in the sweep.
- **Excess against the class sweep.** The measured excess is a median {num(mean_excess_ratio['median'])}× the class-sweep estimate, which fits a full abandoned sweep plus the guard's own overhead.
- **Exceptions.** At H=459, {len(full_sweep)} of 12 instances cost at least a full sweep. The two exceptions are the two instances with more than 7 records per anchor (`{early[0]['id']}`, `{early[1]['id']}`); the budget trips earlier there, after roughly a third of the sweep. Over all 59, the excess is below 0.8× the class-sweep estimate for {len(low)}: these two instances on all five panels (more than 4.5 records per anchor each time), and `kir_r0016:-` at 453 LO ({num(ratio(low_other[0]))}×). Every other fallback costs at least {num(others_min)}× the estimate.
- **Without the guard.** {fs['class_slower_than_ph']} of the 59 would also have been slower than per-haplotype with unguarded class.

### 5.3 Effect on KIR totals

Totals are summed per-haplotype time over summed method time over all 100 reads. The last two columns replace each fallback instance's guarded time with its own per-haplotype or class time in the same repetition. These are reconstructions from measured times, not measurements of another guard.

{impact_table()}

- **Effect.** Late fallback costs KIR's largest panel about {num(i459['total_fallbacks_at_per_haplotype_time'] - i459['total_measured'])} in total speedup ({num(i459['total_measured'])}× against {num(i459['total_fallbacks_at_per_haplotype_time'])}×).
- **Comparison with unguarded class.** On this sample, unguarded class would have done better in total ({num(i459['total_fallbacks_at_class_time'])}×) than the guard, but it has no worst-case protection: it is up to {PAPER['star.maxslowdown']}× slower on the adversarial star panel.
- **Fallback reads.** The 12 fallback reads take {100*i459['fallback_share_of_guarded_time']:.0f}% of the guarded time at H=459.

### 5.4 Recurrence

- **Distinct instances.** The 59 fallbacks are {rec['distinct_components']} distinct instances from {rec['distinct_reads']} reads; no read falls back in both orientations.
- **Recurrence.** {recurrence_text()}.
- **Leave-out panel.** The 453 LO and 459 panels have the same 12 fallback instances.

The full list of the 59 is in the appendix.

## 6. Proposed Figure S1 replacement and text changes

Applied with the user's approval on 2026-10-18; see 6.6. The text below is what was applied, for the clipped figure of 6.1.

### 6.1 Figure

![Proposed Figure S1](preview/fig-s1-proposal.png)

- **Preview source.** The preview ([PDF](preview/fig-s1-proposal.pdf), [TeX](preview/fig-s1-body.tex)) uses the paper's plot styles and the per-read data in `preview/b100_scatter_*.dat`.
- **Axes.** As pre-specified, it plots guarded speedup on a log axis against records per anchor on a linear axis, at the three largest panels, with KIR fallback reads as open markers.
- **Clipped reads.** The x axis ends at {CLIP} so that the 300 reads do not crowd into a fifth of the width. Two KIR fallback reads lie beyond it, at {num(clipped[0]['records'], 2)} and {num(clipped[1]['records'], 2)} records per anchor. They are drawn as triangles at the edge and named in the caption.
- **Full-range alternative.** It is the [full-range PNG](preview/fig-s1-fullrange.png), page 2 of the [preview PDF](preview/fig-s1-proposal.pdf) ([TeX](preview/fig-s1-body-fullrange.tex)); the two outlying reads leave the other 298 in a fifth of the width. Choosing the clipped axis is a presentation choice, not a change to the plan's content.

### 6.2 Caption (proposed, `supplement/experiments-extra.tex`)

```latex
\\caption{{Per-read guarded speedup over the \\perhap engine, both with expanded output, at the largest panels, against compact records per anchor; each read combines its two instances. There are \\data{{mech.n}} reads. Open squares are LRC/KIR reads on which the guard fell back; the two drawn as triangles at the right edge have \\data{{mech.clip.a}} and \\data{{mech.clip.b}} records per anchor. Spearman correlations are \\data{{mech.kir.rec}} on LRC/KIR, \\data{{mech.mhc.rec}} on the MHC, \\data{{mech.chr20.rec}} on chr20 and \\data{{mech.pooled.rec}} pooled. The analysis was specified before it was run.}}\\label{{fig:supp-scatter}}
```

### 6.3 Section 6.5 (main text)

Replace the last two sentences of the second paragraph:

> These counters order the loci as the speedups do. In our original 40-instance benchmark they were computed before any timing and predicted this order, and within loci the speedup fell as records per anchor grew (Supplementary Figure~S1).

with:

```latex
These counters order the loci as the speedups do; in our original 40-instance benchmark they were computed before any timing and predicted this order. They also track the speedup read by read within each locus. In an analysis specified before it was run, the guarded speedup at about 460 haplotypes has Spearman correlation \\data{{mech.kir.rec}} with records per anchor on LRC/KIR, \\data{{mech.mhc.rec}} on the MHC and \\data{{mech.chr20.rec}} on chr20 (Supplementary Figure~S1). The guard falls back on the LRC/KIR instances with the most records, and usually only after most of the class work, so these instances take \\data{{mech.fb.h459.lo}}--\\data{{mech.fb.h459.hi}} times the \\perhap time at 459 haplotypes. Had each cost only the \\perhap time, the LRC/KIR total would be \\xdata{{mech.kir.total.ph}} instead of \\xdata{{b100.kir.h459.total.guard}}.
```

This adds about three lines; the body ended at about mid-page 10 on 2026-10-18. The fallback sentences could move to Section 6.6 instead; they are optional.

### 6.4 Supplement text

In `supplement/experiments-extra.tex`, delete "\\Cref{{fig:supp-scatter}} relates class speedups to record counts in that benchmark." from the Table S6 paragraph, and add after the first paragraph:

```latex
\\Cref{{fig:supp-scatter}} relates each read's guarded speedup to its records per anchor. The analysis was specified before it was run: Spearman correlations of the per-read speedup with records and with bit-sliced decisions per anchor on every panel, with bootstrap intervals over reads and Holm-adjusted permutation tests. Records per anchor track the speedup on the MHC and LRC/KIR at every panel size, with correlations from \\data{{mech.mhckir.rec.weak}} to \\data{{mech.mhckir.rec.strong}}, and on chr20 from 240 haplotypes (\\data{{mech.chr20big.rec.weak}} to \\data{{mech.chr20big.rec.strong}}). At 60 and 120 haplotypes on chr20, where most reads have about one record per anchor, no association is detectable. With the number of anchors held fixed, the correlations at the largest panels are \\data{{mech.kir.rec.partialN}} on LRC/KIR, \\data{{mech.mhc.rec.partialN}} on the MHC and \\data{{mech.chr20.rec.partialN}} on chr20. Bit-sliced decisions track the speedup on the MHC and LRC/KIR but not on chr20. The \\data{{b100.fallbacks}} fallbacks come from \\data{{mech.fb.distinct}} instances, \\data{{mech.fb.allpanels}} of which fall back on all five LRC/KIR panels; records per anchor separate them from the other instances with an area under the ROC curve of \\data{{mech.fb.auroc.lo}}--\\data{{mech.fb.auroc.hi}}. On \\data{{mech.fb.output}} of them the guard refused only the final output charge, after the whole sweep.
```

The paper README layout row for `fig-scatter.tex` and the supplement row also need updating. Table S6 stays the original 40-instance benchmark. The old keys `scatter.spearman` and `scatter.n` can stay unused.

### 6.5 Generator keys (proposed `make_paper_data.py` additions)

The generator would recompute each ρ and fallback statistic from the raw benchmark rows, with its own `spearman()`, and assert agreement with `results.json`. Values with the generator's rounding:

| Key | Value | Source |
| --- | ---: | --- |
| `mech.n` | {P['P2/guard/x_rec']['n']} | reads at the three largest panels |
| `mech.kir.rec`, `mech.mhc.rec`, `mech.chr20.rec` | {disp(big('kir')['rho'])}, {disp(big('mhc')['rho'])}, {disp(big('chr20')['rho'])} | P1, largest panels |
| `mech.pooled.rec` | {disp(p2['rho'])} | P2 |
| `mech.mhckir.rec.weak`, `.strong` | {disp(max(t(P, 'guard', k, 'x_rec')['rho'] for k in mhckir))}, {disp(min(t(P, 'guard', k, 'x_rec')['rho'] for k in mhckir))} | P1, MHC and KIR panels |
| `mech.chr20big.rec.weak`, `.strong` | {disp(max(t(P, 'guard', k, 'x_rec')['rho'] for k in chr20_big))}, {disp(min(t(P, 'guard', k, 'x_rec')['rho'] for k in chr20_big))} | P1, chr20 at H ≥ 240 |
| `mech.<locus>.rec.partialN` | {disp(partial_rec_N['kir']['rho'])}, {disp(partial_rec_N['mhc']['rho'])}, {disp(partial_rec_N['chr20']['rho'])} | E2 (KIR, MHC, chr20) |
| `mech.clip.a`, `mech.clip.b` | {half_up(clipped[0]['records'], 1)}, {half_up(clipped[1]['records'], 1)} | KIR reads beyond the x axis |
| `mech.fb.h459.lo`, `.hi` | {half_up(fb_ratio['kir-hprc459'][0])}, {half_up(fb_ratio['kir-hprc459'][1])} | F1 at KIR 459 |
| `mech.kir.total.ph` | {half_up(i459['total_fallbacks_at_per_haplotype_time'])} | F3 (iii) at KIR 459 |
| `mech.fb.distinct`, `mech.fb.allpanels`, `mech.fb.output` | {rec['distinct_components']}, {rec['components_by_number_of_panels']['5']}, {fs['stage']['output']} | F4, F1 |
| `mech.fb.auroc.lo`, `.hi` | {half_up(min(aur.values()))}, {half_up(max(aur.values()))} | F2, records per anchor |

The values are `$-$` followed by the magnitude, so the paper prints a true minus sign; the generator forbids backslashes in values.

### 6.6 Applied to the paper

- **Generator.** `mechanism()` in `paper/scripts/make_paper_data.py` writes these {len(EDIT['new_keys'])} keys ({EDIT['generator_numbers_before']} → {EDIT['generator_numbers']} numbers) and `b100_scatter_<locus>.dat` and `b100_scatter_clipped.dat`.
  - It recomputes every value from the raw rows that `bench100()` already replays against `summary.json`.
  - It checks each value against `results.json`, and checks the plan's hash against `PLAN.sha256`.
  - It asserts the Holm-tested verdicts that the text relies on, and that the figure's `xmax` equals its clip value.
  - All {EDIT['generator_numbers_before']} earlier keys kept identical values.
- **Text.** The figure, caption, Section 6.5 and supplement paragraph are as above, and the old Figure S1 sentence was removed from the Table S6 paragraph.
  - **Beyond §6.** One provenance sentence and the analysis command were added to the reproducibility appendix.
- **Checks.** `make -C paper check` passed. The body ends on page {EDIT['body_end_page']} of 10, with no undefined references and no overfull boxes. The {EDIT['todo_markers']} TODO markers are the existing administrative ones, and the warning set is unchanged from the previous build. The archive audit passed after the edit.
- **Evidence.** It is in [paper_edit/](paper_edit/):
  - [pre-edit snapshot](paper_edit/before/) and its hashes ([before.sha256](paper_edit/before.sha256));
  - [source diff](paper_edit/manuscript.diff);
  - [build log](paper_edit/build.log);
  - [renders](paper_edit/render/) of main pages 8–10 and supplement pages 9, 14 and 16;
  - [verification.json](paper_edit/verification.json).

## 7. Deviations, additions and limits

- **Deviations from PLAN.md.** None in the pre-specified statistics.
- **Additions beyond the plan, all descriptive.**
  - The share of guarded time spent on fallback reads (5.3).
  - The excess-to-prediction ratio (5.2).
  - The post hoc between-run diagnostic (3.3).
  - The clipped x axis of the figure (presentation).
  - F3 reports exp(median log speedup) as the median speedup of fallback and other reads.
- **Associations only.** These are associations on one benchmark run, one host and these 100 reads per locus. The counters are deterministic; the times are not. Between-run drift (3.3) is larger than the within-run split-half noise suggests.
- **Reads and graphs.** KIR and chr20 reads are simulated. All graphs are SV-level. Panels of one locus share reads, so the 14 panels are not independent replicates.
- **Monte Carlo floor.** With {R['permutations']:,} permutations, "p = {pv(32*min_p)}" after Holm means "at the floor", not an exact value.
- **Exploratory results generate hypotheses.** This includes the work-ratio result.
- **Scope of the fallback account.** It explains KIR's smaller aggregate gain only in part: without fallback reads KIR's total is {num(i459['total_without_fallback_reads'])}×, still below MHC ({totals_other['mhc']:.2f}×) and chr20 ({totals_other['chr20']:.2f}×).

## 8. Reproduce

```sh
python3.13 -B experiments/mechanism_100reads/analysis.py --checks-only
python3.13 -B experiments/mechanism_100reads/analysis.py --out <new path>/results.json
python3.13 -B experiments/mechanism_100reads/posthoc.py --out <new path>/posthoc.json
python3.13 -B experiments/mechanism_100reads/write_report.py
tectonic -X compile --outdir experiments/mechanism_100reads/preview experiments/mechanism_100reads/preview/fig-s1-proposal.tex
```

`analysis.py` refuses to run if `PLAN.md` no longer matches its frozen hash. The full run takes about {R['seconds']:.0f} s with the standard library only. `write_report.py` reads `results.json` and `posthoc.json` from this folder. The first command writes nothing.

| File | Content |
| --- | --- |
| [PLAN.md](PLAN.md), [PLAN.sha256](PLAN.sha256) | Pre-specified plan and its frozen hash |
| [analysis.py](analysis.py), [analysis.log](analysis.log) | The analysis and its single full run |
| [results.json](results.json) | Every statistic, check and per-read figure row, with input hashes |
| [posthoc.py](posthoc.py), [posthoc.json](posthoc.json) | Post hoc between-run diagnostic |
| [write_report.py](write_report.py) | Generates this report and the preview TeX |
| [preview/](preview/) | Proposed Figure S1: data, TeX, PDF and PNG |
| [archive_audit.json](archive_audit.json) | Archive preservation audit after the analysis |
| [paper_edit/](paper_edit/) | Approved paper edit: snapshot, diff, build log, renders, verification |

## Appendix: the 59 fallback instances

"Refused at" is "output" when the refused charge equals the output charge `N + W·records + K`. The per-haplotype time is the median of the five timed repetitions.

{fallback_rows_table()}
"""
(HERE/'REPORT.md').write_text(report)
print('wrote REPORT.md and preview bodies; results.json', hashlib.sha256((HERE/'results.json').read_bytes()).hexdigest()[:12])
