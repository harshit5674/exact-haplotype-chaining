"""Pre-specified association analysis of the saved 100-read benchmark; see PLAN.md.

Reads saved timing and validation rows only. Runs no timing and changes no input.
Standard library only, assertions enabled. A full run refuses to overwrite its output.

  python3.13 -B experiments/mechanism_100reads/analysis.py --checks-only
  python3.13 -B experiments/mechanism_100reads/analysis.py --out experiments/mechanism_100reads/results.json
"""
import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
from decimal import ROUND_HALF_UP, Decimal
import hashlib
import json
import math
import operator
from pathlib import Path
import platform
import random
import re
import statistics as st
import sys
import time

if not __debug__:
    sys.exit('analysis.py: run with assertions enabled (no -O or PYTHONOPTIMIZE)')

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EXP = ROOT/'experiments'
BENCH = EXP/'benchmark_100reads'
NUMBERS = ROOT/'paper/generated/numbers.tex'
PANELS = (('mhc', 'hprc60', 60), ('mhc', 'hprc120', 120), ('mhc', 'hprc240', 240), ('mhc', 'hprc463', 463),
          ('kir', 'hprc60', 60), ('kir', 'hprc120', 120), ('kir', 'hprc240', 240),
          ('kir', 'hprc453_leaveout', 453), ('kir', 'hprc459', 459),
          ('chr20', 'hprc60', 60), ('chr20', 'hprc120', 120), ('chr20', 'hprc240', 240),
          ('chr20', 'hprc457_leaveout', 457), ('chr20', 'hprc463', 463))
LOCI = ('mhc', 'kir', 'chr20')
LARGEST = {'mhc': 'mhc-hprc463', 'kir': 'kir-hprc459', 'chr20': 'chr20-hprc463'}
MODES = ('per_haplotype', 'auto', 'compact', 'class')
OUTCOMES = {'guard': 'auto', 'compact': 'compact', 'class': 'class'}  # outcome -> mode
PREDICTORS = ('x_rec', 'x_sl')
EXPECTED_FALLBACKS = {'kir-hprc60': 20, 'kir-hprc120': 8, 'kir-hprc240': 7,
                      'kir-hprc453_leaveout': 12, 'kir-hprc459': 12}
OLD = {'mhc': (EXP/'hprc_mhc_panel/benchmark_v2.jsonl', 463),
       'kir': (EXP/'second_locus/benchmark_kir.jsonl', 459),
       'chr20': (EXP/'second_locus/benchmark_chr20.jsonl', 463)}
OLD_POOLED, OLD_GAMMA = '-0.78', 10000
RESAMPLES, PERMUTATIONS = 2000, 10000
HALVES = ((0, 2, 4), (1, 3))  # repetitions {1,3,5} and {2,4}

inputs = {}


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def load_json(p):
    inputs[str(p.relative_to(ROOT))] = sha(p)
    return json.loads(p.read_text())


def load_jsonl(p):
    inputs[str(p.relative_to(ROOT))] = sha(p)
    return [json.loads(line) for line in p.read_text().splitlines() if line.strip()]


def close(a, b):
    return math.isclose(a, b, rel_tol=1e-12, abs_tol=0.0)


def f2(x):
    # The paper generator's rounding: half up from the decimal repr.
    return str(Decimal(repr(x)).quantize(Decimal('0.01'), rounding=ROUND_HALF_UP))


# ------------------------------------------------------------------------------------ statistics
def ranks(xs):
    """Average ranks from 0, ties averaged (as spearman() in paper/scripts/make_paper_data.py)."""
    order = sorted(range(len(xs)), key=xs.__getitem__)
    r = [0.0]*len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j)/2
        i = j + 1
    return r


def centered(r):
    m = math.fsum(r)/len(r)
    d = [v - m for v in r]
    return d, math.fsum(v*v for v in d)


def pearson(a, b):
    da, va = centered(a)
    db, vb = centered(b)
    if va == 0 or vb == 0:
        return None  # a constant variable: undefined
    return math.fsum(map(operator.mul, da, db))/math.sqrt(va*vb)


def spearman(x, y):
    return pearson(ranks(x), ranks(y))


def partial(r_yx, r_yz, r_xz):
    if None in (r_yx, r_yz, r_xz) or abs(r_yz) >= 1 or abs(r_xz) >= 1:
        return None
    return (r_yx - r_yz*r_xz)/math.sqrt((1 - r_yz**2)*(1 - r_xz**2))


def quantile(values, fraction):
    # Nearest rank on the sorted list, as experiments/benchmark_100reads/summarize.py.
    ordered = sorted(values)
    return ordered[round(fraction*(len(ordered) - 1))]


def spread(values):
    return dict(min=min(values), q25=quantile(values, .25), median=st.median(values),
                q75=quantile(values, .75), max=max(values))


def rng_for(label):
    return random.Random('mechanism_100reads/' + label)


def bootstrap(groups, stat, label):
    """Percentile interval; rows are resampled with replacement within each group.

    groups: list of tuples of equal-length columns, one tuple per stratum (locus)."""
    rng = rng_for('boot/' + label)
    values, discarded = [], 0
    for _ in range(RESAMPLES):
        sample = []
        for cols in groups:
            n = len(cols[0])
            idx = [rng.randrange(n) for _ in range(n)]
            sample.append([[c[i] for i in idx] for c in cols])
        v = stat(sample)
        if v is None:
            discarded += 1
        else:
            values.append(v)
    return [quantile(values, .025), quantile(values, .975)], discarded


def permutation(groups, label, pooled):
    """Two-sided Monte Carlo permutation test of a rank correlation.

    pooled: one correlation over all rows, the outcome permuted across all rows. Otherwise the
    mean of the per-group correlations, the outcome permuted within each group.
    groups: list of (predictor, outcome) column pairs."""
    rng = rng_for('perm/' + label)
    parts = [([v for g in groups for v in g[0]], [v for g in groups for v in g[1]])] if pooled else groups
    prepared = []
    for x, y in parts:
        dx, vx = centered(ranks(x))
        dy, vy = centered(ranks(y))
        assert vx > 0 and vy > 0
        prepared.append((dx, dy, math.sqrt(vx*vy)))

    def stat(ys):
        return st.fmean(sum(map(operator.mul, dx, dy))/den for (dx, _, den), dy in zip(prepared, ys))
    observed = stat([dy for _, dy, _ in prepared])
    work = [list(dy) for _, dy, _ in prepared]
    hits = 0
    for _ in range(PERMUTATIONS):
        for w in work:
            rng.shuffle(w)
        if abs(stat(work)) >= abs(observed) - 1e-12:
            hits += 1
    return observed, (1 + hits)/(1 + PERMUTATIONS)


def holm(ps):
    """Holm step-down adjusted p-values for a dict label -> p."""
    items = sorted(ps.items(), key=lambda kv: kv[1])
    m, running, adjusted = len(items), 0.0, {}
    for k, (label, p) in enumerate(items):
        running = max(running, min(1.0, (m - k)*p))
        adjusted[label] = running
    return adjusted


def auroc(a, b):
    """P(random a > random b) + P(tie)/2: Mann-Whitney U / (n_a n_b)."""
    u = sum((x > y) + 0.5*(x == y) for x in a for y in b)
    return u/(len(a)*len(b))


def magnitude(r):
    a = abs(r)
    return 'negligible' if a < .10 else 'weak' if a < .30 else 'moderate' if a < .50 else 'strong'


def pooled_rho(sample):
    return spearman([v for g in sample for v in g[0]], [v for g in sample for v in g[1]])


def within_rho(sample):
    per = [spearman(g[0], g[1]) for g in sample]
    return None if None in per else st.fmean(per)


def single_rho(sample):
    return spearman(sample[0][0], sample[0][1])


def corr(groups, label, kind, test=True):
    """kind: 'single' (one group), 'pooled' or 'within' (mean of per-group correlations)."""
    stat = {'single': single_rho, 'pooled': pooled_rho, 'within': within_rho}[kind]
    rho = stat(groups)
    out = dict(n=sum(len(g[0]) for g in groups), rho=rho)
    if kind == 'within':
        out['per_group'] = [spearman(g[0], g[1]) for g in groups]
    if rho is None:
        return dict(out, ci95=None, discarded=None, p=None)
    out['ci95'], out['discarded'] = bootstrap(groups, stat, label)
    if test:
        observed, out['p'] = permutation(groups, label, pooled=kind != 'within')
        assert abs(observed - rho) < 1e-9, (label, observed, rho)
    return out


def partial_corr(groups, label):
    """Partial Spearman of y with x given z; groups: list of (y, x, z); pooled over groups."""
    def stat(sample):
        y, x, z = ([v for g in sample for v in g[k]] for k in range(3))
        return partial(spearman(y, x), spearman(y, z), spearman(x, z))
    rho = stat(groups)
    ci, discarded = bootstrap(groups, stat, label) if rho is not None else (None, None)
    return dict(n=sum(len(g[0]) for g in groups), rho=rho, ci95=ci, discarded=discarded)


def selftest():
    assert ranks([3, 1, 2]) == [2.0, 0.0, 1.0] and ranks([5, 1, 5, 2]) == [2.5, 0.0, 2.5, 1.0]
    rng = random.Random('mechanism_100reads/selftest')
    for n in (5, 17, 100, 300):
        x, y = rng.sample(range(10*n), n), rng.sample(range(10*n), n)
        rx, ry = ranks(x), ranks(y)
        closed = 1 - 6*sum((a - b)**2 for a, b in zip(rx, ry))/(n*(n*n - 1))
        assert math.isclose(spearman(x, y), closed, abs_tol=1e-12)
    assert math.isclose(spearman([1, 2, 3, 4], [2, 4, 6, 9]), 1) and math.isclose(spearman([1, 2, 3, 4], [9, 3, 2, 1]), -1)
    assert spearman([1, 1, 1], [1, 2, 3]) is None
    # Ties by hand: ranks [0,1.5,1.5,3] and [0,2,1,3]; covariance 4.5, variances 4.5 and 5.
    assert math.isclose(spearman([1, 2, 2, 3], [1, 3, 2, 4]), 4.5/math.sqrt(22.5), abs_tol=1e-12)
    adjusted = holm({'a': .01, 'b': .04, 'c': .03, 'd': .005})
    assert all(math.isclose(adjusted[k], v) for k, v in dict(d=.02, a=.03, c=.06, b=.06).items())
    # The partial formula equals the correlation of rank residuals.
    x = [rng.random() for _ in range(200)]
    z = [a + rng.random() for a in x]
    y = [b + rng.random() for b in z]

    def resid(a, b):
        ma, mb = st.fmean(a), st.fmean(b)
        beta = sum((u - ma)*(v - mb) for u, v in zip(a, b))/sum((v - mb)**2 for v in b)
        return [u - ma - beta*(v - mb) for u, v in zip(a, b)]
    rx, ry, rz = ranks(x), ranks(y), ranks(z)
    assert math.isclose(partial(spearman(y, x), spearman(y, z), spearman(x, z)),
                        pearson(resid(ry, rz), resid(rx, rz)), abs_tol=1e-12)
    groups = [(x[:50], y[:50])]
    assert bootstrap(groups, single_rho, 'selftest') == bootstrap(groups, single_rho, 'selftest')
    assert bootstrap(groups, single_rho, 'selftest') != bootstrap(groups, single_rho, 'selftest2')
    observed, p = permutation([(list(range(12)), list(range(12)))], 'selftest', pooled=True)
    assert math.isclose(observed, 1) and p == 1/(1 + PERMUTATIONS)
    assert auroc([3, 4], [1, 2]) == 1 and auroc([1], [1]) == .5 and auroc([1, 2], [2, 3]) == .125
    return 'passed'


# ------------------------------------------------------------------------------------ data
def paper_numbers():
    inputs[str(NUMBERS.relative_to(ROOT))] = sha(NUMBERS)
    found = re.findall(r'\\dataset\{([^}]*)\}\{([^}]*)\}', NUMBERS.read_text())
    return dict(found)


def component(c, v):
    g = c['guard']
    return dict(id=c['id'], N=c['N'], K=c['K'], W=g['W'], records=c['compact_records'],
                sliced=c['class_stats']['pairs_sliced'], pairs=c['class_stats']['pairs'],
                fallback=g['fallback'], guard=g, reps=c['kernel_seconds'],
                t={m: st.median(c['kernel_seconds'][m]) for m in MODES},
                v_class_kernel=v['stats']['class']['kernel_seconds'],
                v_class_export=v['stats']['class']['export_seconds'],
                v_ph_kernel=v['stats']['per_haplotype']['kernel_seconds'])


def load_panel(locus, name, H, summary_panel, reps, numbers, checks):
    key = f'{locus}-{name}'
    trows = load_jsonl(BENCH/'timings'/f'{key}.jsonl')
    vrows = {r['id']: r for r in load_jsonl(BENCH/'validation'/f'{key}.jsonl')}
    tmeta, vmeta = load_json(BENCH/'timings'/f'{key}.json'), load_json(BENCH/'validation'/f'{key}.json')
    assert tmeta['status'] == vmeta['status'] == 'passed'
    assert tmeta['H'] == vmeta['H'] == summary_panel['H'] == H
    assert tmeta['input_hashes'] == vmeta['input_hashes']
    assert len(trows) == len(vrows) == 200 and {r['id'] for r in trows} == set(vrows)
    groups = defaultdict(list)
    for r in trows:
        v = vrows[r['id']]
        assert (r['N'], r['K']) == (v['N'], v['K'])
        assert r['guard'] == v['stats']['auto'] and r['class_stats'] == v['stats']['class']
        assert r['compact_records'] == v['stats']['compact']['compact_records']
        assert r['guard']['auto_fallback_available'] and r['guard']['auto_work'] <= r['guard']['auto_budget']
        assert r['guard']['W'] == r['class_stats']['W'] and r['guard']['H'] == H
        for m in MODES:
            assert len(r['kernel_seconds'][m]) == reps and min(r['kernel_seconds'][m]) > 0
        groups[r['read']].append(r)
    assert len(groups) == 100 and all(len(c) == 2 for c in groups.values())
    srows = summary_panel['read_rows']
    assert [s['read'] for s in srows] == list(groups)
    reads = []
    for s in srows:
        comps = groups[s['read']]
        assert len({c['prior_primary'] for c in comps}) == 1
        # Sum both orientations within each repetition, then take the median (summarize.py).
        per_rep = {m: [sum(c['kernel_seconds'][m][j] for c in comps) for j in range(reps)] for m in MODES}
        t = {m: st.median(per_rep[m]) for m in MODES}
        assert all(close(t[m], s['kernel_seconds'][m]) for m in MODES), (key, s['read'])
        N, K = sum(c['N'] for c in comps), sum(c['K'] for c in comps)
        assert (N, K) == (s['N'], s['K']) and N > 0
        fallback = any(c['guard']['fallback'] for c in comps)
        assert fallback == s['fallback'] and s['prior_primary'] == comps[0]['prior_primary']
        records = sum(c['compact_records'] for c in comps)
        sliced = sum(c['class_stats']['pairs_sliced'] for c in comps)
        work = sum(c['guard']['auto_work'] for c in comps)
        estimate = sum(c['guard']['auto_estimate'] for c in comps)
        reads.append(dict(read=s['read'], prior=s['prior_primary'], fallback=fallback, N=N, K=K,
                          records=records, sliced=sliced, x_rec=records/N, x_sl=sliced/N, s=K/N,
                          w=None if fallback else work/estimate, t=t, per_rep=per_rep,
                          y={o: math.log(t['per_haplotype']/t[m]) for o, m in OUTCOMES.items()},
                          comps=[component(c, vrows[c['id']]) for c in comps]))
    # Panel summaries replay summary.json.
    for o, m in OUTCOMES.items():
        ref = summary_panel['kernel'][m]
        ratios = [r['t']['per_haplotype']/r['t'][m] for r in reads]
        total = sum(r['t']['per_haplotype'] for r in reads)/sum(r['t'][m] for r in reads)
        assert close(st.median(ratios), ref['median_speedup']) and close(total, ref['total_speedup']), (key, m)
        assert sum(x > 1 for x in ratios) == ref['wins']
    fallbacks = sum(c['fallback'] for r in reads for c in r['comps'])
    assert fallbacks == summary_panel['fallback_components'] == EXPECTED_FALLBACKS.get(key, 0), key
    # Predictor definitions reproduce the paper's Table S2 keys.
    rec_med, sl_med = st.median(r['x_rec'] for r in reads), st.median(r['x_sl'] for r in reads)
    paper_rec, paper_sl = numbers[f'b100.{locus}.h{H}.records'], numbers[f'b100.{locus}.h{H}.sliced']
    assert f2(rec_med) == paper_rec and f'{sl_med:.3f}' == paper_sl, (key, rec_med, paper_rec, sl_med, paper_sl)
    checks.append(dict(panel=key, reads=100, components=200, fallback_components=fallbacks,
                       fallback_reads=sum(r['fallback'] for r in reads),
                       median_records_per_anchor=rec_med, paper_records=paper_rec,
                       median_sliced_per_anchor=sl_med, paper_sliced=paper_sl,
                       summary_replayed=True))
    return dict(locus=locus, name=name, key=key, H=H, reads=reads)


def load_all():
    checks = []
    tman, vman = load_json(BENCH/'timings'/'manifest.json'), load_json(BENCH/'validation'/'manifest.json')
    assert tman['status'] == vman['status'] == 'passed'
    assert (tman['phase'], vman['phase']) == ('benchmark', 'validate')
    assert tman['source_hashes'] == vman['source_hashes'] and tman['library'] == vman['library']
    assert tman['gamma'] == 10000 and tman['weight_scale'] == 200 and tman['auto_factor'] == 1.0
    assert tman['pair_cache'] == 'cold' and tman['threads'] == 1 and tman['assertions'] and vman['assertions']
    assert {(p['locus'], p['panel']) for p in tman['panels']} == {(a, b) for a, b, _ in PANELS}
    reps = tman['repetitions']
    assert reps == 5
    summary = load_json(BENCH/'summary.json')
    assert summary['status'] == 'passed' and summary['source_unchanged']
    vsum = load_json(BENCH/'validation_summary.json')
    assert vsum['status'] == 'passed' and vsum['components'] == summary['components'] == 2800
    by_panel = {(p['locus'], p['panel']): p for p in summary['panels']}
    numbers = paper_numbers()
    panels = {}
    for locus, name, H in PANELS:
        p = load_panel(locus, name, H, by_panel[(locus, name)], reps, numbers, checks)
        panels[p['key']] = p
    assert sum(EXPECTED_FALLBACKS.values()) == 59
    meta = dict(repetitions=reps, library=tman['library'], gamma=tman['gamma'],
                timing_boundary='native kernel', source_hashes_equal=True)
    return panels, checks, meta, numbers


def old_bridge(panels, checks_only):
    """S4a: the original 40-instance benchmark at the largest panels, as the paper's Figure S1."""
    groups, mismatches, matched = [], 0, 0
    for locus in LOCI:
        path, H = OLD[locus]
        rows = [r for r in load_jsonl(path) if r['H'] == H and r['gamma'] == OLD_GAMMA]
        assert len(rows) == 40
        recs = [r['stats']['class_compact']['output_records']/r['N'] for r in rows]
        speed = [r['kernel_median']['per_haplotype']/r['kernel_median']['class_compact'] for r in rows]
        groups.append((recs, speed))
        new = {c['id']: c for r in panels[LARGEST[locus]]['reads'] for c in r['comps']}
        for r in rows:
            c = new[r['id']]
            matched += 1
            mismatches += (r['N'], r['stats']['class_compact']['output_records']) != (c['N'], c['records'])
    pooled = pooled_rho(groups)
    assert f2(pooled) == OLD_POOLED, pooled
    out = dict(pooled=pooled, n=120, components_compared=matched, record_count_mismatches=mismatches)
    if not checks_only:
        out['within'] = {locus: spearman(*g) for locus, g in zip(LOCI, groups)}
    return out


# ------------------------------------------------------------------------------------ analyses
def columns(reads, x, outcome):
    return ([r[x] for r in reads], [r['y'][outcome] for r in reads])


def family(panels, outcome):
    tests = {}
    for key, p in panels.items():
        for x in PREDICTORS:
            label = f'P1/{outcome}/{key}/{x}'
            tests[label] = corr([columns(p['reads'], x, outcome)], label, 'single')
    for x in PREDICTORS:
        groups = [columns(panels[LARGEST[locus]]['reads'], x, outcome) for locus in LOCI]
        tests[f'P2/{outcome}/{x}'] = corr(groups, f'P2/{outcome}/{x}', 'pooled')
        tests[f'P3/{outcome}/{x}'] = corr(groups, f'P3/{outcome}/{x}', 'within')
    defined = {k: v['p'] for k, v in tests.items() if v['p'] is not None}
    adjusted = holm(defined)
    for k, v in tests.items():
        if v['rho'] is None:
            v.update(p_holm=None, magnitude=None, verdict='undefined')
            continue
        v['p_holm'] = adjusted[k]
        v['magnitude'] = magnitude(v['rho'])
        v['verdict'] = ('tracks' if v['rho'] < 0 and v['p_holm'] < .05 else
                        'opposite to prediction' if v['rho'] > 0 and v['p_holm'] < .05 else
                        'no detectable association')
    return dict(tests=tests, holm_family_size=len(defined))


def decisions(fam, outcome):
    rec = {locus: fam['tests'][f'P1/{outcome}/{LARGEST[locus]}/x_rec'] for locus in LOCI}
    tracking = [locus for locus in LOCI if rec[locus]['verdict'] == 'tracks']
    sentence = ('supported at all three largest panels' if len(tracking) == 3 else
                'restrict to ' + ', '.join(tracking) if tracking else 'remove or replace')
    across = {}
    for x in PREDICTORS:
        p2, p3 = fam['tests'][f'P2/{outcome}/{x}'], fam['tests'][f'P3/{outcome}/{x}']
        s2, s3 = p2['p_holm'] is not None and p2['p_holm'] < .05, p3['p_holm'] is not None and p3['p_holm'] < .05
        across[x] = ('pooled association reflects the between-locus difference' if s2 and not s3 else
                     'pooled and within-locus associations both significant' if s2 and s3 else
                     'within-locus significant, pooled not' if s3 else 'neither significant')
    return dict(largest_panel_tracking_x_rec=tracking, section_6_5_sentence=sentence, across_loci=across)


def no_fallback_kir(panels):
    out = {}
    for key, p in panels.items():
        if p['locus'] != 'kir':
            continue
        reads = [r for r in p['reads'] if not r['fallback']]
        for x in PREDICTORS:
            label = f'S3/guard/{key}/{x}'
            out[label] = corr([columns(reads, x, 'guard')], label, 'single')
    return out


def bridge(panels, old):
    def comp_groups(selector):
        groups, excluded = [], 0
        for locus in LOCI:
            comps = [c for r in panels[LARGEST[locus]]['reads'] if selector(r) for c in r['comps']]
            excluded += sum(c['N'] == 0 for c in comps)
            comps = [c for c in comps if c['N'] > 0]
            groups.append(([c['records']/c['N'] for c in comps],
                           [c['t']['per_haplotype']/c['t']['compact'] for c in comps]))
        return groups, excluded

    def summarize(groups, excluded):
        return dict(n=sum(len(g[0]) for g in groups), excluded_zero_anchor=excluded, pooled=pooled_rho(groups),
                    within={locus: spearman(*g) for locus, g in zip(LOCI, groups)})
    out = {'a_old_benchmark_120_instances': old}
    out['b_new_timings_earlier_20_reads'] = summarize(*comp_groups(lambda r: r['prior']))
    assert out['b_new_timings_earlier_20_reads']['n'] == 120
    out['c_new_timings_all_components'] = summarize(*comp_groups(lambda r: True))
    for name, outcome in (('d_read_level_compact', 'compact'), ('e_read_level_guarded', 'guard')):
        groups = [columns(panels[LARGEST[locus]]['reads'], 'x_rec', outcome) for locus in LOCI]
        out[name] = summarize(groups, 0)
    return out


def reliability(panels):
    out = {}
    for key, p in panels.items():
        row = dict(H=p['H'])
        for o in ('guard', 'compact'):
            m = OUTCOMES[o]
            halves = [[math.log(st.median(r['per_rep']['per_haplotype'][j] for j in h) /
                                st.median(r['per_rep'][m][j] for j in h)) for r in p['reads']] for h in HALVES]
            row[f'split_half_{o}'] = spearman(*halves)
        for name, vals in (('speedup_guard', [r['t']['per_haplotype']/r['t']['auto'] for r in p['reads']]),
                           ('speedup_compact', [r['t']['per_haplotype']/r['t']['compact'] for r in p['reads']]),
                           ('speedup_class', [r['t']['per_haplotype']/r['t']['class'] for r in p['reads']]),
                           ('x_rec', [r['x_rec'] for r in p['reads']]),
                           ('x_sl', [r['x_sl'] for r in p['reads']])):
            row[name] = spread(vals)
        row['fallback_reads'] = sum(r['fallback'] for r in p['reads'])
        out[key] = row
    return out


def exploratory(panels):
    def var(r, name):
        return r['y']['guard'] if name == 'y_guard' else r['y']['class'] if name == 'y_class' else r[name]
    sets = {locus: panels[LARGEST[locus]]['reads'] for locus in LOCI}
    out = dict(E1={}, E2={}, E3={})
    pairs = (('y_guard', 'N'), ('y_guard', 's'), ('x_rec', 'x_sl'), ('x_rec', 's'), ('x_sl', 's'),
             ('x_rec', 'N'), ('x_sl', 'N'))
    for a, b in pairs:
        for locus in LOCI:
            label = f'E1/{locus}/{a}~{b}'
            out['E1'][label] = corr([([var(r, b) for r in sets[locus]], [var(r, a) for r in sets[locus]])],
                                    label, 'single', test=False)
        label = f'E1/pooled/{a}~{b}'
        out['E1'][label] = corr([([var(r, b) for r in sets[locus]], [var(r, a) for r in sets[locus]]) for locus in LOCI],
                                label, 'pooled', test=False)
    for x in PREDICTORS:
        for z in ('N', 's'):
            for scope in LOCI + ('pooled',):
                label = f'E2/{scope}/y_guard~{x}|{z}'
                loci = LOCI if scope == 'pooled' else (scope,)
                out['E2'][label] = partial_corr([([var(r, 'y_guard') for r in sets[locus]], [var(r, x) for r in sets[locus]],
                                                  [var(r, z) for r in sets[locus]]) for locus in loci], label)
    for y in ('y_guard', 'y_class'):
        for scope in LOCI + ('pooled',):
            label = f'E3/{scope}/{y}~w'
            loci = LOCI if scope == 'pooled' else (scope,)
            groups = [([r['w'] for r in sets[locus] if not r['fallback']],
                       [var(r, y) for r in sets[locus] if not r['fallback']]) for locus in loci]
            out['E3'][label] = corr(groups, label, 'pooled' if scope == 'pooled' else 'single', test=False)
    return out


def fallback_analysis(panels):
    rows, by_panel, impact = [], {}, {}
    for key, p in panels.items():
        if p['locus'] != 'kir':
            continue
        comps = [(r, c) for r in p['reads'] for c in r['comps']]
        assert not any(c['fallback'] for _, c in comps if c['N'] == 0)
        for r, c in comps:
            if not c['fallback']:
                continue
            g, t = c['guard'], c['t']
            output_charge = c['N'] + c['W']*c['records'] + c['K']
            rows.append(dict(panel=key, H=p['H'], id=c['id'], read=r['read'], earlier_20=r['prior'],
                             N=c['N'], K=c['K'], states_per_anchor=c['K']/c['N'], W=c['W'],
                             records_per_anchor=c['records']/c['N'], sliced_per_anchor=c['sliced']/c['N'],
                             pairs_per_anchor=c['pairs']/c['N'], estimate=g['auto_estimate'],
                             budget=g['auto_budget'], admitted=g['auto_work'],
                             admitted_over_budget=g['auto_work']/g['auto_budget'],
                             refused=g['auto_denied_work'], output_charge=output_charge,
                             stage='output' if g['auto_denied_work'] == output_charge else 'sweep',
                             per_haplotype_ms=1000*t['per_haplotype'],
                             guarded_over_ph=t['auto']/t['per_haplotype'],
                             class_over_ph=t['class']/t['per_haplotype'],
                             compact_over_ph=t['compact']/t['per_haplotype'],
                             excess=(t['auto'] - t['per_haplotype'])/t['per_haplotype'],
                             predicted_excess_validation=(c['v_class_kernel'] - c['v_class_export'])/c['v_ph_kernel']))
        nonzero = [(r, c) for r, c in comps if c['N'] > 0]
        fb = [c for _, c in nonzero if c['fallback']]
        other = [c for _, c in nonzero if not c['fallback']]
        variables = {'records_per_anchor': lambda c: c['records']/c['N'], 'sliced_per_anchor': lambda c: c['sliced']/c['N'],
                     'states_per_anchor': lambda c: c['K']/c['N'], 'N': lambda c: c['N'],
                     'class_over_ph': lambda c: c['t']['class']/c['t']['per_haplotype']}
        compare = {name: dict(fallback=spread([f(c) for c in fb]), other=spread([f(c) for c in other]),
                              auroc=auroc([f(c) for c in fb], [f(c) for c in other]))
                   for name, f in variables.items()}
        slow = lambda c: c['t']['class'] > c['t']['per_haplotype']
        compare['two_by_two'] = dict(fallback_class_slower=sum(map(slow, fb)), fallback_class_faster=sum(not slow(c) for c in fb),
                                     other_class_slower=sum(map(slow, other)), other_class_faster=sum(not slow(c) for c in other))
        compare['zero_anchor_components'] = len(comps) - len(nonzero)
        by_panel[key] = compare

        reads = p['reads']
        reps = len(reads[0]['per_rep']['auto'])

        def total(selected, pick):
            num = sum(r['t']['per_haplotype'] for r in selected)
            den = sum(st.median(sum(pick(c, j) for c in r['comps']) for j in range(reps)) for r in selected)
            return num/den
        measured = total(reads, lambda c, j: c['reps']['auto'][j])
        assert close(measured, sum(r['t']['per_haplotype'] for r in reads)/sum(r['t']['auto'] for r in reads))
        fb_reads = [r for r in reads if r['fallback']]
        other_reads = [r for r in reads if not r['fallback']]
        impact[key] = dict(
            fallback_reads=len(fb_reads), fallback_components=len(fb),
            median_speedup_fallback_reads=math.exp(st.median(r['y']['guard'] for r in fb_reads)) if fb_reads else None,
            median_speedup_other_reads=math.exp(st.median(r['y']['guard'] for r in other_reads)),
            fallback_share_of_guarded_time=sum(r['t']['auto'] for r in fb_reads)/sum(r['t']['auto'] for r in reads),
            total_measured=measured,
            total_without_fallback_reads=total(other_reads, lambda c, j: c['reps']['auto'][j]),
            total_fallbacks_at_per_haplotype_time=total(reads, lambda c, j: c['reps']['per_haplotype' if c['fallback'] else 'auto'][j]),
            total_fallbacks_at_class_time=total(reads, lambda c, j: c['reps']['class' if c['fallback'] else 'auto'][j]))
    assert len(rows) == 59
    summary = dict(components=len(rows), stage=dict(Counter(r['stage'] for r in rows)),
                   class_slower_than_ph=sum(r['class_over_ph'] > 1 for r in rows))
    for name in ('guarded_over_ph', 'class_over_ph', 'compact_over_ph', 'excess', 'predicted_excess_validation',
                 'admitted_over_budget', 'records_per_anchor', 'sliced_per_anchor', 'states_per_anchor', 'N'):
        summary[name] = spread([r[name] for r in rows])
    summary['excess_over_predicted'] = spread([r['excess']/r['predicted_excess_validation'] for r in rows])
    ids = defaultdict(list)
    for r in rows:
        ids[r['id']].append(r['panel'])
    read_panels = defaultdict(set)
    for r in rows:
        read_panels[r['read']].add(r['panel'])
    recurrence = dict(distinct_components=len(ids), distinct_reads=len(read_panels),
                      components_by_number_of_panels=dict(sorted(Counter(map(len, ids.values())).items())),
                      reads_by_number_of_panels=dict(sorted(Counter(map(len, read_panels.values())).items())),
                      h453_h459_shared=len({i for i, ps in ids.items() if {'kir-hprc453_leaveout', 'kir-hprc459'} <= set(ps)}))
    return dict(F1=rows, F1_summary=summary, F2=by_panel, F3=impact, F4=recurrence)


def figure_data(panels):
    out = {}
    (HERE/'preview').mkdir(exist_ok=True)
    for locus in LOCI:
        reads = sorted(panels[LARGEST[locus]]['reads'], key=lambda r: (r['x_rec'], r['read']))
        rows = [dict(read=r['read'], records=r['x_rec'], sliced=r['x_sl'], guard=math.exp(r['y']['guard']),
                     compact=math.exp(r['y']['compact']), cls=math.exp(r['y']['class']), fallback=int(r['fallback']))
                for r in reads]
        with open(HERE/'preview'/f'b100_scatter_{locus}.dat', 'w') as fh:
            fh.write('records sliced guard compact class fallback\n')
            for r in rows:
                fh.write(f"{r['records']:.4f} {r['sliced']:.4f} {r['guard']:.4f} {r['compact']:.4f} {r['cls']:.4f} {r['fallback']}\n")
        out[locus] = rows
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path)
    ap.add_argument('--checks-only', action='store_true')
    args = ap.parse_args()
    assert args.checks_only != bool(args.out), 'give exactly one of --out and --checks-only'
    if args.out:
        assert not args.out.exists(), f'{args.out} exists; choose a new path'
    started, t0 = datetime.now(timezone.utc).isoformat(), time.perf_counter()
    plan_hash = sha(HERE/'PLAN.md')
    frozen = (HERE/'PLAN.sha256').read_text().split()[0]
    assert plan_hash == frozen, 'PLAN.md changed after it was frozen; record an amendment first'
    tests_status = selftest()
    panels, checks, meta, _ = load_all()
    old = old_bridge(panels, args.checks_only)
    if args.checks_only:
        for c in checks:
            print(c['panel'], 'fallback components', c['fallback_components'], 'records', f"{c['median_records_per_anchor']:.4f}",
                  c['paper_records'], 'sliced', f"{c['median_sliced_per_anchor']:.4f}", c['paper_sliced'])
        print('self-test', tests_status, '| old pooled rho', f"{old['pooled']:.4f}",
              '| record-count mismatches', old['record_count_mismatches'], 'of', old['components_compared'])
        print('checks passed; no association computed')
        return
    families = {o: family(panels, o) for o in OUTCOMES}
    result = dict(
        status='passed', created_utc=started, command=sys.argv, python=sys.version, platform=platform.platform(),
        analysis_sha256=sha(Path(__file__)), plan_sha256=plan_hash, plan_frozen=(HERE/'PLAN.sha256').read_text().strip(),
        resamples=RESAMPLES, permutations=PERMUTATIONS,
        rng='random.Random("mechanism_100reads/" + boot|perm + "/" + label)', selftest=tests_status,
        benchmark=meta, checks=checks,
        primary=families['guard'], decisions=decisions(families['guard'], 'guard'),
        secondary=dict(S1_compact=families['compact'], S1_decisions=decisions(families['compact'], 'compact'),
                       S2_class=families['class'], S2_decisions=decisions(families['class'], 'class'),
                       S3_kir_without_fallback_reads=no_fallback_kir(panels),
                       S4_bridge=bridge(panels, old), S5_reliability=reliability(panels)),
        exploratory=exploratory(panels), fallback=fallback_analysis(panels), figure=figure_data(panels))
    result['inputs_sha256'] = dict(sorted(inputs.items()))
    result['seconds'] = time.perf_counter() - t0
    args.out.write_text(json.dumps(result, indent=1) + '\n')
    print('wrote', args.out, f"in {result['seconds']:.0f} s")


if __name__ == '__main__':
    main()
