"""Realism check: simulated CHM13 MHC reads against the real CHM13 PacBio MHC reads.

Both read sets are seeded by the same pinned Minichain on the same MHC panel graphs
(experiments/hprc_mhc_panel/data/hprc463 and hprc60): the real reads' instances are those of
the earlier experiment, the simulated ones are data/mhc_sim/<panel>/instances.jsonl (100
reads from the CHM13 MHC walk, simulate_reads.py, the same model as the new loci).
Per read, the major orientation is the instance with more anchors (the true strand is known
for simulated reads and is checked to agree). Compared per read: length, anchors (major,
minor, per kb), states per anchor, the share of major-orientation anchors on CHM13's own path,
anchors per distinct read end position (multi-mapping), and, from counters_mhc*.jsonl, the
class engine's per-anchor counters on the major-orientation instances. Reported as
[p25, median, p75] per read set, the ratio of medians, and the two-sample
Kolmogorov-Smirnov distance D (no p-values are claimed; the reads are not a random
sample of the same distribution: real reads were length-filtered, simulated ones are
lognormal around 18 kb).
Output: realism_check.json.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import locus_common as lc  # noqa: E402

COUNTER_KEYS = ('pairs_sliced', 'sliced_thresholds', 'sliced_minima', 'output_records', 'published_records',
                'remote_candidates', 'prefix_records')


def quartiles(xs):
    s = sorted(xs)
    q = lambda f: s[min(len(s) - 1, int(f * (len(s) - 1) + 0.5))]  # noqa: E731
    return [round(q(0.25), 4), round(q(0.5), 4), round(q(0.75), 4)]


def ks(a, b):
    """Two-sample Kolmogorov-Smirnov distance sup |F_a - F_b|."""
    a, b = sorted(a), sorted(b)
    i = j = 0
    d = 0.0
    while i < len(a) and j < len(b):
        x = min(a[i], b[j])
        while i < len(a) and a[i] <= x:
            i += 1
        while j < len(b) and b[j] <= x:
            j += 1
        d = max(d, abs(i / len(a) - j / len(b)))
    return round(d, 3)


def per_read(panel_dir: Path, truth=None):
    graph = json.loads((panel_dir / 'graph.json').read_text())
    chm13 = graph['haplotypes'].index(next(h for h in graph['haplotypes'] if h.startswith('CHM13#0')))
    on_chm13 = set(graph['paths'][chm13])
    deg = [0] * len(graph['lengths'])
    for p in graph['paths']:
        for v in p:
            deg[v] += 1
    reads = {}
    with (panel_dir / 'instances.jsonl').open() as f:
        for s in f:
            c = json.loads(s)
            reads.setdefault(c['read'], {})[c['strand']] = c
    out = []
    for name, by in reads.items():
        major = max('+-', key=lambda s: (len(by[s]['anchors']), s == '+'))
        minor = '-' if major == '+' else '+'
        A = by[major]['anchors']
        qlen = by[major]['qlen']
        row = dict(read=name, qlen=qlen, major=major, anchors_major=len(A), anchors_minor=len(by[minor]['anchors']),
                   anchors_major_per_kb=1000 * len(A) / qlen,
                   anchors_total_per_kb=1000 * (len(A) + len(by[minor]['anchors'])) / qlen,
                   minor_share=len(by[minor]['anchors']) / max(1, len(A) + len(by[minor]['anchors'])),
                   states_per_anchor=sum(deg[a[0]] for a in A) / max(1, len(A)),
                   on_chm13_share=sum(a[0] in on_chm13 for a in A) / max(1, len(A)),
                   anchors_per_end_position=len(A) / max(1, len({a[4] for a in A})))
        if truth is not None:
            row['true_strand'] = truth[name]
        out.append(row)
    return out


def counters(path: Path, panel: str, majors: dict):
    rows = [json.loads(s) for s in path.read_text().splitlines()]
    sel = [r for r in rows if r['panel'] == panel and r['gamma'] == lc.GAMMA and r['id'] == f"{r['id'].rsplit(':', 1)[0]}:{majors[r['id'].rsplit(':', 1)[0]]}"]
    return {k: [r[k] / r['N'] for r in sel] for k in COUNTER_KEYS} | {
        'share_' + k: [r[k] / max(1, r['pairs']) for r in sel] for k in ('pairs_same_vertex', 'pairs_loose', 'pairs_range', 'pairs_sliced')}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, default=HERE / 'realism_check.json')
    args = ap.parse_args()
    if args.out.exists():
        raise SystemExit('refusing to overwrite ' + str(args.out))
    truth = {}
    for line in (lc.DATA / 'mhc_sim/reads.tsv').read_text().splitlines()[1:]:
        a = line.split('\t')
        truth[a[0]] = a[5]
    result = {}
    for panel in ('hprc463', 'hprc60'):
        real = per_read(lc.MHC / 'data' / panel)
        sim = per_read(lc.DATA / 'mhc_sim' / panel, truth)
        cmp = {}
        for key in ('qlen', 'anchors_major', 'anchors_minor', 'anchors_major_per_kb', 'anchors_total_per_kb', 'minor_share',
                    'states_per_anchor', 'on_chm13_share', 'anchors_per_end_position'):
            a, b = [r[key] for r in real], [r[key] for r in sim]
            qa, qb = quartiles(a), quartiles(b)
            cmp[key] = dict(real=qa, simulated=qb, median_ratio_sim_over_real=round(qb[1] / qa[1], 3) if qa[1] else None,
                            ks_distance=ks(a, b))
        cc = {}
        ra = counters(HERE / 'counters_mhc.jsonl', panel, {r['read']: r['major'] for r in real})
        sa = counters(HERE / 'counters_mhc_sim.jsonl', panel, {r['read']: r['major'] for r in sim})
        for key in ra:
            qa, qb = quartiles(ra[key]), quartiles(sa[key])
            cc[key] = dict(real=qa, simulated=qb, ks_distance=ks(ra[key], sa[key]))
        result[panel] = dict(reads_real=len(real), reads_simulated=len(sim),
                             simulated_major_equals_true_strand=sum(r['major'] == r['true_strand'] for r in sim),
                             per_read=cmp, class_counters_major_orientation=cc)
        print(panel, json.dumps({k: (v['real'][1], v['simulated'][1], v['ks_distance']) for k, v in cmp.items()}))
        print(panel, json.dumps({k: (v['real'][1], v['simulated'][1], v['ks_distance']) for k, v in cc.items()}))
    lc.dump(args.out, dict(script='experiments/second_locus/realism_check.py', gamma=lc.GAMMA, scale=lc.SCALE,
                           note='quartiles are [p25, median, p75] over reads; counters use the major orientation instance at gamma 10000',
                           panels=result))


if __name__ == '__main__':
    main()
