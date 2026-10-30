"""Auto mode (work-budget fallback) on every real panel and locus: does it ever fall back?

For each panel's 40 primary instances (weights x200, gamma 10000) this runs
GapGraph.solve_packed(mode='auto') at the default factor, checks C, B, R against mode
'per_haplotype', and records whether the fallback ran and the ratio of class work (words) to
the per-haplotype estimate (lane operations). Panels: the real-read MHC panels
(experiments/hprc_mhc_panel/data/hprc{60,120,240,463}) and the simulated-read KIR and chr20
panels of experiments/second_locus (all levels). Writes auto_real.json.
"""
from __future__ import annotations

import json
import statistics
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'research/current/new'))
import gap_chaining as gc  # noqa: E402
import gap_native as gn  # noqa: E402

GAMMA, SCALE = 10000, 200


def panels():
    mhc = ROOT / 'experiments/hprc_mhc_panel/data'
    for name in ('hprc60', 'hprc120', 'hprc240', 'hprc463'):
        yield 'mhc', name, mhc / name
    for locus in ('kir', 'chr20'):
        base = ROOT / 'experiments/second_locus/data' / locus
        for d in sorted(p for p in base.iterdir() if p.is_dir() and (p / 'instances.jsonl').exists()):
            yield locus, d.name, d


def load(d: Path):
    graph = json.loads((d / 'graph.json').read_text())
    cases = [json.loads(line) for line in open(d / 'instances.jsonl')]
    return graph, cases


def primary(cases, locus):
    """The 40 primary instances: the MHC panels' fixed list (panel_common.primary_ids), else the first 20
    reads in file order, both strands (experiments/second_locus)."""
    if locus == 'mhc':
        sys.path.insert(0, str(ROOT / 'experiments/hprc_mhc_panel'))
        import panel_common as pc
        want = pc.primary_ids()
        byid = {c['id']: c for c in cases}
        return [byid[i] for i in want]
    reads = []
    for c in cases:
        r = c['id'].rsplit(':', 1)[0]
        if r not in reads:
            reads.append(r)
    keep = set(reads[:20])
    return [c for c in cases if c['id'].rsplit(':', 1)[0] in keep]


def main():
    out = dict(gamma=GAMMA, scale=SCALE, auto_factor=gn.AUTO_FACTOR, library=gn.build()._name, panels=[])
    t0 = time.perf_counter()
    for locus, name, d in panels():
        graph, cases = load(d)
        p = gc.Graph(graph['lengths'], [tuple(e) for e in graph['edges']], graph['paths']).prepare()
        G = gn.GapGraph(p)
        sel = primary(cases, locus)
        ratios, fallbacks, states = [], 0, 0
        for c in sel:
            a = [gc.Anchor(*x[:5], x[5] * SCALE) for x in c['anchors']]
            au = G.solve_packed(a, GAMMA, mode='auto')
            ph = G.solve_packed(a, GAMMA, mode='per_haplotype')
            assert list(au.C) == list(ph.C) and list(au.B) == list(ph.B) and list(au.R) == list(ph.R), (name, c['id'])
            fallbacks += au.stats['fallback']
            ratios.append(au.stats['auto_work'] / max(1.0, au.stats['auto_estimate']))
            states += len(au.C)
        ratios.sort()
        row = dict(locus=locus, panel=name, H=p.H, instances=len(sel), states=states, fallbacks=fallbacks,
                   ratio_min=ratios[0], ratio_median=statistics.median(ratios), ratio_max=ratios[-1])
        out['panels'].append(row)
        print(json.dumps(row), flush=True)
    out['seconds'] = time.perf_counter() - t0
    out['all_identical'] = True
    out['total_fallbacks'] = sum(r['fallbacks'] for r in out['panels'])
    (HERE / 'auto_real.json').write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(dict(total_fallbacks=out['total_fallbacks'], max_ratio=max(r['ratio_max'] for r in out['panels']))))


if __name__ == '__main__':
    main()
