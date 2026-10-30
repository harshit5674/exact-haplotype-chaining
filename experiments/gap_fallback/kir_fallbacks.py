"""On the KIR instances where auto mode fell back: class, per-haplotype and auto kernel times (5 paired reps)."""
import json, statistics, sys, random
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'research/current/new')); sys.path.insert(0, str(HERE))
import gap_chaining as gc, gap_native as gn
import auto_real as ar
rows = []
rng = random.Random(7)
for locus, name, d in ar.panels():
    if locus != 'kir':
        continue
    graph, cases = ar.load(d)
    G = gn.GapGraph(gc.Graph(graph['lengths'], [tuple(e) for e in graph['edges']], graph['paths']).prepare())
    for c in ar.primary(cases, locus):
        a = [gc.Anchor(*x[:5], x[5] * ar.SCALE) for x in c['anchors']]
        b = G.batch(a, ar.GAMMA)
        if not G.solve_packed(b, mode='auto').stats['fallback']:
            continue
        t = {m: [] for m in ('class', 'per_haplotype', 'auto')}
        for _ in range(5):
            modes = list(t); rng.shuffle(modes)
            for m in modes:
                t[m].append(G.solve_packed(b, mode=m).stats['kernel_seconds'] * 1e3)
        row = dict(panel=name, id=c['id'], N=b.n, **{m + '_ms': statistics.median(v) for m, v in t.items()})
        rows.append(row); print(json.dumps(row), flush=True)
(HERE / 'kir_fallbacks.json').write_text(json.dumps(rows, indent=1) + '\n')
