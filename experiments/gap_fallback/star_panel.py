"""Adversarial star panel (experiments/proof_review/pair_cost.py): class vs per-haplotype vs auto.

Label h runs 0 -> h+1 -> sink; one private anchor per label with distinct weights, then k
anchors on the sink, so contested merges have |L| = |R| = H. Outputs of all three modes must be
identical; records kernel times (median of 5 shuffled paired repetitions), class work, the
per-haplotype estimate and whether auto fell back. Writes star_panel.json.
"""
import json, random, statistics, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'research/current/new'))
import gap_chaining as gc, gap_native as gn


def star(m, k):
    sink = m + 1
    lengths = [1] + [1] * m + [4]
    edges = [(0, h + 1) for h in range(m)] + [(h + 1, sink) for h in range(m)]
    paths = [[0, h + 1, sink] for h in range(m)]
    anchors = [gc.Anchor(h + 1, 1, 1, h + 1, h + 1, 10 * m + h) for h in range(m)]
    for t in range(k):
        anchors.append(gc.Anchor(sink, 1 + t % 4, 1 + t % 4, m + 1 + t, m + 1 + t, 1))
    return gc.Graph(lengths, edges, paths), anchors, 10 ** 6


rows, rng = [], random.Random(11)
for m in (8, 32, 64, 128, 256, 512, 1024):
    for k in (16, 64):
        g, a, gamma = star(m, k)
        G = gn.GapGraph(g.prepare()); b = G.batch(a, gamma)
        ref = G.solve_packed(b, mode='per_haplotype')
        t = {x: [] for x in ('class', 'per_haplotype', 'auto')}
        auto_meta = None
        for _ in range(5):
            order = list(t); rng.shuffle(order)
            for x in order:
                r = G.solve_packed(b, mode=x)
                assert list(r.C) == list(ref.C) and list(r.B) == list(ref.B) and list(r.R) == list(ref.R), (m, k, x)
                t[x].append(r.stats['kernel_seconds'] * 1e3)
                if x == 'auto':
                    auto_meta = r.stats
        row = dict(H=m, k=k, N=b.n, **{x + '_ms': statistics.median(v) for x, v in t.items()},
                   fallback=auto_meta['fallback'], work_at_stop=auto_meta['auto_work'],
                   estimate=auto_meta['auto_estimate'], budget=auto_meta['auto_budget'])
        row['class_over_ph'] = row['class_ms'] / row['per_haplotype_ms']; row['auto_over_ph'] = row['auto_ms'] / row['per_haplotype_ms']
        rows.append(row); print(json.dumps(row), flush=True)
(HERE / 'star_panel.json').write_text(json.dumps(dict(library=gn.build()._name, rows=rows), indent=1) + '\n')
