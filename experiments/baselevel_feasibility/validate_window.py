"""Full independent validation of a base-level window panel with the memory-light validator (2026-10-20).

Every instance of windows/<name>/panel: the per-haplotype, class (expanded), guarded and compact (expanded)
outputs are each checked by research/current/new/gap_validation_light.py (graph-major recurrence with its own
last-vertex and distance computation; every C, B, R and every pointer, including acyclicity). The light
validator was cross-checked against gap_validation.py in validator_light/. Nothing is timed.

  python3.13 -B experiments/baselevel_feasibility/validate_window.py --window windows/w31250k_250kb
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'research/current/new'))
import gap_chaining as gc  # noqa: E402
import gap_native as gn  # noqa: E402
import gap_validation_light as gvl  # noqa: E402

GAMMA, SCALE = 10000, 200


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--window', type=Path, required=True)
    args = ap.parse_args()
    win = args.window if args.window.is_absolute() else HERE / args.window
    out = win / 'validation.json'
    assert not out.exists(), 'refuse to overwrite'
    panel = win / 'panel'
    graph = json.loads((panel / 'graph.json').read_text())
    t0 = time.perf_counter()
    p = gc.Graph(graph['lengths'], [tuple(e) for e in graph['edges']], graph['paths']).prepare()
    G = gn.GapGraph(p)
    validator = gvl.Validator(graph)
    index_seconds = time.perf_counter() - t0
    instances = states = checks = 0
    for line in (panel / 'instances.jsonl').open():
        c = json.loads(line)
        A = [gc.Anchor(*x[:5], x[5] * SCALE) for x in c['anchors']]
        batch = G.batch(A, GAMMA)
        vb = validator.batch(A, GAMMA)
        for mode in ('per_haplotype', 'class', 'auto', 'compact'):
            z = G.solve_compact(batch, warm=False).expand() if mode == 'compact' else G.solve_packed(batch, mode=mode, warm=False)
            vb.check(z)
            checks += 1
        vb.close()
        instances += 1
        states += batch.K
    validator.close()
    res = dict(panel=str(panel.relative_to(ROOT)), instances=instances, states=states, engine_outputs_checked=checks,
               validator='research/current/new/gap_validation_light.py', status='passed',
               index_seconds=index_seconds, seconds=time.perf_counter() - t0)
    out.write_text(json.dumps(res, indent=1) + '\n')
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
