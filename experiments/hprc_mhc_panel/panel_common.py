"""Loaders for the real HPRC MHC panels built by build_panel.py and convert.py.

Each panel lives in data/<name>/ (graph.json, instances.jsonl). The 40 primary instances are
the same (read, strand) identifiers as the MHC pilot's primary set
(experiments/mhc_pilot/final_hybrid.jsonl), so results line up with the 60-haplotype work.
The independent oracle is real_check_mosaic.py's graph-major per-haplotype `mine`, executed
from its own source text on the panel's graph (its own offsets, last2reach and distances).
"""
from __future__ import annotations
import bisect
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
NEW = ROOT / 'research/current/new'
if str(NEW) not in sys.path:
    sys.path.insert(0, str(NEW))
import gap_chaining as gc  # noqa: E402

SCALE, GAMMA = 200, 10000
PANELS = ('hprc60', 'hprc120', 'hprc240', 'hprc463')


def primary_ids() -> list[str]:
    rows = (ROOT / 'experiments/mhc_pilot/final_hybrid.jsonl').read_text().splitlines()
    ids = [json.loads(s)['id'] for s in rows]
    assert len(ids) == 40
    return ids


def load_panel(name: str):
    d = HERE / 'data' / name
    graph = json.loads((d / 'graph.json').read_text())
    cases = {json.loads(s)['id']: json.loads(s) for s in (d / 'instances.jsonl').read_text().splitlines()}
    return graph, cases


def primary_cases(name: str) -> list[dict]:
    _, cases = load_panel(name)
    return [cases[i] for i in primary_ids()]


def prepared(name: str) -> 'gc.Prepared':
    graph, _ = load_panel(name)
    return gc.Graph(graph['lengths'], [tuple(e) for e in graph['edges']], graph['paths']).prepare()


def anchors(case: dict, scale: int = SCALE) -> list['gc.Anchor']:
    return [gc.Anchor(*x[:5], x[5] * scale) for x in case['anchors']]


def independent_mine(name: str, gamma: int) -> dict:
    """Namespace with real_check_mosaic.py's index and `mine`, built on this panel's graph."""
    graph, _ = load_panel(name)
    path = ROOT / 'experiments/gap_solution_review/real_check_mosaic.py'
    text = path.read_text()
    body = '# ---- independent index' + text.split('# ---- independent index', 1)[1]
    body = body.split('prep = gc.Graph(lengths, edges, paths).prepare()', 1)[0]
    ns = dict(time=time, bisect=bisect, lengths=list(graph['lengths']), edges=[tuple(e) for e in graph['edges']],
              paths=[list(p) for p in graph['paths']], V=len(graph['lengths']), H=len(graph['paths']), GAMMA=gamma)
    exec(compile(body, str(path), 'exec'), ns)
    return ns
