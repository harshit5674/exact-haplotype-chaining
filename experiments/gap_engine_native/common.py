"""Shared loaders for the gap-engine experiments.

Real data: the MHC graph and the 40 primary reads of experiments/mhc_pilot, with
Minichain scaling (anchor weight x 200) and gamma = 10000 unless stated otherwise.
Mosaic haplotypes come from experiments/class_chaining/scaling.py's generator
(seed 6060, switch rate 0.02), imported unchanged: the first H-60 mosaics of one
900-mosaic stream, exactly as in that experiment and in real_check_mosaic.py.

Independent oracles are loaded from experiments/gap_solution_review without editing
those files: the graph-major per-haplotype implementation `mine` of
real_check_mosaic.py, and the oracle/generator of independent_check.py.
"""
from __future__ import annotations
import importlib.util
import json
from pathlib import Path
import random
import sys

ROOT = Path(__file__).resolve().parents[2]
NEW = ROOT / 'research/current/new'
PILOT = ROOT / 'experiments/mhc_pilot'
REVIEW = ROOT / 'experiments/gap_solution_review'
if str(NEW) not in sys.path:
    sys.path.insert(0, str(NEW))
import gap_chaining as gc  # noqa: E402  (working copy; cached before any oracle import)

SCALE, GAMMA = 200, 10000
MOSAIC_SEED, MOSAIC_RATE = 6060, 0.02
H_LEVELS = (60, 120, 240, 480, 960)


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_graph() -> dict:
    return json.loads((PILOT / 'data/graph.json').read_text())


def primary_cases() -> list[dict]:
    """The 40 primary reads, in instances.jsonl order."""
    primary = {json.loads(s)['id'] for s in (PILOT / 'final_hybrid.jsonl').read_text().splitlines()}
    cases = [json.loads(s) for s in (PILOT / 'data/instances.jsonl').read_text().splitlines()]
    out = [c for c in cases if c['id'] in primary]
    assert len(out) == 40
    return out


_MOSAICS = None


def panel(H: int) -> list[list[int]]:
    """The 60 real paths followed by the first H-60 mosaics of the seed-6060 stream."""
    global _MOSAICS
    real = [list(p) for p in load_graph()['paths']]
    if H < len(real) or H > 960:
        raise ValueError('H must lie in [60, 960]')
    if _MOSAICS is None:
        scaling = load_module('cc_scaling', ROOT / 'experiments/class_chaining/scaling.py')
        _MOSAICS, _ = scaling.mosaics(real, 960 - len(real), MOSAIC_RATE, random.Random(MOSAIC_SEED))
    return real + [list(p) for p in _MOSAICS[:H - len(real)]]


def prepared(H: int) -> 'gc.Prepared':
    g = load_graph()
    return gc.Graph(g['lengths'], [tuple(e) for e in g['edges']], panel(H)).prepare()


def anchors(case: dict, scale: int = SCALE) -> list['gc.Anchor']:
    return [gc.Anchor(*x[:5], x[5] * scale) for x in case['anchors']]


def independent_mine(H: int) -> dict:
    """Namespace holding real_check_mosaic.py's independent per-haplotype `mine(A)`.

    Its index (offsets, last2reach, shortest distances) is its own. Set
    ns['GAMMA'] before calling ns['mine'] to change the penalty.
    """
    path = REVIEW / 'real_check_mosaic.py'
    text = path.read_text()
    head = text.split('prep = gc.Graph(lengths, edges, paths).prepare()')[0]
    saved_argv, saved_path = sys.argv, list(sys.path)
    ns = {'__file__': str(path), '__name__': 'real_check_mosaic_independent'}
    try:
        sys.argv = [str(path), '1', str(H)]
        exec(compile(head, str(path), 'exec'), ns)
    finally:
        sys.argv, sys.path[:] = saved_argv, saved_path
    assert ns['H'] == H and ns['gc'] is gc
    assert [list(p) for p in ns['paths']] == panel(H), 'mosaic stream differs from scaling.py'
    return ns


def independent_oracle() -> dict:
    """Namespace with independent_check.py's own index, oracle(...) and gen(rng)."""
    path = REVIEW / 'independent_check.py'
    head = path.read_text().split('rng = random.Random(777)')[0]
    saved_path = list(sys.path)
    ns = {'__file__': str(path), '__name__': 'independent_check_oracle'}
    try:
        exec(compile(head, str(path), 'exec'), ns)
    finally:
        sys.path[:] = saved_path
    assert ns['gc'] is gc
    return ns
