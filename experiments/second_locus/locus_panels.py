"""Loaders for the new-locus panels (build_locus.py, simulate_reads.py, convert_locus.py).

Mirrors experiments/hprc_mhc_panel/panel_common.py with the locus as a parameter. Each panel
lives in data/<locus>/<panel>/ (graph.json, instances.jsonl). The 40 primary instances are
the first 20 simulated reads in both orientations (simulation.json `primary_ids`), as the
MHC pilot's primary set was its first 20 sampled reads. The independent oracle is
real_check_mosaic.py's graph-major per-haplotype `mine`, executed from its own source text on
the panel's graph (its own offsets, last2reach and distances), exactly as panel_common does.
"""
from __future__ import annotations
import bisect
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import locus_common as lc  # noqa: E402

if str(lc.NEW) not in sys.path:
    sys.path.insert(0, str(lc.NEW))
import gap_chaining as gc  # noqa: E402

SCALE, GAMMA = lc.SCALE, lc.GAMMA


def panels(locus: str) -> dict:
    """{'all': name, 'subsets': [...], 'leaveout': name} from data/<locus>/panels.json."""
    meta = json.loads((lc.DATA / locus / 'panels.json').read_text())
    names = [p['panel'] for p in meta['panels']]
    return dict(all=names[0], subsets=[n for n in names[1:] if not n.endswith('_leaveout')],
                leaveout=[n for n in names if n.endswith('_leaveout')][0], meta=meta)


def primary_ids(locus: str) -> list[str]:
    ids = json.loads((lc.DATA / locus / 'simulation.json').read_text())['primary_ids']
    assert len(ids) == 40
    return ids


def read_groups(locus: str) -> dict:
    """read id -> source group ('absent' or 'present') and source haplotype."""
    rows = (lc.DATA / locus / 'reads.tsv').read_text().splitlines()
    keys = rows[0].split('\t')
    out = {}
    for line in rows[1:]:
        r = dict(zip(keys, line.split('\t')))
        out[r['read']] = dict(group=r['group'], source=r['source'], strand=r['strand'])
    return out


def load_panel(locus: str, name: str):
    d = lc.DATA / locus / name
    graph = json.loads((d / 'graph.json').read_text())
    cases = {}
    with (d / 'instances.jsonl').open() as f:
        for s in f:
            c = json.loads(s)
            cases[c['id']] = c
    return graph, cases


def primary_cases(locus: str, name: str) -> list[dict]:
    _, cases = load_panel(locus, name)
    return [cases[i] for i in primary_ids(locus)]


def prepared(locus: str, name: str) -> 'gc.Prepared':
    graph, _ = load_panel(locus, name)
    return gc.Graph(graph['lengths'], [tuple(e) for e in graph['edges']], graph['paths']).prepare()


def anchors(case: dict, scale: int = SCALE) -> list['gc.Anchor']:
    return [gc.Anchor(*x[:5], x[5] * scale) for x in case['anchors']]


def independent_mine(locus: str, name: str, gamma: int) -> dict:
    """Namespace with real_check_mosaic.py's index and `mine`, built on this panel's graph."""
    graph, _ = load_panel(locus, name)
    path = lc.ROOT / 'experiments/gap_solution_review/real_check_mosaic.py'
    text = path.read_text()
    body = '# ---- independent index' + text.split('# ---- independent index', 1)[1]
    body = body.split('prep = gc.Graph(lengths, edges, paths).prepare()', 1)[0]
    ns = dict(time=time, bisect=bisect, lengths=list(graph['lengths']), edges=[tuple(e) for e in graph['edges']],
              paths=[list(p) for p in graph['paths']], V=len(graph['lengths']), H=len(graph['paths']), GAMMA=gamma)
    exec(compile(body, str(path), 'exec'), ns)
    return ns


def native_library_ready() -> str:
    """Path of the cached native library for the current gap_native.cpp; refuse to compile.

    gap_native.build() compiles into results/native_build/ when its cache entry is missing;
    this experiment must not write outside its folder, so it requires an existing build.
    """
    import hashlib
    import os
    source = lc.NEW / 'gap_native.cpp'
    compiler = os.environ.get('CXX', 'clang++')
    flags = ['-std=c++20', '-O3', '-shared', '-fPIC']
    identity = source.read_bytes() + repr((compiler, flags)).encode()
    suffix = '.dylib' if sys.platform == 'darwin' else '.so'
    target = lc.ROOT / 'results/native_build' / ('gap-' + hashlib.sha256(identity).hexdigest()[:16] + suffix)
    if not target.exists():
        raise SystemExit(f'native library {target.name} is not built; run gap_native.build() once from the shared '
                         'workflow (e.g. scripts/project.py verify --suite gap) before this experiment')
    return str(target)
