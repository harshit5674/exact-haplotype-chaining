"""Differential harness: every solver configuration against the independent brute force.

For one instance (lengths, edges, paths, anchors, gamma) it runs
  spec        gap_chaining.OriginSolver (the external specification)
  oracle      gap_chaining.direct_oracle
  scalar      gap_chaining.scalar_fenwick
  fp_*        gap_fastpath.FastPathSolver in seven configurations (ablations, prompt plan,
              full plan, floor cuts, grouped remote without cut)
  nat_*       gap_native: class compact (+expand), class packed cold/warm, per-haplotype
              lanes with dense_cut in {default, 0 (always dense), 257 (always sparse)}
and requires, for every output:
  * C, B, R equal to brute_chains (explicit chain enumeration when small enough, else the
    explicit chain-step DAG; both are cross-checked when both run);
  * gap_chaining.verify_pointers passes (package checker), and
  * brute_chains.trace_and_score passes (independent traceback under the chain semantics);
and, for fast paths and the native compact engine, record identity with the specification
(per-anchor records; for fast paths also every Fenwick node's contents).
Failures are returned with the full input.
"""
from __future__ import annotations
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'research/current/new'))
sys.path.insert(0, str(HERE))

import gap_chaining as gc  # noqa: E402
import gap_fastpath as fp  # noqa: E402
import brute_chains as bc  # noqa: E402

try:
    import gap_native as gn  # noqa: E402
except Exception as exc:  # pragma: no cover
    gn = None
    NATIVE_IMPORT_ERROR = repr(exc)

FP_CONFIGS = {
    'fp_none': dict(tests='none', remote='labels', cut='none'),
    'fp_i': dict(tests='i', remote='labels', cut='none'),
    'fp_ii': dict(tests='ii', remote='labels', cut='none'),
    'fp_prompt': dict(tests='iii', remote='labels', cut='none'),
    'fp_full': dict(tests='iii', remote='grouped', cut='range'),
    'fp_zero': dict(tests='none', remote='labels', cut='zero'),
    'fp_grouped': dict(tests='iii', remote='grouped', cut='none'),
}
ENUM_CAP = 60000


def node_sets(solver):
    return [sorted((r.coef, r.origin, r.mask, r.witness.kind, r.witness.j, r.witness.g, r.witness.w)
                   for r in node) for node in solver.nodes]


def wit_tuple(w):
    if w.kind == 'restart':
        return ('restart', -1, -1, -1)
    if w.kind == 'local':
        return ('local', w.j, -1, -1)
    if w.kind == 'remote':
        return ('remote', w.j, w.g, w.w)
    return ('same', w.j, w.g, -1)


def flat(C):
    return {(i, h): v for i, row in enumerate(C) for h, v in row.items()}


def flat_ptr(P):
    return {(i, h): wit_tuple(w) for i, row in enumerate(P) for h, w in row.items()}


class GraphCache:
    """Per-graph objects reused across many anchor sets (Prepared, native graph, facts)."""

    def __init__(self, lengths, edges, paths, prepared_cls=None):
        self.key = (tuple(lengths), tuple(map(tuple, edges)), tuple(map(tuple, paths)))
        self.graph = gc.Graph(list(lengths), [tuple(e) for e in edges], [list(p) for p in paths])
        self.p = self.graph.prepare() if prepared_cls is None else prepared_cls(self.graph)
        self.facts = bc.Facts(lengths, edges, paths)
        self.native = None
        self.native_error = None
        if gn is not None:
            try:
                self.native = gn.GapGraph(self.p)
            except gn.NativeRangeError as e:
                self.native_error = 'NativeRangeError: ' + str(e)
        self.fp_bounds = fp.VertexBounds(self.p)
        self.fp_groups = fp.RemoteGroups(self.p)


def brute(gcache, anchors, gamma, stats):
    facts = gcache.facts
    table = bc.step_table(facts, anchors, gamma)
    dag = bc.chain_dag_optimum(facts, anchors, gamma, table)
    en = bc.enumerate_chains(facts, anchors, gamma, cap=ENUM_CAP)
    if en is not None:
        stats['enumerated'] += 1
        stats['chains'] += en[4]
        if (en[0], en[1], en[2], en[3]) != dag:
            return None, ('brute self-check: enumeration and chain DAG disagree', en[:4], dag), table
        return en[:4], None, table
    stats['dag_only'] += 1
    return dag, None, table


def check_instance(gcache, anchors_t, gamma, stats, native=True, fp_configs=FP_CONFIGS,
                   pointer_solvers=None):
    """Run everything on one instance; return a list of failure descriptions."""
    failures = []
    anchors = [gc.Anchor(*a) for a in anchors_t]
    p = gcache.p
    brute_values, err, table = brute(gcache, anchors_t, gamma, stats)
    if err:
        return [dict(solver='brute', kind='self-check', detail=str(err)[:3000])]
    bC, bCnot, bB, bR = brute_values
    outputs = {}

    # Specification.
    spec = gc.OriginSolver(p, anchors, gamma)
    rs = spec.run()
    C, ptr = rs.expand()
    outputs['spec'] = (flat(C), rs.B, rs.R, flat_ptr(ptr), C, ptr)
    spec_nodes = node_sets(spec)
    spec_records = rs.records
    # Package references.
    oC, oB, oR, optr = gc.direct_oracle(p, anchors, gamma)
    outputs['oracle'] = (flat(oC), oB, oR, flat_ptr(optr), oC, optr)
    sC, sB, sR, sptr = gc.scalar_fenwick(p, anchors, gamma)
    outputs['scalar'] = (flat(sC), sB, sR, flat_ptr(sptr), sC, sptr)
    # Fast paths: record identity with the specification.
    for name, cfg in fp_configs.items():
        kw = dict(cfg, bounds=gcache.fp_bounds)
        if cfg['remote'] == 'grouped':
            kw['groups'] = gcache.fp_groups
        fs = fp.FastPathSolver(p, anchors, gamma, **kw)
        rf = fs.run()
        if rf.records != spec_records or rf.B != rs.B or rf.R != rs.R:
            failures.append(dict(solver=name, kind='records differ from spec'))
        if node_sets(fs) != spec_nodes:
            failures.append(dict(solver=name, kind='node contents differ from spec'))
        fC, fptr = rf.expand()
        outputs[name] = (flat(fC), rf.B, rf.R, flat_ptr(fptr), fC, fptr)
        for k in ('pairs_level4', 'sliced_thresholds', 'sliced_minima', 'remote_sliced',
                  'remote_pruned', 'cut_records', 'min_pruned'):
            stats['fp_' + k] += rf.cost.get(k, 0)
    # Native engines.
    if native and gcache.native is not None:
        G = gcache.native
        try:
            batch = G.batch(anchors, gamma)
        except gn.NativeRangeError:
            stats['native_rejected'] += 1
            batch = None
        if batch is not None:
            comp = G.solve_compact(batch)
            if comp.records() != spec_records or list(comp.B) != rs.B or comp.R_values() != rs.R:
                failures.append(dict(solver='nat_compact', kind='records differ from spec'))
            runs = [('nat_compact_expand', comp.expand()),
                    ('nat_class', G.solve_packed(batch, mode='class')),
                    ('nat_class_warm', G.solve_packed(batch, mode='class', warm=True))]
            for cut in (gn.DENSE_CUT, 0, 257):
                try:
                    runs.append(('nat_ph_%d' % cut, G.solve_packed(batch, mode='per_haplotype', dense_cut=cut)))
                except gn.NativeRangeError:
                    stats['native_ph_rejected'] += 1
                    break
            for name, res in runs:
                nC, nptr, nB, nR = res.as_dicts(anchors)
                outputs[name] = (flat(nC), nB, nR, flat_ptr(nptr), nC, nptr)
            for k in ('pairs_sliced', 'sliced_thresholds', 'sliced_minima', 'remote_sliced', 'cut_records'):
                stats['nat_' + k] += comp.stats.get(k, 0)
    elif native and gcache.native is None:
        stats['native_graph_rejected'] += 1

    # Compare every output with the brute force and check pointers twice.
    for name, (fC, fB, fR, fptr, rawC, rawptr) in outputs.items():
        if fC != bC:
            bad = sorted(k for k in set(fC) | set(bC) if fC.get(k) != bC.get(k))[:6]
            failures.append(dict(solver=name, kind='C differs from brute force',
                                 detail=[(k, fC.get(k), bC.get(k)) for k in bad]))
        if list(fB) != bB:
            failures.append(dict(solver=name, kind='B differs', detail=[list(fB), bB]))
        if list(fR) != bR:
            failures.append(dict(solver=name, kind='R differs', detail=[list(fR), bR]))
        if pointer_solvers is not None and name not in pointer_solvers:
            continue
        try:
            counts = gc.verify_pointers(p, anchors, gamma, rawC, list(fB), rawptr)
            for k, v in counts.items():
                stats['ptr_%s' % k] += v
        except AssertionError as e:
            failures.append(dict(solver=name, kind='verify_pointers failed', detail=repr(e)[:500]))
        probs = bc.trace_and_score(gcache.facts, anchors_t, gamma, fC, fptr, table)
        if probs:
            failures.append(dict(solver=name, kind='independent traceback failed', detail=str(probs[:5])))
        stats['outputs_checked'] += 1
    stats['instances'] += 1
    stats['states'] += len(bC)
    return failures


def failure_record(lengths, edges, paths, anchors, gamma, failures, family):
    return dict(family=family, lengths=list(lengths), edges=[list(e) for e in edges],
                paths=[list(q) for q in paths], anchors=[list(a) for a in anchors], gamma=gamma,
                failures=failures)
