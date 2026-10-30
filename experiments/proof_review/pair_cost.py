"""Cost probe: is the fast-path / native merge linear in list length? (PROOF.md section 14)

Star panel: label h runs 0 -> (private vertex h+1) -> sink. One private anchor per label
(distinct weights, so every label keeps its own origin), then 4 anchors on the sink. Each sink
anchor's state list has H singleton records with distinct origins, and consecutive sink
publications meet in shared Fenwick nodes, so contested merges have |L| = |R| = H.

Counted (no wall-clock timing):
  per contested merge of FastPathSolver (full plan), via a counting subclass defined HERE:
    |L|, |R|, records in + out (the unit of Lemma 1's O(W q (|L|+|R|)) merge) and the
    W-word mask operations of the pairwise loop (gap_fastpath.py:300-325, the same loop
    structure as gap_native.cpp:503-518);
  native gap_native compact: overlapping pairs decided (stats 'pairs') and contested merges;
  model  gap_fastpath.predicted_word_ops: class compact word operations against per-haplotype
         Fenwick steps (the unit model of the 2026-10-04 Stage A gate; not a timing).
Note: on this construction Lambda (all records read) also grows quadratically, because each
of the H private anchors reads O(H) earlier records in prefix filtering. The per-merge counts
isolate the merge loop. Writes results/pair_cost.json.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import harness as H  # noqa: E402
gc, fp, gn = H.gc, H.fp, H.gn


class CountingFastPath(fp.FastPathSolver):
    """FastPathSolver that logs (|L|, |R|, records out, mask ops) for every contested merge."""

    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        self.merge_log = []

    def merge(self, left, right):
        m0, c0 = self.cost['mask'], self.cost['merges_contested']
        ans = super().merge(left, right)
        if self.cost['merges_contested'] > c0:
            self.merge_log.append((len(left), len(right), len(ans), self.cost['mask'] - m0))
        return ans


def star(m, k):
    sink = m + 1
    lengths = [1] + [1] * m + [4]
    edges = [(0, h + 1) for h in range(m)] + [(h + 1, sink) for h in range(m)]
    paths = [[0, h + 1, sink] for h in range(m)]
    anchors = [gc.Anchor(h + 1, 1, 1, h + 1, h + 1, 10 * m + h) for h in range(m)]
    for t in range(k):
        anchors.append(gc.Anchor(sink, 1 + t % 4, 1 + t % 4, m + 1 + t, m + 1 + t, 1))
    return gc.Graph(lengths, edges, paths), anchors, 10 ** 6


def main():
    rows = []
    for m in (8, 16, 32, 64, 128, 256, 512):
        g, a, gamma = star(m, 4)
        p = g.prepare()
        rs = gc.OriginSolver(p, a, gamma).run()
        fs = CountingFastPath(p, a, gamma, tests='iii', remote='grouped', cut='range')
        rf = fs.run()
        assert rf.records == rs.records
        log = fs.merge_log
        recs = sum(l + r + o for l, r, o, _ in log)
        masks = sum(x for _, _, _, x in log)
        W = p.W
        K = sum(p.A[z.v].bit_count() for z in a)
        pred = fp.predicted_word_ops(rf.cost, W, K)
        row = dict(H=m, anchors=len(a), W=W, K=K, contested_merges=len(log),
                   largest_contested_merge=max(log, key=lambda t: t[3])[:2] if log else None,
                   contested_records_in_out=recs, contested_mask_ops=masks,
                   mask_ops_per_record=round(masks / recs, 2) if recs else None,
                   all_lambda_records=rf.stats['lambda_records'],
                   predicted_class_compact_word_ops=pred['class_compact'],
                   predicted_per_haplotype_steps=pred['per_haplotype'],
                   predicted_ratio_class_over_per_haplotype=round(pred['class_compact'] / pred['per_haplotype'], 2))
        G = gn.GapGraph(p)
        comp = G.solve_compact(a, gamma)
        assert comp.records() == rs.records
        row['native_pairs'] = int(comp.stats['pairs'])
        row['native_contested_merges'] = int(comp.stats['contested_merges'])
        rows.append(row)
        print(row)
    (HERE / 'results').mkdir(exist_ok=True)
    (HERE / 'results/pair_cost.json').write_text(json.dumps(rows, indent=1) + '\n')


if __name__ == '__main__':
    main()
