"""Effect of each candidate Minichain deviation on its own DP, from DP dumps of the same reads.

Compares a base DP-dump run (the pinned binary's behaviour) with variant runs of the
instrumented build:
  --truedist  order-independent shortest-distance relaxation (fixes MPC_index line 792)
  --gate      RMQ default (INT64_MIN,-1): the remote value is applied to a path's state only
              when that path has a local predecessor (the gate the source appears to intend)
For each (read, component) it reports whether the full state table (hash), B, max(w,R) or
the best score changed, and per read the change in the best chain score over both components.
"""
from __future__ import annotations

import argparse
import collections
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mcdp  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--base', required=True, type=Path)
    ap.add_argument('--variant', required=True, type=Path)
    ap.add_argument('--out', required=True, type=Path)
    args = ap.parse_args()
    base = mcdp.read_dp_dump(args.base / 'dump.tsv')
    var = mcdp.read_dp_dump(args.variant / 'dump.tsv')
    assert set(base) == set(var), 'runs cover different (read, component) pairs'
    counts = collections.Counter()
    per_read = collections.defaultdict(lambda: [None, None])
    examples = []
    for key, b in base.items():
        v = var[key]
        assert b['ha'] == v['ha'] and b['N'] == v['N'], ('different inputs', key)
        counts['pairs'] += 1
        for h in ('hc', 'hb', 'hp'):
            counts[h + '_changed'] += b[h] != v[h]
        counts['best_changed'] += b['best'] != v['best']
        if b['best'] != v['best'] and len(examples) < 50:
            examples.append(dict(read=key[0], cid=key[1], base=b['best'], variant=v['best']))
        r = per_read[key[0]]
        r[0] = b['best'] if r[0] is None else max(r[0], b['best'])
        r[1] = v['best'] if r[1] is None else max(r[1], v['best'])
    deltas = [vb - bb for bb, vb in per_read.values()]
    changed = [d for d in deltas if d]
    summary = dict(base=str(args.base), variant=str(args.variant), counts=dict(counts), reads=len(deltas),
                   reads_best_changed=len(changed), reads_best_increased=sum(d > 0 for d in changed),
                   reads_best_decreased=sum(d < 0 for d in changed),
                   delta_min=min(changed) if changed else 0, delta_max=max(changed) if changed else 0,
                   examples=examples)
    args.out.write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps({k: v for k, v in summary.items() if k != 'examples'}, indent=1))


if __name__ == '__main__':
    main()
