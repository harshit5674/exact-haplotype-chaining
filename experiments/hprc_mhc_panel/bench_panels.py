"""Paired native benchmark on the real HPRC MHC panels (60, 120, 240 and 463 haplotypes).

Reuses experiments/gap_engine_native/benchmark.py unchanged (imported, not edited): the same
modes (per_haplotype, class_expanded, class_compact, class_compact_warm), 5 repetitions in
a freshly shuffled order per read, gc.collect() before every call, kernel and packed times,
and per-read exactness: class compact records equal to the Python specification's and
identical C, B, R across all modes. Inputs are the 40 primary (read, strand) instances,
Minichain scaling (weights x 200) and gamma 10000 and 20000.
"""
from __future__ import annotations
import argparse
import importlib.util
import json
import os
from pathlib import Path
import platform
import random
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panel_common as pc  # noqa: E402
import gap_native as gn  # noqa: E402

spec = importlib.util.spec_from_file_location('gap_benchmark', pc.ROOT / 'experiments/gap_engine_native/benchmark.py')
bm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bm)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--panels', default=','.join(pc.PANELS))
    ap.add_argument('--gamma', type=lambda s: [int(x) for x in s.split(',')], default=[10000, 20000])
    ap.add_argument('--repeat', type=int, default=5)
    ap.add_argument('--seed', type=int, default=9090)
    args = ap.parse_args()
    meta_path = args.out.with_suffix('.meta.json')
    if args.out.exists() or meta_path.exists():
        raise SystemExit('refusing to overwrite')
    lib = gn.build()._name
    rng = random.Random(args.seed)
    load_before, started, prep = os.getloadavg(), time.time(), {}
    with args.out.open('w') as out:
        for name in args.panels.split(','):
            t = time.perf_counter()
            p = pc.prepared(name)
            G = gn.GapGraph(p)
            prep[name] = dict(H=p.H, V=p.V, W=G.W, b=G.b, python_prepared_seconds=time.perf_counter() - t,
                              native_index_seconds=G.native_index_seconds, remote_groups=G.remote_groups,
                              remote_entries=G.remote_entries)
            print(name, prep[name], flush=True)
            cases = pc.primary_cases(name)
            for gamma in args.gamma:
                bm._ANCHORS.clear()  # anchor cache is keyed by instance id, which panels share
                bm.level(G, p, p.H, gamma, cases, args, rng, out)
            del G
    compiler = subprocess.run(['clang++', '--version'], capture_output=True, text=True).stdout.splitlines()[0]
    meta = dict(script='experiments/hprc_mhc_panel/bench_panels.py', argv=sys.argv[1:], python=sys.version,
                platform=platform.platform(), machine=platform.machine(), compiler=compiler, library=lib,
                panels=prep, load_average_before=load_before, load_average_after=os.getloadavg(),
                wall_seconds=round(time.time() - started, 1), modes=bm.MODES, repeat=args.repeat, seed=args.seed)
    meta_path.write_text(json.dumps(meta, indent=1) + '\n')
    print(json.dumps(meta), flush=True)


if __name__ == '__main__':
    main()
