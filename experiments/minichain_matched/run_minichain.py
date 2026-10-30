"""Run the instrumented Minichain build (vendor/minichain-dp) once, recording everything needed
to reproduce the run: command, environment switches, binary/input hashes, host load, timings.

Every run writes into a NEW directory (refuses to reuse one: the dumps append). Mapping runs use
`-x lr -t 1` with OMP_NUM_THREADS=1 and no base-level alignment (-c), so Minichain's time is
dominated by graphUtils::Chaining, which RECOMB_DP_TIME measures directly.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BINARY = HERE / 'vendor/minichain-dp/minichain'
GRAPH = ROOT / 'experiments/mhc_pilot/data/MHC-CHM13.0.gfa.gz'
READS = ROOT / 'experiments/mhc_pilot/data/MHC_CHM13_PacBio_filt.fq.gz'


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True, type=Path, help='new run directory')
    ap.add_argument('--reads', type=Path, default=READS)
    ap.add_argument('--graph', type=Path, default=GRAPH)
    ap.add_argument('--R', type=int, default=10000, help='recombination penalty (-R)')
    ap.add_argument('--scale', type=int, default=None, help='-s scale factor (default: Minichain default 200)')
    ap.add_argument('--seeds', action='store_true', help='export seeds and skip chaining (RECOMB_DUMP_ANCHORS)')
    ap.add_argument('--index', action='store_true', help='dump the component index (RECOMB_DP_INDEX)')
    ap.add_argument('--dump', action='store_true', help='dump DP hashes (RECOMB_DP_DUMP)')
    ap.add_argument('--full', type=int, default=0, help='full DP tables for the first n (read, component) pairs')
    ap.add_argument('--time', action='store_true', help='per-read Chaining()/DP timings (RECOMB_DP_TIME)')
    ap.add_argument('--gate', action='store_true', help='RECOMB_DP_GATE=1')
    ap.add_argument('--truedist', action='store_true', help='RECOMB_DP_TRUEDIST=1')
    ap.add_argument('--threads', type=int, default=1)
    ap.add_argument('--extra', nargs='*', default=[], help='extra Minichain arguments')
    args = ap.parse_args()

    out = args.out
    if out.exists():
        raise SystemExit(f'{out} exists; every run needs a new directory')
    out.mkdir(parents=True)
    env = dict(os.environ, OMP_NUM_THREADS=str(args.threads))
    for k in [k for k in env if k.startswith('RECOMB_')]:
        del env[k]
    switches = {}
    if args.seeds:
        switches['RECOMB_DUMP_ANCHORS'] = str(out / 'seeds.tsv')
    if args.index:
        switches['RECOMB_DP_INDEX'] = str(out / 'index.tsv')
    if args.dump:
        switches['RECOMB_DP_DUMP'] = str(out / 'dump.tsv')
    if args.full:
        switches['RECOMB_DP_FULL'] = str(args.full)
    if args.time:
        switches['RECOMB_DP_TIME'] = str(out / 'time.tsv')
    if args.gate:
        switches['RECOMB_DP_GATE'] = '1'
    if args.truedist:
        switches['RECOMB_DP_TRUEDIST'] = '1'
    env.update(switches)
    cmd = [str(BINARY), '-x', 'lr', '-t', str(args.threads), f'-R{args.R}']
    if args.scale is not None:
        cmd.append(f'-s{args.scale}')
    cmd += args.extra + [str(args.graph), str(args.reads)]
    meta = dict(command=cmd, switches=switches, OMP_NUM_THREADS=args.threads,
                binary_sha256=sha256(BINARY), graph_sha256=sha256(args.graph), reads_sha256=sha256(args.reads),
                host=platform.platform(), machine=platform.machine(), load_before=os.getloadavg(),
                started=time.strftime('%Y-%m-%dT%H:%M:%S%z'))
    t0 = time.perf_counter()
    with open(out / 'minichain.gaf', 'wb') as so, open(out / 'minichain.log', 'wb') as se:
        proc = subprocess.run(cmd, env=env, stdout=so, stderr=se)
    meta.update(returncode=proc.returncode, wall_seconds=time.perf_counter() - t0,
                load_after=os.getloadavg(), finished=time.strftime('%Y-%m-%dT%H:%M:%S%z'))
    (out / 'run.json').write_text(json.dumps(meta, indent=1) + '\n')
    if proc.returncode:
        raise SystemExit(f'minichain failed with {proc.returncode}; see {out}/minichain.log')
    print(json.dumps(dict(out=str(out), wall_seconds=round(meta['wall_seconds'], 1))))


if __name__ == '__main__':
    main()
