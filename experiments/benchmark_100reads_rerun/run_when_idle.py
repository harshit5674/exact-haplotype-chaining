"""Rerun the paper's 100-read guarded benchmark on an idle host (2026-10-20).

The first headline timings ran on a shared laptop
at one-minute load averages of 2-7. This wrapper waits until the host is idle, then runs the unchanged
experiments/benchmark_100reads/run.py timing phase against its existing validation (whose source hashes
must still match), and its unchanged summarize.py. It records the host load and swap before every panel
group and after the run. Nothing else changes: same inputs, engines, repetitions protocol and statistics.

Idle means: one-minute load average below --max-load, and no other project timing job running, for --settle
consecutive minutes. The benchmark itself
is single-threaded. Run it only when no other project benchmark is running.

  python3.13 -B experiments/benchmark_100reads_rerun/run_when_idle.py --reps 5 \
      --out experiments/benchmark_100reads_rerun/timings
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
B100 = ROOT / 'experiments/benchmark_100reads'


def host():
    swap = subprocess.run(['sysctl', '-n', 'vm.swapusage'], capture_output=True, text=True).stdout.strip()
    return dict(time=datetime.now(timezone.utc).isoformat(), load=os.getloadavg(), swap=swap)


OTHER_TIMING = ('benchmark_100reads/run.py', 'minichain_integrated/run.py --phase timing', 'minichain-engine/minichain',
                'bench_windows.py --phase benchmark')


def other_timing_running() -> bool:
    out = subprocess.run(['ps', '-Ao', 'args'], capture_output=True, text=True).stdout
    return any(pat in line for line in out.splitlines() for pat in OTHER_TIMING)


def wait_idle(max_load: float, settle: int, log: list):
    quiet = 0
    while quiet < settle:
        quiet = quiet + 1 if (os.getloadavg()[0] < max_load and not other_timing_running()) else 0
        if quiet < settle:
            time.sleep(60)
    log.append(dict(event='idle', **host()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--reps', type=int, default=5)
    ap.add_argument('--max-load', type=float, default=1.5)
    ap.add_argument('--settle', type=int, default=5, help='consecutive idle minutes before starting')
    args = ap.parse_args()
    assert not args.out.exists(), 'refuse to reuse an output directory'
    log = [dict(event='waiting', **host())]
    wait_idle(args.max_load, args.settle, log)
    cmd = ['python3.13', '-B', str(B100 / 'run.py'), '--phase', 'benchmark', '--reps', str(args.reps),
           '--validation', str(B100 / 'validation'), '--out', str(args.out)]
    log.append(dict(event='start', command=cmd, **host()))
    rc = subprocess.run(cmd, cwd=ROOT).returncode
    log.append(dict(event='finished', returncode=rc, **host()))
    if rc == 0:
        summary = args.out.parent / (args.out.name + '-summary.json')
        cmd2 = ['python3.13', '-B', str(B100 / 'summarize.py'), '--run', str(args.out), '--validation',
                str(B100 / 'validation'), '--out', str(summary)]
        log.append(dict(event='summarize', command=cmd2, returncode=subprocess.run(cmd2, cwd=ROOT).returncode))
    (args.out.parent / (args.out.name + '-host.json')).write_text(json.dumps(log, indent=1) + '\n')
    raise SystemExit(rc)


if __name__ == '__main__':
    main()
