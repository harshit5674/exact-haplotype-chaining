#!/bin/zsh
# Start the second idle-host headline rerun only after the base-level window timing has finished (2026-10-21).
# bench_windows.py waits while any run_when_idle.py process exists, and run_when_idle.py waits while
# bench_windows.py runs, so queueing both at once made each wait for the other (about 11 hours on 2026-10-21).
# This script's own command line names neither, so neither guard sees it while it waits.
#
#   zsh experiments/benchmark_100reads_rerun/queue_after_windows.sh > experiments/benchmark_100reads_rerun/run2.log 2>&1
cd "${0:A:h}/../.." || exit 1
while pgrep -f 'bench_windows.py --phase' > /dev/null; do sleep 60; done
exec python3.13 -B experiments/benchmark_100reads_rerun/run_when_idle.py --reps 5 \
    --out experiments/benchmark_100reads_rerun/timings2
