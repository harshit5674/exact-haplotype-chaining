#!/bin/sh
# Matched timing campaign. Run from the repository root on an otherwise idle machine; every step
# is sequential and single-threaded. Requires the seed-export and DP-dump runs (see README.md):
#   data/all_seeds, data/dp_R10000, data/<panel>_seeds, data/<panel>_dp_R10000.
# Writes data/time_* (Minichain's own timers) and results/bench_*, results/timing_*.
set -eu
E=experiments/minichain_matched
PY=python3.13
REPS="1 2 3"

# 1. Authors' MHC graph (60 haplotypes), all 2,104 reads, R = 10000.
for r in $REPS; do
  $PY $E/run_minichain.py --out $E/data/time_mhc_R10000_rep$r --time --R 10000
done
$PY $E/bench_matched.py --seeds $E/data/all_seeds --dump $E/data/dp_R10000 --R 10000 --repeat 3 \
  --out $E/results/bench_mhc2104_R10000.jsonl
$PY $E/summarize_timing.py --minichain $E/data/time_mhc_R10000_rep1 $E/data/time_mhc_R10000_rep2 \
  $E/data/time_mhc_R10000_rep3 --ours $E/results/bench_mhc2104_R10000.jsonl --label mhc2104_R10000 \
  --out $E/results/timing_mhc2104_R10000.json

# 2. Real HPRC release 2 panels, the 100 sampled reads, R = 10000.
for p in hprc60 hprc120 hprc240 hprc463; do
  for r in $REPS; do
    $PY $E/run_minichain.py --out $E/data/time_${p}_R10000_rep$r --time --R 10000 \
      --graph experiments/hprc_mhc_panel/data/$p/graph.gfa.gz --reads experiments/mhc_pilot/data/sample100.fq
  done
  $PY $E/bench_matched.py --seeds $E/data/${p}_seeds --dump $E/data/${p}_dp_R10000 --R 10000 --repeat 3 \
    --out $E/results/bench_${p}_R10000.jsonl
  $PY $E/summarize_timing.py --minichain $E/data/time_${p}_R10000_rep1 $E/data/time_${p}_R10000_rep2 \
    $E/data/time_${p}_R10000_rep3 --ours $E/results/bench_${p}_R10000.jsonl --label ${p}_R10000 \
    --out $E/results/timing_${p}_R10000.json
done

# 3. Context: the pinned binary end to end (no instrumentation), all 2,104 reads.
/usr/bin/time -l env OMP_NUM_THREADS=1 experiments/mhc_pilot/vendor/minichain/minichain -x lr -t 1 -R10000 \
  experiments/mhc_pilot/data/MHC-CHM13.0.gfa.gz experiments/mhc_pilot/data/MHC_CHM13_PacBio_filt.fq.gz \
  > $E/data/pinned_end_to_end_R10000.gaf 2> $E/data/pinned_end_to_end_R10000.log
