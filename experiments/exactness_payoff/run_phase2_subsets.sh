#!/bin/sh
# Phase 2 on random 16-haplotype subsets (no CHM13) of both 60-haplotype panels, sequentially (one worker, run alongside the two-worker Phase 3/4 jobs).
set -e
cd "$(dirname "$0")/../.."
for base in hprc60 mhc60; do
  for r in 0 1 2; do
    ds="${base}_s16_r${r}"
    python3.13 experiments/exactness_payoff/run_phase2.py --dataset "$ds" \
      --out "experiments/exactness_payoff/phase2_${ds}.jsonl" --workers 1 \
      > "experiments/exactness_payoff/phase2_${ds}.log" 2>&1
  done
done
