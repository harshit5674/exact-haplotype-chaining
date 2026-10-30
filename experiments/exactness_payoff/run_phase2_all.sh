#!/bin/sh
# Phase 2 on all four 60/59-haplotype datasets, sequentially (3 workers each).
set -e
cd "$(dirname "$0")/../.."
for ds in hprc60 hprc60_noCHM13 mhc60 mhc60_noCHM13; do
  python3.13 experiments/exactness_payoff/run_phase2.py --dataset "$ds" \
    --out "experiments/exactness_payoff/phase2_${ds}.jsonl" --workers 3 \
    > "experiments/exactness_payoff/phase2_${ds}.log" 2>&1
done
