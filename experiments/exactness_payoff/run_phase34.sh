#!/bin/sh
# Phases 3 and 4, sequentially. WORKERS (default 3) sets the pool size of each evaluation step.
set -e
cd "$(dirname "$0")/../.."
E=experiments/exactness_payoff
WK=${WORKERS:-3}
mkdir -p $E/sim

# Phase 3: 300 simulated reads from the 60-haplotype panel (80 single-haplotype controls,
# 110 two-haplotype and 110 three-haplotype mosaics), anchors from pinned Minichain.
python3.13 $E/simulate_reads.py --panel hprc60 --prefix $E/sim/hprc60 --seed 20261006 --counts 80,110,110
python3.13 $E/export_anchors.py --panel hprc60 --prefix $E/sim/hprc60
python3.13 $E/run_phase3.py --prefix $E/sim/hprc60 --panel hprc60 --rhos 2,10,56 --workers $WK \
  --out $E/phase3_hprc60.jsonl > $E/phase3_hprc60.log 2>&1

# Phase 3, informative set: 200 mosaics (100 with two and 100 with three haplotypes) whose own vertex
# walk (vertices >= 17 bp) needs at least one haplotype switch with this panel.
python3.13 $E/simulate_reads.py --panel hprc60 --prefix $E/sim/hprc60_inf --seed 20261008 --counts 0,100,100 \
  --min-walk-switches 1
python3.13 $E/export_anchors.py --panel hprc60 --prefix $E/sim/hprc60_inf
python3.13 $E/run_phase3.py --prefix $E/sim/hprc60_inf --panel hprc60 --rhos 2,10,56 --workers $WK \
  --out $E/phase3_hprc60_inf.jsonl > $E/phase3_hprc60_inf.log 2>&1

# Phase 4: 120 simulated reads from the 463-haplotype panel, chained on the full panel and on
# random 16- and 64-haplotype subsets (three replicates each).
python3.13 $E/simulate_reads.py --panel hprc463 --prefix $E/sim/hprc463 --seed 20261007 --counts 40,40,40
python3.13 $E/export_anchors.py --panel hprc463 --prefix $E/sim/hprc463
python3.13 $E/make_subsets.py --out $E/phase4_subsets.json --seed 4630
python3.13 $E/run_phase3.py --prefix $E/sim/hprc463 --panel hprc463 --rhos 2,56 --workers $WK \
  --out $E/phase4_full.jsonl > $E/phase4_full.log 2>&1
for s in $(python3.13 -c "import json; print(' '.join(x['name'] for x in json.load(open('$E/phase4_subsets.json'))['subsets']))"); do
  haps=$(python3.13 -c "import json; d={x['name']: x for x in json.load(open('$E/phase4_subsets.json'))['subsets']}; print(','.join(map(str, d['$s']['haplotypes'])))")
  python3.13 $E/run_phase3.py --prefix $E/sim/hprc463 --panel hprc463 --rhos 2,56 --workers $WK --subset "$haps" \
    --out $E/phase4_$s.jsonl > $E/phase4_$s.log 2>&1
done
