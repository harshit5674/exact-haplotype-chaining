# Pre-specified mechanism analysis of the 100-read benchmark

Written 2026-10-18, before any association between a speedup and a work counter was computed in this session. `analysis.py` records this file's SHA-256 in `results.json`. Amendments, if any, are appended at the end with their date and reason; nothing above them is edited after the first full run.

## Question

1. Within each locus and panel, does a read's guarded speedup over the per-haplotype engine track the class engine's work counters?
2. Across loci at the largest panels, does it track them, and is any pooled association more than the between-locus difference?
3. What are the 59 KIR fallback components, and what do abandoned work plus fallback cost relative to the per-haplotype engine?

Section 6.5 of the paper currently says, from the original 40-instance benchmark, that "within loci the speedup fell as records per anchor grew (Supplementary Figure S1)". Question 1 tests that sentence on the 100-read guarded benchmark. The H=60 stage breakdown (Table S5) is not redone.

## What had been seen before this plan

- Panel-level medians, totals, wins and fallback counts in [benchmark_100reads/REPORT.md](../benchmark_100reads/REPORT.md); the per-read median counters in `paper/generated/numbers.tex`.
- Figure S1's pooled Spearman correlation of −0.78 (120 instances, compact speedup, original benchmark).
- The component-level counters and guard work ratios in the 2026-10-18 advisory entry of `docs/SESSION_LOG.md`.
- The field layout of the timing and validation rows, and the guard fields of one fallback component (the first fallback row of `validation/kir-hprc459.jsonl`), read to learn the field semantics.
- No per-read or per-component association between any speedup and any counter.

## Data

Read-only inputs; `results.json` records the SHA-256 of each.
- `experiments/benchmark_100reads/timings/`: `manifest.json`, and `<panel>.json` and `<panel>.jsonl` for the 14 panels.
- `experiments/benchmark_100reads/validation/`: the same files.
- `experiments/benchmark_100reads/summary.json` and `validation_summary.json`.
- Bridge to the old figure only (S4a): `experiments/hprc_mhc_panel/benchmark_v2.jsonl`, `experiments/second_locus/benchmark_kir.jsonl`, `experiments/second_locus/benchmark_chr20.jsonl`.
- Cross-checks only: `paper/generated/numbers.tex`.

No timing is run. No read is filtered, except where a sensitivity analysis below says so.

### Integrity checks (the run stops on any failure)

1. Both manifests passed, and their source hashes and native library agree.
2. Each panel has 200 components, 100 reads with two each. Timing and validation ids, `N` and `K` agree. The guard, class and compact statistics in timing rows equal the validation rows. Every mode has 5 positive repetitions. Every guarded call had its fallback available and stayed within budget.
3. Per-read times replay `summary.json` `read_rows` to relative 1e-12. Panel guarded, compact and class medians, totals and win counts replay `summary.json`.
4. Per-panel medians of the two predictors reproduce the paper keys `b100.<locus>.h<H>.records` (2 decimals) and `.sliced` (3 decimals).
5. There are 59 fallback components, all KIR: 20, 8, 7, 12 and 12 at H = 60, 120, 240, 453 (leave-out) and 459.
6. S4a reproduces the old pooled correlation, −0.78 to 2 decimals.

## Unit and variables

The unit is the **read**, with its two orientation components `c` summed.

For each mode `m` (per-haplotype `ph`, guarded `auto`, `compact`, `class`), `t_m` is the median over the 5 repetitions `j` of `sum_c kernel_seconds[m][j]`. This is the native-kernel boundary of `summarize.py`.

**Outcomes**
- `y_guard = ln(t_ph / t_auto)` is primary: same expanded output, guarded.
- `y_compact = ln(t_ph / t_compact)` is secondary. It is unguarded and omits expansion, so it compares different output representations.
- `y_class = ln(t_ph / t_class)` is secondary. It has the same expanded output but no guard, so it separates the guard's decision from the class engine's cost.

**Predictors** (paper Table S2 definitions, `bench_panel()` in `paper/scripts/make_paper_data.py`)
- Records per anchor: `x_rec = sum_c compact_records / sum_c N`.
- Bit-sliced decisions per anchor: `x_sl = sum_c class_stats.pairs_sliced / sum_c N`.

Both are deterministic counts from the validation run's unguarded class and compact calls. The guarded call does not save class counters when its attempt is refused, so these unguarded counts of the same deterministic algorithm also describe fallback components.

**Context, exploratory only**
- Anchors `N = sum_c N`.
- States per anchor `s = sum_c K / N`.
- Fallback flag: true if any component fell back.
- Accounting work ratio `w = sum_c auto_work / sum_c auto_estimate`, defined only for reads without a fallback, because a refused attempt's admitted work is truncated.

**Predicted direction:** negative for both predictors. More records or bit-sliced decisions per anchor mean more class work per anchor, so a lower speedup.

## Statistics

- **Correlation.** Spearman's ρ is the Pearson correlation of average ranks, with ties averaged, as `spearman()` in the paper generator. The logarithm does not change ρ; it is used for dispersion and plotting.
- **Interval.** Percentile bootstrap: 2,000 resamples of reads with replacement, ranks recomputed in every resample, and nearest-rank 2.5% and 97.5% quantiles as in `summarize.py`. A resample in which either variable is constant is discarded and counted. Pooled and within-locus summaries resample 100 reads within each locus.
- **p-value.** Two-sided Monte Carlo permutation: 10,000 permutations of the outcome against the predictor, `p = (1 + #{|ρ*| ≥ |ρ|}) / 10,001`. The pooled test permutes all 300 reads. The within-locus summary permutes within each locus.
- **Seeds.** Each statistic uses `random.Random('mechanism_100reads/' + label)`, with a unique label.
- **Multiplicity.** Holm step-down within each family below. Holm is valid under the dependence that shared reads across panel sizes create.

### Primary family: guarded outcome, Holm over 32 tests

- **P1.** In each of the 14 panels, ρ(`y_guard`, `x_rec`) and ρ(`y_guard`, `x_sl`): 28 tests.
- **P2.** Pooled at the largest panels, MHC 463, KIR 459 and chr20 463 (300 reads): one ρ per predictor, 2 tests.
- **P3.** Within-locus summary at those panels: the mean of the three within-locus ρ, one per predictor, 2 tests.

Read sets repeat across panel sizes, so panels are not independent and are never pooled across H.

### Decision rules, fixed now

- **Magnitude.** |ρ| < 0.10 is negligible, 0.10–0.29 weak, 0.30–0.49 moderate and ≥ 0.50 strong.
- **Tracking.** A panel tracks a counter if ρ < 0 and the Holm-adjusted p < 0.05. If ρ > 0 and the adjusted p < 0.05, the association is opposite to the prediction. Otherwise there is no detectable association, which is not evidence of no association.
- **The Section 6.5 sentence** ("within loci the speedup fell as records per anchor grew") is supported at about 460 haplotypes only if ρ(`y_guard`, `x_rec`) tracks at all three largest panels. If it tracks at some, the sentence must name those loci. If at none, the sentence must be removed or replaced.
- **Across loci.** If P2 is Holm-significant but P3 is not, the pooled association reflects the difference between loci, not a per-read mechanism, and will be described that way. If both are significant, both are reported with their sizes.
- **Scope.** All results are conditional on these 100 reads, these SV-level graphs and this host. They are associations, not causal effects.

### Secondary analyses

Each of S1 and S2 is its own Holm family of 32 tests.
- **S1.** P1–P3 with `y_compact`.
- **S2.** P1–P3 with `y_class`.
- **S3.** On the five KIR panels, reads without a fallback component: ρ(`y_guard`, `x_rec`) and ρ(`y_guard`, `x_sl`), with intervals and raw permutation p. Descriptive, not adjusted.
- **S4.** Bridge to the old Figure S1. It uses records per anchor and compact speedup. A component's records per anchor is `compact_records / N`; its compact speedup is `median_j ph / median_j compact`. Components with `N = 0` are excluded and counted. For each item below, the pooled ρ and the three within-locus ρ are reported; this is descriptive, with no tests.
  - (a) Old benchmark, old definitions (`class_compact` `output_records / N`, ratio of `kernel_median`), 120 instances at the largest panels. Old record counts are compared with the new `compact_records` of the same components, and mismatches are counted.
  - (b) New timings, the same 120 components (the earlier 20 reads).
  - (c) New timings, all components at the largest panels.
  - (d) New timings, per read, compact (S1 at the pooled level).
  - (e) New timings, per read, guarded (P2).
- **S5. Outcome reliability.** In each panel, the Spearman correlation between per-read `y_guard` from repetitions {1, 3, 5} and from {2, 4}, and the same for `y_compact`. Also the per-panel median, quartiles and range of each speedup and predictor.

### Exploratory analyses

These give point estimates and bootstrap intervals, with no tests, at the three largest panels and pooled.
- **E1.** ρ(`y_guard`, `N`) and ρ(`y_guard`, `s`); between predictors, ρ(`x_rec`, `x_sl`), ρ(`x_rec`, `s`), ρ(`x_sl`, `s`), ρ(`x_rec`, `N`) and ρ(`x_sl`, `N`).
- **E2.** Partial Spearman ρ(`y_guard`, `x` | `z`) = (ρ_yx − ρ_yz ρ_xz) / sqrt((1 − ρ_yz²)(1 − ρ_xz²)), for `x` in {`x_rec`, `x_sl`} and `z` in {`N`, `s`}.
- **E3.** ρ(`y_guard`, `w`) and ρ(`y_class`, `w`), on reads without a fallback.

## Fallback characterization (KIR, 59 components)

- **F1. One row per fallback component**:
  - panel, H, id, read, earlier-benchmark flag;
  - `N`, `K`, `s`, mask words `W`, `x_rec`, `x_sl`, pair decisions per anchor;
  - estimate `E`, budget `B`, admitted work `A`, `A/B`, refused charge `D`;
  - refusal stage: "output" if `D == N + W·compact_records + K` exactly, which is the output charge of `class_solve`, else "sweep";
  - per-haplotype time (median of 5 repetitions);
  - timed guarded/per-haplotype, class/per-haplotype and compact/per-haplotype ratios, each a ratio of component medians;
  - excess (guarded − per-haplotype) / per-haplotype;
  - the single validation run's (class kernel − class export) / per-haplotype kernel, which predicts the excess when the whole sweep was abandoned.
- **F2. Fallback against other components, per KIR panel.** Medians and quartiles of `x_rec`, `x_sl`, `s`, `N` and class/per-haplotype time. The AUROC of each for separating fallback components: Mann–Whitney U / (n1·n0), ties ½. A 2×2 count of fallback against "class slower than per-haplotype" by timed medians.
- **F3. Read-level impact, per KIR panel.**
  - Number of fallback reads.
  - Median `y_guard` of fallback reads and of other reads.
  - Total-readset guarded speedup: (i) as measured; (ii) without fallback reads; (iii) with each fallback component's guarded time replaced by its own per-haplotype time in the same repetition; (iv) replaced by its class time.
  - (iii) and (iv) are reconstructions from measured times, not measurements of another guard.
- **F4. Recurrence.** How often fallback components and reads recur across the KIR panel sizes.

## Proposed Figure S1 replacement (content fixed now)

- **Plot.** Per-read guarded speedup (log scale, dotted line at 1) against records per anchor (linear), at MHC 463, KIR 459 and chr20 463. There are 300 points, coloured by locus as in Figure 1. KIR reads with a fallback get open markers.
- **Caption.** It reports P2 and the three largest-panel P1 values for `x_rec`.
- **Section 6.5 text.** It follows the decision rules above.

The figure is proposed whatever the results; its text depends on them. Axis ranges are cosmetic and are set after the data are seen. The paper is edited only with the user's approval.

## Procedure

1. Write this file first.
2. Test `analysis.py` before the first full run. Its statistics functions are tested on synthetic data: ranks, Spearman against the tie-free closed form, bootstrap determinism, Holm and the partial formula. The integrity checks and the predictor and time replays are run with the association code disabled (`--checks-only`).
3. One full run writes `results.json`. Any later code fix is recorded in `REPORT.md` with its effect on the results.
4. `REPORT.md` takes its numbers from `results.json` through `write_report.py`.
