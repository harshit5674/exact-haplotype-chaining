# Exact haplotype-aware chaining with origin-normalized classes

Code for the RECOMB 2027 submission "Exact haplotype-aware chaining with origin-normalized classes". It computes CGJ's gap-sensitive haplotype-aware chaining
(the recurrence of Chandra, Gibney and Jain's Algorithm 2) exactly for every (anchor, haplotype) state,
sharing work across haplotypes with origin-normalized records.

Author: Harshit Verma, Yale University.

## Contents

| Path | What |
| --- | --- |
| `research/current/new/gap_chaining.py` | Python specification of the origin-class algorithm, with the direct recurrence oracle |
| `research/current/new/gap_fastpath.py` | Specification of the certified scalar shortcuts |
| `research/current/new/gap_native.cpp`, `gap_native.py` | Native engines (per-haplotype, class, guarded) and their Python adapter; the adapter builds the library with clang++ on first use |
| `research/current/new/gap_validation*.{cpp,py}` | Independent graph-major validators (validation only, never timed) |
| `research/current/new/test_*.py` | Differential and regression tests |
| `experiments/proof_review/` | Exhaustive and adversarial small-instance tests: every solver against an independent brute-force enumeration of chains |
| `experiments/*/` | Scripts of every experiment in the paper (inputs are downloaded, see below) |
| `paper/scripts/make_paper_data.py` | Generates every number in the paper from saved results |

## Requirements

Python 3.10 or later (no third-party packages), a C++20 compiler (clang++ by default; set `CXX` to change),
and for the experiments: GCC (to build Minichain), Rust/cargo (gbz-base), and the public inputs listed below.
Keep Python assertions enabled (no `-O`).

## Quick check

```sh
python3 verify_release.py
```

runs reduced versions of the specification, native-engine, budget and validator test suites (a few minutes).
Run it once before the experiments: it also builds the native library that the benchmark scripts load.
The full-size suites of the paper's runs are

```sh
python3 -B research/current/new/test_gap_chaining.py --cases 5000 --arithmetic 2000 --out spec.json
python3 -B research/current/new/test_gap_native.py --cases 2600 --wide 200 --bigkey 100 --out native.json
```

and `experiments/proof_review/` holds the exhaustive and adversarial small-instance tests; each script's
docstring gives its usage.

## Third-party software (fetched, not included)

The scripts expect four builds, each from a fresh checkout of a pinned commit:

| Build directory | Source | Patches, in order |
| --- | --- | --- |
| `experiments/mhc_pilot/vendor/minichain` | Minichain v1.3, commit `d21e6941572c9152ec20668be47c63c672bf8177` (https://github.com/at-cg/minichain) | `experiments/mhc_pilot/minichain-export.patch` (seed export and macOS integer types) |
| `experiments/minichain_matched/vendor/minichain-dp` | the same | the export patch, then `experiments/minichain_matched/minichain-dpdump.patch` (DP-table dumps; apply with `-p3`) |
| `experiments/minichain_integrated/vendor/minichain-engine` | the same | the two above, then `experiments/minichain_integrated/minichain-engine.patch` (the drop-in chainer) |
| `experiments/baselevel_feasibility/vendor/gbz-base` | gbz-base, commit `2b48da95bfca5b1fef717a4b6006233b7535a385` (https://github.com/jltsiren/gbz-base) | none |

```sh
git clone https://github.com/at-cg/minichain experiments/mhc_pilot/vendor/minichain
git -C experiments/mhc_pilot/vendor/minichain checkout d21e6941572c9152ec20668be47c63c672bf8177
patch -d experiments/mhc_pilot/vendor/minichain -p1 < experiments/mhc_pilot/minichain-export.patch
mkdir -p experiments/minichain_matched/vendor experiments/minichain_integrated/vendor
cp -R experiments/mhc_pilot/vendor/minichain experiments/minichain_matched/vendor/minichain-dp
patch -d experiments/minichain_matched/vendor/minichain-dp -p3 < experiments/minichain_matched/minichain-dpdump.patch
cp -R experiments/minichain_matched/vendor/minichain-dp experiments/minichain_integrated/vendor/minichain-engine
patch -d experiments/minichain_integrated/vendor/minichain-engine -p1 < experiments/minichain_integrated/minichain-engine.patch
```

Build each Minichain copy with GCC and Minichain's own flags, for example
`make -C experiments/mhc_pilot/vendor/minichain -j4 CXX=g++-15`. On recent macOS, GCC 15.1 may need
`MACOSX_DEPLOYMENT_TARGET` set to the running macOS version. Build gbz-base with `cargo build --release` in its
directory.

## Inputs (public)

- Minichain's MHC benchmark, Zenodo record 10665350: `MHC-CHM13.0.gfa.gz`, `MHC_CHM13_PacBio_filt.fq.gz`.
- Li and Hickey's 464-haplotype minigraph graph and calls, Zenodo record 16728828 (v4.1): `CHM13-464.gfa.gz`,
  `CHM13-464.call.bed.gz`.
- HPRC release 2 Minigraph-Cactus graph v2.1, GBZ-base database `hprc-v2.1-mc-chm13.gbz.db`
  (s3://human-pangenomics/pangenomes/freeze/release2/minigraph-cactus/v2.1/hprc-v2.1-mc-chm13/).

## Reproducing the paper

Each experiment directory's script documents its command line in its module docstring. The paper's
Supplement D lists the commands in order. Every result script refuses to overwrite an existing output.
The MHC anchors come from the seed export, run between the two preparation steps on a new, empty dump path
(the exporter appends):

```sh
python3 experiments/mhc_pilot/prepare.py sample
RECOMB_DUMP_ANCHORS="$PWD/experiments/mhc_pilot/data/anchors.tsv" OMP_NUM_THREADS=1 \
  experiments/mhc_pilot/vendor/minichain/minichain -x lr -t 1 \
  experiments/mhc_pilot/data/MHC-CHM13.0.gfa.gz experiments/mhc_pilot/data/sample100.fq
python3 experiments/mhc_pilot/prepare.py convert
```

The identity with Minichain's DP tables (Section 6.2 of the paper) uses the DP-dump build. On Minichain's MHC
benchmark, then on each HPRC MHC panel:

```sh
M=experiments/minichain_matched
python3 $M/run_minichain.py --out $M/data/all_seeds --seeds --index
python3 $M/run_minichain.py --out $M/data/dp_R10000 --dump --index --R 10000
python3 $M/run_minichain.py --out $M/data/dp_R10000_truedist --dump --index --truedist --R 10000
python3 $M/run_minichain.py --out $M/data/dp_R10000_gate --dump --gate --R 10000
python3 $M/compare.py --run $M/data/dp_R10000 --seeds $M/data/all_seeds --distance minichain \
  --out $M/results/identity_mhc2104_R10000_minichain.jsonl
python3 $M/compare.py --run $M/data/dp_R10000_truedist --seeds $M/data/all_seeds --distance true \
  --out $M/results/identity_mhc2104_R10000_truedist_vs_alg2.jsonl
python3 $M/variants.py --base $M/data/dp_R10000 --variant $M/data/dp_R10000_truedist --out $M/results/variant_truedist.json
python3 $M/variants.py --base $M/data/dp_R10000 --variant $M/data/dp_R10000_gate --out $M/results/variant_gate.json
for P in hprc60 hprc120 hprc240 hprc463; do
  GRAPH=experiments/hprc_mhc_panel/data/$P/graph.gfa.gz; READS=experiments/mhc_pilot/data/sample100.fq
  python3 $M/run_minichain.py --out $M/data/seeds_$P --seeds --index --graph $GRAPH --reads $READS
  python3 $M/run_minichain.py --out $M/data/dp_R10000_$P --dump --index --R 10000 --graph $GRAPH --reads $READS
  python3 $M/compare.py --run $M/data/dp_R10000_$P --seeds $M/data/seeds_$P --distance minichain \
    --out $M/results/identity_R10000_minichain_$P.jsonl
done
```

`paper/scripts/make_paper_data.py` shows how every number in the paper is derived from the saved results.
Those result files (about 47 MB) are not part of this repository.

## License

See `LICENSE`.
