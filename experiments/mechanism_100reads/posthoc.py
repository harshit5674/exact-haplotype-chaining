"""Post hoc diagnostic, NOT in PLAN.md; written after results.json was seen.

Question: the same 120 components (earlier 20 reads, largest panels) have identical record
counts in the old 40-instance benchmark and the new 100-read benchmark (S4a), yet their pooled
Spearman correlation of records per anchor with compact speedup is -0.78 with the old timings
and -0.86 with the new ones (S4b). How well do the old and new speedups of the same components
agree, and which mode's time changed? Reads saved files only; runs no timing.

  python3.13 -B experiments/mechanism_100reads/posthoc.py --out experiments/mechanism_100reads/posthoc.json
"""
import argparse
import hashlib
import json
from pathlib import Path
import statistics as st
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import analysis as an  # the pre-specified analysis's loaders and rank correlation

if not __debug__:
    sys.exit('posthoc.py: run with assertions enabled')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    assert not args.out.exists()
    out = dict(status='passed', note='post hoc diagnostic, not pre-specified', loci={})
    for locus in an.LOCI:
        path, H = an.OLD[locus]
        old = {r['id']: r for r in an.load_jsonl(path) if r['H'] == H and r['gamma'] == an.OLD_GAMMA}
        new = {r['id']: r for r in an.load_jsonl(an.BENCH/'timings'/f'{an.LARGEST[locus]}.jsonl')}
        ids = sorted(old)
        assert len(ids) == 40 and all(new[i]['prior_primary'] for i in ids)
        assert all(old[i]['stats']['class_compact']['output_records'] == new[i]['compact_records'] for i in ids)

        def med(r, m):
            return st.median(r['kernel_seconds'][m])
        old_speed = [old[i]['kernel_median']['per_haplotype']/old[i]['kernel_median']['class_compact'] for i in ids]
        new_speed = [med(new[i], 'per_haplotype')/med(new[i], 'compact') for i in ids]
        out['loci'][locus] = dict(
            components=len(ids), H=H,
            spearman_old_vs_new_compact_speedup=an.spearman(old_speed, new_speed),
            median_new_over_old_compact_speedup=st.median(b/a for a, b in zip(old_speed, new_speed)),
            median_new_over_old_per_haplotype_time=st.median(med(new[i], 'per_haplotype')/old[i]['kernel_median']['per_haplotype'] for i in ids),
            median_new_over_old_compact_time=st.median(med(new[i], 'compact')/old[i]['kernel_median']['class_compact'] for i in ids),
            spread_new_over_old_compact_speedup=an.spread([b/a for a, b in zip(old_speed, new_speed)]))
    out['inputs_sha256'] = dict(sorted(an.inputs.items()))
    out['posthoc_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.out.write_text(json.dumps(out, indent=1) + '\n')
    for locus, v in out['loci'].items():
        print(locus, {k: (round(x, 3) if isinstance(x, float) else x) for k, x in v.items() if k != 'spread_new_over_old_compact_speedup'},
              {k: round(x, 3) for k, x in v['spread_new_over_old_compact_speedup'].items()})


if __name__ == '__main__':
    main()
