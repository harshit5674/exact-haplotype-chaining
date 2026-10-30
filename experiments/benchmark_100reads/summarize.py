"""Read-level summaries from paired timings; never modifies the raw runs."""
import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import random
import statistics as st

MODES = ('class', 'auto', 'compact')


def quantile(values, fraction):
    ordered = sorted(values)
    return ordered[round(fraction*(len(ordered)-1))]


def summarize_reads(reads, boundary, bootstrap_seed):
    answer = {}
    for mode in MODES:
        pairs = [(r[boundary]['per_haplotype'], r[boundary][mode]) for r in reads]
        ratios = [a/b for a,b in pairs]
        rng = random.Random(bootstrap_seed)
        medians, totals = [], []
        for _ in range(2000):
            selected = [pairs[rng.randrange(len(pairs))] for _ in pairs]
            medians.append(st.median(a/b for a,b in selected))
            totals.append(sum(a for a,b in selected)/sum(b for a,b in selected))
        answer[mode] = dict(
            reads=len(pairs), median_speedup=st.median(ratios),
            median_ci95=[quantile(medians,.025), quantile(medians,.975)],
            total_baseline_seconds=sum(a for a,b in pairs),
            total_method_seconds=sum(b for a,b in pairs),
            total_speedup=sum(a for a,b in pairs)/sum(b for a,b in pairs),
            total_ci95=[quantile(totals,.025), quantile(totals,.975)],
            p10_speedup=quantile(ratios,.1), p90_speedup=quantile(ratios,.9),
            min_speedup=min(ratios), max_speedup=max(ratios),
            max_slowdown=1/min(ratios), wins=sum(x>1 for x in ratios))
    return answer


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run',type=Path,required=True)
    ap.add_argument('--validation',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args = ap.parse_args()
    assert not args.out.exists()
    manifest = json.loads((args.run/'manifest.json').read_text())
    valid = json.loads((args.validation/'manifest.json').read_text())
    assert manifest['status'] == valid['status'] == 'passed'
    assert manifest['source_hashes'] == valid['source_hashes']
    panels = []
    for panel in manifest['panels']:
        key = panel['locus']+'-'+panel['panel']
        rows = [json.loads(line) for line in (args.run/(key+'.jsonl')).read_text().splitlines()]
        assert len(rows) == 200 and len({r['id'] for r in rows}) == 200
        assert sum(r['K'] for r in rows) == panel['states']
        groups = defaultdict(list)
        for row in rows: groups[row['read']].append(row)
        assert len(groups) == 100 and all(len(v)==2 for v in groups.values())
        reads = []
        for name, components in groups.items():
            assert len({c['prior_primary'] for c in components}) == 1
            r = dict(read=name, prior_primary=components[0]['prior_primary'],
                     fallback=any(c['guard']['fallback'] for c in components),
                     N=sum(c['N'] for c in components), K=sum(c['K'] for c in components))
            for boundary in ('kernel_seconds', 'packed_seconds'):
                # Sum both orientations within each repetition, then take the
                # median over repetitions. Reads, not orientations/repetitions,
                # are the units resampled in the conditional bootstrap.
                r[boundary] = {mode: st.median(
                    sum(c[boundary][mode][j] for c in components)
                    for j in range(manifest['repetitions']))
                    for mode in ('per_haplotype', *MODES)}
            reads.append(r)
        seed = 20261017+sum(key.encode())
        result = dict(locus=panel['locus'], panel=panel['panel'], H=panel['H'],
                      reads=100, components=len(rows), states=sum(r['K'] for r in reads),
                      anchors=sum(r['N'] for r in reads),
                      fallback_components=sum(c['guard']['fallback'] for c in rows),
                      fallback_reads=sum(r['fallback'] for r in reads),
                      kernel=summarize_reads(reads,'kernel_seconds',seed),
                      packed=summarize_reads(reads,'packed_seconds',seed),
                      subgroups={})
        for label, selected in [('previous_20',[r for r in reads if r['prior_primary']]),
                                 ('additional_80',[r for r in reads if not r['prior_primary']])]:
            assert len(selected) == (20 if label=='previous_20' else 80)
            result['subgroups'][label] = summarize_reads(selected,'kernel_seconds',seed)
        result['component_median_auto_speedup'] = st.median(
            st.median(c['kernel_seconds']['per_haplotype'])/st.median(c['kernel_seconds']['auto']) for c in rows)
        result['median_records_per_anchor'] = st.median(c['compact_records']/max(1,c['N']) for c in rows)
        result['median_state_record_ratio'] = st.median(c['K']/max(1,c['compact_records']) for c in rows)
        result['read_rows'] = reads
        panels.append(result)
    output = dict(status='passed', panels=panels,
                  states=sum(p['states'] for p in panels), components=sum(p['components'] for p in panels),
                  source_manifest=str(args.run/'manifest.json'),
                  validation_manifest=str(args.validation/'manifest.json'),
                  analysis_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in sorted(args.run.glob('*.json*'))},
                  methodology='5 paired shuffled repetitions, cold pair caches; median of the per-repetition sum of both orientation times per read; paired whole-read bootstrap 2000 resamples; 95% percentile intervals conditional on each fixed panel/read sample, not independent biological replications across panel sizes',
                  source_unchanged=manifest['source_hashes']==valid['source_hashes'])
    args.out.write_text(json.dumps(output,indent=2)+'\n')
    for p in panels:
        a=p['kernel']['auto']; b=p['packed']['auto']
        print(p['locus'],p['H'],p['panel'],f"kernel {a['median_speedup']:.3f} {a['median_ci95']}, total {a['total_speedup']:.3f}, packed {b['median_speedup']:.3f}, wins {a['wins']}, fallbacks {p['fallback_components']}")


if __name__ == '__main__': main()
