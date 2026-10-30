"""Generate this experiment's report from its saved summary, never the paper."""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main():
    s = json.loads((HERE/'summary.json').read_text())
    assert s['status'] == 'passed'
    out = HERE/'REPORT.md'
    assert not out.exists(), 'Preserve prior reports'
    panels = sorted(s['panels'],key=lambda p: (['mhc','kir','chr20'].index(p['locus']),p['H']))
    big = [next(p for p in panels if p['locus']==l and p['H']==h)
           for l,h in [('mhc',463),('kir',459),('chr20',463)]]
    small = [p for p in panels if p['H']==60]
    fallback = sum(p['fallback_components'] for p in panels)
    kir_total_ci = big[1]['kernel']['auto']['total_ci95']
    lines = [
        '# All-100-read benchmark — 2026-10-17', '',
        '**Completed:** all 100 saved reads, both orientations, on all 14 existing panels. '
        'All correctness checks and all five paired timing repetitions passed. '
        'The manuscript and production solvers were not changed.', '',
        '## Main result', '',
        'At the largest primary panels, repaired auto mode has median **read-level native-kernel** '
        'speedups of '+', '.join(f"**{p['kernel']['auto']['median_speedup']:.2f}× {p['locus'].upper()}**" for p in big)+
        ' over the unchanged native per-haplotype baseline. Total-readset speedups are '+
        ', '.join(f"{p['kernel']['auto']['total_speedup']:.2f}×" for p in big)+' in the same order. '
        'The expanded sample supports the large-panel benefit, but KIR has a modest aggregate gain '
        'and some slower reads.', '',
        f"KIR's total-readset speedup has a conditional bootstrap interval of "
        f"[{kir_total_ci[0]:.2f}, {kir_total_ci[1]:.2f}]×. "
        'It includes parity, so these data do not establish a reliable aggregate KIR gain.', '',
        f"At H=60, median guarded speedups range from {min(p['kernel']['auto']['median_speedup'] for p in small):.2f}× "
        f"to {max(p['kernel']['auto']['median_speedup'] for p in small):.2f}×. "
        'The solver remains slower at the median on all three small panels. '
        'No uniform practical speedup is established.', '',
        '## Correctness and coverage', '',
        f"- **{s['components']:,} panel/orientation components; {s['states']:,} valid states**, counted once per input. "
        'These reuse the same 100-read sets across panel sizes and are not independent biological replicates.',
        '- All `C`, `B`, `R` values and every traceback passed the independent graph-major recurrence and '
        'transition/score/acyclicity checks for per-haplotype, class, repaired auto and expanded compact outputs.',
        '- All 560 earlier primary components reproduced their saved exact value signatures; source input hashes agree.',
        '- 224 overlapping bounded-anchor checks, drawn from the additional 80 reads, passed the direct Python '
        'recurrence and pointer checker at four penalties. The new native verifier was separately checked on '
        '301 in-domain fixed/random inputs and rejected all 1,160 injected corruptions.',
        '- Every warmup and timed output reproduced its mode-specific validated signature, including traceback '
        'or exact compact records. All auto calls had a legal fallback and respected the work budget.',
        f'- Pre-timing validation records **{fallback} fallback components**, all in KIR. Fallback counts below '
        'are from that validation pass; the timed rows retain these diagnostics.',
        '- The existing full gap suite, smoke suite and preservation audit passed. Production algorithms were unchanged.', '',
        'Evidence: [validation summary](validation_summary.json), [validation manifest](validation/manifest.json), '
        '[timing manifest](timings/manifest.json), [complete numerical summary](summary.json), '
        '[independent verifier design](../baselevel_mhc_pilot/VALIDATOR.md).', '',
        '## Compatible expanded-output performance', '',
        'Speedup >1 means faster. Each read combines both orientations within a repetition, then uses the '
        'median of five repetitions. The headline is the median of 100 paired read ratios. '
        'Intervals use 2,000 paired whole-read bootstrap resamples; they are conditional on this fixed panel '
        'and sample. “Total” is the ratio of summed per-read median times. “Packed” also includes adapter '
        'output allocation/copy work. Win counts are out of 100 reads. LO denotes an existing leaveout panel.', '',
        '| Locus | H | Guarded kernel median [95% CI] | Total speedup | Packed median | Read wins | Fallback components / 200 |',
        '| --- | ---: | ---: | ---: | ---: | ---: | ---: |',
    ]
    for p in panels:
        a=p['kernel']['auto'];lo,hi=a['median_ci95'];label=p['locus'].upper()+(' LO' if 'leaveout' in p['panel'] else '')
        lines.append(f"| {label} | {p['H']} | {a['median_speedup']:.2f}× [{lo:.2f}, {hi:.2f}] | "
                     f"{a['total_speedup']:.2f}× | {p['packed']['auto']['median_speedup']:.2f}× | {a['wins']} | {p['fallback_components']} |")
    lines += ['', 'The largest primary panels have the following tails and subgroup results. '
              'The previous 20 and additional 80 are workload subsets, not independent external validation cohorts.', '',
              '| Locus | Slowest read: guarded / baseline time | Previous 20 median speedup | Additional 80 median speedup | Component median speedup |',
              '| --- | ---: | ---: | ---: | ---: |']
    for p in big:
        lines.append(f"| {p['locus'].upper()} | {p['kernel']['auto']['max_slowdown']:.2f}× | "
                     f"{p['subgroups']['previous_20']['auto']['median_speedup']:.2f}× | "
                     f"{p['subgroups']['additional_80']['auto']['median_speedup']:.2f}× | {p['component_median_auto_speedup']:.2f}× |")
    lines += ['', 'The earlier repaired-guard headline used medians over individual orientation components. '
              'That statistic differs from the read-level statistic here; differences between those headlines '
              'should not all be attributed to adding 80 reads. Both measures are saved.', '',
              '## Other measured output modes', '',
              'These are native-kernel median read speedups over the same expanded per-haplotype baseline. '
              '**Compact mode is unguarded and omits expansion**, so its column compares different output representations. '
              'Unguarded class mode returns the same expanded output but has no runtime safeguard.', '',
              '| Locus | H | Class expanded | Class compact |', '| --- | ---: | ---: | ---: |']
    for p in panels:
        label=p['locus'].upper()+(' LO' if 'leaveout' in p['panel'] else '')
        lines.append(f"| {label} | {p['H']} | {p['kernel']['class']['median_speedup']:.2f}× | {p['kernel']['compact']['median_speedup']:.2f}× |")
    lines += ['', '## Timing boundaries and limits', '',
              'Prepared graphs and packed anchor batches are shared between paired modes. Native timings include '
              'read setup, initialization, sweep and export. Packed-call timings additionally include output allocations '
              'and adapter work. Neither includes graph preparation, anchor packing, seed generation, checking, or '
              'compact expansion. Preparation/packing costs, raw repetitions, source/input hashes, compiler/library '
              'versions and workstation load are retained in the run files. Pair caches start cold; mode and component '
              'orders are shuffled. All agent-launched validation/test jobs finished before timing. This is one shared '
              'workstation, not a dedicated benchmark host.', '',
              'MHC uses real CHM13 HiFi reads; KIR and chr20 use the existing simulated HiFi-like reads. '
              'All these timed panels remain SV-level graphs from the existing construction process. The '
              '[separate base-level MICB pilot](../baselevel_mhc_pilot/REPORT.md) establishes a small release-1 '
              'input/correctness feasibility result, not a release-2 full-MHC speedup. The expanded benchmark therefore '
              'strengthens the sample-size evidence while leaving broad graph generalization and end-to-end mapping '
              'performance open. It does not certify the theoretical proofs or RECOMB acceptance.', '',
              'Protocol and commands: [README.md](README.md). The report is generated by `write_report.py` from '
              '`summary.json`; no manuscript generator or paper artifact was run or edited.', '']
    out.write_text('\n'.join(lines))


if __name__ == '__main__': main()
