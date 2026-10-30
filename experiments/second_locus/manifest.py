"""Write run_manifest.json: environment, commands and SHA-256 of sources, inputs and outputs."""
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SOURCES = ['research/current/new/gap_native.cpp', 'research/current/new/gap_native.py',
           'research/current/new/gap_chaining.py',
           'experiments/gap_engine_native/benchmark.py', 'experiments/gap_engine_native/summarize.py',
           'experiments/gap_engine_native/common.py',
           'experiments/gap_solution_review/real_check_mosaic.py',
           'experiments/hprc_mhc_panel/build_panel.py', 'experiments/hprc_mhc_panel/convert.py',
           'experiments/hprc_mhc_panel/bench_panels.py', 'experiments/hprc_mhc_panel/panel_common.py',
           'experiments/hprc_mhc_panel/data/hprc463/graph.gfa.gz', 'experiments/hprc_mhc_panel/data/hprc463/instances.jsonl',
           'experiments/hprc_mhc_panel/data/hprc60/graph.gfa.gz', 'experiments/hprc_mhc_panel/data/hprc60/instances.jsonl',
           'experiments/mhc_pilot/data/sample100.fq', 'experiments/mhc_pilot/final_hybrid.jsonl',
           'experiments/mhc_pilot/minichain-export.patch', 'experiments/mhc_pilot/vendor/minichain/minichain']
INPUTS = ['experiments/hprc_mhc_panel/data/CHM13-464.gfa.gz', 'experiments/hprc_mhc_panel/data/CHM13-464.call.bed.gz']
COMMANDS = [
    'python3.13 -B experiments/second_locus/scan_calls.py',
    'python3.13 -B experiments/second_locus/survey_loci.py',
    'python3.13 -B experiments/second_locus/build_locus.py',
    'python3.13 -B experiments/second_locus/simulate_reads.py',
    'python3.13 -B experiments/second_locus/convert_locus.py kir hprc459 hprc60 hprc120 hprc240 hprc453_leaveout',
    'python3.13 -B experiments/second_locus/convert_locus.py chr20 hprc463 hprc60 hprc120 hprc240 hprc457_leaveout',
    'python3.13 -B experiments/second_locus/convert_locus.py mhc_sim hprc463 hprc60',
    'python3.13 -B experiments/second_locus/verify_locus.py kir --out experiments/second_locus/verify_kir.json   '
    '(attempt 1 stopped by hand after 3 passing instances, host swapping: verify_kir.aborted1.log; attempt 2 stopped '
    'after 25 passing instances of its first level to respect the 2-process limit: verify_kir.log)',
    'python3.13 -B experiments/second_locus/verify_locus.py kir --levels subsets,gamma20000,leaveout '
    '--out experiments/second_locus/verify_kir_levels2to5.json   (stopped after 21 passing instances, same reason: '
    'verify_kir_levels2to5.log)',
    'python3.13 -B experiments/second_locus/verify_locus.py kir --out experiments/second_locus/verify_kir.json '
    '--resume experiments/second_locus/verify_kir.log experiments/second_locus/verify_kir_levels2to5.log   '
    '(log verify_kir_resume.log; reuses the 46 rows printed in full by the stopped runs, computes the rest)',
    'python3.13 -B experiments/second_locus/verify_locus.py chr20 --out experiments/second_locus/verify_chr20.json',
    'python3.13 -B experiments/second_locus/counters.py',
    'python3.13 -B experiments/second_locus/counters.py --summarize experiments/second_locus/counters_summary.json',
    'python3.13 -B experiments/second_locus/chain_switches.py',
    'python3.13 -B experiments/second_locus/realism_check.py',
    'python3.13 -B experiments/second_locus/analyze.py',
    'python3.13 -B experiments/second_locus/manifest.py',
]


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1 << 22), b''):
            h.update(block)
    return h.hexdigest()


def main():
    out = HERE / 'run_manifest.json'
    assert not out.exists(), 'refusing to overwrite'
    compiler = subprocess.run(['clang++', '--version'], capture_output=True, text=True).stdout.splitlines()[0]
    files = sorted(p for p in HERE.rglob('*') if p.is_file() and p.name != 'run_manifest.json'
                   and '__pycache__' not in p.parts)
    manifest = dict(
        python=sys.version, platform=platform.platform(), machine=platform.machine(), compiler=compiler,
        commands=COMMANDS, sources={s: sha(ROOT / s) for s in SOURCES},
        inputs={s: dict(bytes=(ROOT / s).stat().st_size, checksums='MD5 checked against the published value and SHA-256 '
                        'recorded at build time in data/<locus>/panels.json; MD5 also in data/scan/scan.json') for s in INPUTS},
        experiment_files={str(p.relative_to(HERE)): dict(bytes=p.stat().st_size, sha256=sha(p)) for p in files},
        notes=('Write-up: RESULTS.md (the harness blocks REPORT.md for subagents). No timing benchmark was run (shared, loaded host); bench_locus.py is prepared for an '
               'idle machine. Verification and counters ran single-process with python3.13 -B (no bytecode written outside '
               'this folder) against the existing native library in results/native_build (not rebuilt).'))
    out.write_text(json.dumps(manifest, indent=1) + '\n')
    print('wrote', out, len(files), 'files')


if __name__ == '__main__':
    main()
