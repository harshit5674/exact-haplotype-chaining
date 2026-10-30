"""Write run_manifest.json: environment, commands and SHA-256 of sources, inputs and outputs."""
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SOURCES = ['research/current/new/gap_fastpath.py', 'research/current/new/gap_native.cpp',
           'research/current/new/gap_native.py', 'research/current/new/test_gap_native.py',
           'research/current/new/gap_chaining.py', 'research/current/new/test_gap_chaining.py',
           'experiments/gap_solution_review/real_check_mosaic.py', 'experiments/gap_solution_review/independent_check.py',
           'experiments/class_chaining/scaling.py', 'experiments/mhc_pilot/final_hybrid.jsonl',
           'experiments/mhc_pilot/data/graph.json', 'experiments/mhc_pilot/data/instances.jsonl',
           'scripts/project.py']
COMMANDS = [
    'python3.13 experiments/gap_engine_native/stageA_verify.py --out experiments/gap_engine_native/stageA_verification.json '
    '--random 3000 --wide 400 --independent 3000 --real60 40 --real960 6 --gamma20000 6 --workers 4',
    'python3.13 experiments/gap_engine_native/stageA_cost.py --out experiments/gap_engine_native/stageA_cost.jsonl --workers 5',
    'python3.13 experiments/gap_engine_native/verify_native.py --out experiments/gap_engine_native/native_verification.json',
    'python3.13 experiments/gap_engine_native/benchmark.py --out experiments/gap_engine_native/benchmark.jsonl',
    'python3.13 experiments/gap_engine_native/summarize.py experiments/gap_engine_native/benchmark.jsonl '
    '--out experiments/gap_engine_native/benchmark_summary.json',
    'python3.13 scripts/project.py verify --suite gap   (includes research/current/new/test_gap_native.py)',
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    out = HERE / 'run_manifest.json'
    assert not out.exists(), 'refusing to overwrite'
    compiler = subprocess.run(['clang++', '--version'], capture_output=True, text=True).stdout.splitlines()[0]
    files = sorted(p for p in HERE.iterdir() if p.is_file() and p.name != 'run_manifest.json' and p.suffix != '.pyc')
    manifest = dict(
        python=sys.version, platform=platform.platform(), machine=platform.machine(), compiler=compiler,
        compile_flags=['-std=c++20', '-O3', '-shared', '-fPIC'], commands=COMMANDS,
        sources={s: sha(ROOT / s) for s in SOURCES},
        experiment_files={p.name: sha(p) for p in files},
        notes=('Timings are single-threaded, sequential and paired on one Apple M4 host with background '
               'system load; see benchmark.meta.json for load averages. Python timings in Stage A are not '
               'native evidence; Stage A reports operation counts only.'))
    out.write_text(json.dumps(manifest, indent=1) + '\n')
    print('wrote', out)


if __name__ == '__main__':
    main()
