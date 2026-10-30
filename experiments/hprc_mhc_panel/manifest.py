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
           'research/current/new/test_gap_native.py', 'research/current/new/gap_chaining.py',
           'experiments/gap_engine_native/benchmark.py', 'experiments/gap_engine_native/summarize.py',
           'experiments/gap_solution_review/real_check_mosaic.py', 'experiments/mhc_pilot/data/MHC-CHM13.0.gfa.gz',
           'experiments/mhc_pilot/data/sample100.fq', 'experiments/mhc_pilot/final_hybrid.jsonl',
           'experiments/mhc_pilot/minichain-export.patch', 'experiments/mhc_pilot/vendor/minichain/minichain']
COMMANDS = [
    'curl -o experiments/hprc_mhc_panel/data/CHM13-464.gfa.gz https://zenodo.org/records/16728828/files/CHM13-464.gfa.gz',
    'curl -o experiments/hprc_mhc_panel/data/CHM13-464.call.bed.gz https://zenodo.org/records/16728828/files/CHM13-464.call.bed.gz',
    'curl -o experiments/hprc_mhc_panel/data/CHM13-464.call.vcf.gz https://zenodo.org/records/16728828/files/CHM13-464.call.vcf.gz (inspected, not used)',
    'curl -o experiments/hprc_mhc_panel/data/hprc-v2.1-mc-grch38.sv.gfa.gz https://s3-us-west-2.amazonaws.com/human-pangenomics/'
    'pangenomes/freeze/release2/minigraph-cactus/v2.1/hprc-v2.1-mc-grch38/hprc-v2.1-mc-grch38.sv.gfa.gz (checked: no walks; not used)',
    'python3.13 experiments/hprc_mhc_panel/build_panel.py',
    'python3.13 experiments/hprc_mhc_panel/build_leaveout.py',
    'python3.13 experiments/hprc_mhc_panel/convert.py hprc60 hprc120 hprc240 hprc463 && ... convert.py hprc462_noCHM13',
    'python3.13 experiments/hprc_mhc_panel/bench_panels.py --out experiments/hprc_mhc_panel/benchmark.jsonl   (always-dense baseline, v1)',
    'python3.13 experiments/hprc_mhc_panel/bench_panels.py --panels hprc462_noCHM13 --out experiments/hprc_mhc_panel/benchmark_leaveout.jsonl   (v1)',
    'python3.13 experiments/hprc_mhc_panel/verify_panels.py --out experiments/hprc_mhc_panel/panel_verification.json   (v1 binary)',
    'python3.13 experiments/hprc_mhc_panel/calibrate_dense.py --out experiments/hprc_mhc_panel/calibration_dense.json',
    'python3.13 experiments/hprc_mhc_panel/bench_panels.py --panels hprc60,hprc120,hprc240,hprc463,hprc462_noCHM13 '
    '--out experiments/hprc_mhc_panel/benchmark_v2.jsonl --seed 9191   (strengthened baseline, primary)',
    'python3.13 experiments/hprc_mhc_panel/verify_panels.py --full 10 --subset 3 --gamma20000 3 --leaveout 3 '
    '--out experiments/hprc_mhc_panel/panel_verification_v2.json',
    'python3.13 experiments/gap_engine_native/summarize.py <rows> --out <summary>   (for each benchmark file)',
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
        compile_flags=['-std=c++20', '-O3', '-shared', '-fPIC'], commands=COMMANDS,
        sources={s: sha(ROOT / s) for s in SOURCES},
        experiment_files={str(p.relative_to(HERE)): dict(bytes=p.stat().st_size, sha256=sha(p)) for p in files},
        notes=('Benchmarks are single-threaded, sequential and paired on one Apple M4 host with background system '
               'load (see the *.meta.json files). benchmark.jsonl and benchmark_leaveout.jsonl used the always-dense '
               'per-haplotype baseline (library gap-19525d3b6d895e9a); benchmark_v2.jsonl used the strengthened '
               'baseline (hybrid lanes, grouped remote queries; see calibration_dense.json).'))
    out.write_text(json.dumps(manifest, indent=1) + '\n')
    print('wrote', out, len(files), 'files')


if __name__ == '__main__':
    main()
