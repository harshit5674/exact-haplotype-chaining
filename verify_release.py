#!/usr/bin/env python3
"""Reduced test suites of the released code (a few minutes). Writes results to a temporary directory."""
import subprocess, sys, tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
NEW = HERE / 'research/current/new'
if not __debug__:
    sys.exit('run with assertions enabled (no -O)')
with tempfile.TemporaryDirectory() as tmp:
    steps = [
        ('specification', ['test_gap_chaining.py', '--cases', '300', '--arithmetic', '200']),
        ('native engines', ['test_gap_native.py', '--cases', '300', '--wide', '30', '--bigkey', '20', '--midgraph', '30']),
        ('guard budget', ['test_gap_budget.py']),
        ('validators', ['test_gap_validation_light.py']),
    ]
    for k, (name, args) in enumerate(steps):
        out = Path(tmp) / f'step{k}.json'
        print(f'== {name}', flush=True)
        subprocess.run([sys.executable, '-B', *args, '--out', str(out)], cwd=NEW, check=True)
    print('all release checks passed')
