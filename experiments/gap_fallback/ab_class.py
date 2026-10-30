"""A/B check that the patch did not slow the class modes: one process per implementation.

--impl old imports the pre-patch gap_native saved in prepatch/ (it reuses the cached library);
--impl new imports research/current/new. Runs class compact and class expanded on the 40
primary hprc463 instances, 3 shuffled repetitions, and prints per-instance median kernel times
as JSON. The driver (--drive) runs old, new, new, old and compares per-instance means.
"""
import argparse, json, random, statistics, subprocess, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
OLD = str(HERE / 'prepatch')  # byte-identical pre-patch sources; they build the cached gap-a63a2efa... library


def measure(impl):
    sys.path.insert(0, str(ROOT / 'research/current/new'))
    if impl == 'old':
        sys.path.insert(0, OLD)
    sys.path.insert(0, str(ROOT / 'experiments/hprc_mhc_panel'))
    import gap_native as gn, panel_common as pc
    assert (impl == 'old') == gn.__file__.startswith(OLD), gn.__file__
    G = gn.GapGraph(pc.prepared('hprc463'))
    rng = random.Random(5); out = {}
    for case in pc.primary_cases('hprc463'):
        b = G.batch(pc.anchors(case), 10000)
        t = {'compact': [], 'expanded': []}
        for _ in range(3):
            for m in rng.sample(list(t), 2):
                r = G.solve_compact(b) if m == 'compact' else G.solve_packed(b, mode='class')
                t[m].append(r.stats['kernel_seconds'])
        out[case['id']] = {m: statistics.median(v) for m, v in t.items()}
    print(json.dumps(dict(impl=impl, library=gn.build()._name, times=out)))


def drive():
    runs = []
    for impl in ('old', 'new', 'new', 'old'):
        p = subprocess.run([sys.executable, '-B', __file__, '--impl', impl], capture_output=True, text=True, check=True)
        runs.append(json.loads(p.stdout.strip().splitlines()[-1]))
    ids = list(runs[0]['times'])
    res = {}
    for m in ('compact', 'expanded'):
        ratios = []
        for i in ids:
            old = statistics.mean(r['times'][i][m] for r in runs if r['impl'] == 'old')
            new = statistics.mean(r['times'][i][m] for r in runs if r['impl'] == 'new')
            ratios.append(new / old)
        ratios.sort()
        res[m] = dict(median_new_over_old=statistics.median(ratios), min=ratios[0], max=ratios[-1], n=len(ratios))
    out = dict(libraries=sorted({r['library'] for r in runs}), result=res, runs=runs)
    (HERE / 'ab_class.json').write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(dict(libraries=out['libraries'], result=res), indent=1))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--impl'); ap.add_argument('--drive', action='store_true')
    a = ap.parse_args()
    drive() if a.drive else measure(a.impl)
