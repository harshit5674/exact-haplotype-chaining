"""Fresh validation/paired timings for the repaired guard; never overwrites old results.

python3.13 -B experiments/gap_safeguard_repair/run.py --phase validate --out <new.json>
python3.13 -B experiments/gap_safeguard_repair/run.py --phase benchmark --out <new.json>
Both phases use the saved primary cases at weights x200, gamma=10000. Benchmarking
compares old/new implementations in one process, in shuffled order, on identical batches.
Native kernel timings include read setup, class initialization, sweep and packed export.
Prepared-graph construction is excluded and shared. Every timed output is checked after
the native timer stops: C/B/R must match the per-haplotype baseline and pointers must
match the unchanged class solver or baseline, according to the branch that completed.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import random
import statistics
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'research/current/new'))
sys.path.insert(0, str(ROOT/'experiments/gap_fallback'))
import gap_chaining as gc
import gap_native as new
import auto_real as data
from test_gap_budget import collision_case
from test_gap_native import star_case

spec = importlib.util.spec_from_file_location('gap_native_before_guard_repair', HERE/'prepatch/gap_native.py')
old = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = old
spec.loader.exec_module(old)


def digest_file(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def summary(values):
    v = sorted(values)
    return dict(n=len(v), median=statistics.median(v), min=v[0], max=v[-1],
                p10=v[int(.1*(len(v)-1))], p90=v[int(.9*(len(v)-1))])


def validate(out, checkpoint):
    for locus, panel, directory in data.panels():
        graph, cases = data.load(directory)
        p = gc.Graph(graph['lengths'], [tuple(e) for e in graph['edges']], graph['paths']).prepare()
        G = new.GapGraph(p); previous = old.GapGraph(p)
        rows = []
        for c in data.primary(cases, locus):
            a = [gc.Anchor(*x[:5], x[5]*data.SCALE) for x in c['anchors']]
            b = G.batch(a, data.GAMMA)
            ref = G.solve_packed(b, mode='per_haplotype')
            z = G.solve_packed(b, mode='auto')
            assert z.arrays_digest() == ref.arrays_digest(), (locus, panel, c['id'])
            ptr_ref = ref if z.stats['fallback'] else previous.solve_packed(a, data.GAMMA, mode='class')
            assert z.back == ptr_ref.back, (locus, panel, c['id'], 'pointers')
            assert z.stats['auto_work'] <= z.stats['auto_budget']
            rows.append(dict(id=c['id'], N=b.n, K=b.K, digest=z.arrays_digest(),
                             pointer_digest=hashlib.sha256(z.back.tobytes()).hexdigest(),
                             guard={k:v for k,v in z.stats.items() if k.startswith('auto_') or k=='fallback'}))
        row = dict(locus=locus, panel=panel, H=p.H, instances=len(rows), states=sum(r['K'] for r in rows),
                   fallbacks=sum(r['guard']['fallback'] for r in rows), rows=rows,
                   input_hashes={f:digest_file(directory/f) for f in ('graph.json','instances.jsonl')})
        out['panels'].append(row); checkpoint()
        print(json.dumps({k:v for k,v in row.items() if k not in ('rows','input_hashes')}), flush=True)


def measure(p, anchors, gamma, reps, rng):
    graphs = {'old':old.GapGraph(p), 'new':new.GapGraph(p)}
    batches = {k:g.batch(anchors,gamma) for k,g in graphs.items()}
    b = batches['new']
    ref = graphs['new'].solve_packed(b, mode='per_haplotype')
    expected = ref.arrays_digest()
    class_ref = graphs['old'].solve_packed(batches['old'], mode='class')
    assert class_ref.arrays_digest() == expected
    def check(k,z):
        checked = z.expand() if k.endswith('compact') else z
        assert checked.arrays_digest() == expected, k
        ptr_ref = ref if k == 'per_haplotype' or z.stats.get('fallback',False) else class_ref
        assert checked.back == ptr_ref.back, (k,'pointers')
    actions = {
        'old_auto': lambda: graphs['old'].solve_packed(batches['old'], mode='auto'),
        'new_auto': lambda: graphs['new'].solve_packed(b, mode='auto'),
        'old_class': lambda: graphs['old'].solve_packed(batches['old'], mode='class'),
        'new_class': lambda: graphs['new'].solve_packed(b, mode='class'),
        'per_haplotype': lambda: graphs['new'].solve_packed(b, mode='per_haplotype'),
        'old_compact': lambda: graphs['old'].solve_compact(batches['old']),
        'new_compact': lambda: graphs['new'].solve_compact(b),
    }
    times = {k:[] for k in actions}; packed = {k:[] for k in actions}; guards = {}
    # Untimed warm-up for code/data pages. Pair caches still start cold on every call.
    for k in rng.sample(list(actions), len(actions)):
        z = actions[k]()
        check(k,z)
    for _ in range(reps):
        for k in rng.sample(list(actions), len(actions)):
            z = actions[k]()
            times[k].append(z.stats['kernel_seconds'])
            packed[k].append(z.stats['packed_seconds'])
            if k.endswith('auto'):
                guards[k] = {f:v for f,v in z.stats.items() if f.startswith('auto_') or f=='fallback'}
                if k == 'new_auto': assert z.stats['auto_work'] <= z.stats['auto_budget']
            check(k,z)
    med = {k:statistics.median(v) for k,v in times.items()}
    return dict(N=b.n, H=p.H, K=b.K, digest=expected, kernel_seconds=times, packed_seconds=packed,
                medians=med, guards=guards,
                ratios=dict(new_auto_over_old=med['new_auto']/med['old_auto'],
                            new_class_over_old=med['new_class']/med['old_class'],
                            new_compact_over_old=med['new_compact']/med['old_compact'],
                            new_auto_over_class=med['new_auto']/med['new_class'],
                            new_auto_over_ph=med['new_auto']/med['per_haplotype'],
                            old_auto_over_ph=med['old_auto']/med['per_haplotype'],
                            ph_over_new_auto=med['per_haplotype']/med['new_auto']))


def benchmark(out, checkpoint, reps, limit):
    rng = random.Random(20261014)
    groups = {}
    for locus, panel, directory in data.panels():
        groups.setdefault(locus, []).append((panel, directory))
    for locus, panels in groups.items():
        loaded = [(panel,d,data.load(d)) for panel,d in panels]
        loaded.sort(key=lambda x:len(x[2][0]['paths']))
        for panel, directory, (graph,cases) in (loaded[0], loaded[-1]):
            p = gc.Graph(graph['lengths'], [tuple(e) for e in graph['edges']], graph['paths']).prepare()
            sel = data.primary(cases,locus)
            if limit: sel=sel[:limit]
            rows=[]; load=os.getloadavg()
            for c in sel:
                a=[gc.Anchor(*x[:5],x[5]*data.SCALE) for x in c['anchors']]
                rows.append(dict(id=c['id'], **measure(p,a,data.GAMMA,reps,rng)))
            row=dict(locus=locus,panel=panel,H=p.H,instances=len(rows),load_before=load,load_after=os.getloadavg(),
                     ratios={k:summary([r['ratios'][k] for r in rows]) for k in rows[0]['ratios']},
                     fallbacks={k:sum(r['guards'][k]['fallback'] for r in rows) for k in ('old_auto','new_auto')},rows=rows,
                     input_hashes={f:digest_file(directory/f) for f in ('graph.json','instances.jsonl')})
            out['panels'].append(row); checkpoint()
            print(json.dumps({k:v for k,v in row.items() if k not in ('rows','input_hashes')}),flush=True)
    for family, params in [('collision',m) for m in (128,1024,2048,4096,8192)]+[('star',m) for m in (8,64,256,1024)]:
        g,a,gamma=collision_case(params) if family=='collision' else star_case(params,64)
        row=dict(family=family,size=params,**measure(g.prepare(),a,gamma,reps,rng))
        out['adversarial'].append(row); checkpoint()
        print(json.dumps({k:v for k,v in row.items() if k not in ('kernel_seconds','packed_seconds','digest')}),flush=True)


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--phase',choices=('validate','benchmark'),required=True)
    ap.add_argument('--out',type=Path,required=True); ap.add_argument('--reps',type=int,default=5)
    ap.add_argument('--limit',type=int,default=0); args=ap.parse_args()
    assert not args.out.exists(), 'Use a fresh result path'
    started=time.time()
    out=dict(phase=args.phase,command=sys.argv,started=started,platform=platform.platform(),python=sys.version,
             load_start=os.getloadavg(),libraries={'old':old.build()._name,'new':new.build()._name},
             source_hashes={str(p.relative_to(ROOT)):digest_file(p) for p in
                           [ROOT/'research/current/new/gap_native.cpp',ROOT/'research/current/new/gap_native.py',
                            HERE/'prepatch/gap_native.cpp',HERE/'run.py']},
             repetitions=args.reps,limit=args.limit,panels=[],adversarial=[],status='running')
    def checkpoint(): args.out.write_text(json.dumps(out,indent=2)+'\n')
    checkpoint()
    try:
        validate(out,checkpoint) if args.phase=='validate' else benchmark(out,checkpoint,args.reps,args.limit)
        out['status']='all passed'
    except BaseException as e:
        out['status']='failed'; out['error']=repr(e); raise
    finally:
        out['seconds']=time.time()-started; out['load_end']=os.getloadavg(); checkpoint()


if __name__=='__main__': main()
