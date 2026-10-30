"""Same-input basic-objective pilot. No upstream mapping-speed claim.

Timed modes use cached graph indexes, include validation and all per-query setup,
and compute all valid states. The published reconstruction lacks witness output;
it is supplementary, not a like-for-like full-output baseline.
"""
import argparse, collections, gc, hashlib, json, math, platform, random, resource
import statistics, sys, time
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'research/current/new'))
import chaining_final as cf
import permutation_engine as pe
import chaining as original
import continuation
from chaining import Anchor, Instance, NEG

def load():
    g=json.loads((HERE/'data/graph.json').read_text())
    cases=[json.loads(s) for s in (HERE/'data/instances.jsonl').read_text().splitlines()]
    base=Instance(g['lengths'],[tuple(e) for e in g['edges']],g['paths'],[],100)
    t=time.perf_counter();idx=cf.index_cover(base,base.paths);seconds=time.perf_counter()-t
    return base,cases,idx,seconds

def flat(ins,inc,index,backend=pe.FenwickPrefix):
    """Incidence/owner sweep, no decomposition. Reuses only owner geometry."""
    ins.validate(); A=ins.anchors;n=len(A); pi,leaves=cf.endpoint_orders(ins)
    graph_leaves=[0]*(2*n)
    for rank,e in zip(pi,leaves):graph_leaves[rank]=e
    owner,updates,urank,cuts,stats=cf.compile_owner_geometry(ins,ins.paths,graph_leaves,index)
    counts=[0]*len(ins.paths);qr=[{}for _ in A];ur=[{}for _ in A]
    for e in graph_leaves:
        i=e//2
        for h in inc[i]:
            if e%2:counts[h]+=1;ur[i][h]=counts[h]
            else:qr[i][h]=counts[h]
    local=[backend(k) if k else None for k in counts];glob=[backend(len(u)) for u in updates]
    C=[{}for _ in A];back=[{}for _ in A];B=[NEG]*n;R=[NEG]*n;who=[None]*n
    for e in leaves:
        i=e//2
        if not e%2:
            r=pe.BOTTOM
            for s,t in enumerate(cuts[i]):
                if t:r=pe.maximum(r,glob[s].query(t))
            R[i]=r[0]
            for h in inc[i]:
                best=pe.maximum((0,None),local[h].query(qr[i][h]));best=pe.maximum(best,(r[0]-ins.gamma,r[1]))
                C[i][h]=A[i].weight+best[0];back[i][h]=best[1]
                if C[i][h]>B[i]:B[i]=C[i][h];who[i]=(i,h)
        else:
            s=owner[i];glob[s].update(urank[i],B[i],who[i])
            for h in inc[i]:local[h].update(ur[i][h],C[i][h],(i,h))
    stats.update(N=n,K=sum(map(len,inc)))
    return dict(C=C,B=B,R=R,back=back,stats=stats)

def canonical_digest(out):
    # Exact integer results and membership keys; pointers may differ on ties.
    obj=[[[[h,v]for h,v in sorted(row.items())]for row in out['C']],out['B'],out['R']]
    return hashlib.sha256(json.dumps(obj,separators=(',',':')).encode()).hexdigest()

def reachability(base):
    # Independent graph DFS, no path-cover reachability index.
    adj=[[]for _ in base.lengths]
    for u,v in base.edges:adj[u].append(v)
    reach=[]
    for start in range(len(adj)):
        seen={start};stack=[start]
        while stack:
            for v in adj[stack.pop()]:
                if v not in seen:seen.add(v);stack.append(v)
        reach.append(seen)
    return reach

def check_pointers(ins,out,inc,reach):
    total=0
    for i,row in enumerate(out['C']):
        assert set(row)==set(inc[i]);assert out['B'][i]==max(row.values())
        a=ins.anchors[i]
        for h,value in row.items():
            ptr=out['back'][i][h]
            if ptr is None:assert value==a.weight
            else:
                j,g=ptr;b=ins.anchors[j]
                assert b.d<a.c and a.v in reach[b.v] and (a.v!=b.v or b.y<a.x)
                assert g in inc[j] and value==out['C'][j][g]+a.weight-ins.gamma*(g!=h)
            total+=1
    # Strict query advancement makes these local witness checks an inductive
    # certificate for every full traceback, without reconstructing K chains.
    return total

def direct(ins,inc,reach):
    """Literal original transition over every eligible predecessor state."""
    A=ins.anchors;C=[{}for _ in A];B=[NEG]*len(A);R=B.copy()
    for i in sorted(range(len(A)),key=lambda i:A[i].c):
        a=A[i];pred=[j for j,b in enumerate(A) if b.d<a.c and a.v in reach[b.v] and (b.v!=a.v or b.y<a.x)]
        R[i]=max((B[j]for j in pred),default=NEG)
        for h in inc[i]:
            best=0
            for j in pred:
                for g,v in C[j].items():best=max(best,v-ins.gamma*(g!=h))
            C[i][h]=a.weight+best
        B[i]=max(C[i].values())
    return dict(C=C,B=B,R=R)

def run_mode(mode,ins,inc,index):
    if mode.startswith('flat'):return flat(ins,inc,index,pe.FenwickPrefix if mode=='flat_fenwick' else cf.SkylinePrefix)
    if mode=='published':
        # Reuse graph preprocessing without changing the archived solver.
        old=original.preprocess;original.preprocess=lambda obj:index
        try:rows=original.published_sweep(ins,structure=pe.FenwickPrefix)
        finally:original.preprocess=old
        return dict(C=[{h:row[h]for h in hs}for row,hs in zip(rows,inc)],stats={})
    policy=mode.removeprefix('final_')
    return cf.solve(ins,inc=inc,cover_index=index,quotient_policy=policy,owner_policy=policy)

def timed(args):
    base,cases,index,prep=load();rng=random.Random(741);rng.shuffle(cases)
    if args.reads:
        selected=set(json.loads((HERE/'input_manifest.json').read_text())['read_names'][:args.reads])
        cases=[c for c in cases if c['read'] in selected]
    cases=cases[:args.limit] if args.limit else cases
    outpath=Path(args.out);assert not outpath.exists();outpath.parent.mkdir(parents=True,exist_ok=True)
    reach=reachability(base)
    with outpath.open('w') as f:
        for num,case in enumerate(cases):
            A=[Anchor(*a)for a in case['anchors']];ins=Instance(base.lengths,base.edges,base.paths,A,args.gamma)
            t=time.perf_counter();inc=[index[0][a.v]for a in A];membership_seconds=time.perf_counter()-t
            timings=[];digests=[];states=0
            for rep in range(args.repeat):
                gc.collect();t=time.perf_counter();out=run_mode(args.mode,ins,inc,index);timings.append(time.perf_counter()-t)
                if args.mode!='published':
                    digests.append(canonical_digest(out))
                    if rep==0:states=check_pointers(ins,out,inc,reach)
                else:
                    digests.append(hashlib.sha256(json.dumps([sorted(r.items())for r in out['C']],separators=(',',':')).encode()).hexdigest())
                stats=out['stats'];del out
            assert len(set(digests))==1
            row=dict(id=case['id'],mode=args.mode,gamma=args.gamma,N=len(A),H=len(base.paths),K=sum(map(len,inc)),
                times=timings,seconds=statistics.median(timings),digest=digests[0],witness_states=states,stats=stats,
                graph_index_seconds=prep,membership_seconds=membership_seconds,
                peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
            f.write(json.dumps(row)+'\n');f.flush()
            if num%10==0:print(args.mode,num+1,'/',len(cases),'N',len(A),'seconds',round(row['seconds'],3),flush=True)

def verify(args):
    base,cases,index,_=load();reach=reachability(base);rng=random.Random(931)
    checks=[]
    # Small independently evaluated subinstances retain original graph geometry.
    for case in cases[::10]:
        anchors=rng.sample(case['anchors'],min(24,len(case['anchors'])))
        for gamma in (0,17,100):
            ins=Instance(base.lengths,base.edges,base.paths,[Anchor(*a)for a in anchors],gamma)
            inc=[index[0][a.v]for a in ins.anchors];truth=direct(ins,inc,reach)
            for mode in ('flat_fenwick','flat_skyline','final_hybrid','final_run','final_prefix'):
                out=run_mode(mode,ins,inc,index)
                assert all(out[k]==truth[k]for k in ('C','B','R')),(case['id'],gamma,mode)
                check_pointers(ins,out,inc,reach)
            checks.append(dict(id=case['id'],N=len(anchors),gamma=gamma,K=sum(map(len,inc))))
    Path(args.out).write_text(json.dumps(dict(status='passed',instances=len(checks),modes=5,checks=checks),indent=2))
    print('Independent direct-recurrence checks passed:',len(checks),'instances x 5 modes')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--mode',default='final_hybrid');p.add_argument('--out',required=True)
    p.add_argument('--limit',type=int,default=0);p.add_argument('--reads',type=int,default=0);p.add_argument('--repeat',type=int,default=1);p.add_argument('--gamma',type=int,default=100)
    a=p.parse_args();(verify if a.mode=='verify' else timed)(a)
