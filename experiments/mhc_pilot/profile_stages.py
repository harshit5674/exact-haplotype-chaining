"""Diagnostic stage timings, deliberately separate from main benchmark timings."""
import collections, gc, json, time
from pathlib import Path
import benchmark as b

def run(case,base,index):
    ins=b.Instance(base.lengths,base.edges,base.paths,[b.Anchor(*a)for a in case['anchors']],100)
    inc=[index[0][a.v]for a in ins.anchors]
    times=collections.defaultdict(float);saved={};capture={};owner_depth=0
    def instrument(name):
        old=getattr(b.cf,name);saved[name]=old
        def wrapped(*args,**kw):
            nonlocal owner_depth
            if name=='owner_schedule_engines':owner_depth+=1
            t=time.perf_counter();result=old(*args,**kw);elapsed=time.perf_counter()-t
            if name!='ScalarEngine' or owner_depth==0:times[name]+=elapsed
            if name=='owner_schedule_engines':owner_depth-=1
            if name=='restrict_all':
                hist=collections.Counter();allarity=q_arity=0
                for h,spec in result[0].items():
                    for u,kids in enumerate(spec.children):
                        allarity+=len(kids)
                        if spec.kind[u]=='Q':
                            ranks=spec.ranks[u];r=len(ranks);d=1+sum(x>y for x,y in zip(ranks,ranks[1:]));hist[r,d]+=1;q_arity+=r
                capture.update(quotients=[dict(arity=r,runs=d,count=n)for (r,d),n in sorted(hist.items())],
                    q_arity=q_arity,total_arity=allarity)
            return result
        setattr(b.cf,name,wrapped)
    for name in ('endpoint_orders','decompose','restrict_all','ScalarEngine','compile_owner_geometry','owner_schedule_engines'):instrument(name)
    gc.collect();t=time.perf_counter()
    try:out=b.cf.solve(ins,inc=inc,cover_index=index)
    finally:
        total=time.perf_counter()-t
        for name,old in saved.items():setattr(b.cf,name,old)
    times['evaluation_and_other_wrapper']=total-sum(times.values())
    return dict(id=case['id'],N=len(ins.anchors),K=sum(map(len,inc)),total_seconds=total,
                stage_seconds=dict(times),structure=capture,stats=out['stats'],digest=b.canonical_digest(out))

if __name__=='__main__':
    base,cases,index,prep=b.load()
    selected=set(json.loads((b.HERE/'input_manifest.json').read_text())['read_names'][:20])
    cases=sorted([c for c in cases if c['read']in selected],key=lambda c:len(c['anchors']))
    chosen=[cases[i]for i in (9,19,29,37)]
    rows=[]
    for c in chosen:
        row=run(c,base,index);rows.append(row);print(row['N'],row['stage_seconds'],flush=True)
    (b.HERE/'stage_profile.json').write_text(json.dumps(dict(graph_index_seconds=prep,profiles=rows),indent=2))
