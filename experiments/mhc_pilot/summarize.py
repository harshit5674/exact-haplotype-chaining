"""Check all matching results and derive the report's numerical evidence."""
import hashlib, json, statistics
from pathlib import Path
import benchmark as b

P=Path(__file__).resolve().parent
def read(name):return [json.loads(s)for s in (P/name).read_text().splitlines()]
def dist(xs):
    xs=sorted(xs)
    return dict(min=min(xs),median=statistics.median(xs),p90=xs[int(.9*(len(xs)-1))],max=max(xs))
def indexed(name):return {r['id']:r for r in read(name)}

members=json.loads((P/'membership_profile.json').read_text())
fen=indexed('flat_fenwick.jsonl');sky=indexed('flat_skyline.jsonl');hy=indexed('final_hybrid.jsonl')
assert len(fen)==len(sky)==len(hy)==40 and fen.keys()==sky.keys()==hy.keys()
assert all(fen[k]['digest']==sky[k]['digest']==hy[k]['digest']for k in hy)
summary=dict(population=dict(reads=100,orientations=200,anchors=sum(r['N']for r in members),
    stats={k:dist([r[k]for r in members])for k in ('N','K','incidence_density','support_ratio','support_classes','p')}),
    timed=dict(reads=20,orientations=40,anchors=sum(r['N']for r in hy.values()),states=sum(r['K']for r in hy.values()),
    all_state_digests_match=True,stats={k:dist([r['stats'][k]for r in hy.values()])for k in
        ('label_max_quotient_runs','label_max_quotient_arity','label_phi_selected','p','Q')},
    phi_per_K=dist([r['stats']['label_phi_selected']/r['K']for r in hy.values()]),
    selected_run_quotients=sum(r['stats']['label_run_dc_quotients']for r in hy.values()),
    all_quotients=sum(r['stats']['label_quotient_nodes']for r in hy.values()),
    method_metrics={name:dict(total_seconds=sum(r['seconds']for r in rows.values()),
        seconds=dist([r['seconds']for r in rows.values()]),peak_rss_bytes=max(r['peak_rss_bytes']for r in rows.values()),
        graph_index_seconds=next(iter(rows.values()))['graph_index_seconds'])for name,rows in [('flat_fenwick',fen),('flat_skyline',sky),('final_hybrid',hy)]},
    hybrid_over_fenwick=dist([hy[k]['seconds']/fen[k]['seconds']for k in hy]),
    hybrid_over_skyline=dist([hy[k]['seconds']/sky[k]['seconds']for k in hy]),
    hybrid_wins_fenwick=sum(hy[k]['seconds']<fen[k]['seconds']for k in hy),
    hybrid_wins_skyline=sum(hy[k]['seconds']<sky[k]['seconds']for k in hy)))
rh=indexed('repeat_hybrid.jsonl');rf=indexed('repeat_fenwick.jsonl');rs=indexed('repeat_skyline.jsonl')
assert len(rh)==len(rf)==len(rs)==4
assert all(rh[k]['digest']==rf[k]['digest']==rs[k]['digest']==hy[k]['digest']for k in rh)
summary['repeated']=dict(instances=4,repeats_per_method=3,
    hybrid_over_fenwick=dist([rh[k]['seconds']/rf[k]['seconds']for k in rh]),
    hybrid_over_skyline=dist([rh[k]['seconds']/rs[k]['seconds']for k in rh]))
summary['ablations']={}
for name in ('run','prefix'):
    rows=read('ablation_'+name+'.jsonl');assert len(rows)==4
    assert all(r['digest']==hy[r['id']]['digest']for r in rows)
    summary['ablations'][name]=dict(instances=4,all_states_match=True,
        seconds=[r['seconds']for r in rows],over_repeated_hybrid=dist([r['seconds']/rh[r['id']]['seconds']for r in rows]))
base,cases,index,_=b.load();byid={c['id']:c for c in cases};pub=read('published.jsonl')
diversity=[]
for r in pub:
    c=byid[r['id']];ins=b.Instance(base.lengths,base.edges,base.paths,[b.Anchor(*a)for a in c['anchors']],100)
    inc=[index[0][a.v]for a in ins.anchors];out=b.flat(ins,inc,index)
    assert b.canonical_digest(out)==hy[r['id']]['digest']
    digest=hashlib.sha256(json.dumps([sorted(x.items())for x in out['C']],separators=(',',':')).encode()).hexdigest()
    assert digest==r['digest']
    distinct=[len(set(row.values()))for row in out['C']]
    diversity.append(dict(id=r['id'],N=len(distinct),K=sum(map(len,out['C'])),
        distinct_scores_total=sum(distinct),distinct_per_anchor=dist(distinct),
        all_equal_rows=sum(x==1 for x in distinct),
        note='Descriptive numeric equality only; different eligibility histories cannot be merged without a proof.'))
summary['published_reconstruction']=dict(instances=len(pub),all_C_match=True,seconds=[r['seconds']for r in pub],
    note='Supplementary C-only comparison, reconstructed source Algorithm 1, not upstream Minichain; no witness/R output and cached validation bypass.')
summary['numeric_diversity_diagnostic']=diversity
(P/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
