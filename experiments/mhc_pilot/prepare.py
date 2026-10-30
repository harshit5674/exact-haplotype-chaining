"""Reproducible public-data sampling and exact seed conversion for the MHC pilot."""
import argparse, collections, gzip, hashlib, heapq, json, random, re, subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / 'data'
SEED = 20261004
RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')

def graph():
    seqs, edges, paths, names = {}, [], [], []
    with gzip.open(DATA/'MHC-CHM13.0.gfa.gz', 'rt') as f:
        for line in f:
            a = line.rstrip().split('\t')
            if a[0] == 'S': seqs[a[1]] = a[2]
            elif a[0] == 'L':
                assert a[2] == a[4] == '+' and a[5] == '0M'
                edges.append((a[1], a[3]))
            elif a[0] == 'W':
                assert '<' not in a[6]
                paths.append(re.findall(r'>([^><]+)', a[6])); names.append('|'.join(a[1:4]))
    assert set(seqs) == {v for p in paths for v in p}
    indeg = dict.fromkeys(seqs, 0); adj = collections.defaultdict(list)
    for u,v in edges: adj[u].append(v); indeg[v] += 1
    ready = [v for v in seqs if not indeg[v]]; heapq.heapify(ready); order=[]
    while ready:
        u=heapq.heappop(ready);order.append(u)
        for v in adj[u]:
            indeg[v]-=1
            if not indeg[v]:heapq.heappush(ready,v)
    assert len(order)==len(seqs), 'cycle'
    ids={v:i for i,v in enumerate(order)}
    return dict(lengths=[len(seqs[v]) for v in order], edges=[[ids[u],ids[v]] for u,v in edges],
                paths=[[ids[v] for v in p] for p in paths], names=order, haplotypes=names), seqs, ids

def sample():
    rng=random.Random(SEED); reservoir=[]; total=0
    with gzip.open(DATA/'MHC_CHM13_PacBio_filt.fq.gz','rt') as f:
        while True:
            h=f.readline()
            if not h:break
            rec=(h,f.readline(),f.readline(),f.readline())
            assert h.startswith('@') and rec[2].startswith('+') and len(rec[1])==len(rec[3])
            if total<100:reservoir.append((total,rec))
            else:
                j=rng.randrange(total+1)
                if j<100:reservoir[j]=(total,rec)
            total+=1
    reservoir.sort()
    with (DATA/'sample100.fq').open('w') as f:
        for _,r in reservoir:f.writelines(r)
    g,_,_=graph(); (DATA/'graph.json').write_text(json.dumps(g))
    files={p.name:dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
        md5=hashlib.md5(p.read_bytes()).hexdigest()) for p in DATA.glob('*.gz')}
    assert files['MHC-CHM13.0.gfa.gz']['md5']=='96c9a32af21f741237aef8353b0e1f78'
    assert files['MHC_CHM13_PacBio_filt.fq.gz']['md5']=='fe92bba8f586209583e264eddd0c15b8'
    meta=dict(dataset='https://zenodo.org/records/10665350', seed=SEED, population_reads=total,
        sample_size=100, sample_indices=[i for i,_ in reservoir], read_names=[r[0].split()[0][1:] for _,r in reservoir],
        files=files, graph_vertices=len(g['lengths']), graph_edges=len(g['edges']), H=len(g['paths']),
        minichain_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=HERE/'vendor/minichain',text=True).strip())
    (HERE/'input_manifest.json').write_text(json.dumps(meta,indent=2)); print(json.dumps({k:v for k,v in meta.items() if k not in ('read_names','sample_indices','files')},indent=2))

def convert():
    g,seqs,ids=graph(); reads={}; reverse_seqs={}
    with (DATA/'sample100.fq').open() as f:
        while h:=f.readline():
            reads[h.split()[0][1:]]=f.readline().strip();f.readline();f.readline()
    cases=[]; seen=set(); seeds=0; current=None
    with (DATA/'anchors.tsv').open() as f:
        for line in f:
            a=line.rstrip().split('\t')
            if a[0]=='R':
                if current:
                    assert sum(len(x['anchors']) for x in current)==expected
                    cases.extend(current)
                name,qlen,expected=a[1],int(a[2]),int(a[3]);assert name not in seen;seen.add(name)
                assert len(reads[name])==qlen
                current=[dict(id=name+':'+strand,read=name,strand=strand,qlen=qlen,anchors=[]) for strand in ('+','-')]
            else:
                assert a[0]=='A'; seg,strand=a[1],int(a[2]);y,d,span=map(int,a[3:]);x=y-span+1;c=d-span+1
                if strand and seg not in reverse_seqs:reverse_seqs[seg]=seqs[seg].translate(RC)[::-1]
                ref=seqs[seg] if strand==0 else reverse_seqs[seg]
                assert 0<=x<=y<len(ref) and 0<=c<=d<qlen
                assert ref[x:y+1].upper()==reads[name][c:d+1].upper(), (name,seg,x,y,c,d)
                if strand:
                    x,y=len(ref)-1-y,len(ref)-1-x;c,d=qlen-1-d,qlen-1-c
                current[strand]['anchors'].append([ids[seg],x+1,y+1,c+1,d+1,span]);seeds+=1
        if current:
            assert sum(len(x['anchors']) for x in current)==expected
            cases.extend(current)
    assert seen==set(reads)
    with (DATA/'instances.jsonl').open('w') as f:
        for c in cases:f.write(json.dumps(c)+'\n')
    meta=dict(reads=len(seen),orientation_instances=len(cases),nonempty_instances=sum(bool(c['anchors']) for c in cases),
              exact_sequence_checks=seeds,anchor_weight='minimizer span',coordinates='1-based closed',
              reverse_policy='reverse-complement both graph path and query to the forward graph',
              filter='unchanged Minichain lr seed collection; no further anchor filtering')
    (HERE/'conversion.json').write_text(json.dumps(meta,indent=2));print(json.dumps(meta,indent=2))

def profile():
    g,_,_=graph();mem=[[]for _ in g['lengths']]
    for h,path in enumerate(g['paths']):
        for v in path:mem[v].append(h)
    rows=[]
    for line in (DATA/'instances.jsonl').read_text().splitlines():
        case=json.loads(line);aa=case['anchors'];verts={a[0]for a in aa};groups={}
        for h,path in enumerate(g['paths']):
            sig=tuple(v for v in path if v in verts)
            if sig:groups.setdefault(sig,[]).append(h)
        counts=collections.Counter(a[0]for a in aa)
        kh=sum(sum(counts[v]for v in sig)for sig in groups)
        k=sum(len(mem[a[0]])for a in aa);H=len(g['paths'])
        rows.append(dict(id=case['id'],N=len(aa),K=k,H=H,Khat=kh,support_classes=len(groups),
            incidence_density=k/(H*len(aa)) if aa else None,support_ratio=kh/k if k else None,
            p=len({mem[a[0]][0]for a in aa})))
    (HERE/'membership_profile.json').write_text(json.dumps(rows,indent=2));print('Profiled',len(rows),'orientation instances')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['sample','convert','profile']);a=p.parse_args()
    {'sample':sample,'convert':convert,'profile':profile}[a.action]()
