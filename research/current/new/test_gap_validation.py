"""Cross-check the independent native validator, including rejected corruptions."""
import argparse,json,sys,time
from pathlib import Path
import gap_chaining as gc
import gap_native as gn
import gap_validation as gv
from test_gap_chaining import random_case,fixed_cases

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    if args.out.exists():raise SystemExit('refusing overwrite')
    cases=[case for _,case in fixed_cases()]+[random_case(10000+s,wide=True) for s in range(350)]
    passed=states=skipped=mutations=0;t=time.perf_counter()
    for graph,A,gamma in cases:
        p=graph.prepare()
        try:
            G=gn.GapGraph(p);b=G.batch(A,gamma)
        except (gn.NativeRangeError,OverflowError):skipped+=1;continue
        if not b.max_key<1<<38:skipped+=1;continue
        validator=gv.Validator(dict(lengths=graph.lengths,edges=graph.edges,paths=graph.paths))
        vb=validator.batch(A,gamma)
        expected=gc.direct_oracle(p,A,gamma)
        for mode in ('per_haplotype','class','auto'):
            result=G.solve_packed(b,mode=mode);vb.check(result)
            C,ptr,B,R=result.as_dicts(A)
            assert (C,B,R)==expected[:3]
            gc.verify_pointers(p,A,gamma,C,B,ptr)
        if b.K:
            for field in ('C','B','R','back'):
                arr=getattr(result,field);old=arr[0]
                arr[0]=len(result.C) if field=='back' else old+1
                try:vb.check(result)
                except AssertionError:mutations+=1
                else:raise AssertionError(('corruption accepted',field))
                finally:arr[0]=old
        vb.close();validator.close();passed+=1;states+=b.K
    args.out.write_text(json.dumps(dict(passed=passed,states=states,skipped_outside_native_domain=skipped,
        corruptions_rejected=mutations,seconds=time.perf_counter()-t,argv=sys.argv),indent=2)+'\n')
    print(args.out.read_text())
if __name__=='__main__':main()
