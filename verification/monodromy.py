"""Standalone exact source-cover and monodromy verifier.

Usage: python3 monodromy.py --data .. --loop all
Requires Python 3.9+ and SymPy. Reads only math_core.json and
monodromy_tubes.json.gz from the supplied public data directory.
No historic receipts, paths, compiled binaries, or proposed permutations
are used as mathematical premises.
"""
import argparse, gzip, hashlib, json, signal, time
from fractions import Fraction as F
from pathlib import Path
import sympy as S
from cover import check as check_cover
from monodromy_intervals import Box, g_at, rouche_at, one_segment, point
from monodromy_geometry import check_manifest, verify_path, identify_psl

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def run(directory,chosen):
    started=time.monotonic()
    core_path=directory/'math_core.json';tube_path=directory/'monodromy_tubes.json.gz'
    core=json.loads(core_path.read_text())
    with gzip.open(tube_path,'rt') as f:trace=json.load(f)
    cover=check_cover(core)
    manifest=trace['manifest'];k=S.Rational(core['cover']['kappa'])
    assert list(map(S.Rational,manifest['B_coefficients_ascending']))==[-2*k,-3*k,S.Integer(0),S.Integer(1)]
    isolation=check_manifest(manifest,k)
    assert trace['fixed_point_scale']=='10^100' and trace['base_disc_scaled_radius']=='10^80'
    assert point(trace['basepoint'])==point(manifest['basepoint'])
    c=core['cover'];assert list(map(F,c['D']))==[3,0,1]
    data={'P':list(map(Box.rational,c['P0'])),'R':list(map(Box.rational,c['R0']))}
    D=list(map(Box.rational,c['D']))
    base_roots=list(map(Box.pair,trace['base_centers']));assert len(base_roots)==14
    base_t=Box.pair(trace['basepoint']);base_radius=10**80
    for center in base_roots:rouche_at(g_at(data,D,base_t),center,base_radius)
    assert all((a-b).lower()>2*base_radius for i,a in enumerate(base_roots) for b in base_roots[i+1:])
    names=['b1','b2','b3'];assert set(trace['loops'])==set(names)
    names=names if chosen=='all' else [chosen]
    reports={};permutations=[]
    for name in names:
        loop_start=time.monotonic()
        frames=trace['loops'][name]['frames'];verify_path(name,frames,manifest)
        assert len(frames)-1=={'b1':183,'b2':179,'b3':144}[name]
        maximum=F(0);last_radii=None;maximum_segment_seconds=0
        for j,(start,end) in enumerate(zip(frames,frames[1:])):
            seg_start=time.monotonic()
            try:radii,ratio=one_segment(data,D,start,end)
            except AssertionError as e:raise AssertionError((name,j,e.args)) from e
            duration=time.monotonic()-seg_start
            assert duration<35
            maximum_segment_seconds=max(maximum_segment_seconds,duration)
            maximum=max(maximum,ratio)
            if j==0:
                for i,center in enumerate(map(Box.pair,start['centers'])):
                    assert (center-base_roots[i]).upper()+base_radius<radii[i]
            last_radii=radii
        assert last_radii is not None
        permutation=[]
        for i,center in enumerate(map(Box.pair,frames[-1]['centers'])):
            matches=[j for j,b in enumerate(base_roots) if (center-b).upper()+base_radius<last_radii[i]]
            assert len(matches)==1;permutation.append(matches[0])
        assert sorted(permutation)==list(range(14))
        permutations.append(tuple(permutation))
        reports[name]={'certified_segments':len(frames)-1,'permutation':permutation,'max_uniform_rouche_ratio':str(maximum),'maximum_segment_seconds':maximum_segment_seconds,'seconds':time.monotonic()-loop_start}
    group=identify_psl(permutations) if chosen=='all' else None
    return {'status':'PASS_EXACT_COVER_AND_GEOMETRIC_PSL2_F13' if group else 'PASS_EXACT_COVER_AND_SINGLE_MONODROMY_LOOP','input_sha256':{'math_core.json':digest(core_path),'monodromy_tubes.json.gz':digest(tube_path)},'source_sha256':{name:digest(Path(__file__).parent/name) for name in ['monodromy.py','monodromy_intervals.py','monodromy_geometry.py','cover.py']},'cover':cover,'branch_root_isolation':isolation,'loops':reports,'total_segments':sum(r['certified_segments'] for r in reports.values()),'finite_group':group,'arithmetic_specialization_and_SL2_lift_checked':False,'seconds':time.monotonic()-started}

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data',type=Path,default=Path(__file__).resolve().parent.parent)
    ap.add_argument('--loop',choices=['all','b1','b2','b3'],default='all')
    args=ap.parse_args()
    # Each individual invocation is deliberately bounded; on slower machines
    # select one loop. For all-loop certification the aggregate must also fit.
    if hasattr(signal,'alarm'):signal.alarm(35)
    print(json.dumps(run(args.data,args.loop),indent=2))
if __name__=='__main__':main()
