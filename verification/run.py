"""Source-only exact verification. Run one bounded phase per invocation.

Python 3.11+, GMP and a C++17 compiler; no archive binaries, CAS or network.
All mathematical inputs are relative to this source's parent directory.
The work directory is newly generated, append-only computational output.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,gzip,hashlib,json,os,shlex,shutil,signal,subprocess,sys,time
HERE=Path(__file__).resolve().parent; DATA=HERE.parent
sys.set_int_max_str_digits(0)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):
    p=Path(p);return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
def put(p,obj):
    raw=(json.dumps(obj,separators=(',',':'))+'\n').encode()
    if p.suffix=='.gz':raw=gzip.compress(raw,compresslevel=1,mtime=0)
    with p.open('xb') as f:f.write(raw)
def Eblock(v):
    qq=list(map(Q,v));d=1
    import math
    for q in qq:d=math.lcm(d,q.denominator)
    return {'denominator':str(d),'numerator_coefficients_a_ascending':[str(q.numerator*(d//q.denominator)) for q in qq]}
def E_to_V(e):return {'denominator':e['denominator'],'numerator_rows_x_then_y':[[n]+['0']*12 for n in e['numerator_coefficients_a_ascending']]}
def emitV(v):
    assert len(v['numerator_rows_x_then_y'])==14 and all(len(r)==13 for r in v['numerator_rows_x_then_y'])
    return v['denominator']+' '+' '.join(n for row in v['numerator_rows_x_then_y'] for n in row)+'\n'
def emitE(e):return e['denominator']+' '+' '.join(e['numerator_coefficients_a_ascending'])+'\n'
def execute(cmd,packet=None,timeout=35):
    child=subprocess.Popen(cmd,stdin=subprocess.PIPE if packet is not None else subprocess.DEVNULL,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,start_new_session=True)
    try:out,err=child.communicate(packet,timeout=timeout)
    except subprocess.TimeoutExpired:
        os.killpg(child.pid,signal.SIGKILL);out,err=child.communicate();raise RuntimeError('resource cap reached; NOT PASS: '+err[:400])
    assert child.returncode==0,(child.returncode,err[:3000]);return out
def main():
    assert __debug__ and not os.environ.get('PYTHONOPTIMIZE'),'Assertions must remain enabled'
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('phase');ap.add_argument('--work',type=Path,default=HERE/'work');ap.add_argument('--cxx',default='c++');ap.add_argument('--gmp-prefix',type=Path);a=ap.parse_args()
    work=a.work.resolve();work.mkdir(parents=True,exist_ok=True);phase=a.phase;began=time.monotonic();dest=work/(phase+'.json.gz');assert not dest.exists(),'Append-only: select a new work directory or an unexecuted phase'
    if phase=='list':
        print('build, g, spin0..spin6, spin-conclusion, norm01..norm13, h27..h00, encoding, conclusion');return
    if phase=='build':
        compiler=shutil.which(a.cxx);assert compiler;flags=[];pkg=shutil.which('pkg-config')
        if pkg:
            q=subprocess.run([pkg,'--cflags','--libs','gmpxx'],capture_output=True,text=True,timeout=5)
            if q.returncode==0:flags=shlex.split(q.stdout)
        if not flags:flags=['-lgmpxx','-lgmp']
        if a.gmp_prefix:flags=['-I'+str(a.gmp_prefix/'include'),'-L'+str(a.gmp_prefix/'lib')]+flags
        entries=[]
        for name in ['norm13','norm28']:
            binary=work/name;assert not binary.exists();cmd=[compiler,'-O2','-std=c++17',str(HERE/(name+'.cpp')),'-o',str(binary),*flags];execute(cmd)
            entries.append({'source':name+'.cpp','source_sha256':sha(HERE/(name+'.cpp')),'binary':name,'binary_sha256':sha(binary),'compile_command':cmd})
        out={'status':'PASS_FRESH_SOURCE_BUILD','entries':entries,'header_sha256':sha(HERE/'exact_V.hpp'),'runner_sha256':sha(Path(__file__)),'seconds':time.monotonic()-began};put(dest,out);print(json.dumps({'phase':phase,'status':out['status'],'seconds':out['seconds']}));return
    inputs={};deps={}
    def read_data(name):
        p=DATA/name;inputs[name]=sha(p);return load(p)
    def previous(name):
        p=work/(name+'.json.gz');obj=load(p);deps[name]=sha(p)
        assert obj['runner_sha256']==sha(Path(__file__))
        for n,h in obj.get('inputs_sha256',{}).items():assert sha(DATA/n)==h
        for n,h in obj.get('dependencies_sha256',{}).items():assert sha(work/(n+'.json.gz'))==h
        assert obj['result']['status'].startswith('PASS');return obj['result']
    core=read_data('math_core.json');g=core['g_Q_monic'];rad=core['r'];beta=[Eblock(core['beta_x']),Eblock(core['beta_y'])]
    if phase in ['g','spin-conclusion'] or phase.startswith('spin'):
        import elementary
        local=read_data('local_spin.json.gz')
        if phase=='g':result=elementary.check_g(core,local)
        elif phase=='spin-conclusion':result=elementary.spin_conclusion(core,local,[previous('spin'+str(i)) for i in range(7)])
        else:
            k=int(phase[4:]);assert 0<=k<7;result=elementary.spin_slice(core,local,k)
        inputs['verification/elementary.py']=sha(HERE/'elementary.py')
    elif phase=='conclusion':
        checked={name:previous(name)['status'] for name in ['g','spin-conclusion']+[f'norm{i:02}' for i in range(1,14)]+[f'h{i:02}' for i in range(27,-1,-1)]+['encoding']}
        result={'status':'PASS_SOURCE_ONLY_PUBLISHED_ARITHMETIC_SPINE','checked':checked,'meaning':'Recomputed g real/irreducible, fixed local Spin data, Norm13(u)=beta, P28(beta)=0, and primitive degree56 encoding.','not_reproved_here':['geometric/arithmetic PSL cover theorem','all characteristic-zero root labels','trace-isometry origin beyond exact U check','integral overorder and p-adic norm-separation identification of delta','cubic Norm_K/F(delta)=u and canonical M embedding','Spin central-extension theorem and stabilizer-degree theorem'],'full_inverse_Galois_proof_is_paper_plus_these_separately_identified_premises':True}
    else:
        build=load(work/'build.json.gz');assert build['runner_sha256']==sha(Path(__file__)) and build['header_sha256']==sha(HERE/'exact_V.hpp');deps['build']=sha(work/'build.json.gz')
        for b in build['entries']:assert sha(HERE/b['source'])==b['source_sha256'] and sha(work/b['binary'])==b['binary_sha256']
        if phase.startswith('norm'):
            k=int(phase[4:]);assert 1<=k<=13;desc=read_data('norm13_input.json.gz');u=[desc['A_in_V'],desc['B_in_canonical_w_basis']]
            prev=previous(f'norm{k-1:02}') if k>1 else None;power=prev['power_u_k'] if prev else u;history=prev['history'] if prev else []
            packet=str(k)+'\n'+' '.join(g)+'\n'+' '.join(rad)+'\n'+''.join(emitV(v) for v in u+power+[E_to_V(e) for e in beta])
            for item in history:packet+=''.join(emitV(v) for v in item['trace_u_k']+item['elementary_e_k'])
            result=json.loads(execute([str(work/'norm13')],packet));assert result['stage']==k
            history.append({key:result[key] for key in ['trace_u_k','elementary_e_k']});result['history']=history
        elif phase.startswith('h') or phase=='encoding':
            p=DATA/'f56_Z.txt.gz';inputs[p.name]=sha(p);z=gzip.decompress(p.read_bytes()).decode('ascii').splitlines();assert len(z)==57 and all(z[i]=='0' for i in range(1,57,2));D=z[56];assert int(D)>0
            if phase=='encoding':
                result=json.loads(execute([str(work/'norm28')],'content\n57\n'+'\n'.join(z)+'\n'));result.update({'degree':56,'monic_by_definition':True,'positive_scale':True,'odd_coefficients_zero':True,'coefficients':57})
            else:
                k=int(phase[1:]);assert 0<=k<=27
                if k==27:pair=[{'denominator':'1','numerator_coefficients_a_ascending':['1']+['0']*13},{'denominator':'1','numerator_coefficients_a_ascending':['0']*14}]
                else:pair=previous(f'h{k+1:02}')['residual_pair']
                packet='horner\n'+str(k)+'\n'+' '.join(g)+'\n'+' '.join(rad)+'\n'+''.join(emitE(e) for e in beta+pair)+z[2*k]+'/'+D+'\n'
                result=json.loads(execute([str(work/'norm28')],packet));assert result['coefficient_index']==k
        else:raise ValueError('Unknown phase; use list')
    record={'phase':phase,'result':result,'inputs_sha256':inputs,'dependencies_sha256':deps,'runner_sha256':sha(Path(__file__)),'seconds':time.monotonic()-began,'native_cap_seconds':35,'wrapper_cap_seconds':42,'legacy_archive_imported':False,'historical_binary_hash_required':False}
    put(dest,record);print(json.dumps({'phase':phase,'status':result['status'],'seconds':record['seconds'],'receipt_sha256':sha(dest)}))
if __name__=='__main__':
    signal.alarm(42)
    try:main()
    except BaseException as exc:
        # A failed run never creates a PASS checkpoint. Preserve an additive log.
        try:
            if '--work' in sys.argv:w=Path(sys.argv[sys.argv.index('--work')+1])
            else:w=HERE/'work'
            if w.exists():put(w/('failure-'+str(time.time_ns())+'.json'),{'status':'FAIL_NOT_PASS','phase':sys.argv[1:2],'exception':repr(exc)})
        except Exception:pass
        raise
