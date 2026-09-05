"""Summarize an actual fresh-directory test without publishing machine paths.

This records tests already executed by run.py. It is not a substitute for
executing the arithmetic. No historical archive is read or imported.
"""
from pathlib import Path
import argparse,gzip,hashlib,json,platform,sys
HERE=Path(__file__).resolve().parent;sys.set_int_max_str_digits(0)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(gzip.decompress(p.read_bytes()))
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--tested-root',type=Path,required=True);ap.add_argument('--work-name',default='fresh-work');ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    root=a.tested_root.resolve();work=root/a.work_name;assert root!=HERE.parent and not a.output.exists()
    names=['build','g']+['spin'+str(i) for i in range(7)]+['spin-conclusion']+[f'norm{i:02}' for i in range(1,14)]+[f'h{i:02}' for i in range(27,-1,-1)]+['encoding','conclusion']
    sources=['run.py','elementary.py','exact_V.hpp','norm13.cpp','norm28.cpp','CMakeLists.txt'];source_hashes={n:sha(root/'verification'/n) for n in sources}
    assert all(sha(HERE/n)==h for n,h in source_hashes.items()),'Changed source after test'
    rows=[];inputs={};native_seconds=[]
    for name in names:
        p=work/(name+'.json.gz');d=load(p);assert d['runner_sha256']==source_hashes['run.py']
        if name=='build':assert d['status']=='PASS_FRESH_SOURCE_BUILD';status=d['status']
        else:
            result=d['result'];assert result['status'].startswith('PASS');status=result['status']
            for n,h in d['inputs_sha256'].items():assert sha(root/n)==h;inputs[n]=h
            for n,h in d['dependencies_sha256'].items():assert sha(work/(n+'.json.gz'))==h
            if 'seconds' in result:native_seconds.append(result['seconds'])
            if 'native_seconds' in result:native_seconds.append(result['native_seconds'])
        rows.append({'phase':name,'status':status,'checkpoint_sha256':sha(p),'checkpoint_bytes':p.stat().st_size,'wall_seconds':d['seconds']})
    assert all(sha(HERE.parent/n)==h for n,h in inputs.items()),'Changed public data after test'
    assert load(work/'norm13.json.gz')['result']['status']=='PASS_RECOMPUTED_NORM13_EQUALS_PUBLISHED_BETA'
    assert load(work/'h00.json.gz')['result']['zero_coordinates']==28
    assert load(work/'encoding.json.gz')['result']['content']==1
    assert not list(work.glob('failure-*.json')),'A failure must be reported separately, not hidden'
    out={'schema':'SL2F13-fresh-directory-source-only-test-v1','status':'PASS_FRESH_DIRECTORY_SELECTED_ARITHMETIC_SPINE','tested_source_sha256':source_hashes,'build_interface_tested':'run.py build; the alternative CMake configuration was not executed','tested_input_sha256':inputs,'phases':rows,'phase_count':len(rows),'native_max_seconds':max(native_seconds),'phase_max_wall_seconds':max(r['wall_seconds'] for r in rows),'sum_phase_wall_seconds':sum(r['wall_seconds'] for r in rows),'platform':{'system':platform.system(),'machine':platform.machine(),'python':platform.python_version()},'legacy_archive_used':False,'precompiled_archive_executables_used':False,'expected_intermediate_beta_or_Horner_values_loaded':False,'fresh_directory_test':True,'independent_second_machine_test':False,'runtime_checkpoints_shipped_as_required_inputs':False,'phase_results':{'g_real_roots':14,'g_irreducible_over_Q':True,'trace_U_exact':True,'Spin_integer_determinants':1092,'Spin_coset_norms':28,'new_norm13_steps':13,'new_Horner_steps':28,'P28_beta_exact_zero_coordinates':28,'published_Z_content':1},'remaining_proof_layers':'See the paper and generic exact field/cover reference; this test alone is not full inverse-Galois certification.'}
    with a.output.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
    print(json.dumps({'status':out['status'],'sha256':sha(a.output),'phases':len(rows),'maximum_seconds':out['phase_max_wall_seconds']}))
if __name__=='__main__':main()
