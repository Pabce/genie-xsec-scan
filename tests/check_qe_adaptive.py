#!/usr/bin/env python3
"""Runtime smoke/regression checks; use a fresh output directory. Requires GENIE."""
import argparse,csv,json,math,pathlib,subprocess,time
p=argparse.ArgumentParser();p.add_argument('--scanner',type=pathlib.Path,required=True);p.add_argument('--output',type=pathlib.Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
base=[str(a.scanner),'--mode','tune','--tune','G18_10a_02_11a','--event-generator-list','EMQE','--probe','11','--target','1000200400','--observable','d2','--diff','Eprime,costheta_l','--fixed','E=0.15','--fixed','costheta_l=0','--scan','Eprime:0.12:0.12:1','--xsec-unit','nb','--components','--fold','auto','--qel-fold-density','adaptive-theta']
def run(name,extra,ok):
 print('START',name,flush=True)
 out=a.output/(name+'.csv');start=time.monotonic()
 with (a.output/(name+'.log')).open('w') as log:r=subprocess.run(base+extra+['--output',str(out)],stdout=log,stderr=subprocess.STDOUT)
 rows=list(csv.DictReader(l for l in out.open() if not l.startswith('#'))) if out.exists() else []
 if ok:assert r.returncode==0 and rows and all(x['status']=='ok' for x in rows)
 else:assert r.returncode!=0 or any(x['status']!='ok' for x in rows)
 print('PASS',name,round(time.monotonic()-start,2),'seconds',flush=True)
 return rows,time.monotonic()-start
normal,elapsed=run('normal',[],True);tight,_=run('tight',['--qel-rel-tol','0.00025'],True)
automatic,_=run('automatic',['--qel-fold-density','auto'],True)
for x,y in zip(normal,automatic):assert x['value']==y['value']
for x,y in zip(normal,tight):assert math.isclose(float(x['value']),float(y['value']),rel_tol=.005,abs_tol=1e-6)
norm_tight,_=run('norm_tight',['--fold-norm-rel-tol','0.0001'],True)
for x,y in zip(tight,norm_tight):assert math.isclose(float(x['value']),float(y['value']),rel_tol=.005,abs_tol=1e-6)
raw,_=run('raw',['--fold-normalization','raw'],True)
assert any(abs(float(x['value'])/float(y['value'])-1)>.02 for x,y in zip(normal,raw) if float(y['value'])>0)
failed,_=run('budget',['--qel-max-eval','15'],False);assert any('nonconverged' in x['message'] for x in failed)
saved=base[:]
base=[x.replace('E=0.15','E=0.24').replace('costheta_l=0','costheta_l=0.5').replace('Eprime:0.12:0.12:1','Eprime:0.19:0.19:1').replace('1000200400','1000060120') for x in base]
carbon,_=run('carbon',[],True)
assert any(x['component']=='total' and float(x['value'])>0 for x in carbon)
carbon_tight,_=run('carbon_tight',['--qel-rel-tol','0.00025'],True)
for x,y in zip(carbon,carbon_tight):assert math.isclose(float(x['value']),float(y['value']),rel_tol=.002,abs_tol=1e-6)
base=saved
run('unsupported_diff',['--diff','W,Q2'],False)
run('invalid_tolerance',['--qel-rel-tol','nan'],False)
run('infinite_tolerance',['--qel-abs-tol','inf'],False)
run('invalid_norm_tolerance',['--fold-norm-rel-tol','nan'],False)
run('invalid_norm_mode',['--fold-normalization','unknown'],False)
run('legacy_requires_raw',['--qel-fold-density','exact-theta'],False)
mec,_=run('mec',['--event-generator-list','EMMEC','--initial-state-fold-samples','8192'],True)
assert any(x['component'].startswith('MEC:') and float(x['value'])>0 for x in mec)
assert '# fold_normalization,event' in (a.output/'mec.csv').read_text()
base=[x.replace('E=0.15','E=0.12').replace('Eprime:0.12:0.12:1','Eprime:0.10:0.10:1') for x in saved]
threshold_mec,_=run('threshold_mec',['--event-generator-list','EMMEC','--initial-state-fold-samples','8192'],True)
assert any(x['component']=='total' and float(x['value'])>0 for x in threshold_mec)
base=saved
assert any(x['status']=='nonconverged' for x in failed)
(a.output/'receipt.json').write_text(json.dumps({'status':'passed','pilot_seconds':elapsed,'checks':['normal_vs_tight','norm_tolerance_stability','raw_differs','MEC_nonzero','120MeV_MEC_converges','budget_not_success','carbon_normal_vs_tight','unsupported_diff_not_success','invalid_tolerance_rejected','legacy_requires_raw']},indent=2))
print('adaptive runtime tests passed')
