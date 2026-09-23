"""NNLL small-R acceptance, endpoint, domain guards and optional binary regression."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import json, math, os, subprocess
from test_support import ROOT, CARD, PDF

OUT=ROOT/'build/test-results/small-r'
OUT.mkdir(parents=True,exist_ok=True)
BASE=['-in',str(CARD),'-scheme','anew','-cross-section']

def run(name,args,exe=None,error=None):
    target=OUT/(name+'.dat')
    base=['-cross-section','-pdf_name',PDF] if '-resum-only' in args else BASE
    proc=subprocess.run([str(exe or ROOT/'jetvheto'),*base,*args,'-out',str(target)],
                        cwd=ROOT,text=True,capture_output=True,timeout=180)
    (OUT/(name+'.log')).write_text(proc.stdout+proc.stderr)
    if error:
        assert proc.returncode!=0 and error in proc.stdout+proc.stderr,(name,proc.stdout,proc.stderr)
        return
    assert proc.returncode==0,(name,proc.stdout,proc.stderr)
    rows=[[float(v) for v in line.split()] for line in target.read_text().splitlines()
          if line.strip() and not line.startswith('#')]
    assert rows and all(math.isfinite(v) for row in rows for v in row)
    return rows

def positive(case):
    fo,profile=case
    args=['-order','2','-fixed_order',str(fo)]
    if profile: args+=['-new_modlog','-xM','0.5']
    name=f'fo{fo}-profile{int(profile)}'
    off=run(name+'-off',args)
    on=run(name+'-on',args+['-small-r'])
    assert on!=off
    if profile:
        for a,b in zip(on,off):
            if a[0]>62.5:
                assert a==b,('endpoint flag dependence',a,b)
                assert math.isclose(a[1],a[3],rel_tol=2e-9),('fixed-order endpoint',a)
    # Exercise R0 (not a claim of prediction validity for the changed R0).
    varied=run(name+'-R0',args+['-small-r','-R0','0.8'])
    assert varied!=on
    return name,{'off':next(r[1] for r in off if r[0]==30),
                 'on':next(r[1] for r in on if r[0]==30)}

with ThreadPoolExecutor(max_workers=2) as pool:
    results=dict(pool.map(positive,[(fo,p) for fo in (1,2,3) for p in (False,True)]))
for order in (0,1,3):
    run(f'reject-order-{order}',['-order',str(order),'-small-r'],error='only at NNLL')
for name,args,error in [
    ('ln2z',['-small-r-ln2z'],'not supported'),
    ('ln2z-with-small-r',['-small-r','-small-r-ln2z'],'not supported'),
    ('negative-R',['-resum-only','-small-r','-R','-0.4'],'0 < R <= R0'),
    ('bad-R0',['-small-r','-R0','0.1'],'0 < R <= R0'),
    ('nan-R0',['-small-r','-R0','NaN'],'nonfinite'),
    ('pole',['-small-r','-R0','1e30'],'Landau pole'),
    ('algorithm',['-small-r','-jet-algorithm','kt','-fixed-order-algorithm','kt'],'requires ptj and antikt'),
]:
    run('reject-'+name,['-order','2','-fixed_order','2',*args],error=error)

baseline=os.environ.get('SMALL_R_BASELINE')
if baseline:
    def regression(case):
        order,fo,profile=case
        args=['-order',str(order),'-fixed_order',str(fo)]
        if profile:args+=['-new_modlog','-xM','0.5']
        name=f'regression-res{order}-fo{fo}-p{int(profile)}'
        old=run(name+'-old',args,exe=baseline)
        new=run(name+'-new',args)
        assert old==new,('flag-off regression',name)
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(regression,[(o,f,p) for o in range(4) for f in (1,2,3) for p in (False,True)]))
    print('All 24 flag-off accuracy/logarithm combinations unchanged against baseline binary')
(OUT/'summary.json').write_text(json.dumps(results,indent=2)+'\n')
print('18 NNLL small-R/off/R0 runs and 10 rejection cases passed; numerical values are test fixtures, not predictions')
