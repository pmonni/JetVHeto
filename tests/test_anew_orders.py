"""All supported accuracy combinations, both logarithms, and removed aliases."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import subprocess,math,json,os
from test_support import CARD, NROWS
JET=Path(__file__).resolve().parents[1]
OUT=Path(os.environ.get('ANEW_TEST_OUTPUT',JET/'build/test-results/anew'));OUT.mkdir(parents=True,exist_ok=True)
def case(args):
    order,fo,profile=args;name=f'res{order}-fo{fo}-profile{int(profile)}'
    cmd=[str(JET/'jetvheto'),'-in',str(CARD),'-out',name+'.dat','-order',str(order),'-fixed_order',str(fo),'-scheme','anew','-cross-section','-n3ll-data',str(JET/'data/TMDs_ptj')]
    if profile:cmd+=['-new_modlog','-xM','0.5']
    p=subprocess.run(cmd,cwd=OUT,text=True,capture_output=True,timeout=180)
    (OUT/(name+'.log')).write_text(p.stdout+p.stderr)
    assert p.returncode==0,(name,p.stdout[-2000:],p.stderr[-2000:])
    data=(OUT/(name+'.dat')).read_text()
    rows=[[float(x) for x in l.split()] for l in data.splitlines() if l.strip() and not l.startswith('#')]
    assert len(rows)==NROWS and all(math.isfinite(v) for r in rows for v in r),name
    assert '# matching scheme = anew' in data
    if profile and not (order==3 and fo==1):
        assert all(math.isclose(r[1],r[3],rel_tol=2e-9,abs_tol=1e-12) for r in rows if r[0]>=62.5),('matching cutoff',name)
    return name,next(r[1]*1000 for r in rows if r[0]==30)
with ThreadPoolExecutor(max_workers=2) as pool:
    result=dict(pool.map(case,[(o,f,p) for o in range(4) for f in range(1,4) for p in [False,True]]))
for scheme in ['draft','a','b','c','d','moda','modRa','R','loga','cunx']:
    p=subprocess.run([str(JET/'jetvheto'),'-in',str(CARD),'-out','removed-scheme.dat','-scheme',scheme],cwd=OUT,text=True,capture_output=True)
    assert p.returncode!=0 and 'Only matching scheme anew' in p.stderr+p.stdout,scheme
(OUT/'matrix.json').write_text(json.dumps(result,indent=2)+'\n')
print('24 accuracy/logarithm combinations passed; 10 removed schemes rejected')
