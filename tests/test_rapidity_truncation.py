"""Opt-in radius truncation: routing, matching expansion and cutoff regressions."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json, math, os, subprocess

from test_support import CARD, NROWS
JET=Path(__file__).resolve().parents[1]
OUT=Path(os.environ.get('RAPIDITY_TEST_OUTPUT',JET/'build/test-results/rapidity'))
OUT.mkdir(parents=True,exist_ok=True)
FLAG='-truncate-rapidity-rge-at-n3ll'

def case(args):
    profile,fo=args
    results=[]
    for truncate in (False,True):
        name=f'profile{int(profile)}-fo{fo}-truncate{int(truncate)}'
        cmd=[str(JET/'jetvheto'),'-in',str(CARD),'-out',name+'.dat',
             '-order','3','-fixed_order',str(fo),
             '-cross-section',
             '-n3ll-data',str(JET/'data/TMDs_ptj')]
        if profile:cmd+=['-new_modlog','-xM','0.5']
        if truncate:cmd+=[FLAG]
        p=subprocess.run(cmd,cwd=OUT,capture_output=True,text=True,timeout=180)
        (OUT/(name+'.log')).write_text(p.stdout+p.stderr)
        assert p.returncode==0,(name,p.stderr)
        text=(OUT/(name+'.dat')).read_text()
        assert ('strict N3LL rapidity RGE' if truncate else 'RadISH inclusive products') in text
        rows=[[float(v) for v in line.split()] for line in text.splitlines() if line.strip() and not line.startswith('#')]
        assert all(math.isfinite(v) for row in rows for v in row)
        results.append(rows)
    native,strict=results
    assert len(native)==len(strict)==NROWS
    assert any(x[2]!=y[2] for x,y in zip(native,strict) if x[0]<40),'Flag did not reach exponent'
    for x,y in zip(native,strict):
        assert x[0]==y[0]
        assert x[3:]==y[3:],'Fixed order changed'
        if profile and x[0]>62.5: assert x==y,'Cutoff not preserved'
    return f'profile{int(profile)}-fo{fo}',{'native_30':next(x[1] for x in native if x[0]==30),
                                           'strict_30':next(x[1] for x in strict if x[0]==30)}

with ThreadPoolExecutor(max_workers=2) as pool:
    summary=dict(pool.map(case,[(p,f) for p in (False,True) for f in (1,2,3)]))
for order in (0,1,2):
    p=subprocess.run([str(JET/'jetvheto'),'-in',str(CARD),'-out','invalid.dat','-order',str(order),FLAG],cwd=OUT,capture_output=True,text=True)
    assert p.returncode!=0 and FLAG+' requires -order 3' in p.stdout+p.stderr
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print('12 runs passed: both logarithms and all 3 matching orders; unchanged fixed order and cutoff; 3 invalid orders rejected')
