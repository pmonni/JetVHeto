"""CLI integration regressions. Recorded numbers are not independent physics benchmarks."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import math
import subprocess
import tempfile

from test_support import CARD, PDF
ROOT=Path(__file__).resolve().parents[1]
BASE=['-resum-only','-pdf_name',PDF,'-nbins','3','-ptmin','20','-ptmax','40',
      '-order','3','-cross-section']

def run_case(case):
    name,args,expected_error=case
    with tempfile.TemporaryDirectory() as tmp:
        output=Path(tmp)/'result.dat'
        if name=='n3ll-long-output':
            output=Path(tmp)/('x'*120)/('y'*120)/'result.dat'
            output.parent.mkdir(parents=True)
        result=subprocess.run([str(ROOT/'jetvheto'),*args,'-out',str(output)],cwd=ROOT,
                              text=True,capture_output=True,timeout=180)
        if expected_error:
            assert result.returncode!=0,(name,'unexpected success')
            assert expected_error in result.stdout+result.stderr,(name,result.stderr)
            return name,None
        assert result.returncode==0,(name,result.stdout,result.stderr)
        content=output.read_text()
        rows=[[float(x) for x in line.split()] for line in content.splitlines()
              if line.strip() and not line.startswith('#')]
        assert rows and all(math.isfinite(x) for row in rows for x in row),(name,rows)
        if name.startswith('n3ll-'):
            assert 'Unprimed N3LL resummation' in content
            assert all(row[1]==0 and row[3]==0 and row[2]>0 for row in rows)
        return name,rows

cases=[('n3ll-'+proc+'-'+alg,BASE+['-proc',proc,'-jet-algorithm',alg],None)
       for proc in ['H','DY'] for alg in ['antikt','kt','ca']]
cases += [
 ('n3ll-long-output',BASE,None),
 ('missing-input',[],'Missing fixed-order input'),
 ('removed-draft',BASE+['-new_modlog','-scheme','draft','-xM','0.5'],'Only matching scheme anew'),
 ('profile-anew',BASE+['-new_modlog','-scheme','anew','-xM','1.0'],None),
 ('bad-xM',BASE+['-new_modlog','-xM','1.5'],'0.5 <= xM <= 1'),
 ('unused-xM',BASE+['-xM','0.5'],'-xM requires'),
 ('unused-p',BASE+['-new_modlog','-p','5'],'-p is not used'),
 ('bad-profile-Q',BASE+['-new_modlog','-Q','1'],'Q >= xM*M/4'),
 ('nan-profile',BASE+['-new_modlog','-xM','NaN'],'nonfinite xM'),
 ('bad-radius',BASE+['-R','0.02'],'radius outside'),
 ('bad-algorithm',BASE+['-jet-algorithm','unknown'],'algorithm must'),
 ('bad-small-r',BASE+['-small-r'],'not supported'),
 ('bad-scheme',BASE+['-scheme','b'],'Only matching scheme anew'),
 ('removed-gate',BASE+['-n3ll-experimental'],'Unrecognized command-line option'),
 ('missing-data',BASE+['-n3ll-data','/nonexistent-jetvheto-data'],'missing boundary data'),
 ('nnll',['-resum-only','-pdf_name',PDF,'-nbins','5','-ptmin','10','-ptmax','60'],None),
 ('nnll-explicit',['-resum-only','-pdf_name',PDF,'-nbins','5','-ptmin','10','-ptmax','60','-scheme','anew'],None),
 ('matched',['-in',str(CARD),
             '-order','3','-fixed_order','3','-cross-section','-ptmax','60'],None),
]
with ThreadPoolExecutor(max_workers=3) as pool:
    results=dict(pool.map(run_case,cases))
assert results['nnll']==results['nnll-explicit'],'Default is not anew'
matched=next(row for row in results['matched'] if row[0]==30)
assert results['n3ll-H-antikt']!=results['n3ll-H-kt'],'Missing algorithm dependence'
print(f'{len(cases)} CLI cases passed; default equals explicit anew; matched value at 30 GeV: {matched[1]*1000:.7f} pb')
