"""Unit validation and source-unit conversion for the optional figure runner."""
from pathlib import Path
import importlib.util,tempfile,unittest

script=Path(__file__).resolve().parents[1]/'scripts/reproduce_figures.py'
spec=importlib.util.spec_from_file_location('runner',script)
runner=importlib.util.module_from_spec(spec);spec.loader.exec_module(runner)

CARD='''proc = H
collider = pp
sqrts = 13600
M = 125
pdf_name = NNPDF40_nnlo_as_01180.LHgrid
xmur = 0.5
xmuf = 0.5
xQ = 0.5
xsct_lo = 16000
delta_xsct_nlo = 21000
delta_xsct_nnlo = 9000
delta_xsct_nnnlo = 1600
10 -20000 -3000 -100
30 -11000 -7000 200
125 -100 -10 1
'''

class FigureInputs(unittest.TestCase):
    def convert(self,text,units='fb'):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)/'card.fxd';p.write_text(text)
            return runner.convert_card(p,(.5,.5,.5,.5),units)
    def test_units(self):
        result,x=self.convert(CARD)
        parameters={l.split('=',1)[0].strip():l.split('=',1)[1].strip() for l in result.splitlines() if '=' in l}
        self.assertEqual(float(parameters['xsct_lo']),.016)
        self.assertEqual(float(parameters['xmuR']),.5)
        self.assertEqual(parameters['pdf_name'],'NNPDF40_nnlo_as_01180')
        self.assertEqual(list(x),[10,30,125])
        rows=[list(map(float,l.split())) for l in result.splitlines() if l.strip() and not l.startswith('#') and '=' not in l]
        self.assertAlmostEqual(rows[1][1],-.011)
        nb,_=self.convert(result,'nb')
        self.assertEqual(nb,result)
    def test_rejections(self):
        for modified in [CARD.replace('13600','13000'),CARD.replace('xmur = 0.5','xmur = 1.0'),
                         CARD.replace('30 -11000','25 -11000'),CARD+'R = 0.7\n',
                         CARD.replace('16000','nan'),CARD+'M = 125\n']:
            with self.subTest(card=modified),self.assertRaises(ValueError):self.convert(modified)

if __name__=='__main__':unittest.main()
