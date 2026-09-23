"""Independent numerical evaluation of the supplied Wolfram expressions.

The oracle does not use the importer's colour-polynomial parser. It tokenizes
implicit multiplication and evaluates an allowlisted Python AST, never eval().
"""
import ast
import importlib.util
import io
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import tokenize
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('import_rad', ROOT/'scripts/import_rad.py')
imp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(imp)
SOURCE = ROOT/'data/RapidityAnomalousDimension.wl'


def numeric(expression, values):
    s = expression.replace('`', '').replace('*^', 'e').replace('^', '**')
    s = s.replace(r'\[Pi]', 'pi').replace(r'\[Sqrt]', 'sqrt')
    s = s.replace('[', '(').replace(']', ')')
    ts = [(t.type, t.string) for t in tokenize.generate_tokens(io.StringIO(s).readline)
          if t.type not in (tokenize.NEWLINE, tokenize.ENDMARKER)]
    out = []
    previous = None
    for typ, value in ts:
        starts = typ in (tokenize.NAME, tokenize.NUMBER) or value == '('
        if previous:
            ptyp, pvalue = previous
            ends = ptyp in (tokenize.NAME, tokenize.NUMBER) or pvalue == ')'
            call = pvalue in ('Log', 'Zeta', 'sqrt') and value == '('
            if ends and starts and not call: out.append('*')
        out.append(value)
        previous = typ, value
    def visit(node):
        if isinstance(node, ast.Constant) and type(node.value) in (int,float): return node.value
        if isinstance(node, ast.Name): return values[node.id]
        if isinstance(node, ast.UnaryOp):
            if isinstance(node.op, ast.USub): return -visit(node.operand)
            if isinstance(node.op, ast.UAdd): return visit(node.operand)
        if isinstance(node, ast.BinOp):
            a,b = visit(node.left),visit(node.right)
            if isinstance(node.op,ast.Add): return a+b
            if isinstance(node.op,ast.Sub): return a-b
            if isinstance(node.op,ast.Mult): return a*b
            if isinstance(node.op,ast.Div): return a/b
            if isinstance(node.op,ast.Pow): return a**b
        if isinstance(node, ast.Call) and isinstance(node.func,ast.Name) and len(node.args)==1:
            arg=visit(node.args[0])
            if node.func.id=='Log': return math.log(arg)
            if node.func.id=='sqrt': return math.sqrt(arg)
            if node.func.id=='Zeta': return {3:1.2020569031595942854,5:1.0369277551433699263}[arg]
        raise ValueError(f'Unsupported AST: {ast.dump(node)}')
    return visit(ast.parse(' '.join(out),mode='eval').body)


class TestRAD(unittest.TestCase):
    def test_regeneration(self):
        with tempfile.TemporaryDirectory() as tmp:
            generated=Path(tmp)/'rad_grid.inc'
            subprocess.run([sys.executable,str(ROOT/'scripts/import_rad.py'),str(SOURCE),str(generated)],check=True)
            self.assertEqual(generated.read_bytes(),(ROOT/'src/rad_grid.inc').read_bytes())

    def test_parser_rejects_malformed_polynomials(self):
        basis=[dict(CA=1),dict(CR=1)]
        for expression in ['1 CA + 2 CA','1 CA + 2 CR + 3','1 CA + 2 CR junk',
                           '1 CA + 2 CR; Print[1]','1e999 CA + 2 CR']:
            with self.subTest(expression=expression), self.assertRaises(ValueError):
                imp.polynomial(expression,basis)

    def test_input_rejections(self):
        for mode in ['low','high','nan','bad-algorithm']:
            result=subprocess.run([str(ROOT/'tests/test_rad'),mode],capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('RAD:',result.stderr)

    def test_all_grid_points_against_original_expressions(self):
        lines=SOURCE.read_text().splitlines()
        expressions={}
        for key in [r'\[Gamma]RSV[2]',r'\[Gamma]RSV[3]',r'\[Gamma]\[Nu][2,R_]',
                    r'\[CapitalDelta]\[Gamma]Analytic[3,R_]']:
            line=next(l for l in lines if l.startswith(key))
            expressions[key]=line.split('=',1)[1].rstrip(';')
        tables={}
        for alg,key in [('antikt','AKT'),('kt','KT'),('ca','CA')]:
            line=next(l for l in lines if l.split('=',1)[0].strip()==r'\[CapitalDelta]\[Gamma]'+key)
            # Only the balanced-list splitter is shared with the converter.
            rows=imp.split_top(line.split('=',1)[1].strip().rstrip(';')[1:-1])
            tables[alg]=[imp.split_top(row.strip()[1:-1]) for row in rows]
        result=subprocess.run([str(ROOT/'tests/test_rad'),'dump'],capture_output=True,text=True,check=True)
        count=0
        for line in result.stdout.splitlines():
            alg,*fields=line.split()
            r,cr,nf,v,err,g2,rsv2,rsv3=map(float,fields)
            values=dict(CA=3,CF=4/3,CR=cr,nf=nf,TR=0.5,R=r,pi=math.pi)
            row=tables[alg][round(r/0.05)-1]
            expected=numeric(row[1],values)+numeric(expressions[r'\[CapitalDelta]\[Gamma]Analytic[3,R_]'],values)
            pairs=[(v,expected),(err,numeric(row[2],values)),
                   (g2,numeric(expressions[r'\[Gamma]\[Nu][2,R_]'],values)),
                   (rsv2,numeric(expressions[r'\[Gamma]RSV[2]'],values)),
                   (rsv3,numeric(expressions[r'\[Gamma]RSV[3]'],values))]
            for actual,expected in pairs:
                self.assertTrue(math.isclose(actual,expected,rel_tol=3e-12,abs_tol=2e-8),
                                (alg,r,cr,nf,actual,expected))
            count+=1
        self.assertEqual(count,360)


if __name__=='__main__': unittest.main()
