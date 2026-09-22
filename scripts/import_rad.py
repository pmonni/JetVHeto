"""Convert the authors' Wolfram RAD tables to portable Fortran includes.

No Wolfram evaluation is used: only the documented numeric colour polynomials
are accepted. The analytic CCI contribution is added separately at runtime.
"""
import argparse
import hashlib
import math
import re
from pathlib import Path

def split_top(s):
    out, start, stack = [], 0, []
    for i, c in enumerate(s):
        if c in '{[(': stack.append(c)
        elif c in '}])':
            if not stack or stack.pop() != dict(zip('}])', '{[('))[c]:
                raise ValueError('Unbalanced Wolfram expression')
        elif c == ',' and not stack:
            out.append(s[start:i]); start = i + 1
    if stack: raise ValueError('Unbalanced Wolfram expression')
    out.append(s[start:])
    return out

def polynomial(s, powers):
    s = s.replace('`', '').replace('*^', 'e')
    pattern = re.compile(r'\s*([+-]?)\s*((?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?)\s+((?:[A-Za-z]+(?:\^\d+)?\s*)+)')
    result = {tuple(sorted(p.items())): 0.0 for p in powers}
    pos, seen = 0, set()
    while pos < len(s.rstrip()):
        match = pattern.match(s, pos)
        if match is None: raise ValueError(f'Unparsed polynomial at: {s[pos:]}')
        sign, val, monomial = match.groups()
        if pos and not sign: raise ValueError('Missing term separator')
        pos = match.end()
        p = {}
        for name, power in re.findall(r'([A-Za-z]+)(?:\^(\d+))?', monomial):
            if name in p: raise ValueError(f'Duplicate symbol: {name}')
            p[name] = int(power or 1)
        key = tuple(sorted(p.items()))
        if key not in result: raise ValueError(f'Unexpected monomial: {monomial}')
        if key in seen: raise ValueError(f'Duplicate monomial: {monomial}')
        seen.add(key)
        result[key] = float(sign + val)
        if not math.isfinite(result[key]): raise ValueError('Nonfinite coefficient')
    if len(seen) != len(powers): raise ValueError('Incomplete colour polynomial')
    return [result[tuple(sorted(p.items()))] for p in powers]

def read_tables(path):
    src = path.read_text()
    basis = [dict(CA=2,CR=1), dict(CA=1,CR=2), dict(CR=3),
             dict(CA=1,CR=1,nf=1,TR=1), dict(CF=1,CR=1,nf=1,TR=1),
             dict(CR=2,nf=1,TR=1), dict(CR=1,nf=2,TR=2)]
    tables = {}
    for alg in ['AKT','KT','CA']:
        line = next(l for l in src.splitlines() if l.startswith(r'\[CapitalDelta]\[Gamma]'+alg+'=') or l.startswith(r'\[CapitalDelta]\[Gamma]'+alg+' '))
        body = line.split('=',1)[1].strip().rstrip(';')[1:-1]
        rows=[]
        for row in split_top(body):
            fields=split_top(row.strip()[1:-1])
            if len(fields)!=4: raise ValueError('Expected radius, value, error, relative error')
            radius=float(fields[0].replace('`',''))
            values=polynomial(fields[1],basis)
            err=fields[2].strip()
            if not err.startswith(r'\[Sqrt]('): raise ValueError(err)
            variances=polynomial(err[len(r'\[Sqrt]('):-1],[{k:2*v for k,v in p.items()} for p in basis])
            if not math.isfinite(radius) or radius<=0: raise ValueError('Invalid radius')
            if any(v<0 for v in variances): raise ValueError('Negative variance')
            rows.append((radius,values,variances))
        tables[alg]=rows
    radii=[row[0] for row in tables['AKT']]
    if len(radii)<4 or any(a>=b for a,b in zip(radii,radii[1:])):
        raise ValueError('Grid must contain at least four strictly increasing radii')
    if any([row[0] for row in rows]!=radii for rows in tables.values()):
        raise ValueError('Algorithm radius grids differ')
    return tables

def main():
    p=argparse.ArgumentParser(); p.add_argument('input',type=Path); p.add_argument('output',type=Path)
    a=p.parse_args(); tables=read_tables(a.input)
    lines=['! Generated from RapidityAnomalousDimension.wl', '! SHA256 '+hashlib.sha256(a.input.read_bytes()).hexdigest()]
    sizes={len(rows) for rows in tables.values()}
    if len(sizes)!=1: raise ValueError('Algorithm grids differ in length')
    n=sizes.pop(); lines.append(f'integer, parameter :: rad_grid_size={n}')
    def array(name, shape, vals):
        dims='('+','.join(str(d) for d in shape)+')'
        lines.append(f'real(dp), parameter :: {name}{dims} = reshape([ &')
        lines.extend('  '+', '.join(f'{v:.17e}_dp' for v in vals[i:i+3])+', &' for i in range(0,len(vals),3))
        lines[-1]=lines[-1].replace(', &',' &')
        lines.append(f'], {str(list(shape)).replace("(","[").replace(")","]")})')
    # Emit dimensions explicitly to keep the generated include compiler independent.
    array('rad_radii',(n,),[r[0] for r in tables['AKT']])
    array('rad_values',(7,n,3),[x for rows in tables.values() for r in rows for x in r[1]])
    array('rad_variances',(7,n,3),[x for rows in tables.values() for r in rows for x in r[2]])
    a.output.write_text('\n'.join(lines)+'\n')

if __name__=='__main__': main()
