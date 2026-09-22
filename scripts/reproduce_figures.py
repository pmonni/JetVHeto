#!/usr/bin/env python3
"""Run the draft's 40-prediction scan and plot Figures 2 and 3.

Requires NumPy, Matplotlib, a built JetVHeto executable, its numerical
dependencies, NNPDF40_nnlo_as_01180, and the external NNLOJET study inputs.
No unpublished input files or personal paths are embedded in this script.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse
import hashlib
import json
import math
import os
import shlex
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
# muF/M, muR/M, initial Q/M, matching scale/M; central is first.
SCALES=[(.5,.5,.5,.5),(1.,1.,.5,.5),(.25,.25,.5,.5),
        (.5,1.,.5,.5),(1.,.5,.5,.5),(.5,.25,.5,.5),(.25,.5,.5,.5),
        (.5,.5,1.,.5),(.5,.5,.25,.5),(.5,.5,.5,1.)]
ACCURACIES=[(1,1),(2,2),(3,2),(3,3)]  # fixed order, resummed order

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def convert_card(path, scales, units):
    """Validate draft settings and convert cross sections, never momenta, to nb."""
    import numpy as np
    f,r,q,_=scales
    parameters={};rows=[]
    for line in path.read_text().splitlines():
        line=line.split('#',1)[0].strip()
        if not line:continue
        if '=' in line:
            key,value=map(str.strip,line.split('=',1))
            key={'xmur':'xmuR','xmuf':'xmuF'}.get(key,key)
            if key in parameters:raise ValueError(f'Duplicate parameter {key}: {path}')
            parameters[key]=value
        else:
            rows.append([float(v) for v in line.split()])
    numeric=dict(M=125,sqrts=13600,xmuF=f,xmuR=r,xQ=q)
    for key,value in numeric.items():
        if key not in parameters or not math.isclose(float(parameters[key]),value,rel_tol=0,abs_tol=1e-10):
            raise ValueError(f'{path}: expected {key}={value}')
    for key,value in [('proc','H'),('collider','pp')]:
        if parameters.get(key)!=value:raise ValueError(f'{path}: expected {key}={value}')
    pdf=parameters.get('pdf_name','').removesuffix('.LHgrid')
    if pdf!='NNPDF40_nnlo_as_01180':raise ValueError(f'{path}: wrong PDF for draft figures')
    parameters['pdf_name']=pdf
    if 'R' in parameters and float(parameters['R'])!=.4:raise ValueError(f'{path}: expected R=0.4')
    parameters['R']='0.4'
    for key in ('xsct_lo','delta_xsct_nlo','delta_xsct_nnlo','delta_xsct_nnnlo'):
        if key not in parameters or not math.isfinite(float(parameters[key])):
            raise ValueError(f'{path}: missing/nonfinite {key}')
    if float(parameters['xsct_lo'])<=0:raise ValueError(f'{path}: nonpositive Born normalization')
    array=np.asarray(rows,dtype=float)
    if array.ndim!=2 or array.shape[1]!=4 or not np.isfinite(array).all():
        raise ValueError(f'{path}: expected finite pT plus three sigmabar columns')
    if not np.all(np.diff(array[:,0])>0) or array[0,0]>10 or array[-1,0]<125 or 30 not in array[:,0]:
        raise ValueError(f'{path}: require increasing pT grid, including 30, covering 10 through 125 GeV')
    factor=1e-6 if units=='fb' else 1.
    converted=['# Converted to nb; momenta remain GeV. Original saved under inputs/.']
    for key,value in parameters.items():
        if 'xsct' in key:value=format(float(value)*factor,'.17g')
        converted.append(f'{key} = {value}')
    array[:,1:]*=factor
    converted.extend(' '.join(format(v,'.17g') for v in row) for row in array)
    return '\n'.join(converted)+'\n',array[:,0]

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input-dir',type=Path,help='Directory containing the nine NNLOJET ptj1_xmuR_*_xmuF_*_xQ_*.fxd files')
    parser.add_argument('--output-dir',type=Path,required=True,help='New or empty destination; no existing results are overwritten')
    parser.add_argument('--input-units',choices=['fb','nb'],default='fb',help='Units of external fixed-order inputs (default fb)')
    parser.add_argument('--executable',type=Path,default=ROOT/'jetvheto')
    parser.add_argument('--data-dir',type=Path,default=ROOT/'data/TMDs_ptj')
    parser.add_argument('--workers',type=int,default=2,help='Independent concurrent processes (default 2)')
    parser.add_argument('--timeout',type=int,default=600,help='Maximum seconds per run')
    parser.add_argument('--truncate-rapidity-rge-at-n3ll',action='store_true',help='Use strict radius evolution for N3LL curves; default retains higher-order coupling terms')
    parser.add_argument('--small-r-nnll',action='store_true',help='Enable small-R resummation only for the NNLL curves, with R0=1')
    parser.add_argument('--plots-only',action='store_true',help='Regenerate plots/tables from the validated manifest in output-dir; no runs')
    args=parser.parse_args(argv)
    if args.workers<1 or args.timeout<1:parser.error('workers and timeout must be positive')
    try:
        import numpy as np
        import matplotlib
    except ImportError as error:
        parser.error(f'{error}; install numpy and matplotlib in the Python environment used to run this script')
    out=args.output_dir.resolve()
    if args.plots_only:
        manifest=json.loads((out/'run_manifest.json').read_text())
        if args.truncate_rapidity_rge_at_n3ll and not manifest['truncate_rapidity_rge_at_n3ll']:
            parser.error('plots-only cannot change the prescription of existing runs')
        if args.small_r_nnll and not manifest.get('small_r_nnll',False):
            parser.error('plots-only cannot enable small-R for existing flag-off runs')
        from paper_figures_plot import plot_figures
        plot_figures(out,manifest['truncate_rapidity_rge_at_n3ll'])
        return
    if args.input_dir is None:parser.error('--input-dir is required unless --plots-only is selected')
    exe=args.executable.resolve();data=args.data_dir.resolve();inputs=args.input_dir.resolve()
    if not exe.is_file() or not os.access(exe,os.X_OK):parser.error('Executable missing/not executable; build JetVHeto first')
    if not data.is_dir() or not (data/'BoundaryConditions').is_dir():parser.error('Missing TMD boundary-data directory')
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        parser.error('Output directory must be empty/new; use a different directory or --plots-only')
    # Preflight every card before creating the output directory or launching runs.
    cards={};grid=None
    for scales in SCALES:
        f,r,q,_=scales
        name=f'ptj1_xmuR_{r}_xmuF_{f}_xQ_{q}.fxd'
        source=inputs/name
        if not source.is_file():parser.error(f'Missing input: {source}')
        converted,x=convert_card(source,scales,args.input_units)
        if grid is not None:np.testing.assert_array_equal(x,grid)
        grid=x;cards[name]=(source,converted)
    out.mkdir(parents=True,exist_ok=True)
    for directory in ('inputs','cards-nb','runs','logs'):(out/directory).mkdir()
    shutil.copy2(exe,out/'jetvheto-used')
    shutil.copytree(data,out/'data-used')
    inputs_manifest={}
    for name,(source,converted) in cards.items():
        shutil.copy2(source,out/'inputs'/name)
        (out/'cards-nb'/name).write_text(converted)
        inputs_manifest[name]=dict(source=str(source),sha256=sha(out/'inputs'/name),card_sha256=sha(out/'cards-nb'/name))
    jobs=[]
    for scales in SCALES:
        f,r,q,m=scales
        card=out/'cards-nb'/f'ptj1_xmuR_{r}_xmuF_{f}_xQ_{q}.fxd'
        for fo,res in ACCURACIES:
            stem=(f'ptj_NNPDF40_nnlo_as_01180_LHC13600_R_0.4_fo{fo}_res{res}'
                  f'_xmuF_{f}_xmuR_{r}_xQ_{q}_xM_{m}_newmodlog_anew_exp_jetveto')
            output=out/'runs'/(stem+'.dat')
            cmd=[str(out/'jetvheto-used'),'-in',str(card),'-out',output.name,
                 '-order',str(res),'-fixed_order',str(fo),'-scheme','anew',
                 '-new_modlog','-xM',str(m),'-cross-section','-n3ll-data',str(out/'data-used')]
            if res==3 and args.truncate_rapidity_rge_at_n3ll:cmd.append('-truncate-rapidity-rge-at-n3ll')
            if res==2 and args.small_r_nnll:cmd+=['-small-r','-R0','1.0']
            jobs.append(dict(command=cmd,cwd=str(out/'runs'),output=str(output),res=res,fo=fo,scales=scales))
    manifest=dict(truncate_rapidity_rge_at_n3ll=args.truncate_rapidity_rge_at_n3ll,
                  small_r_nnll=args.small_r_nnll,
                  input_units=args.input_units,output_units='nb',plot_units='pb',inputs=inputs_manifest,
                  executable_sha256=sha(out/'jetvheto-used'),
                  data_sha256={str(p.relative_to(out/'data-used')):sha(p) for p in sorted((out/'data-used').rglob('*')) if p.is_file()},jobs=jobs)
    (out/'planned_runs.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (out/'commands.sh').write_text('#!/bin/sh\nset -eu\ncd '+shlex.quote(str(out/'runs'))+'\n'+'\n'.join(shlex.join(j['command']) for j in jobs)+'\n')
    def run(job):
        p=subprocess.run(job['command'],cwd=job['cwd'],capture_output=True,text=True,timeout=args.timeout)
        log=out/'logs'/(Path(job['output']).stem+'.log')
        log.write_text(p.stdout+p.stderr)
        if p.returncode:raise RuntimeError(f'Run failed; see {log}')
        path=Path(job['output']);a=np.loadtxt(path)
        if a.shape!=(len(grid),4) or not np.isfinite(a).all():raise ValueError(f'Invalid output: {path}')
        np.testing.assert_array_equal(a[:,0],grid)
        mask=a[:,0]>=125*job['scales'][3]
        np.testing.assert_allclose(a[mask,1],a[mask,3],rtol=2e-9,atol=1e-12)
        if job['res']==3:
            if '# radius exponent = ' not in path.read_text():raise ValueError(f'Wrong radius prescription: {path}')
        job['sha256']=sha(path)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for i,future in enumerate(as_completed([pool.submit(run,j) for j in jobs]),1):
            future.result();print(f'{i}/40 predictions completed',flush=True)
    if sha(out/'jetvheto-used')!=manifest['executable_sha256']:raise RuntimeError('Executable changed during run')
    (out/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    from paper_figures_plot import plot_figures
    plot_figures(out,args.truncate_rapidity_rge_at_n3ll)
    print('30-GeV cross sections and uncertainties: '+str(out/'cross_sections_30GeV.csv'))

if __name__=='__main__':
    try:main()
    except (ValueError,RuntimeError,FileNotFoundError,subprocess.TimeoutExpired) as error:
        print(f'Error: {error}',file=sys.stderr)
        sys.exit(1)
