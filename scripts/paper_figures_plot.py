"""Plotting helper for reproduce_figures.py; no input data are bundled here."""

def plot_figures(output_directory, truncate):
    from pathlib import Path
    import os
    import hashlib
    import json
    import csv

    OUT = Path(output_directory).resolve()
    os.environ.setdefault('MPLCONFIGDIR', str(OUT/'mpl-cache'))
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.ticker import MultipleLocator
    from matplotlib.patches import Patch

    BASE = OUT
    run_manifest=json.loads((BASE/'run_manifest.json').read_text())
    RUN_JOBS=run_manifest['jobs']
    SMALL_R=run_manifest.get('small_r_nnll',False)
    # Tuple order: muF/M, muR/M, initial Q/M, matching scale/M.
    SCALES = [(0.5,0.5,0.5,0.5), (1.,1.,0.5,0.5), (.25,.25,.5,.5),
              (.5,1.,.5,.5), (1.,.5,.5,.5), (.5,.25,.5,.5), (.25,.5,.5,.5),
              (.5,.5,1.,.5), (.5,.5,.25,.5), (.5,.5,.5,1.)]
    MANIFEST = []

    def load(fo, res):
        arrays = []
        for f,r,q,m in SCALES:
            path = BASE / 'runs' / (f'ptj_NNPDF40_nnlo_as_01180_LHC13600_R_0.4_fo{fo}_res{res}'
                           f'_xmuF_{f}_xmuR_{r}_xQ_{q}_xM_{m}_newmodlog_anew_exp_jetveto.dat')
            raw = path.read_bytes()
            header = raw.decode().splitlines()[0]
            assert '-scheme anew' in header and '-new_modlog' in header
            assert f'-order {res}' in header and f'-fixed_order {fo}' in header
            if res==3:
                assert '# radius exponent = ' in raw.decode()
            job=next(j for j in RUN_JOBS if j['output']==str(path))
            assert ('-small-r' in job['command']) == (SMALL_R and res==2)
            assert hashlib.sha256(raw).hexdigest()==job['sha256']
            if res==3:assert ('-truncate-rapidity-rge-at-n3ll' in job['command']) == truncate
            data = np.loadtxt(path)
            assert data.ndim == 2 and data.shape[1] == 4 and np.isfinite(data).all()
            if arrays:
                np.testing.assert_array_equal(data[:,0],arrays[0][:,0])
            arrays.append(data)
            MANIFEST.append({'path':str(path), 'sha256':hashlib.sha256(raw).hexdigest(),
                             'fixed_order':fo,'resummed_order':res,
                             'muF_over_M':f,'muR_over_M':r,'Q_over_M':q,'muM_over_M':m,
                             'command':job['command']})
        return np.stack(arrays)

    RAW = {name:load(fo,res) for name,fo,res in
           [('NLL_NLO',1,1),('NNLL_NNLO',2,2),('NNLL_N3LO',3,2),('N3LL_N3LO',3,3)]}
    X = RAW['N3LL_N3LO'][0,:,0]

    def envelope(a):
        return a[0],a.min(axis=0),a.max(axis=0)

    CURVES = {name:envelope(data[:,:,1]*1000) for name,data in RAW.items()}
    for name,source in [('NLO','NLL_NLO'),('NNLO','NNLL_NNLO'),('N3LO','N3LL_N3LO')]:
        # FO bands have seven unique muR/muF combinations, no Q or muM dependence.
        arr = RAW[source]
        for row in arr[7:]:
            np.testing.assert_allclose(row[:,3],arr[0,:,3],rtol=0,atol=1e-12)
        CURVES[name] = envelope(arr[:7,:,3]*1000)

    index30 = np.where(X==30)[0].item()
    with (OUT/'cross_sections_30GeV.csv').open('w',newline='') as handle:
        writer=csv.writer(handle)
        writer.writerow(['order','fixed_order_pb','fixed_down_percent','fixed_up_percent',
                         'matched_pb','matched_down_percent','matched_up_percent'])
        for n,fo,matched,target in [(1,'NLO','NLL_NLO',24.32),
                                    (2,'NNLO','NNLL_NNLO',28.98),(3,'N3LO','N3LL_N3LO',30.12)]:
            values=[]
            for name in [fo,matched]:
                c,lo,hi=(a[index30] for a in CURVES[name])
                values.extend([c,100*(lo/c-1),100*(hi/c-1)])
            # Fresh predictions are not forced to reproduce archived table values.
            writer.writerow([n,*values])
    with (OUT/'plot_data.csv').open('w',newline='') as handle:
        writer=csv.writer(handle)
        writer.writerow(['pT_GeV','prediction','central_pb','lower_pb','upper_pb'])
        for name,band in CURVES.items():
            for x,c,lo,hi in zip(X,*band):
                writer.writerow([x,name,c,lo,hi])
    (OUT/'input_manifest.json').write_text(json.dumps(MANIFEST,indent=2)+'\n')

    plt.rcParams.update({'font.family':'serif','font.serif':['DejaVu Serif'],
                         'mathtext.fontset':'dejavuserif','font.size':10,
                         'axes.labelsize':11,'legend.fontsize':9,
                         'axes.linewidth':.7,'lines.linewidth':1.25,
                         'hatch.linewidth':.45,'pdf.fonttype':42,'ps.fonttype':42})
    STYLES = {
     'N3LL_N3LO':(r'N$^3$LL+N$^3$LO','#2879b9',r'\\\\','--'),
     'NNLL_N3LO':(r'NNLL+N$^3$LO','#8cbcc8','////','--'),
     'NNLL_NNLO':('NNLL+NNLO','#43a83b','///','-'),
     'NLL_NLO':('NLL+NLO','#9844a7',r'\\','-.'),
     'N3LO':(r'N$^3$LO','#e33a36','....','-.')}
    if SMALL_R:
        for name in ('NNLL_NNLO','NNLL_N3LO'):
            label,color,hatch,ls=STYLES[name]
            STYLES[name]=(label+' + small-R',color,hatch,ls)

    def axes(fig,spec):
        gs=spec.subgridspec(2,1,height_ratios=[2,1],hspace=.055)
        upper=fig.add_subplot(gs[0]); lower=fig.add_subplot(gs[1],sharex=upper)
        for ax in (upper,lower):
            ax.set_xlim(10,100)
            ax.xaxis.set_major_locator(MultipleLocator(10))
            ax.grid(True,color='#b6b6b6',linewidth=.45)
            ax.set_axisbelow(True)
            ax.tick_params(direction='in',top=True,right=True)
        upper.tick_params(labelbottom=False)
        upper.set_ylim(0,52)
        upper.yaxis.set_major_locator(MultipleLocator(10))
        upper.set_ylabel(r'$\sigma(p_T^{\rm veto})$ [pb]')
        lower.set_xlabel(r'$p_T^{\rm veto}$ [GeV]')
        return upper,lower

    def panel(fig,spec,names,reference,ratio_range,legend_location):
        up,down=axes(fig,spec)
        denominator=CURVES[reference][0]
        for name in names:
            label,color,hatch,ls=STYLES[name]
            central,low,high=CURVES[name]
            for ax,norm in [(up,1),(down,denominator)]:
                ax.fill_between(X,low/norm,high/norm,facecolor='none',edgecolor=color,
                                linewidth=.55,hatch=hatch)
                ax.plot(X,central/norm,color=color,linestyle=ls,linewidth=1.4)
        down.set_ylim(*ratio_range)
        down.yaxis.set_major_locator(MultipleLocator(.2 if ratio_range[0]<.7 else .1))
        ratio_label='Ratio to '+STYLES[reference][0]
        if SMALL_R and reference.startswith('NNLL_'):
            ratio_label=ratio_label.removesuffix(' + small-R')+'\n(small-R on)'
        down.set_ylabel(ratio_label,fontsize=10)
        down.axhline(1,color='#555555',linewidth=.45,zorder=0)
        handles=[Patch(facecolor='none',edgecolor=STYLES[n][1],hatch=STYLES[n][2],
                       label=STYLES[n][0]) for n in names]
        anchor={'bbox_to_anchor':(.99,.40)} if 'NLL_NLO' in names else {}
        up.legend(handles=handles,loc=legend_location,frameon=False,handlelength=1.7,
                  **anchor,
                  borderpad=.2,labelspacing=.25)
        up.text(.22,.08,'JetVHeto+NNLOJET \nNNPDF4.0 (NNLO)\n'
                r'13.6 TeV, $pp\to H+X$, $R=0.4$, anti-$k_t$'+'\n'
                r'uncertainties: $\mu_R,\mu_F,\mu_L,\mu_M$ variations',
                transform=up.transAxes,fontsize=8.7,linespacing=1.55)

    fig2=plt.figure(figsize=(6.2,5.3))
    outer=fig2.add_gridspec(1,1,left=.14,right=.97,bottom=.11,top=.98)
    panel(fig2,outer[0],['NNLL_NNLO','NNLL_N3LO','N3LO'],'NNLL_N3LO',(.8,1.2),'center right')
    fig3=plt.figure(figsize=(12.0,5.3))
    outer=fig3.add_gridspec(1,2,left=.075,right=.985,bottom=.11,top=.98,wspace=.19)
    panel(fig3,outer[0],['N3LL_N3LO','N3LO'],'N3LL_N3LO',(.8,1.2),'center right')
    panel(fig3,outer[1],['N3LL_N3LO','NNLL_NNLO','NLL_NLO'],'N3LL_N3LO',(.6,1.2),'center right')
    with PdfPages(OUT/'figures_2_3.pdf',metadata={'Title':'JetVHeto: Figures 2 and 3 ('+('strict N3LL' if truncate else 'RadISH products')+')'}) as pdf:
        for number,fig in [(2,fig2),(3,fig3)]:
            fig.savefig(OUT/f'figure_{number}.pdf')
            fig.savefig(OUT/f'figure_{number}.svg')
            fig.savefig(OUT/f'figure_{number}.png',dpi=180)
            pdf.savefig(fig)
            plt.close(fig)
    print('Figures 2 and 3 written to '+str(OUT), flush=True)
