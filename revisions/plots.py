"""Revision plots that reuse Figure 2 and Figure 4's original plotting code."""
from .common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap, BoundaryNorm
from .manuscript_style import figure2, figure4


def save(fig, out, name, tight=False):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    for suffix in ['png', 'pdf', 'svg']:
        fig.savefig(out/f'{name}.{suffix}', transparent=True, bbox_inches='tight' if tight else None)
    plt.close(fig)


def _bottom_legend(fig, handles, labels, ncol=5, anchor=.015):
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(.5, anchor),
               ncol=ncol, frameon=True, fancybox=False, edgecolor='black',
               facecolor='white', handlelength=2.2, columnspacing=1.8, borderpad=.8)


def _expand_shared_limits(axes, lower=None):
    """Keep manuscript limits unless the new data would be clipped."""
    old = [ax.get_ylim() for ax in axes]
    values = np.array([ax.dataLim.intervaly for ax in axes]).ravel()
    values = values[np.isfinite(values)]
    lo, hi = min(v[0] for v in old), max(v[1] for v in old)
    if len(values):
        pad = .08 * max(np.ptp(values), hi-lo)
        if values.min() < lo: lo = values.min()-pad
        if values.max() > hi: hi = values.max()+pad
    if lower is not None: lo = lower
    for ax in axes: ax.set_ylim(lo, hi)


def plot_panel(directory, mode):
    directory = Path(directory)
    d = pd.read_csv(directory/'traces.csv', low_memory=False)
    metadata = json.loads((directory/'metadata.json').read_text())
    expected = metadata['n'] * (4 if mode=='concave' else len(metadata['variants']))
    if len(d)!=expected or d.duplicated(['source_row', 'method']).any():
        raise ValueError('Simulation output is incomplete or has duplicate trace-method pairs')
    groups = {'SI_concave': (METHODS, LABELS)} if mode=='concave' else {
        'SI_derpsalsa_p': (['default','p_0.001','p_0.1'], ['Reference (p = 0.01)','p = 0.001','p = 0.1']),
        'SI_derpsalsa_k': (['default','k_0.5x','k_2x'], ['Reference','k = 0.5×','k = 2×']),
        'SI_derpsalsa_numerics': (['default','diff_order_1','diff_order_3','max_iter_100','tol_0.0001'],
                                ['Reference','Order = 1','Order = 3','100 iter.',r'Tol. = $10^{-4}$']),
        'SI_derpsalsa_smoothing': (['default','smooth_half_window_0','num_smooths_8','num_smooths_32'],
                                 ['Reference','Half-window = 0','Passes = 8','Passes = 32']),
    }
    for name, (methods, labels) in groups.items():
        pairs = [(m,l) for m,l in zip(methods,labels) if m in d.method.unique()]
        if not pairs: continue
        methods, labels = map(list, zip(*pairs))
        ns = figure2(methods, labels)
        if mode=='concave':
            original = ns['make_baseline_examples']
            def inverted_examples():
                x, signal, noisy = original()
                for values in noisy.values():
                    bg = values['baseline'].copy()
                    values['baseline'] = -bg
                    values['signal_plus_baseline'] -= 2*bg
                    values['noisy'] -= 2*bg
                return x, signal, noisy
            ns['make_baseline_examples'] = inverted_examples
        grouped = d[d.method.isin(methods)].groupby(['method','psnr_db','baseline']).bias
        stats = pd.concat({'mean_err': grouped.mean().unstack('baseline'),
                           'std_err': grouped.std().unstack('baseline')}, axis=1)
        fig = ns['plot_fig_2'](stats_df=stats)
        if mode=='derpsalsa':
            # Open reference markers reveal coincident variants inside the blue
            # ring; a dashed foreground line reveals their solid lines beneath.
            reference_label = labels[methods.index('default')]
            for ax in fig.axes:
                for line in ax.lines:
                    if line.get_label() == reference_label:
                        line.set_zorder(10)
                        line.set_linestyle('--')
                        line.set_markerfacecolor('none')
                        line.set_markeredgecolor(ns['METHOD_STYLE']['default']['color'])
                        line.set_markeredgewidth(1.5)
                        line.set_markersize(line.get_markersize()+3)
        # Original axes order: five examples, bias bar, bias curve, SD bar, SD curve.
        _expand_shared_limits([fig.axes[5], fig.axes[6]])
        _expand_shared_limits([fig.axes[7], fig.axes[8]], lower=0)
        if mode=='concave':
            _expand_shared_limits(fig.axes[:5])
            for ax in fig.axes[:5]: ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
        for legend in list(fig.legends): legend.remove()
        handles, legend_labels = ns['make_method_legend_handles']()
        if mode=='derpsalsa':
            reference = handles[methods.index('default')]
            reference.set_linestyle('--')
            reference.set_markerfacecolor('none')
            reference.set_markeredgecolor(ns['METHOD_STYLE']['default']['color'])
            reference.set_markeredgewidth(1.5)
            reference.set_markersize(reference.get_markersize()+3)
        handles.append(Line2D([0],[0],color='black',lw=1.,ls=(0,(4,2))))
        _bottom_legend(fig,handles,legend_labels+['Ground Truth'],ncol=len(methods)+1)
        save(fig, directory, name, tight=True)


def plot_atp(raw, result, summary, out, known=False):
    if not known: raise ValueError('Figure 4-style ATP plotting requires concentration mapping')
    ns = figure4()
    aliases = dict(zip(METHODS, ['aswift','poly_linear','plaxco_gauss','justine_ewd']))
    methods = dict(zip(aliases.values(), LABELS))
    freqs = sorted(raw.hz.unique())
    if len(freqs)!=3: raise ValueError('ATP layout expects three frequencies')
    fig, axes = plt.subplots(3, 3, figsize=(10.2, 9.2), sharex=False)
    for row, hz in enumerate(freqs):
        primary = raw[(raw.hz==hz) & raw.num.between(5,39)].copy()
        # Final scan of each initial-series concentration, matching Figure 4's examples.
        selected = primary.sort_values('num').groupby('concentration_M', sort=True).tail(1).sort_values('concentration_M')
        values = selected.concentration_M.to_numpy()*1e6
        cmap = ListedColormap(plt.colormaps['turbo'](np.linspace(0,1,len(values))))
        norm = BoundaryNorm(np.arange(len(values)+1)-.5, cmap.N)
        ax = axes[row,0]
        for j,r in enumerate(selected.itertuples()):
            ax.plot(r.voltage,r.current,color=cmap(j),linewidth=1.,alpha=.95)
        sm = matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap)
        sm.set_array([])
        ns['_add_inset_colorbar'](fig,ax,sm,title='',tick_values=np.arange(len(values)),
                                 tick_labels=[ns['_value_tick_label'](v) for v in values])
        ax.set_ylabel(fr'Current, $\mu A$ ({hz:g} Hz)')
        ax.set_xlabel(r'Potential, $V$')

        # Call original Figure 4 row renderers with the already-normalized ATP
        # stage summaries. This keeps the exact markers, offsets, opacity and axes.
        def summary_adapter(_df, *, freq_hz, **_kwargs):
            stats = {}
            for method, alias in aliases.items():
                stage = active_summary[(active_summary.hz==freq_hz)&(active_summary.method==method)].sort_values('concentration_M')
                x = stage.concentration_M.to_numpy()*1e6
                y, sd = stage['mean'].to_numpy(), stage.sd.to_numpy()
                stats[alias] = dict(x=x, xplot=np.where(x>0,x*10**ns['METHOD_STYLE'][alias]['offset_decades'],x),
                    mean=y,std=sd,keep=np.isfinite(y))
            return x, stats
        ns['_compute_gain_stats'] = summary_adapter
        active_summary = summary[summary.stage.str.startswith(('baseline','dose'))]
        ns['_plot_gain_row'](axes[row,1],df=pd.DataFrame(),freq_hz=hz,methods=methods,
                             show_uncertainty=False,filter_bad_linear_means=False)
        ns['_plot_gain_std_row'](axes[row,2],df=pd.DataFrame(),freq_hz=hz,methods=methods,
                                 filter_bad_linear_means=False,std_log_scale=True)
        # Rechallenge remains separate, with the original circular markers and a
        # dashed line. Washout crosses and buffer-test scans are not displayed.
        for method in [m for m in METHODS if m!='aswift']+['aswift']:
            alias=aliases[method];st=ns['METHOD_STYLE'][alias]
            s=summary[(summary.hz==hz)&(summary.method==method)&summary.stage.str.startswith('rechallenge')].sort_values('concentration_M')
            x=s.concentration_M.to_numpy()*1e6*10**st['offset_decades']
            for ax,y in [(axes[row,1],s['mean']), (axes[row,2],s.sd)]:
                ax.plot(x,y,color=st['color'],linestyle='--',linewidth=st['linewidth'],
                    marker=st['marker'],markersize=st['markersize'],markerfacecolor=st['color'],
                    markeredgecolor='black',markeredgewidth=.45,
                    alpha=1. if method=='aswift' else .55,zorder=st.get('zorder',3)+1)
    handles=[]
    for alias,label in methods.items():
        st=ns['METHOD_STYLE'][alias]
        handles.append(Line2D([0],[0],color=st['color'],linestyle=st['linestyle'],linewidth=st['linewidth'],
            marker=st['marker'],markersize=st['markersize'],markerfacecolor=st['color'],
            markeredgecolor='black',markeredgewidth=.45))
    _bottom_legend(fig,handles,LABELS,ncol=4,anchor=.02)
    fig.subplots_adjust(top=.96,bottom=.14,left=.08,right=.98,wspace=.40,hspace=.40)
    save(fig,out,'SI_ATP')


def plot_stress(directory):
    from .stress import make_signals, TAIL_RATIOS
    directory=Path(directory)
    meta=json.loads((directory/'metadata.json').read_text())
    summary=pd.read_csv(directory/'summary.csv')
    if len(summary)!=20 or not (summary.n==meta['replicates']).all():
        raise ValueError('Incomplete stress-test summary')
    ns=figure2()
    x,peak,backgrounds,signals=make_signals(meta['replicates'],meta['seed'],meta['snr_db'],meta['relative_peak_width'])
    fig=plt.figure(figsize=(12,5.8))
    outer=fig.add_gridspec(2,2,height_ratios=[.95,1.5],hspace=.5,wspace=.30)
    top=outer[0,:].subgridspec(1,5,wspace=.40)
    for i,ratio in enumerate(TAIL_RATIOS):
        ax=fig.add_subplot(top[0,i])
        ax.scatter(x,signals[i,0],s=4,color='black',alpha=.9,linewidths=0)
        ax.plot(x,backgrounds[i],color='black',ls='--',lw=1.,alpha=.8)
        ax.set_title(f'Tail / peak = {ratio:g}')
        ax.set_xticks([-.4,-.2,0.])
        ax.set_ylim(-.2,ratio*1.08+.1)
        if i==0:ax.set_ylabel(r'Current, $\mu A$')
        if i==2:ax.set_xlabel(r'Voltage, $V$',labelpad=0)
    bias_ax=fig.add_subplot(outer[1,0]);sd_ax=fig.add_subplot(outer[1,1])
    bias_ax.axhline(0,color='black',lw=1.,ls=(0,(4,2)),zorder=1)
    for method in METHODS:
        sub=summary[summary.method==method].sort_values('log2_tail_ratio')
        # SD is shown separately; no scaled uncertainty ribbon on mean bias.
        for ax,column in [(bias_ax,'mean_bias'),(sd_ax,'sd')]:
            ns['plot_method_series'](ax,sub.log2_tail_ratio.to_numpy(),sub[column].to_numpy(),
                                     np.zeros(len(sub)),method,use_ribbon=False)
    for ax in [bias_ax,sd_ax]:
        ax.set_xticks(np.arange(5));ax.set_xlim(-.2,4.2)
        ax.set_xlabel(r'$\log_2$(polynomial tail height / peak height)')
    bias_ax.set_ylabel(r'Bias, $\mu A$');sd_ax.set_ylabel('Std. Dev.')
    sd_ax.set_yscale('log')
    handles,labels=ns['make_method_legend_handles']()
    handles.append(Line2D([0],[0],color='black',lw=1.,ls=(0,(4,2))))
    _bottom_legend(fig,handles,labels+['Ground Truth'])
    fig.subplots_adjust(bottom=.20)
    save(fig,directory,'SI_polynomial_stress',tight=True)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser()
    p.add_argument('mode',choices=['derpsalsa','concave'])
    p.add_argument('directory',type=Path)
    a=p.parse_args();plot_panel(a.directory,a.mode)
