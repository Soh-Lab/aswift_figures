"""Build figure-specific audit tables and shared preprocessing diagnostics."""
from .common import *
import argparse


def table_row(name,df,peak_col='peak'):
    peaks=df[peak_col]
    median=peaks/df.median_absolute_residual
    mean=peaks/df.mean_absolute_residual
    idx=median.idxmin()
    return dict(figure=name,n=len(df),missing_saved_peaks=int(peaks.isna().sum()),
        refit_rejections=int(df.status.str.startswith('rejected').sum()),
        refit_other_errors=int((df.status=='fit_error').sum()),
        min_height_over_median=median.min(),min_height_over_mean=mean.min(),
        minimum_source_row=df.loc[idx,'source_row'],
        maximum_saved_refit_difference=float(np.max(np.abs(peaks-df.refit_peak))))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',type=Path,default=ROOT/'revisions/results')
    p.add_argument('--dox-methods',type=Path,required=True)
    p.add_argument('--kana-folder',type=Path,required=True)
    p.add_argument('--structured-folder',type=Path,required=True)
    a=p.parse_args();rows=[]
    d=pd.read_csv(a.results/'fig2/traces.csv');rows.append(table_row('Figure 2',d))
    d=pd.read_csv(a.results/'dox/traces.csv');saved=pd.read_csv(a.dox_methods)
    # Figure 4 defaults: first 45 acquisitions, 10 and 200 Hz; absent channels are not plotted.
    selected=[]
    for _,r in d.iterrows():
        col=f'signal-{r.hz}hz-ch{r.channel}-aswift'
        if r.hz in [10,200] and r.num<45 and col in saved:
            r=r.copy();r['saved_peak']=saved.loc[int(r.num),col];selected.append(r)
    selected=pd.DataFrame(selected)
    selected.to_csv(a.results/'figure4_membership.csv',index=False)
    rows.append(table_row('Figure 4',selected,'saved_peak'))
    d=pd.read_csv(a.results/'kana/traces.csv')
    path=a.kana_folder/(a.kana_folder.name+'_methods_results_10mad.csv');saved=pd.read_csv(path)
    d['saved_peak']=[saved.loc[int(r.num),f'signal-200hz-ch{r.channel}-aswift'] for r in d.itertuples()]
    d.to_csv(a.results/'figure5_membership.csv',index=False)
    for channel,name in [(0,'Figure 5 blood'),(1,'Figure 5 ISF')]:rows.append(table_row(name,d[d.channel==channel],'saved_peak'))
    atp=pd.read_csv(a.results/'atp_audit/traces.csv');rows.append(table_row('ATP channel 0',atp))
    if (a.results/'fig3/traces.csv').exists():
        rows.insert(1,table_row('Figure 3 regenerated',pd.read_csv(a.results/'fig3/traces.csv')))
    else:rows.insert(1,dict(figure='Figure 3',n=np.nan,missing_saved_peaks=np.nan,
        note='Exact raw input sim_binding_curve_mg_ewd.feather absent; checked-in binding raw files do not match summary amplitudes'))
    pd.DataFrame(rows).to_csv(a.results/'figure_noise_table.csv',index=False)
    print(pd.DataFrame(rows).to_string(index=False))
    # Timestamp spacing is not scan duration. Keep per-channel/frequency summaries separate.
    intervals=[]
    for path in a.structured_folder.glob('*.csv'):
        raw=pd.read_csv(path)
        if path.name.startswith('ATP'):raw=raw[raw.channel==0]
        for (hz,ch),g in raw.groupby(['hz','channel']):
            dt=np.diff(g.sort_values('num').time.to_numpy())*3600
            positive=dt[dt>0]
            intervals.append(dict(dataset=path.stem,hz=hz,channel=ch,n=len(g),
                minimum_positive_interval_s=positive.min() if len(positive) else np.nan,
                median_positive_interval_s=np.median(positive) if len(positive) else np.nan,
                zero_or_negative_intervals=int((dt<=0).sum())))
    pd.DataFrame(intervals).to_csv(a.results/'acquisition_intervals.csv',index=False)
    # Supplementary Table S2: match shared masking and full-trace membership.
    from peak_extraction.extract_peaks import mask_multichannel_outliers
    raw=structured(a.structured_folder/'kana_invivo_dual_structured.csv')
    proportions={}
    for cutoff in [10,20]:
        vals=[]
        for _,g in raw.groupby(['hz','num']):
            x=g.iloc[0].voltage;xx,_=mask_multichannel_outliers(x,np.stack(g.current),cutoff=cutoff)
            vals.append(1-len(xx)/len(x))
        proportions[cutoff]=np.mean(vals)
    scatter=[]
    for cutoff,suffix in [(None,'inf'),(10,'10'),(20,'20')]:
        path=a.kana_folder/(a.kana_folder.name+f'_methods_results_{suffix}mad.csv')
        saved=pd.read_csv(path)
        for ch in [0,1]:
            nums=raw[raw.channel==ch].num.astype(int).tolist()
            for method in ['aswift','poly_linear','plaxco_gauss','justine_ewd']:
                values=saved.loc[nums,f'signal-200hz-ch{ch}-{method}'].to_numpy()
                reference=np.mean(values[:5]);norm=values/reference
                # Do not bridge failed fits: NaNs invalidate both adjacent differences.
                delta=np.diff(norm);finite=delta[np.isfinite(delta)]
                scatter.append(dict(matrix='blood' if ch==0 else 'ISF',method=method,
                    cutoff='none' if cutoff is None else cutoff,n=len(values),missing=int(np.isnan(values).sum()),
                    fraction_masked=0 if cutoff is None else proportions[cutoff],
                    consecutive_difference_pairs=len(finite),
                    scatter_sd_diff_over_sqrt2=np.std(finite,ddof=1)/np.sqrt(2)))
    pd.DataFrame(scatter).to_csv(a.results/'in_vivo_scatter.csv',index=False)
    save_json(a.results/'tables_metadata.json',provenance([a.dox_methods,
        a.kana_folder/(a.kana_folder.name+'_methods_results_10mad.csv'),
        a.kana_folder/(a.kana_folder.name+'_methods_results_20mad.csv'),
        a.kana_folder/(a.kana_folder.name+'_methods_results_infmad.csv')]))

if __name__=='__main__':main()
