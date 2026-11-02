"""Smoothed measured tail minus fitted baseline at peak, over peak height."""
from .common import *
from peak_extraction import extract_peaks as legacy


def main():
    out=ROOT/'revisions/results/curvature_full'
    fitted=pd.read_csv(out/'traces.csv')
    folder=Path('/Users/Max_1/PycharmProjects/aswift/examples/example_data/structured_csv_example')
    files={'Dox full':'dox_invitro_shouldering_structured.csv','Kana full':'kana_invivo_dual_structured.csv'}
    rows=[]
    for name,file in files.items():
        d=structured(folder/file)
        if name=='Kana full':
            for _,g in d.groupby(['hz','num'],sort=False):
                xx,yy=legacy.mask_multichannel_outliers(g.iloc[0].voltage,np.stack(g.current),cutoff=10)
                for idx,y in zip(g.index,yy):d.at[idx,'voltage']=xx;d.at[idx,'current']=y
        ref=fitted[fitted.dataset==name].set_index('source_row')
        for i,r in enumerate(d.itertuples(),1):
            a=ref.loc[r.source_row]
            row=dict(dataset=name,source_row=r.source_row,num=r.num,hz=r.hz,channel=r.channel,
                baseline_at_peak=a.baseline_at_peak,peak_height=a.peak_height)
            try:
                legacy._y_smooth=None
                bg,_,_=legacy.get_aswift_methods()
                _,window=bg(legacy.aswift_prefilter_current(r.current),r.voltage)
                smooth=legacy._y_smooth
                if smooth is None or not np.isfinite(smooth).all():raise ValueError('No finite smoother output')
                ends=np.argsort(r.voltage)[[0,-1]]
                tail=float(max(smooth[ends]));outside=np.flatnonzero(~window)
                imax=outside[np.argmax(smooth[outside])]
                row.update(smoothed_left_tail=smooth[ends[0]],smoothed_right_tail=smooth[ends[1]],
                    maximum_smoothed_flank=smooth[imax],maximum_flank_voltage=r.voltage[imax],
                    endpoint_ratio=(tail-a.baseline_at_peak)/a.peak_height,
                    maximum_flank_ratio=(smooth[imax]-a.baseline_at_peak)/a.peak_height,status='ok',error='')
            except Exception as exc:row.update(status='error',error=str(exc))
            rows.append(row)
            if i%500==0:print(f'{name}: {i}/{len(d)}',flush=True)
    result=pd.DataFrame(rows);result.to_csv(out/'smoothed_tail_ratios_per_trace.csv',index=False)
    summaries=[]
    for name,g in list(result.groupby('dataset',sort=False))+[('Dox + kana pooled',result)]:
        for metric in ['endpoint_ratio','maximum_flank_ratio']:
            v=g.loc[g.status=='ok',metric]
            summaries.append(dict(dataset=name,metric=metric,n_total=len(g),n_valid=len(v),mean=v.mean(),median=v.median(),p95=v.quantile(.95),p99=v.quantile(.99),maximum=v.max(),n_ge_8=int((v>=8).sum()),percent_ge_8=100*(v>=8).mean()))
    table=pd.DataFrame(summaries);table.to_csv(out/'smoothed_tail_ratio_summary.csv',index=False)
    save_json(out/'smoothed_tail_metadata.json',dict(**provenance([out/'traces.csv']+[folder/f for f in files.values()]),definition='(max smoothed measured current at tails - selected baseline at fitted peak) / fitted peak height. Endpoints and maximum outside selected peak window reported separately. Full-trace robust smoother; no subtraction of selected baseline at tails.',limitation='Smoothed tails may include electrode background, peak tails, and edge effects; this is not independently known true background.'))
    print(table.to_string(index=False))

if __name__=='__main__':main()
