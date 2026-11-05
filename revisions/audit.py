"""Audit every trace, retaining failure reasons and original-row identifiers."""
from .common import *
import argparse
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor


def audit_one(task):
    meta,x,y,saved_peak = task
    from peak_extraction import extract_peaks as legacy
    legacy._y_smooth = None
    try:
        result=legacy.baseline_peak_fit(y,x)
        peak=float(result[0]); state='ok';error=''
    except Exception as exc:
        peak=np.nan;state=status(exc);error=str(exc)
    smooth=legacy._y_smooth
    mean=median=np.nan
    if smooth is not None:
        residual=np.abs(y-smooth)
        mean=float(np.mean(residual));median=float(np.median(residual))
    # For historical simulations, use the saved peak and recompute only its noise reference.
    reference=peak if saved_peak is None else saved_peak
    return dict(**meta, status=state,error=error, refit_peak=peak,peak=reference,
                saved_peak_missing=saved_peak is not None and not np.isfinite(saved_peak),
                mean_absolute_residual=mean,median_absolute_residual=median,
                height_over_mean=reference/mean if mean>0 else np.nan,
                height_over_median=reference/median if median>0 else np.nan)


def summary(df):
    group=[c for c in ['dataset','condition','hz','channel'] if c in df]
    return df.groupby(group,dropna=False).agg(n=('status','size'),
        ok=('status',lambda s:(s=='ok').sum()),
        rejected_noise=('status',lambda s:(s=='rejected_noise').sum()),
        rejected_no_peak=('status',lambda s:(s=='rejected_no_peak').sum()),
        other_errors=('status',lambda s:(s=='fit_error').sum()),
        saved_peak_missing=('saved_peak_missing','sum'),
        minimum_height_over_mean=('height_over_mean','min'),
        minimum_height_over_median=('height_over_median','min')).reset_index()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--dataset',required=True);p.add_argument('--workers',type=int,default=8)
    p.add_argument('--mask-cutoff',type=float);p.add_argument('--channel',type=int)
    p.add_argument('--limit',type=int)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    tasks=[]
    if a.input.suffix=='.feather':
        df=pd.read_feather(a.input)
        for i,r in df.iterrows():
            meta=r.drop(labels=[c for c in df if c=='signal' or c.endswith(('_peak','_bl'))]).to_dict()
            meta.update(dataset=a.dataset,source_row=i,condition='saved_simulation')
            tasks.append((meta,np.linspace(-.4,0,len(r.signal)),np.asarray(r.signal),r.aswift_peak))
    else:
        df=structured(a.input)
        # Mask jointly across all channels before selecting channel, matching manuscript extraction.
        if a.mask_cutoff is not None:
            from peak_extraction.extract_peaks import mask_multichannel_outliers
            for _,g in df.groupby(['hz','num'],sort=False):
                x=g.iloc[0].voltage;ys=np.stack(g.current)
                if not all(np.array_equal(x,v) for v in g.voltage): raise ValueError('Channel grids differ')
                xx,yy=mask_multichannel_outliers(x,ys,cutoff=a.mask_cutoff)
                for (idx,_),y in zip(g.iterrows(),yy):
                    df.at[idx,'voltage']=xx;df.at[idx,'current']=y
        if a.channel is not None: df=df[df.channel==a.channel]
        for _,r in df.iterrows():
            meta={k:r[k] for k in ['source_row','file','hz','num','channel','time','trace_type']}
            meta.update(dataset=a.dataset,condition='raw' if a.mask_cutoff is None else f'{a.mask_cutoff:g}MAD')
            tasks.append((meta,r.voltage,r.current,None))
    if a.limit:tasks=tasks[:a.limit]
    save_json(a.out/'metadata.json',dict(**provenance([a.input]),arguments=vars(a),n=len(tasks),
        noise_definition='Absolute residual to full-trace robust smoother. Both mean and median reported. Historical manuscript fitter; actual rejection uses mean, not median.'))
    rows=[]
    with ProcessPoolExecutor(a.workers,mp_context=mp.get_context('spawn')) as pool:
        for row in pool.map(audit_one,tasks,chunksize=16):
            rows.append(row)
            if len(rows)%500==0: print(f'{a.dataset}: {len(rows)}/{len(tasks)}',flush=True)
    result=pd.DataFrame(rows);result.to_csv(a.out/'traces.csv',index=False)
    if a.input.suffix=='.feather' and a.dataset=='FIG_2' and not a.limit:
        from simulation.statistics import compute_error_stats
        compute_error_stats(df,METHODS).to_feather(a.out/'panel_summary.feather')
    s=summary(result);s.to_csv(a.out/'summary.csv',index=False);print(s.to_string(index=False))
    result.loc[result.groupby(['dataset','condition']).height_over_median.idxmin()].to_csv(a.out/'minimum_traces.csv',index=False)

if __name__=='__main__': main()
