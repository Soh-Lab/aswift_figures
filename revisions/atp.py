"""Channel-zero ATP analysis with all manuscript comparators."""
from .common import *
import argparse
from .timing import timed


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--concentrations',type=Path,help='CSV with hz,num,concentration_M; one row per acquisition')
    p.add_argument('--baseline-n',type=int,default=5)
    p.add_argument('--reuse-fits',action='store_true',help='Replot existing traces.csv with a revised concentration mapping')
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    df=structured(a.input);df=df[df.channel==0].sort_values(['hz','num'])
    if a.concentrations:
        mapping=pd.read_csv(a.concentrations)
        df=df.merge(mapping,on=['hz','num'],validate='one_to_one')
        if len(df)!=len(structured(a.input).query('channel==0')):raise ValueError('Incomplete concentration mapping')
    rows=[]
    existing=pd.read_csv(a.out/'traces.csv') if a.reuse_fits else None
    for r in df.itertuples():
        for method in METHODS:
            result=(existing[(existing.source_row==r.source_row)&(existing.method==method)].iloc[0].to_dict() if existing is not None else timed((method,r.source_row,r.voltage,r.current)))
            result.update(hz=r.hz,num=r.num,channel=0,time=r.time)
            if a.concentrations:
                for col in ['concentration_M','stage','normalization_reference','included']:result[col]=getattr(r,col)
            rows.append(result)
    result=pd.DataFrame(rows)
    for _,idx in result.groupby(['hz','method']).groups.items():
        sub=result.loc[idx].sort_values('num')
        baseline=sub[sub.normalization_reference] if a.concentrations else sub.iloc[:a.baseline_n]
        if len(baseline)==0 or not np.isfinite(baseline.peak).all():raise ValueError('Missing or failed baseline fits')
        reference=baseline.peak.mean()
        if reference<=0:raise ValueError('Invalid normalization reference')
        result.loc[idx,'normalized_peak']=result.loc[idx,'peak']/reference
    result.to_csv(a.out/'traces.csv',index=False)
    # Without metadata, blocks are descriptive consecutive-scan groups, not concentrations or biological replicates.
    result['block']=result.groupby(['hz','method']).cumcount()//a.baseline_n
    grouping=['hz','method','stage','concentration_M'] if a.concentrations else ['hz','method','block']
    if a.concentrations:result_for_summary=result[result.included]
    else:result_for_summary=result
    s=result_for_summary.groupby(grouping).agg(n=('normalized_peak','count'),mean=('normalized_peak','mean'),sd=('normalized_peak','std'),num=('num','mean')).reset_index()
    s.to_csv(a.out/'summary.csv',index=False)
    save_json(a.out/'metadata.json',dict(**provenance([a.input]+([a.concentrations] if a.concentrations else [])),
        baseline_n=a.baseline_n,channel=0,concentrations_known=bool(a.concentrations),
        normalization='Mean scans 5–14 (initial zero ATP); washout excluded from normalization' if a.concentrations else f'Mean first {a.baseline_n} scans per frequency and method; provisional reference',
        uncertainty='SD across repeated scans on one electrode; not independent electrode replicates'))
    from .plots import plot_atp
    plot_atp(df,result,s,a.out,bool(a.concentrations))

if __name__=='__main__':main()
