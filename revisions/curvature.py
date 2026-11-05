"""Experimental fitted-baseline convexity via distance to its convex envelope."""
from .common import *
import argparse


def convex_gap(x, y):
    """Maximum vertical distance above the greatest convex minorant."""
    x, y = np.asarray(x), np.asarray(y)
    order = np.argsort(x); x, y = x[order], y[order]
    if len(x)<3 or np.any(np.diff(x)<=0) or not np.isfinite(y).all():
        raise ValueError('Need >=3 finite points at distinct voltages')
    hull=[]
    for i in range(len(x)):
        while len(hull)>1:
            a,b=hull[-2:]
            if (y[b]-y[a])/(x[b]-x[a]) <= (y[i]-y[b])/(x[i]-x[b]): break
            hull.pop()
        hull.append(i)
    return max(0., float(np.max(y-np.interp(x,x[hull],y[hull]))))


def endpoint_gaps(x, y):
    """Maximum excursions above/below the line through window endpoints."""
    x, y = np.asarray(x), np.asarray(y)
    order = np.argsort(x); x, y = x[order], y[order]
    if len(x)<3 or np.any(np.diff(x)<=0) or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Need >=3 finite points at distinct voltages')
    chord = np.interp(x, x[[0,-1]], y[[0,-1]])
    delta = y-chord
    return max(0., float(delta.max())), max(0., float(-delta.min()))


def summarize_endpoints(d):
    rows=[]
    groups=list(d.groupby('dataset',sort=False))
    groups += [('Figures 4–5 pooled' if 'Figure 4' in set(d.dataset) else 'Dox + kana full datasets',d[d.dataset!='ATP channel 0'])]
    if 'ATP channel 0' in set(d.dataset): groups.append(('All including ATP',d))
    for name,g in groups:
        ok=g[g.status=='ok']
        for tol in [0.,.001,.005,.01,.02]:
            above=ok.endpoint_above_over_peak>tol+1e-12
            below=ok.endpoint_below_over_peak>tol+1e-12
            r=dict(dataset=name,tolerance_fraction=tol,n_total=len(g),n_classified=len(ok),n_failed=len(g)-len(ok))
            for label,mask in [('stays_below_with_tolerance',~above),('dips_below',below),
                               ('below_only',below&~above),('approximately_linear',~above&~below),
                               ('above_only',above&~below),('crosses_both',above&below)]:
                r['n_'+label]=int(mask.sum());r['percent_'+label]=100*mask.mean()
            rows.append(r)
    return pd.DataFrame(rows)


def measure(x,y):
    from peak_extraction import extract_peaks as legacy
    current=legacy.aswift_prefilter_current(y)
    bg_method,bg_eval,pk_method=legacy.get_aswift_methods()
    params,window=bg_method(current,x)
    background=bg_eval(x,params)
    fitted_peak,height,_=pk_method(current,x,background,window)
    xx,bb=np.asarray(x)[window],np.asarray(background)[window]
    if not np.isfinite(height) or height<=0: raise ValueError('Invalid peak height')
    from .tail_ratios import metrics
    if pk_method is not legacy.fit_tikhonov_peak:
        raise ValueError('Tail-ratio peak location currently requires manuscript tikhonov peak fitter')
    selected=np.flatnonzero(window)
    peak_index=selected[np.argmax(np.asarray(fitted_peak)[window])]
    tails=metrics(x,background,window,peak_index,height)
    smooth = np.asarray(legacy._y_smooth)
    ends = np.argsort(x)[[0, -1]]
    tails["smoothed_endpoint_ratio"] = (smooth[ends].max() - background[peak_index]) / height
    tails["smoothed_flank_ratio"] = (smooth[~window].max() - background[peak_index]) / height
    above,below=endpoint_gaps(xx,bb)
    return dict(**tails,endpoint_above=above,endpoint_below=below,endpoint_above_over_peak=above/height,endpoint_below_over_peak=below/height,peak_height=height,window_min_V=xx.min(),window_max_V=xx.max(),window_points=len(xx),
        convex_gap=convex_gap(xx,bb),concave_gap=convex_gap(xx,-bb),
        convex_gap_over_peak=convex_gap(xx,bb)/height,concave_gap_over_peak=convex_gap(xx,-bb)/height)


def summarize(d):
    rows=[]
    groups=list(d.groupby('dataset',sort=False))
    groups += [('Figures 4–5 pooled' if 'Figure 4' in set(d.dataset) else 'Dox + kana full datasets',d[d.dataset!='ATP channel 0'])]
    if 'ATP channel 0' in set(d.dataset): groups.append(('All including ATP',d))
    for name,g in groups:
        ok=g[g.status=='ok']
        for tol in [0.,.001,.005,.01,.02]:
            cv=ok.convex_gap_over_peak<=tol+1e-12
            cc=ok.concave_gap_over_peak<=tol+1e-12
            rows.append(dict(dataset=name,tolerance_fraction=tol,n_total=len(g),n_classified=len(ok),
                n_failed=len(g)-len(ok),n_convex_or_linear=int(cv.sum()),percent_convex_or_linear=100*cv.mean(),
                n_linear=int((cv&cc).sum()),n_convex=int((cv&~cc).sum()),
                n_concave=int((~cv&cc).sum()),n_mixed=int((~cv&~cc).sum())))
    return pd.DataFrame(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--structured-folder',type=Path,required=True)
    p.add_argument('--include-atp',action='store_true',help='Optional ATP channel 0; excluded by default')
    p.add_argument('--scope',choices=['figures','full'],default='figures')
    p.add_argument('--out',type=Path,default=ROOT/'revisions/results/curvature')
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    sources=[('Figure 4','dox_invitro_shouldering_structured.csv'),('Figure 5','kana_invivo_dual_structured.csv'),('ATP channel 0','ATP_badNPAu_structured.csv')]
    if not a.include_atp: sources=[(name,file) for name,file in sources if name!='ATP channel 0']
    rows=[]
    save_json(a.out/'metadata.json',dict(**provenance([a.structured_folder/f for _,f in sources]),
        definition='Maximum baseline excess above greatest convex minorant within native ASWIFT peak window, divided by fitted peak height. Primary threshold 0.01. Concavity uses negative baseline.',
        tail_ratio_definition='Maximum baseline rise above baseline at fitted peak maximum, divided by fitted peak height. Endpoint and all-outside-window maxima reported; unshifted ratios also retained. Negative rises clipped to zero.',
        endpoint_definition='Maximum upward and downward departures from the line joining fitted-baseline values at the first and last selected peak-window points, normalized by fitted peak height. Primary tolerance 0.01; approximate linear means both departures <= tolerance.',
        scope=a.scope, include_atp=a.include_atp, membership='Dox: saved Figure 4 membership for figures scope, all rows for full scope. Kana: all traces with joint 10-MAD mask. ATP excluded by default; optional inclusion uses all channel 0 scans.',
        limitation='Estimated ASWIFT baseline, not independently observed true background. Trace-weighted descriptive percentages, not independent electrode observations.'))
    for name,file in sources:
        d=structured(a.structured_folder/file)
        if name=='Figure 4' and a.scope=='figures':
            selected=pd.read_csv(ROOT/'revisions/results/figure4_membership.csv').source_row
            d=d[d.source_row.isin(selected)]
        elif name=='Figure 5':
            from peak_extraction.extract_peaks import mask_multichannel_outliers
            for _,g in d.groupby(['hz','num'],sort=False):
                x=g.iloc[0].voltage
                xx,yy=mask_multichannel_outliers(x,np.stack(g.current),cutoff=10)
                for idx,y in zip(g.index,yy): d.at[idx,'voltage']=xx;d.at[idx,'current']=y
        elif name=='ATP channel 0': d=d[d.channel==0]
        if a.scope=='full': name={'Figure 4':'Dox full','Figure 5':'Kana full'}.get(name,name)
        for i,r in enumerate(d.itertuples(),1):
            meta=dict(dataset=name,source_row=r.source_row,hz=r.hz,channel=r.channel,num=r.num)
            try: result=dict(**measure(r.voltage,r.current),status='ok',error='')
            except Exception as exc: result=dict(status=status(exc),error=str(exc))
            rows.append(dict(**meta,**result))
            if i%100==0: print(f'{name}: {i}/{len(d)}',flush=True)
        pd.DataFrame(rows).to_csv(a.out/'traces.csv',index=False)
    table=summarize(pd.DataFrame(rows));table.to_csv(a.out/'summary.csv',index=False)
    endpoints=summarize_endpoints(pd.DataFrame(rows));endpoints.to_csv(a.out/'endpoint_summary.csv',index=False)
    from .tail_ratios import summarize as summarize_tails
    summarize_tails(pd.DataFrame(rows)).to_csv(a.out/'tail_ratio_summary.csv',index=False)
    print(endpoints[endpoints.tolerance_fraction==.01].to_string(index=False))

if __name__=='__main__':main()
