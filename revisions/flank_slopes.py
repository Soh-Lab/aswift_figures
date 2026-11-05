"""Raw-data flank slopes, independent of the selected baseline values."""
from .common import *
from scipy.stats import theilslopes
import argparse


def measure(x,y,lo,hi,fraction=.5):
    x,y=np.asarray(x),np.asarray(y)
    width=(hi-lo)*fraction
    if width<=0: raise ValueError('Invalid window width')
    left=(x>=lo-width)&(x<lo);right=(x>hi)&(x<=hi+width)
    if left.sum()<5 or right.sum()<5: raise ValueError('Fewer than five points on a flank')
    # Robust local straight-line slopes: median of pairwise secant slopes.
    sl=float(theilslopes(y[left],x[left]).slope)
    sr=float(theilslopes(y[right],x[right]).slope)
    return dict(left_slope=sl,right_slope=sr,slope_difference=sr-sl,
                right_ge_left=sr>=sl,n_left=int(left.sum()),n_right=int(right.sum()))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--structured-folder',type=Path,required=True)
    p.add_argument('--windows',type=Path,default=ROOT/'revisions/results/curvature/traces.csv')
    p.add_argument('--out',type=Path,default=ROOT/'revisions/results/flank_slopes')
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    window_path=a.windows
    windows=pd.read_csv(window_path)
    files={'Figure 4':'dox_invitro_shouldering_structured.csv','Figure 5':'kana_invivo_dual_structured.csv','ATP channel 0':'ATP_badNPAu_structured.csv'}
    if 'Dox full' in set(windows.dataset):
        files={'Dox full':files['Figure 4'],'Kana full':files['Figure 5'],'ATP channel 0':files['ATP channel 0']}
    files={name:file for name,file in files.items() if name in set(windows.dataset)}
    rows=[]
    for dataset,file in files.items():
        raw=structured(a.structured_folder/file)
        if dataset in ['Figure 5','Kana full']:
            from peak_extraction.extract_peaks import mask_multichannel_outliers
            for _,g in raw.groupby(['hz','num'],sort=False):
                xx,yy=mask_multichannel_outliers(g.iloc[0].voltage,np.stack(g.current),cutoff=10)
                for idx,y in zip(g.index,yy):raw.at[idx,'voltage']=xx;raw.at[idx,'current']=y
        raw=raw.set_index('source_row')
        for w in windows[windows.dataset==dataset].itertuples():
            r=raw.loc[w.source_row]
            for f in [.25,.5,1.]:
                row=dict(dataset=dataset,source_row=w.source_row,hz=w.hz,channel=w.channel,num=w.num,flank_width_fraction=f)
                try:
                    if w.status!='ok':raise ValueError('Peak window unavailable')
                    row.update(measure(r.voltage,r.current,w.window_min_V,w.window_max_V,f),status='ok',error='')
                except Exception as exc:row.update(status='unclassified',error=str(exc))
                rows.append(row)
    d=pd.DataFrame(rows);d.to_csv(a.out/'traces.csv',index=False)
    groups=list(d.groupby('dataset',sort=False))+[('Figures 4–5 pooled' if 'Figure 4' in set(d.dataset) else 'Dox + kana full datasets',d[d.dataset!='ATP channel 0'])]
    if 'ATP channel 0' in set(d.dataset): groups.append(('All including ATP',d))
    summary=[]
    for name,g in groups:
        for f,h in g.groupby('flank_width_fraction'):
            ok=h[h.status=='ok'];n=int(ok.right_ge_left.sum())
            summary.append(dict(dataset=name,flank_width_fraction=f,n_total=len(h),n_classified=len(ok),n_unclassified=len(h)-len(ok),n_right_ge_left=n,percent_right_ge_left=100*n/len(ok) if len(ok) else np.nan,median_slope_difference=ok.slope_difference.median()))
    s=pd.DataFrame(summary);s.to_csv(a.out/'summary.csv',index=False)
    save_json(a.out/'metadata.json',dict(**provenance([window_path]+[a.structured_folder/f for f in files.values()]),
        protocol='Theil-Sen local slopes of raw currents strictly outside previously selected peak window; Figure 5 retains joint 10-MAD mask. Primary flank span = 0.5 window width on each side, sensitivity 0.25 and 1.0. Minimum five observed points per side; no smoothing or fitted-baseline values. Compare right slope >= left slope with zero tolerance.',
        limitation='Window locations remain ASWIFT-derived. Raw flanks may include peak tails and noise. Endpoint slope ordering is a necessary, not sufficient, condition for a globally convex background and cannot exclude hidden interior concavity. Descriptive repeated-trace percentages, not independent electrode replication.'))
    print(s.to_string(index=False))

if __name__=='__main__':main()
