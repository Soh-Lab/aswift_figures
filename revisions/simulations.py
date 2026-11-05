"""Paired derpsalsa sensitivity and inverted-background Figure 2 experiments."""
from .common import *
import argparse
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
from pybaselines import Baseline

# One factor at a time; k is a multiplier of the documented default std(y)/10.
VARIANTS={
 'default':{}, 'p_0.001':{'p':.001}, 'p_0.1':{'p':.1},
 'k_0.5x':{'k_scale':.5}, 'k_2x':{'k_scale':2.},
 'diff_order_1':{'diff_order':1}, 'diff_order_3':{'diff_order':3},
 'max_iter_100':{'max_iter':100}, 'tol_0.0001':{'tol':.0001},
 'smooth_half_window_0':{'smooth_half_window':0}, 'num_smooths_8':{'num_smooths':8},
 'num_smooths_32':{'num_smooths':32},
}

class ConfiguredBaseline(Baseline):
    options={}
    def derpsalsa(self,data,*args,**kwargs):
        opts=self.options.copy()
        scale=opts.pop('k_scale',None)
        if scale is not None: opts['k']=np.std(data)/10*scale
        return super().derpsalsa(data,*args,**kwargs,**opts)


def work(task):
    mode,idx,row,variants=task
    from peak_extraction import extract_peaks as legacy
    from simulation.peak_simulation import BASELINE_BUILDERS
    x=np.linspace(-.4,0,len(row['signal'])); y=np.asarray(row['signal'])
    meta={k:row[k] for k in ['baseline','class','psnr_db','width','amp','location','replicate']}
    meta['source_row']=idx
    output=[]
    if mode=='concave':
        baseline=BASELINE_BUILDERS[row['baseline']](x)[int(row['class'])]*row['amp']
        # Preserve exact peak and noise realization from the saved Figure 2 panel.
        y=y-2*baseline
        sigma=.4*row['width']/2/np.sqrt(-2*np.log(.0111))
        window=np.abs(x-row['location'])<=3*sigma
        curvature=np.gradient(np.gradient(-baseline,x),x)[window]
        meta.update(concave_fraction=float(np.mean(curvature<=1e-8)),
                    fully_concave=bool(np.all(curvature<=1e-8)),
                    arch_height=float(np.max(-baseline[window])-np.min(-baseline[window])))
        for method in METHODS:
            try:
                peak=legacy.baseline_peak_fit(y,x)[0] if method=='aswift' else fit(method,x,y)
                state='ok' if np.isfinite(peak) else 'fit_error';error=''
            except Exception as exc: peak=np.nan;state=status(exc);error=str(exc)
            output.append(dict(**meta,method=method,peak=peak,error=error,status=state,bias=peak-row['amp']))
    else:
        _,smooth,_=legacy.choose_lambda_lcurve(y)
        _,smooth,_=legacy.huber_reweighted_lcurve(y,smooth)
        legacy._y_smooth=smooth
        lower,upper,_=legacy.get_background_range(smooth,rel_height=legacy.config.parameters.peak_prominence)
        indices=np.zeros_like(y,dtype=bool);indices[lower:upper]=True
        original=legacy.Baseline
        try:
            legacy.Baseline=ConfiguredBaseline
            for name in variants:
                ConfiguredBaseline.options=VARIANTS[name]
                try:
                    lam,bg,_=legacy.choose_lambda_area(smooth,x)
                    _,peak,_=legacy.fit_tikhonov_peak(y,x,bg,indices)
                    state='ok';error=''
                except Exception as exc: peak=lam=np.nan;state=status(exc);error=str(exc)
                output.append(dict(**meta,method=name,peak=peak,baseline_lambda=lam,error=error,status=state,bias=peak-row['amp']))
        finally: legacy.Baseline=original
    return output


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=['derpsalsa','concave']);p.add_argument('--workers',type=int,default=8)
    p.add_argument('--input',type=Path,default=ROOT/'simulation/simulation_data/paper_data/panel_random_bl.feather')
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--variants',nargs='+',choices=list(VARIANTS),default=list(VARIANTS))
    p.add_argument('--snr',type=float,nargs='+');p.add_argument('--limit',type=int)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    df=pd.read_feather(a.input)
    if a.snr:df=df[df.psnr_db.isin(a.snr)]
    if a.limit:df=df.iloc[:a.limit]
    save_json(a.out/'metadata.json',dict(**provenance([a.input]),arguments=vars(a),
        variants={k:VARIANTS[k] for k in a.variants},n=len(df),
        protocol='Saved Fig 2 panel, identical noise and peak across variants. Inverted baselines are not assumed concave: per-trace curvature is recorded.'))
    tasks=[(a.mode,i,r.to_dict(),a.variants) for i,r in df.iterrows()]
    rows=[]
    with ProcessPoolExecutor(a.workers,mp_context=mp.get_context('spawn')) as pool:
        for count,result in enumerate(pool.map(work,tasks,chunksize=4),1):
            rows.extend(result)
            if count%250==0:
                print(f'{a.mode}: {count}/{len(tasks)}',flush=True)
                pd.DataFrame(rows).to_csv(a.out/'traces.csv',index=False)
    result=pd.DataFrame(rows);result.to_csv(a.out/'traces.csv',index=False)
    summary=result.groupby(['method','baseline','psnr_db']).agg(n=('status','size'),
        failures=('status',lambda s:(s!='ok').sum()),mean_bias=('bias','mean'),sd_bias=('bias','std')).reset_index()
    summary.to_csv(a.out/'summary.csv',index=False)

if __name__=='__main__':main()
