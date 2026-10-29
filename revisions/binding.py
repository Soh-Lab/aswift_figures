"""Regenerate lost Figure 3 raw data using the approved amplitude grid."""
from .common import *
from .audit import audit_one
import argparse
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor


def work(task):
    i,row=task
    x=np.linspace(-.4,0,len(row['signal']));y=row['signal']
    from peak_extraction import extract_peaks as legacy
    metrics=audit_one((dict(source_row=i,dataset='FIG_3_regenerated',condition='regenerated'),x,y,None))
    row['aswift_peak']=metrics['peak']
    for method in METHODS[1:]:
        try:row[f'{method}_peak']=fit(method,x,y)
        except Exception:row[f'{method}_peak']=np.nan
    return row,metrics


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,default=ROOT/'revisions/results/fig3')
    p.add_argument('--seed',type=int,default=20261028);p.add_argument('--workers',type=int,default=6)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    from simulation.simulation_config import load_config
    config=ROOT/'simulation/simulation_config/binding_curve.toml';load_config(str(config))
    from simulation.peak_simulation import generate_simulation_data
    from simulation.statistics import compute_peak_stats
    from revisions.audit import summary
    np.random.seed(a.seed);d=generate_simulation_data()
    save_json(a.out/'metadata.json',dict(**provenance([config]),seed=a.seed,n=len(d),
        note='New realization, not recovery of lost original Figure 3 raw data.'))
    rows=[];metrics=[]
    with ProcessPoolExecutor(a.workers,mp_context=mp.get_context('spawn')) as pool:
        for count,(r,m) in enumerate(pool.map(work,enumerate(d.to_dict('records')),chunksize=4),1):
            rows.append(r);metrics.append(m)
            if count%250==0:print(f'binding: {count}/{len(d)}',flush=True)
    result=pd.DataFrame(rows);result.to_feather(a.out/'sim_binding_curve_mg_ewd.feather')
    compute_peak_stats(result,METHODS,bad_threshold=np.inf).to_feather(a.out/'binding_summary.feather')
    metrics=pd.DataFrame(metrics);metrics.to_csv(a.out/'traces.csv',index=False)
    summary(metrics).to_csv(a.out/'summary.csv',index=False)

if __name__=='__main__':main()
