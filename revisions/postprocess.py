"""Complete-result checks and simulation interpretation tables."""
from .common import *
import argparse


def main():
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['concave','derpsalsa']);p.add_argument('directory',type=Path)
    a=p.parse_args();meta=json.loads((a.directory/'metadata.json').read_text());d=pd.read_csv(a.directory/'traces.csv',low_memory=False)
    expected=meta['n']*(4 if a.mode=='concave' else len(meta['variants']))
    if len(d)!=expected or d.duplicated(['source_row','method']).any():raise ValueError('Incomplete or duplicate simulation results')
    if a.mode=='concave':
        from simulation.peak_simulation import BASELINE_BUILDERS
        x=np.linspace(-.4,0,400);arrays={k:v(x) for k,v in BASELINE_BUILDERS.items()}
        arches={}
        for _,r in d.drop_duplicates('source_row').iterrows():
            cls=int(r['class'])
            bg=-arrays[r.baseline][cls]*r.amp
            sigma=.4*r.width/2/np.sqrt(-2*np.log(.0111));mask=np.abs(x-r.location)<=3*sigma
            xx=x[mask];bb=bg[mask];chord=np.interp(xx,[xx[0],xx[-1]],[bb[0],bb[-1]])
            arches[r.source_row]=max(0.,float(np.max(bb-chord)))
        d['arch_height']=d.source_row.map(arches)
        d.to_csv(a.directory/'traces.csv',index=False)
        d.groupby(['baseline','method','fully_concave']).agg(n=('bias','size'),
            successful=('bias','count'),mean_bias=('bias','mean'),sd_bias=('bias','std'),
            mean_arch_height=('arch_height','mean')).reset_index().to_csv(a.directory/'concavity_table.csv',index=False)
        d.drop_duplicates('source_row').groupby('baseline').fully_concave.agg(['size','mean']).to_csv(a.directory/'concavity_fractions.csv')
    else:
        d['squared_error']=d.bias**2;d['absolute_error']=d.bias.abs()
        d.groupby('method').agg(n=('status','size'),failures=('status',lambda s:(s!='ok').sum()),
            mean_bias=('bias','mean'),mae=('absolute_error','mean'),rmse=('squared_error',lambda s:np.sqrt(s.mean()))).to_csv(a.directory/'selection_table.csv')
    from .plots import plot_panel
    plot_panel(a.directory,a.mode)

if __name__=='__main__':main()
