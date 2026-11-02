"""Show three largest absolute fitted-baseline tail / peak-height ratios."""
from .common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from peak_extraction import extract_peaks as legacy


def main():
    directory=ROOT/'revisions/results/curvature_full'
    selected=pd.read_csv(directory/'tail_ratios_per_trace.csv').nlargest(3,'maximum_flank_absolute_over_peak')
    folder=Path('/Users/Max_1/PycharmProjects/aswift/examples/example_data/structured_csv_example')
    files={'Dox full':'dox_invitro_shouldering_structured.csv','Kana full':'kana_invivo_dual_structured.csv'}
    datasets={}
    for name in selected.dataset.unique():
        d=structured(folder/files[name])
        if name=='Kana full':
            for _,g in d.groupby(['hz','num'],sort=False):
                xx,yy=legacy.mask_multichannel_outliers(g.iloc[0].voltage,np.stack(g.current),cutoff=10)
                for idx,y in zip(g.index,yy):d.at[idx,'voltage']=xx;d.at[idx,'current']=y
        datasets[name]=d.set_index('source_row')
    plt.rcParams.update({'font.family':'sans-serif','font.size':11,'axes.spines.top':True,'axes.spines.right':True})
    fig,axes=plt.subplots(1,3,figsize=(13.5,4.5),sharey=True)
    for rank,(ax,r) in enumerate(zip(axes,selected.itertuples()),1):
        raw=datasets[r.dataset].loc[r.source_row];x,y=raw.voltage,raw.current
        current=legacy.aswift_prefilter_current(y)
        bgmethod,bgeval,pkmethod=legacy.get_aswift_methods()
        params,window=bgmethod(current,x);b=bgeval(x,params)
        peakfit,h,_=pkmethod(current,x,b,window)
        indices=np.flatnonzero(window);ip=indices[np.argmax(peakfit[window])]
        outside=np.flatnonzero(~window);it=outside[np.argmax(b[outside])]
        assert np.isclose(b[it]/h,r.maximum_flank_absolute_over_peak,rtol=1e-8)
        ax.plot(x,y,color='.25',lw=1,label='Measured current')
        ax.plot(x,b,color='#1f77b4',lw=2,label='ASWIFT baseline')
        ax.axvspan(x[window].min(),x[window].max(),color='#1f77b4',alpha=.09,label='Peak window')
        ax.scatter([x[it]],[b[it]],color='#d62728',s=35,zorder=5,label='Maximum baseline tail')
        ax.annotate('',xy=(x[ip],b[ip]+h),xytext=(x[ip],b[ip]),arrowprops={'arrowstyle':'<->','color':'#2ca02c','lw':1.8})
        ax.set_title(f'{rank}. Kanamycin · scan {r.num}\n{r.hz} Hz · channel {r.channel}')
        ax.text(.04,.05,f'Tail: {b[it]:.3f} µA\nPeak height: {h:.3f} µA\nTail / peak: {b[it]/h:.2f}×',transform=ax.transAxes,fontsize=10,bbox={'facecolor':'white','alpha':.9,'edgecolor':'.8'})
        ax.set_xlabel('Potential (V)')
        ax.set_ylim(bottom=0)
    axes[0].set_ylabel('Current (µA)')
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',ncol=4,frameon=False)
    fig.suptitle('Largest absolute baseline-tail / peak-height ratios',fontsize=14)
    fig.subplots_adjust(bottom=.20,top=.77,wspace=.15)
    for ext in ['png','pdf','svg']:fig.savefig(directory/f'largest_tail_ratios.{ext}',dpi=180,bbox_inches='tight',facecolor='white')
    selected.to_csv(directory/'largest_tail_ratios.csv',index=False)
    plt.close(fig)

if __name__=='__main__':main()
