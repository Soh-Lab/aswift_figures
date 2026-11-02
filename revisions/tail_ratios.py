"""Tail-to-peak ratios with and without removal of baseline offset."""
from .common import *


def metrics(x, baseline, window, peak_index, height):
    x,b=np.asarray(x),np.asarray(baseline)
    window=np.asarray(window,dtype=bool)
    if height<=0 or not np.isfinite(height) or not np.isfinite(b).all():
        raise ValueError('Invalid baseline or peak height')
    order=np.argsort(x);ends=order[[0,-1]]
    center=b[peak_index];outside=~window
    if not outside.any():raise ValueError('No points outside peak window')
    endpoint_rise=max(0.,float(np.max(b[ends])-center))
    flank_rise=max(0.,float(np.max(b[outside])-center))
    return dict(peak_voltage=x[peak_index],baseline_at_peak=center,
        baseline_left_endpoint=b[ends[0]],baseline_right_endpoint=b[ends[1]],
        maximum_baseline_outside_peak=float(np.max(b[outside])),
        endpoint_rise_over_peak=endpoint_rise/height,
        maximum_flank_rise_over_peak=flank_rise/height,
        endpoint_absolute_over_peak=float(np.max(b[ends]))/height,
        maximum_flank_absolute_over_peak=float(np.max(b[outside]))/height)


def summarize(d):
    rows=[]
    for name,g in list(d.groupby('dataset',sort=False))+[('Dox + kana pooled',d)]:
        ok=g[g.status=='ok']
        for metric in ['endpoint_rise_over_peak','maximum_flank_rise_over_peak','endpoint_absolute_over_peak','maximum_flank_absolute_over_peak']:
            v=ok[metric].dropna()
            rows.append(dict(dataset=name,metric=metric,n_total=len(g),n_valid=len(v),n_failed_or_missing=len(g)-len(v),
                mean=v.mean(),median=v.median(),p95=v.quantile(.95),p99=v.quantile(.99),maximum=v.max(),
                n_ge_8=int((v>=8).sum()),percent_ge_8=100*(v>=8).mean()))
    return pd.DataFrame(rows)
