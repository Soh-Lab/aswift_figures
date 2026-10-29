"""Convex quadratic tail-height stress test with paired Gaussian noise."""
from .common import *
import argparse
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor

TAIL_RATIOS = np.array([1., 2., 4., 8., 16.])


def make_signals(replicates=100, seed=20261028, snr_db=30., width=.3):
    x = np.linspace(-.4, 0., 400)
    sigma = .4 * width / 2 / np.sqrt(-2 * np.log(.0111))
    peak = np.exp(-.5 * ((x + .2) / sigma)**2)
    backgrounds = TAIL_RATIOS[:, None] * ((x[None, :] + .2) / .2)**2
    noise = np.random.default_rng(seed).normal(0., 10**(-snr_db/20), (replicates, len(x)))
    signals = backgrounds[:, None, :] + peak[None, None, :] + noise[None, :, :]
    return x, peak, backgrounds, signals


def work(task):
    ratio, replicate, x, y = task
    from peak_extraction import extract_peaks as legacy
    rows = []
    for method in METHODS:
        try:
            height = legacy.baseline_peak_fit(y, x)[0] if method == 'aswift' else fit(method, x, y)
            state, error = ('ok', '') if np.isfinite(height) else ('fit_error', 'Nonfinite peak')
        except Exception as exc:
            height, state, error = np.nan, status(exc), str(exc)
        rows.append(dict(tail_ratio=ratio, log2_tail_ratio=np.log2(ratio), replicate=replicate,
                         method=method, peak=height, bias=height-1., status=state, error=error))
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, default=ROOT/'revisions/results/polynomial_stress')
    p.add_argument('--replicates', type=int, default=100)
    p.add_argument('--seed', type=int, default=20261028)
    p.add_argument('--snr-db', type=float, default=30.)
    p.add_argument('--workers', type=int, default=6)
    p.add_argument('--plot-only', action='store_true')
    a = p.parse_args()
    if a.replicates < 2 or a.workers < 1: p.error('Need at least two replicates and one worker')
    a.out.mkdir(parents=True, exist_ok=True)
    if not a.plot_only:
        x, peak, bg, signals = make_signals(a.replicates, a.seed, a.snr_db)
        save_json(a.out/'metadata.json', dict(**provenance(), seed=a.seed, snr_db=a.snr_db,
            replicates=a.replicates, tail_ratios=TAIL_RATIOS.tolist(), peak_height=1.,
            peak_center_V=-.2, relative_peak_width=.3, x=[-.4, 0., 400],
            noise_pairing='Identical noise vector for each replicate across all five tail ratios',
            baseline='tail_ratio * ((voltage + 0.2) / 0.2)**2'))
        tasks = [(ratio, rep, x, signals[i, rep]) for i, ratio in enumerate(TAIL_RATIOS) for rep in range(a.replicates)]
        rows = []
        with ProcessPoolExecutor(a.workers, mp_context=mp.get_context('spawn')) as pool:
            for i, result in enumerate(pool.map(work, tasks, chunksize=4), 1):
                rows.extend(result)
                if i % 100 == 0: print(f'Stress test: {i}/{len(tasks)} traces', flush=True)
        result = pd.DataFrame(rows)
        result.to_csv(a.out/'traces.csv', index=False)
        result.groupby(['method', 'tail_ratio', 'log2_tail_ratio']).agg(
            n=('status', 'size'), failures=('status', lambda s: (s!='ok').sum()),
            mean_bias=('bias', 'mean'), sd=('peak', 'std')).reset_index().to_csv(a.out/'summary.csv', index=False)
    from .plots import plot_stress
    plot_stress(a.out)

if __name__ == '__main__': main()
