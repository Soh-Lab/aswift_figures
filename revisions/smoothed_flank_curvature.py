"""Mean second differences of smoothed experimental flanks outside saved peak windows."""
from .common import *
from peak_extraction import extract_peaks as legacy
import argparse


def measure(x, smooth, lo, hi, fraction):
    order = np.argsort(x)
    x, smooth = np.asarray(x)[order], np.asarray(smooth)[order]
    span = (hi - lo) * fraction
    masks = [(x >= lo-span) & (x < lo), (x > hi) & (x <= hi+span)]
    differences, derivatives = [], []
    result = {}
    for side, mask in zip(['left', 'right'], masks):
        xx, yy = x[mask], smooth[mask]
        if len(xx) < 5 or np.any(np.diff(xx) <= 0):
            raise ValueError('Need five distinct voltage points on each flank')
        delta = np.diff(yy, n=2)
        # Divided differences account for unequal spacing after artifact masking.
        derivative = 2 * np.diff(np.diff(yy) / np.diff(xx)) / (xx[2:] - xx[:-2])
        result[side + '_mean_second_difference'] = delta.mean()
        result[side + '_n_points'] = len(xx)
        differences.append(delta)
        derivatives.append(derivative)
    result['mean_second_difference'] = np.concatenate(differences).mean()
    result['mean_second_derivative'] = np.concatenate(derivatives).mean()
    result['both_flanks_nonnegative'] = all(a.mean() >= 0 for a in differences)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--structured-folder', type=Path, required=True)
    parser.add_argument('--windows', type=Path, default=ROOT/'revisions/results/experimental_checks/traces.csv')
    parser.add_argument('--out', type=Path, default=ROOT/'revisions/results/smoothed_flank_curvature')
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    windows = pd.read_csv(args.windows)
    files = {'Dox full': 'dox_invitro_shouldering_structured.csv',
             'Kana full': 'kana_invivo_dual_structured.csv',
             'ATP channel 0': 'ATP_badNPAu_structured.csv'}
    rows = []
    for dataset, filename in files.items():
        raw = structured(args.structured_folder/filename)
        if dataset == 'Kana full':
            for _, group in raw.groupby(['hz', 'num']):
                x, ys = legacy.mask_multichannel_outliers(group.iloc[0].voltage, np.stack(group.current), cutoff=10)
                for idx, y in zip(group.index, ys):
                    raw.at[idx, 'voltage'], raw.at[idx, 'current'] = x, y
        raw = raw.set_index('source_row')
        for w in windows[windows.dataset == dataset].itertuples():
            r = raw.loc[w.source_row]
            error = ''
            try:
                if w.status != 'ok':
                    raise ValueError('Saved window unavailable')
                current = legacy.aswift_prefilter_current(r.current)
                _, smooth, _ = legacy.choose_lambda_lcurve(current)
                if legacy.config.parameters.huber_reweight:
                    _, smooth, _ = legacy.huber_reweighted_lcurve(current, smooth)
            except Exception as exc:
                error = str(exc)
            for fraction in [.25, .5, 1.]:
                row = dict(dataset=dataset, source_row=w.source_row, hz=w.hz, channel=w.channel,
                           num=w.num, flank_width_fraction=fraction)
                try:
                    if error:
                        raise ValueError(error)
                    row.update(measure(r.voltage, smooth, w.window_min_V, w.window_max_V, fraction),
                               status='ok', error='')
                except Exception as exc:
                    row.update(status='unclassified', error=str(exc))
                rows.append(row)
        print(dataset, 'complete', flush=True)
    data = pd.DataFrame(rows)
    data.to_csv(args.out/'traces.csv', index=False)
    groups = list(data.groupby('dataset', sort=False))
    groups += [('Blood', data[(data.dataset == 'Kana full') & (data.channel == 0)]),
               ('ISF', data[(data.dataset == 'Kana full') & (data.channel == 1)]),
               ('All experimental', data)]
    groups += [(f'Dox {hz} Hz', g) for hz, g in data[data.dataset == 'Dox full'].groupby('hz')]
    summaries = []
    for name, group in groups:
        for fraction, g in group.groupby('flank_width_fraction'):
            ok = g[g.status == 'ok']
            summaries.append(dict(dataset=name, flank_width_fraction=fraction, n_total=len(g),
                n_classified=len(ok), n_unclassified=len(g)-len(ok),
                n_mean_nonnegative=int((ok.mean_second_difference >= 0).sum()),
                percent_mean_nonnegative=100*(ok.mean_second_difference >= 0).mean(),
                percent_both_flanks_nonnegative=100*ok.both_flanks_nonnegative.astype(bool).mean(),
                percent_spacing_corrected_nonnegative=100*(ok.mean_second_derivative >= 0).mean()))
    summary = pd.DataFrame(summaries)
    summary.to_csv(args.out/'summary.csv', index=False)
    save_json(args.out/'metadata.json', dict(**provenance([args.windows]+[args.structured_folder/f for f in files.values()]),
        protocol='Recompute manuscript full-trace robust smoother. Saved peak windows. Primary flank width 0.5 window width per side; sensitivity 0.25 and 1.0. Five points minimum per side. Second differences stay entirely within each flank; pool all valid triples equally. Zero sign threshold. Also report separate-flank signs and unequal-voltage-spacing corrected second derivatives.',
        limitation='A nonnegative mean second difference does not imply pointwise convexity. Smoothed total current includes peak tails; windows remain ASWIFT-derived. All ATP channel-0 stages included; in vivo joint 10-MAD mask retained.'))
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()
