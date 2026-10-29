# ACS revision results

Completed on branch `acs_rev`. Full per-trace outputs are saved locally under
`revisions/results`; the large generated trace files are ignored by Git and can
be regenerated using the commands in README.md. Summary tables, provenance,
notebooks and PNG figure exports are available for review.

## Rejections and noise ratios

No ASWIFT rejection or other fitting error occurred in the current-code audits
of 10,800 Figure 2 panel traces, 9,720 regenerated Figure 3 traces, all 3,000
structured doxorubicin traces, 412 full in vivo traces after joint-channel
10-MAD masking, or 180 ATP channel-0 traces. Available saved ASWIFT peak heights
for the selected Figure 2/4/5 inputs also contain no missing values.

| Figure/data subset | Traces | Minimum height / median absolute residual | Minimum height / mean absolute residual |
|---|---:|---:|---:|
| Figure 2 verified raw panel | 10,800 | 7.351746 | 5.750582 |
| Figure 3 regenerated | 9,720 | 21.149855 | 18.687019 |
| Figure 4 included ASWIFT heights | 990 | 7.910252 | 3.284569 |
| Figure 5 blood | 206 | 6.280463 | 5.491561 |
| Figure 5 ISF | 206 | 9.861385 | 8.474557 |
| ATP channel 0, including buffer-test scans | 180 | 12.761386 | 6.259338 |

The minimum across these figure inputs is **6.280463** using median absolute
residuals. The actual code rejects below twice the **mean** absolute residual,
not twice the median. Both columns use full-trace robust-smoother residuals.
Minimum row identifiers are in `results/figure_noise_table.csv`; full records
and frequency/channel breakdowns are in each audit directory.

These are not an unqualified reconstruction of all historical processing:

- Figure 2's historical summary does not match its checked-in raw panel. The
  new `FIG_2_revision.ipynb` uses a summary rebuilt from that panel's saved
  heights. Noise was recomputed under current smoothing code; current refitted
  heights differ from saved heights by up to 0.067401 µA.
- Figure 3 is a new realization with seed 20261028 and the user-approved
  27-amplitude grid. The lost historical data cannot be audited retrospectively.
- Figure 4 uses its actual first-45-scan, 10/200 Hz membership and available
  ASWIFT columns for channels 1–11. Saved heights are retained; reconstructed
  noise uses current smoothing code. Maximum saved/refit peak difference within
  this subset is 0.007419 µA.
- Figure 5's selected saved peaks are reproduced exactly. The 23 partial
  acquisitions excluded by the figure are separate from fit rejections.

## Derpsalsa sensitivity

All 12 variants were applied to all 10,800 paired inputs (129,600 fits).
The reference settings had zero failures, mean bias −0.003303 µA and RMSE
0.026849 µA. Difference orders 1 and 3 increased RMSE to 0.060757 and
0.037674 µA, respectively. Doubling maximum iterations made negligible
numerical difference; tighter tolerance and changes to smoothing produced
small differences. Changing p to 0.001 produced one non-positive-definite
solver failure (source row 1621), not a noise rejection. Its denominator is
retained in the failure table.

These results support the current defaults as a reasonable reference. They do
not establish a unique optimum: several alternatives differ by only a small
amount in RMSE, and the experiment is a retrospective sensitivity analysis.
The four SI pages and `selection_table.csv` report the results without
suppressing high-variance conditions.

## Inverted backgrounds

All four methods returned finite fits on all 10,800 inverted-background traces.
ASWIFT mean bias was +0.214602 µA for exponential backgrounds and +0.310693 µA
for polynomial backgrounds, against a true 1 µA peak. The linear control had
mean bias +0.008363 µA. Only 79.26% of the inverted multi-Gaussian traces were
concave throughout the defined true peak window; mean ASWIFT bias in that
subset was +0.366926 µA. See `concavity_table.csv` for the mixed-curvature subset.
This experiment demonstrates systematic overestimation even without rejection.

## ATP and additional corrections

The channel-0 analysis includes all 180 traces at 150/200/250 Hz and all four
methods. The concentration figure excludes the first five buffer-test scans
per frequency and normalizes to the ten initial zero-ATP scans (5–14). Washout
and rechallenge remain separate. The CSV exports include stage-specific sample
counts, normalized means and sample SDs. Repeated scans are one-electrode
measurements, not independent electrode replication.

Added `in_vivo_scatter.csv` and `acquisition_intervals.csv` as supporting tables.
The current masking code uses rolling-median residuals, not ASWIFT-smoother
residuals; revise the letter accordingly. Figure 5's example traces now use an
explicit cutoff of 10, matching its fitted-data mask. The actual prominence
default is 0.5, whereas the checklist says 0.3; this discrepancy remains a
manuscript/settings decision and was not silently changed.

## Timing and validation

The package-based timing pipeline passed a 64-fit smoke test covering all four
methods with one and two workers. All 64 succeeded. This run occurred alongside
other computation and is **not publication timing evidence**. Run the documented
full-data commands on the M5 Pro and i5-6500 to substantiate hardware claims.
Neither the 200 ms assertion nor a named-machine comparison has been confirmed.

Seven focused tests passed. The regenerated Figure 2/3 plots, corrected
Figure 5 plot, ATP plot, concave SI and all four derpsalsa SI pages were rendered
and visually inspected. The regenerated Figure 3 SD axis was expanded to show
the full multi-Gaussian variation instead of clipping it at the historical limit.
No expected-potential feature was added to this repository.

## Symmetric quadratic extremes and plotting update

SI 12 is restored from the committed notebook. ATP now uses Figure 4's exact
plotting helpers in a 3×3 layout, with the four-method legend underneath;
washout crosses are omitted. The concave and derpsalsa figures reuse Figure 2's
original layout, example generator, markers, ribbons and legend. Fit results
were not recomputed for these formatting changes.

The added polynomial stress test used a 1 µA positive peak centered at −0.2 V,
a zero-centered quadratic baseline with tail/peak ratios 1, 2, 4, 8 and 16,
30 dB noise and 100 paired replicates. All 2,000 fits succeeded. ASWIFT mean
bias was respectively −0.007713, −0.017409, −0.042969, −0.173969 and −0.477911 µA;
sample SD was 0.013453, 0.013588, 0.012744, 0.017754 and 0.048072 µA. This test
shows deterioration at the largest tail heights even without fit rejection.
The new SI includes all four methods without clipping their errors.
