> **Current baseline-shape analysis scope:** doxorubicin and kanamycin only
> (3,000 + 412 = 3,412 traces). ATP is excluded from the curvature,
> endpoint-line, and flank-slope tables and notebook. Earlier ATP-inclusive
> counts below are historical. At 1% tolerance, 100% stay below the endpoint
> line within tolerance; 167/3,412 (4.89%) dip more than 1% below it, and
> 3,245/3,412 (95.11%) remain within ±1%. The separate ATP response figure
> is unaffected. Curvature runs now exclude ATP by default; `--include-atp`
> explicitly restores optional channel-0 inclusion. Flank runs follow their
> input window-table membership.

# ACS revision analyses

Run commands from the repository root with Python 3.12. Install the manuscript
requirements and the published ASWIFT package (activate `venv` first if using it):

```bash
python -m pip install -r revisions/requirements.txt
```

The timing adapter calls `aswift.aswift_fit` from **PyPI ASWIFT 1.0.3**, including
smoothing, parameter selection, baseline fitting, peak fitting, and returned
metrics. The three comparators use the manuscript implementations in
`peak_extraction/extract_peaks.py`. Historical audits and sensitivity experiments
use that manuscript implementation, not a replacement algorithm.

The example notebook is
https://github.com/Soh-Lab/aswift/blob/main/examples/02_csv_workflow.ipynb.
Its data archive is
https://github.com/Soh-Lab/aswift/releases/download/v1.0.3/example_data.zip.
Unlike that notebook's demonstration, timing uses **every row**, without channel
or frequency selection. The local doxorubicin CSV contains 3,000 traces.
Set `DATA` to the extracted `structured_csv_example` directory:

```bash
DATA='/Users/Max_1/Downloads/example_data/structured_csv_example'
```

## Timing protocol and Table S1

Close other computational workloads, use AC power, and record the machine label.
Run the same command on each target machine; do not infer old-desktop timings
from the laptop. Do not run benchmarks alongside the simulations below.

```bash
python -m revisions.timing \
  --input "$DATA/dox_invitro_shouldering_structured.csv" \
  --out revisions/results/timing_M5 --hardware-label 'M5 Pro' \
  --workers 1 18 --repeats 3

# On the older computer, with DATA pointing to its local archive:
python -m revisions.timing \
  --input "$DATA/dox_invitro_shouldering_structured.csv" \
  --out revisions/results/timing_i5 --hardware-label 'i5-6500, 8 GB' \
  --workers 1 --repeats 3

# Rebuild the table from an existing run without fitting again:
python -m revisions.timing --out revisions/results/timing_M5 --summarize-only
```

Outputs:

- `timings.csv`: elapsed full-fit time, trace ID, method, worker count, process
  ID, repeat, peak, success/rejection/error for every attempted fit.
- `batches.csv`: measured batch wall time and effective milliseconds per trace.
- `timing_table.csv`: median, 95th percentile, 99th percentile, maximum worker
  latency, and median effective batch time per trace across repeats.
- `metadata.json`: Python/package versions, input SHA256, settings, hardware/OS
  descriptors and measurement scope. The operator-supplied hardware label is
  not automatically verified against a CPU model or RAM capacity.
- `warmup.csv` / `worker_warmup.csv`: warmup observations. Every process warms
  the current method once in its initializer before measured tasks.

Array parsing and imports are outside individual fit timing. The total batch
clock includes startup, warmup, serialization, collection and shutdown, but not
CSV parsing or writing. Thus throughput is conservative for a persistent worker
pool. Individual latency is never divided by the worker count. All attempted
fits remain in the latency distribution, including failures. Numerical thread
limits are set to one before importing numerical libraries. Worker counts refer
to processes. The serial run is one process with numerical thread limits of one.

`--limit 8 --workers 1 2 --repeats 1` is only a smoke test. The checked-in smoke
output was produced while other analyses were active; **it does not substantiate
any runtime claim**, including the draft's unverified 200 ms bound.

`acquisition_intervals.csv` gives within-channel/frequency timestamp spacing,
not integration duration or the interval between every trace across all
channels. Duplicate/coarsely rounded timestamps are counted separately. Do not
use it to claim hardware scan duration without the acquisition protocol.

## Rejection and noise audits

```bash
python -m revisions.audit \
  --input simulation/simulation_data/paper_data/panel_random_bl.feather \
  --dataset FIG_2 --out revisions/results/fig2
python -m revisions.audit \
  --input "$DATA/dox_invitro_shouldering_structured.csv" \
  --dataset FIG_4 --out revisions/results/dox
python -m revisions.audit \
  --input "$DATA/kana_invivo_dual_structured.csv" --mask-cutoff 10 \
  --dataset FIG_5 --out revisions/results/kana
python -m revisions.audit \
  --input "$DATA/ATP_badNPAu_structured.csv" --channel 0 \
  --dataset ATP --out revisions/results/atp_audit

python -m revisions.tables \
  --dox-methods '/Users/Max_1/Soh Lab/echem/20260213_Dox_BC/combined_ch0.csv' \
  --kana-folder '/Users/Max_1/Soh Lab/echem/20250502_Fig5a_BICEP_csv' \
  --structured-folder "$DATA"
```

The tables command requires the historical methods CSVs referenced by Figures
4/5, including the no-cutoff, 10-MAD and 20-MAD in vivo files. Paths are explicit
arguments so another machine can supply its copies.

The actual rejection test is **height < 2 × mean absolute residual**, not median
absolute residual. Both ratios are exported. Noise means the absolute residual
to the full-trace robust smoother; it is not the residual only inside the peak,
the standard deviation, or the scaled MAD. No-peak rejections, noise rejections,
and other solver errors are counted separately. A missing saved peak is not
silently treated as a confirmed rejection.

For the Figure 2 raw simulation panel, use saved peak heights with recomputed
noise. The historical `panel_summary.feather` does not reproduce from the
checked-in raw panel (maximum ASWIFT summary-entry difference approximately
0.0498 µA). The audit regenerates `results/fig2/panel_summary.feather`;
`FIG_2_revision.ipynb` uses it so the revised figure, audit and SI share inputs.
The reported Figure 2 minimum pertains to this verified raw panel, not an
unavailable realization behind the historical summary. For Figure 4, use saved ASWIFT heights for the first 45 acquisitions at
10/200 Hz and the available ASWIFT channels (1–11; channel 0 has no ASWIFT column
in the historical combined file). This gives 990 fitted heights. Figure 5 uses
206 complete acquisitions per channel, with the same joint-channel 10-MAD mask
as extraction. Its structured file already omits the 23 partial acquisitions;
these exclusions are not ASWIFT rejections. The current refit exactly reproduces
the selected saved Figure 5 peaks. Figure 2/4 do not exactly reproduce saved
heights; the maximum difference is disclosed in the audit table. Their noise
ratios are retrospective estimates under the current smoothing implementation,
not recovered historical noise diagnostics.

Figure 5's example-trace masking now explicitly uses cutoff 10 (it previously
fell back to the helper's cutoff-5 default), matching the saved 10-MAD fitted
results. The figure was regenerated.

The masking implementation currently uses a **rolling-median residual with a
global scaled MAD**, shared across channels. The draft's claim that it uses the
ASWIFT smoother does not describe this code. `in_vivo_scatter.csv` reports mean
fraction of voltage points removed and SD of consecutive normalized differences
(ddof=1)/sqrt(2), using the same full-trace membership and first-five-point
normalization for every method. Missing fits do not get bridged by differences.

## Figure 3 regeneration

The original `sim_binding_curve_mg_ewd.feather` was lost. The two older binding
raw files in `paper_data` have different amplitude grids and must not be used to
validate the final Figure 3. The approved 27-amplitude grid is already in
`simulation/simulation_config/binding_curve.toml`.

```bash
python -m revisions.binding --seed 20261028 --workers 6
```

This generates 9,720 traces (4 families × 30 backgrounds × 27 amplitudes × 3
noise replicates) at 30 dB, using the configuration's constant-noise convention.
It writes a new raw feather file, summary, rejection audit and metadata in
`results/fig3`. Open `figures/FIG_3_revision.ipynb` for the revised figure. This is
a new seeded realization, not recovery of the historical noise. Rebuild the
audit tables after it finishes. Historical `FIG_3.ipynb` remains unchanged.

## Derpsalsa parameter selection SI

```bash
python -m revisions.simulations derpsalsa --out revisions/results/derpsalsa --workers 8
python -m revisions.postprocess derpsalsa revisions/results/derpsalsa
```

All 10,800 saved Figure 2 traces are paired across 12 one-factor-at-a-time
settings; the robust smoother is computed once per trace and shared among
variants. Only the derpsalsa options change. The baseline lambda search remains
50 logarithmic values from 0.1 to 1e7, followed by 10 refinement values (the existing grid helper enforces at least a
three-decade span),
with the existing 0.1 weight-support test on both sides of the peak.

Reference settings are the pybaselines 1.2.1 defaults actually called by the
manuscript: p=0.01, k=SD(smoothed input)/10, difference order=2, maximum
iterations=50, tolerance=0.001, smoothing half-window=floor(N/200), 16 smoothing
passes. Alternatives: p=0.001/0.1; k=0.5×/2× default; order=1/3; iterations=100;
tolerance=0.0001; half-window=0; smoothing passes=8/32. These tests provide
retrospective sensitivity evidence, not proof of how settings were originally
chosen, and not an exhaustive joint optimization.

The other settings remain Huber cutoff 2, local MAD window 21 and peak
prominence **0.5**, matching current code. The checklist's 0.3 is a discrepancy
requiring a manuscript/code decision; this branch does not silently switch it.
The baseline grid above also differs from the smoothing/peak-fitting grids.

`selection_table.csv` reports bias, MAE, RMSE and failure counts. The four
`SI_16.png through SI_19.png` pages follow the five-example/top and four-summary/bottom
layout of SI 10–12. No high-variance outcomes are suppressed. Lower panels use
all classes/locations/widths: mean bias and RMS within-family SD, then the same
metrics versus SNR. These SDs describe the panel, not a confidence interval.

## Inverted-background SI

```bash
python -m revisions.simulations concave --out revisions/results/concave --workers 6
python -m revisions.postprocess concave revisions/results/concave
```

The exact stored signal becomes `signal - 2 * true_baseline`: the peak and noise
are unchanged, and the baseline is negated. Baseline reconstruction is checked
against the stored baseline-at-peak values in tests. All four methods are run
on all 10,800 traces. Linear backgrounds remain linear. Some multi-Gaussian
backgrounds have mixed curvature even before inversion, so the full inversion
panel must not be described as uniformly concave.

`concavity_table.csv` separates traces by whether numerical second derivatives
are ≤1e-8 throughout the true ±3σ Gaussian peak window; the threshold is in the
simulation's current/potential² units. `arch_height` is maximum height above the
chord joining the background at that window's endpoints. `concavity_fractions.csv`
reports how much of each family satisfies the test. Error bars/curves report
successful fits; failure counts remain in `summary.csv`.

## ATP channel 0

```bash
python -m revisions.atp --input "$DATA/ATP_badNPAu_structured.csv" \
  --concentrations revisions/atp_concentrations.csv --out revisions/results/atp
```

The supplied concentration mapping is recorded explicitly in
`atp_concentrations.csv`, at 150, 200 and 250 Hz. Scans 0–4 are buffer tests and
remain in the per-trace fit output/audit, but are excluded from the binding
summary. Scans 5–14 define the initial zero-ATP normalization. Stages 15–39 are
the initial dose series; 40–49 are washout; 50–59 are rechallenge. Repeated
concentrations in different stages are never pooled. Error bars are sample SD
across five scans (ten for the pooled initial zero-ATP reference) on **one electrode**, not biological/electrode replication.
Use `--reuse-fits` to update mapping/plots without fitting again.

The plot contains representative voltammograms, normalized concentration
responses, and repeated-scan SD in the Figure 4/5 colors and three-column layout.
Solid curves show initial dosing and dashed curves show rechallenge, both with
circle markers. Washout measurements remain in the output tables but are omitted
from the figure.

## Verification and boundaries

```bash
python -m unittest revisions.test_revisions -v
```

Tests cover input membership, retained errors, default derpsalsa equivalence,
exact baseline reconstruction/noise preservation, and latency-versus-throughput
bookkeeping. Results carry source hashes and trace IDs. Long simulations write
partial trace checkpoints every 250 inputs; only run postprocessing after the
simulation exits successfully. Postprocessing rejects incomplete results.

New notebooks: `SI_16-19.ipynb`, `SI_4.ipynb`, `SI_10.ipynb`,
`SI_revision_tables.ipynb`, `FIG_2_revision.ipynb`, and `FIG_3_revision.ipynb`. The extra in vivo scatter table supports the checklist;
no expected-potential-range feature is implemented here. Hardware-specific final
timing runs remain for the user on the named machines.

## Consistent manuscript formatting and polynomial extremes

SI 12 was restored byte-for-byte from the current committed version. Revision
plots now load only the definition sections of `FIG_2.ipynb` and `FIG_4.ipynb`
through `manuscript_style.py`. This reuses the original fonts, spines, markers,
colors, axis labels, panel order and boxed bottom legends. The concave and
parameter figures use Figure 2's original example generator (backgrounds are
negated for the concave illustration); numerical summaries remain unchanged.
Figure 2's bias ribbons/error bars are displayed with its original 0.35/0.8
visual scale factors for the across-family spread; they are not confidence
intervals. ATP is 3×3, with Figure 4's discrete concentration colorbars and
logarithmic SD axes. Washout crosses are omitted from the plot, but their data
remain in the exported table. Rechallenge is a separate dashed series with the
same circular markers. Zero-ATP normalization is unchanged.

The symmetric quadratic stress experiment is reproducible with:

```bash
python -m revisions.stress --workers 6
# Redraw without fitting:
python -m revisions.stress --plot-only
```

`b(V) = H*((V+0.2)/0.2)^2` has zero baseline at −0.2 V and endpoint heights H
at −0.4 and 0 V. H is 1, 2, 4, 8 or 16 times the 1 µA peak. The Gaussian center
is −0.2 V, the relative width is 0.3 using the existing simulation convention,
and SNR is 30 dB. The seed is 20261028; each of 100 noise replicates is reused
across the five severities and all four methods. No new solver settings are
introduced. All fits, including errors, are counted in the summary denominator.

`figures/SI_3.ipynb` displays five example traces above two
panels showing mean bias and sample SD versus log₂(tail height / peak height).
The upper panels use separate current scales to show the tails in full. Output
files are in `revisions/results/polynomial_stress/` (CSV, metadata, PNG/PDF/SVG).

## Experimental baseline convexity

```bash
python -m revisions.curvature --structured-folder "$DATA"
```

Outputs: `results/curvature/traces.csv`, `summary.csv`, and `metadata.json`.
The analysis refits the manuscript ASWIFT baseline and peak, using the native
selected peak window. Figure 4 uses its saved 990-trace membership; Figure 5
uses the joint 10-MAD preprocessing (412 traces); ATP uses all 180 channel-0
traces, including buffer and washout scans, and is reported separately.

For each fitted baseline within the window, compute its greatest convex
minorant (lower convex envelope) by a monotone-slope hull. The maximum vertical
gap above this envelope is zero for a convex or linear sequence and positive
for a concave or mixed-curvature sequence. Normalize the gap by fitted peak
height. The primary practical convexity threshold is 1% of peak height, chosen
before inspecting the results; sensitivity tables include 0, 0.1%, 0.5%, 1%,
and 2%. The zero threshold includes only a 1e-12 numerical allowance on the
normalized gap. This is an approximate shape criterion, not a significance test
or a pointwise second-derivative sign test. Applying the same calculation to
the negative baseline distinguishes approximately linear (both gaps small),
convex, concave, and mixed shapes. Failures are retained and explicitly counted;
percentages use successfully classified traces as the denominator.

These are **estimated baselines**, not independently measured backgrounds.
ASWIFT itself selects these estimates, so the result cannot independently
establish the true background's convexity or exclude peak-height bias. Percentages
are trace-weighted; repeated scans are not independent electrodes. Use the
Figures 4–5 pooled row for the main experimental manuscript claim, and report
ATP separately unless the intended claim explicitly includes the supplement.

### Endpoint-line comparison (preferred description of selected shapes)

The same command now also writes `results/curvature/endpoint_summary.csv`.
Join the **fitted baseline values** at the lowest and highest selected-window
voltages by a straight line. Compute the maximum upward and downward departures
of the baseline from that line, divided by fitted peak height. The endpoints
are baseline values, not raw-current values or peak heights.

At each tolerance (0, 0.1%, 0.5%, 1%, 2%), report:
- `stays_below_with_tolerance`: no point exceeds the line by more than tolerance.
- `dips_below`: at least one point lies below the line by more than tolerance.
- Four exclusive categories: `below_only`, `approximately_linear`, `above_only`,
  and `crosses_both`, according to which excursions exceed tolerance.

The first two metrics overlap and must not be added. The four exclusive categories
sum to all successfully classified traces. Endpoints themselves lie on the line;
“above_only” therefore describes excursions, not every point being strictly above.
This endpoint comparison is weaker than full convexity: a baseline can have
local concavity and still stay below the endpoint line. Use it to describe the
shapes selected by ASWIFT, not to infer the unobserved true background. The
numerical allowance is 1e-12 of peak height, as in the envelope analysis.

### Independent raw-flank slope diagnostic

```bash
python -m revisions.flank_slopes --structured-folder "$DATA"
```

Uses saved peak-window bounds from `curvature/traces.csv`, but no selected
baseline values. Fits each raw flank with Theil–Sen (median pairwise slopes),
strictly outside each endpoint, over 0.5 times the window width per side;
checks 0.25 and 1.0 as sensitivity. At least five measured points per flank
are required. Figure 5 retains the joint 10-MAD preprocessing. Counts right
slope >= left slope with no practical tolerance, retaining unclassified cases.
Results and provenance are in `results/flank_slopes/`.

This is a diagnostic of the measured trace, not an isolated background:
residual peak tails can reverse the slope ordering even for a known convex
background (covered by a synthetic regression test). Moreover, ordered endpoint
slopes are necessary but not sufficient for full convexity; they cannot rule
out hidden interior concavity. This analysis therefore cannot close the gap
between selected-baseline geometry and true background shape.

### Full experimental files (rather than figure membership)

The original 1,402-trace denominator was Figure 4 (990 plotted doxorubicin
traces) plus Figure 5 (412 kanamycin traces). The full files contain 3,000
doxorubicin and 412 kanamycin traces, totaling 3,412. ATP contains 540 traces
across three channels, of which channel 0 contributes 180, for a full analysis
of 3,592 traces consistent with the requested ATP channel restriction.

```bash
DATA='/Users/Max_1/PycharmProjects/aswift/examples/example_data/structured_csv_example'
python -m revisions.curvature --structured-folder "$DATA" \
  --scope full --out revisions/results/curvature_full
python -m revisions.flank_slopes --structured-folder "$DATA" \
  --windows revisions/results/curvature_full/traces.csv \
  --out revisions/results/flank_slopes_full
```

These outputs remain separate from the original figure-membership tables.
`curvature` defaults to figure membership; pass `--scope full` explicitly for
all doxorubicin traces. Both scopes use all kanamycin traces with the same
joint 10-MAD mask and all ATP channel-0 scans. The source paths and hashes are
recorded independently for each run.

### Experimental tail-to-peak ratios

The full `revisions.curvature` command also produces `tail_ratio_summary.csv`
and per-trace tail metrics in `traces.csv`, using dox/kana only by default.
At the maximum of the fitted background-subtracted peak, obtain baseline value
B_peak. The stress-comparable endpoint ratio is
`max(0, max(B_left_endpoint, B_right_endpoint) - B_peak) / fitted_peak_height`.
This removes arbitrary vertical current offset and recovers 1, 2, 4, 8, 16
for the symmetric quadratic stress examples. A more conservative metric takes
the maximum over **all points outside the selected peak window**, rather than
only the two full-trace endpoints. Both are reported, along with the unshifted
baseline/peak ratios for transparency (the latter depend on current offset and
are not the stress-test parameter). The peak position uses the maximum of the
manuscript Tikhonov fitted peak within its selected window.

Tables include mean, median, 95th/99th percentiles, maximum, and count/fraction
>=8. These compare the scale of fitted baselines, not true hidden backgrounds.
An experimental sloped or asymmetric baseline with a given ratio is not equivalent
to a symmetric quadratic baseline of that ratio. Ratios alone do not establish
peak-extraction accuracy; repeated traces also weight this descriptive average.

### Corrected smoothed measured-tail ratio

Run `python -m revisions.smoothed_tail_ratios` after the full curvature run.
The requested metric uses **smoothed measured current at the tails**, not the
selected fitted baseline at those tails. Subtract the selected baseline at the
fitted peak and divide by fitted peak height. The earlier baseline-to-baseline
rise metric (maximum 0.560) is a different quantity: the selected baseline can
lie well below the measured tails. The corrected output files start with
`smoothed_tail_` in `results/curvature_full/`. Both exact endpoint maxima and
maxima over all points outside the selected peak window are reported. Retain
this distinction because smoothed tail maxima need not occur at an endpoint.
