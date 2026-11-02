> **Current baseline-shape analysis scope:** doxorubicin and kanamycin only
> (3,000 + 412 = 3,412 traces). ATP is excluded from the curvature,
> endpoint-line, and flank-slope tables and notebook. Earlier ATP-inclusive
> counts below are historical. At 1% tolerance, 100% stay below the endpoint
> line within tolerance; 167/3,412 (4.89%) dip more than 1% below it, and
> 3,245/3,412 (95.11%) remain within ±1%. The separate ATP response figure
> is unaffected. Curvature runs now exclude ATP by default; `--include-atp`
> explicitly restores optional channel-0 inclusion. Flank runs follow their
> input window-table membership.

# Draft revision figure captions

## Derpsalsa sensitivity

Sensitivity of ASWIFT peak-height extraction to derpsalsa parameters. The same
10,800 simulated voltammograms used for the Figure 2 panel are analyzed with
one parameter changed at a time. Upper panels illustrate the peak and four
background families. Lower panels show mean peak-height bias and root mean
square within-family standard deviation, summarized over the panel and as a
function of SNR. Background class, peak width and location vary within each
family. The noise realization is fixed across parameter settings. The default
is p=0.01, k=SD of the smoothed input/10, difference order 2, 50 maximum
iterations, tolerance 0.001, smoothing half-window floor(N/200), and 16 smoothing
passes. Individual pages vary p, k, numerical controls, and smoothing controls.
Failed fits are enumerated in the accompanying table and are not assigned zero
error. These are retrospective sensitivity experiments.

## Inverted backgrounds

ASWIFT and comparator performance on the Figure 2 simulation panel after
negating each true background while holding the positive peak and noise fixed.
The panel spans four background families, 30 background shapes per family,
three peak locations, three widths and SNRs from 15 to 60 dB. Linear backgrounds
remain linear. The inverted polynomial and exponential backgrounds are concave;
the multi-Gaussian family can have mixed curvature. The accompanying concavity
table separates traces whose second derivative is nonpositive, within a stated
numerical tolerance, throughout the true ±3σ peak window. Upper panels show
representative inputs; lower panels show mean bias and RMS within-family SD
across the full inversion panel. Statistics include successful fits; failures
are reported separately rather than omitted from the denominator of failure
rates.

## ATP channel 0

Analysis of ATP binding measurements from one electrode (channel 0) at 150,
200 and 250 Hz. Left: representative raw voltammograms. Middle: normalized peak
height versus ATP concentration, using ASWIFT, polynomial, multi-Gaussian and
extend-wavelet analysis. Right: standard deviation of normalized peak heights
across repeated scans. Each method is normalized by its mean over scans 5–14,
the initial zero-ATP measurements. Buffer-test scans 0–4 are excluded from the
binding plot. Solid curves indicate the initial series; dashed curves indicate
the later 100 µM and 1 mM rechallenge, with circles for both series. Washout
measurements are omitted from the plot. Stages are not pooled across the
washout/rechallenge sequence. The right column reports sample standard deviations
over five scans, except the ten-scan initial zero-ATP reference. These repeated scans are not independent electrodes.

## Regenerated Figure 3

The revised Figure 3 uses the approved 27-amplitude binding grid, three noise
replicates and random seed 20261028. The original raw simulation file was lost;
this is a new realization under the recorded configuration. Use the original
Figure 3 caption's model description, with this provenance note retained in the
reproducibility record. Consult the new audit table before claiming a minimum
peak-height-to-noise ratio for this realization.

## Formatting update for the ATP and parameter figures

ATP now uses the Figure 4 format: a 3×3 layout with a shared four-method bottom
legend, concentration colorbars, normalized response and logarithmic SD plots.
Solid circles/lines show the initial series; dashed lines with the same circles
show rechallenge. Washout crosses are omitted, while all stage results remain
in the accompanying data table. The normalization still uses scans 5–14.

The concave/derpsalsa panels now use Figure 2's original illustrative examples
and panel order. Bias ribbons and error bars show the original plot's scaled
across-family spread (0.35 and 0.8 times that spread, respectively), not confidence
intervals. SD panels show RMS within-family SD. The separate per-trace and
failure tables are unchanged.

## Symmetric quadratic tail-height stress test

Effect of increasingly steep convex quadratic backgrounds on peak-height
extraction. A 1 µA Gaussian peak centered at −0.2 V is added to a quadratic
background with its minimum at the same potential and zero current at its
minimum. Background current at both endpoints (−0.4 and 0 V) is 1, 2, 4, 8 or
16 µA. The Gaussian relative width is 0.3 under the manuscript simulation
convention. Upper panels show one paired Gaussian-noise realization at 30 dB;
current-axis scales vary between examples. Lower panels show mean bias and
sample standard deviation over 100 noise replicates as a function of
log₂(background tail height / peak height). Identical noise realizations are
used across severities and methods. All 2,000 fits returned finite values;
large bias therefore reflects inaccurate recovery rather than fit rejection.

Derpsalsa parameter plots use a dashed blue reference with larger open circles,
drawn above the other curves to expose overlap. The same reference settings are
used in all four parameter figures. The quadratic stress-test SD axis is logarithmic.

## Experimental estimated-baseline convexity table

Trace-weighted classification of fitted ASWIFT baselines in the native selected
peak window. Convexity departure is the maximum vertical excess above the lower
convex envelope, normalized by fitted peak height. Thresholds range from strict
convexity (with numerical allowance) to 2%; the prespecified primary threshold
is 1%. The strict, 0.1%, 0.5%, and 1% results for Figures 4–5 are 61.34%, 96.50%,
99.93%, and 100%, respectively. All fits succeeded. These results describe
estimated baselines, not independent measurements of the background.


## Endpoint-line comparison table

Maximum upward and downward departures of fitted ASWIFT baselines from the
line joining their selected peak-window endpoints, normalized by fitted peak
height. At 1% tolerance, all 1,402 Figure 4–5 baselines stay below the line
within tolerance; 167 (11.91%) dip more than 1% below it and 1,235 (88.09%)
remain within ±1%. ATP is reported separately. Tolerance sensitivity is shown
from 0 to 2%. This describes fitted baseline shapes, not independently measured
backgrounds; remaining below one chord does not establish full convexity.

## Experimental baseline tail-rise table

Per-trace maximum fitted-baseline rise outside the selected peak window,
relative to the baseline beneath the fitted peak, divided by fitted peak
height. All 3,412 dox/kana traces are included. The mean is 0.153 and maximum
0.560; no ratio reaches 8. Endpoint-only and outside-window maxima agree in
this dataset. Absolute ratios are reported separately because they include a
vertical current offset: 206 kanamycin traces have absolute ratios >=8.
