# Publication inventory

## Figures

Run figure notebooks from the figures directory. The export directory contains the nine panels of Figure 1, Figures 2–5, and Figures S1–S19. PNG exports are tracked in Git. PDF and SVG exports, where generated, are local files. Revision plotting routines for S3, S4, S11 and S16–S19 write into their analysis results directories. When regenerating these figures, copy the numbered exports to the figure export directory before submission.

| Figure | Content | Notebook in figures/ |
|---|---|---|
| 1 | Method illustration, panels a–i | FIG_1.ipynb |
| 2 | Simulation benchmark across baseline families | FIG_2.ipynb |
| 3 | Simulated binding curves | FIG_3.ipynb |
| 4 | In vitro doxorubicin measurements | FIG_4.ipynb |
| 5 | In vivo kanamycin measurements | FIG_5.ipynb |
| S1 | Huber weighting, Pareto curvature and derpsalsa weighting | SI_1.ipynb |
| S2 | Simulation benchmark with Laplacian outliers | SI_2.ipynb |
| S3 | Quadratic baseline tail-height stress test | SI_3.ipynb |
| S4 | Inverted baseline benchmark | SI_4.ipynb |
| S5 | Full 10 Hz doxorubicin response plots | SI_5.ipynb |
| S6–S10 | Representative fits in order: doxorubicin at 10 and 200 Hz, blood, ISF, and ATP on nanoporous gold | SI_6-10.ipynb |
| S11 | Full ATP titration series | SI_11.ipynb |
| S12 | In vivo masking comparison | SI_12.ipynb |
| S13 | Huber cutoff | SI_13.ipynb |
| S14 | Rolling MAD window | SI_14.ipynb |
| S15 | Relative-height threshold | SI_15.ipynb |
| S16–S19 | Derpsalsa asymmetry, scale, numerical settings and smoothing, respectively | SI_16-19.ipynb |

Figure 2 uses the summary rebuilt from the saved raw panel in results/fig2. Figure 3 uses a regenerated 27-amplitude simulation with seed 20261028. The historical raw Figure 3 realization was unavailable. The superseded renders are available in the branch history.

## Main-text tables

| Table | Content | Source |
|---|---|---|
| 1 | Four simulated baseline families, equations and parameter ranges | simulation/peak_simulation.py, panel shape generators |
| 2 | Figure 2 and Figure 3 simulation grids and counts | simulation/simulation_config/simul_config.toml and binding_curve.toml, results/fig2 and results/fig3 |

## SI tables

[tables_S2_S4.tex](tables_S2_S4.tex) contains the two table environments available in this workspace.

| Table | Content | Source under results/ |
|---|---|---|
| S1 | Glossary of symbols in the main text | Manuscript-owned table |
| S2 | Processing time and acquisition intervals | timing_M5/ and timing_R5/ (three runs per configuration), acquisition_intervals.csv |
| S3 | Fixed ASWIFT settings, rationale and sensitivity comparisons | Manuscript-owned table; see S13–S19 and the manuscript checks |
| S4 | In vivo masking, normalization and point-to-point scatter | in_vivo_scatter.csv, figure5_membership.csv, Figure 5 and S12 notebooks |

The S4 rejection row cites Table S3. The noise and rejection audit is available in results/figure_noise_table.csv, with figure membership and input provenance in results/figure4_membership.csv, results/figure5_membership.csv and results/tables_metadata.json. The workspace contains no separate typeset S1 or S3 table. SI_revision_tables.ipynb is a viewer for the audit, timing and scatter summaries.

## Supporting analyses

| Results directory | Role |
|---|---|
| fig2, fig3 | Main simulation figures and noise audit |
| dox, kana, atp_audit | Experimental rejection and noise audits |
| polynomial_stress, concave, derpsalsa, atp | S3, S4, S16–S19 and S11, respectively |
| timing_M5, timing_R5 | Full timing runs used in S2 |
| timing_smoke | Pipeline check only |
| curvature, curvature_full, flank_slopes, flank_slopes_full | Supporting baseline-shape diagnostics. SI_baseline_convexity.ipynb displays these analyses. Retained for the main-text tail-height and peak-window endpoint claims. These plots are diagnostic examples. |
| experimental_checks | All-dataset endpoint and smoothed-tail audit, including ATP channel 0 |

See [manuscript checks](MANUSCRIPT_CHECKS.md) for evidence, scope and wording corrections.

Raw traces, per-fit results, tests, metadata and simulation configurations support reproduction. Earlier drafts and superseded renders are preserved in Git history.

## Regeneration

Install the root requirements and revisions/requirements.txt. Experimental workflows require the structured CSV data and, for the manuscript plots and scatter table, the original method-output CSVs. Set DATA to the structured CSV directory. Some figure notebooks contain local experimental paths that must be set for the machine running them.

```bash
python -m revisions.timing \
  --input "$DATA/dox_invitro_shouldering_structured.csv" \
  --out revisions/results/timing_M5 --hardware-label 'M5 Pro' \
  --workers 1 18 --repeats 3

# On the older computer, with DATA pointing to its local archive:
python -m revisions.timing \
  --input "$DATA/dox_invitro_shouldering_structured.csv" \
  --out revisions/results/timing_R5 --hardware-label 'Ryzen 5 2400G' \
  --workers 1 --repeats 3

# Rebuild the table from an existing run without fitting again:
python -m revisions.timing --out revisions/results/timing_M5 --summarize-only
```

Set DOX_METHODS to the combined doxorubicin methods CSV and KANA_FOLDER to the directory containing the in vivo methods CSVs.

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
  --dox-methods "$DOX_METHODS" \
  --kana-folder "$KANA_FOLDER" \
  --structured-folder "$DATA"
```

```bash
python -m revisions.binding --seed 20261028 --workers 6
```

```bash
python -m revisions.simulations derpsalsa --out revisions/results/derpsalsa --workers 8
python -m revisions.postprocess derpsalsa revisions/results/derpsalsa
```

```bash
python -m revisions.simulations concave --out revisions/results/concave --workers 6
python -m revisions.postprocess concave revisions/results/concave
```

```bash
python -m revisions.atp --input "$DATA/ATP_badNPAu_structured.csv" \
  --concentrations revisions/atp_concentrations.csv --out revisions/results/atp
```

Regenerate the S3 stress test with `python -m revisions.stress --workers 6`. Use `--plot-only` to redraw from the saved results.

Reproduce the all-dataset endpoint-line and measured-tail checks with the manuscript fitter:

```bash
python -m revisions.curvature --structured-folder "$DATA" --scope full --include-atp --out revisions/results/experimental_checks
```

The per-trace file contains the maximum smoothed endpoint and outside-window current, expressed as rise above the baseline beneath the fitted peak, divided by fitted peak height. The endpoint-line metric uses only fitted-baseline values within the selected window.

The smoothed-flank second-difference check uses the saved experimental peak windows and recomputes the full-trace smoother. Results and flank-width sensitivity are in results/smoothed_flank_curvature.

```bash
python -m revisions.smoothed_flank_curvature --structured-folder "$DATA"
```
