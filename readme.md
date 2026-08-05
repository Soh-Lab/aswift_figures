# aswift_figures

This repository contains the notebooks, simulation scripts, and figure exports used
to generate ASWIFT manuscript figures. Its purpose is transparency and
reproducibility: each figure notebook is kept with the data-processing and
simulation code used to produce it.

The main ASWIFT software repository is available at
[Soh-Lab/aswift](https://github.com/Soh-Lab/aswift).

## Setup
Set up Python 3.12 and install a C++ compiler.

Install requirements:

```bash
pip install -r requirements.txt
```

## Figure Notebooks

Primary figure notebooks live in `figures/figure1.ipynb` through
`figures/figure5.ipynb`. Supporting information notebooks live in
`figures/SI_1.ipynb` onward.

Rendered outputs are collected in `figures/figure_exports/`.

## Peak Extraction

To extract peaks and generate JSON/CSV outputs, edit
`peak_extraction/config/example_config.toml`, then run:

```bash
python -m peak_extraction.app -c peak_extraction/config/example_config.toml --save
```

To visualize individual fits:

```bash
streamlit run analysis/fit_visualization.py "peak_extraction/config/example_config.toml"
```

Use a different config path in the command when needed.

## Simulations and Statistics

Run one simulation config directly:

```bash
python -m simulation.simul_app -c simulation/simulation_config/simul_config.toml --save
```

Run the manuscript simulation and statistics workflow:

```bash
python -m simulation.workflow -c simulation/workflow_config/paper_simulations.toml
```

To recompute only the statistics from existing simulation outputs:

```bash
python -m simulation.workflow -c simulation/workflow_config/paper_simulations.toml --skip-simulations
```

The binding-curve summary uses `kind = "peak"` in the workflow config; all other
configured summaries use peak-height error statistics.

## Citation

If you use ASWIFT in your research, please cite:

> Yates, M., Ji, J., Yee, S., & Soh, H. T. (2026). Robust Regularization
> Enables Automated, Real-Time Square-Wave Voltammetry Signal Quantification.
> *bioRxiv*. [https://doi.org/10.64898/2026.07.25.740173](https://doi.org/10.64898/2026.07.25.740173)

Machine-readable citation metadata is available in the main ASWIFT repository's
[`CITATION.cff`](https://github.com/Soh-Lab/aswift/blob/main/CITATION.cff).
