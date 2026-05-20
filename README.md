# CSB_binary
Reproducible simulation code for Robust Estimation and Inference in Hybrid Controlled Trials for Binary Outcomes. This repository contains the R code used for the primary simulation studies and supplementary sensitivity analyses.

## Overview
The simulation code is organized by scenario. Each folder contains the data-generating code and simulation code for one specific setting.

- `primary/`: primary simulation studies reported in the main text.
- `weaknull_hetero/`: supplementary sensitivity analysis under the weak null scenario with treatment effect heterogeneity.
- `X_dim/`: supplementary sensitivity analysis with higher covariate dimensionality.

## Folder structure

```text
.
├── fun_JL_all.R
├── 3_plot.R
├── primary/
│   ├── 1_simdata_JL.R
│   ├── 2_sim_JL_all.R
│   └── 2_sim_JL_all_mse_gamma.R
├── weaknull_hetero/
│   ├── 1_simdata_JL_weaknull_hetero.R
│   └── 2_sim_JL_all_weaknull_hetero.R
└── X_dim/
    ├── 1_simdata_JL_Xdim.R
    └── 2_sim_JL_all_Xdim.R
```

## File Description

For each simulation scenario:

- Files beginning with `1_` generate the simulated data.
- Files beginning with `2_` run the simulation analysis.
- `fun_JL_all.R` contains shared functions used across all simulation scenarios.
- `3_plot.R` contains plotting code for summarizing and visualizing simulation results.

## How to Run the Code

Before running a simulation scenario, please place the shared files

```text
fun_JL_all.R
3_plot.R
```
in the same folder as the corresponding `1_` and `2_` files.

For example, to run the primary simulation studies, the `primary/` folder should contain:

```text
primary/
├── 1_simdata_JL.R
├── 2_sim_JL_all.R
├── 2_sim_JL_all_mse_gamma.R
├── fun_JL_all.R
└── 3_plot.R
```

Then run the scripts in order:

```r
source("1_simdata_JL.R")
source("2_sim_JL_all.R")
source("3_plot.R")
```
## Notes
- The files `fun_JL_all.R` and `3_plot.R` are intended to be reused across all simulation scenarios. To avoid path issues, we recommend copying these two files into the same folder as the scenario-specific data-generation and simulation scripts before running the code.
- Simulation results and figures can be generated separately within each scenario folder.
- To generate datasets under a given hidden-bias magnitude, please first run the corresponding `1_` data-generation file with `b = 0`, which generates the no-hidden-bias scenario. Then change the hidden-bias parameter to the desired value, for example `b = 6`, and rerun the same data-generation script. This will generate the dataset under hidden bias magnitude 6. Other hidden-bias magnitudes can be generated in the same way by changing the value of `b` and rerunning the script.

