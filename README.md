# FTB_design
R code supporting the simulation study in the manuscript Flexible Threshold Biomarker Design for Biomarker-Guided Clinical Trials.
# A Flexible Threshold Biomarker (FTB) Trial Design for Treatment Intensification and De-intensification

This repository contains R code for the simulation study of a flexible threshold biomarker-based (FTB) design for jointly evaluating non-inferiority and superiority across biomarker-defined subgroups.

## Overview

The design considers two treatment comparisons within a single protocol:

* **Non-inferiority**: de-intensified treatment (**D**) versus control (**C**) in lower-risk patients
* **Superiority**: intensified treatment (**I**) versus control (**C**) in higher-risk patients

Rather than relying on a single fixed biomarker cutoff, the FTB design uses a two-stage procedure:

1. **Stage 1: Bootstrap-based threshold selection**
   Candidate thresholds are screened using Wald-type statistics, and a bootstrap-based global test is used to control family-wise error across the initial threshold search.

2. **Stage 2: Fixed-order sequential expansion with gatekeeping**
   Starting from the selected threshold, the subgroup is expanded one biomarker level at a time using smoothed odds ratios and one-sided likelihood ratio tests.

The simulation study evaluates type I error and operating characteristics under several manuscript scenarios.

## Repository structure

```text
FTB_design/
├── README.md
├── R/
│   ├── 00_packages.R
│   ├── 10_NI_tests.R
│   ├── 20_Sup_tests.R
│   ├── 30_NI_smoothing.R
│   ├── 40_Sup_smoothing.R
│   ├── 50_NI_bootstrap_stage1_stage2.R
│   ├── 60_Sup_bootstrap_stage1_stage2.R
│   ├── 70_scenario_definitions.R
│   ├── 80_simulation_functions.R
│   ├── 90_scenario_results.R
│   └── 95_figures.R
├── scripts/
│   ├── run_scenarioC_type1.R
│   ├── run_scenarioD_benchmark.R
│   ├── run_scenarioD_targeted_power.R
│   ├── run_scenarioE_targeted_power.R
│   └── run_scenarioF_targeted_power.R
├── results/
└── figures/
```

## Description of files

### Core testing and smoothing functions

* **`10_NI_tests.R`**
  Functions for non-inferiority testing, including Wald-type Z statistics and one-degree-of-freedom likelihood ratio tests.

* **`20_Sup_tests.R`**
  Functions for superiority testing, including Wald-type Z statistics and one-degree-of-freedom likelihood ratio tests.

* **`30_NI_smoothing.R`**
  LOESS smoothing of bin-wise odds ratios for the non-inferiority comparison (**D vs C**).

* **`40_Sup_smoothing.R`**
  LOESS smoothing of bin-wise odds ratios for the superiority comparison (**I vs C**).

### Stage 1 / Stage 2 procedures

* **`50_NI_bootstrap_stage1_stage2.R`**
  Two-stage procedure for non-inferiority threshold selection and expansion.

* **`60_Sup_bootstrap_stage1_stage2.R`**
  Two-stage procedure for superiority threshold selection and expansion.

### Simulation setup

* **`70_scenario_definitions.R`**
  Defines the manuscript scenarios and their underlying bin-specific event probabilities.

* **`80_simulation_functions.R`**
  Contains the simulation functions used to generate trial data and apply the FTB procedure.

### Results and plotting

* **`90_scenario_results.R`**
  Summarizes simulation results and prepares threshold-selection and cumulative-inclusion data.

* **`95_figures.R`**
  Generates the figures for threshold selection frequencies and cumulative inclusion percentages.

## Manuscript scenarios

The simulations consider four manuscript scenarios:

* **Scenario C**
  Global null scenario used to assess **type I error**.
  There is no non-inferiority region and no superiority region.

* **Scenario D**
  Alternative scenario used to assess the flexible design under a setting where a fixed cutoff may perform reasonably well.

* **Scenario E**
  Alternative scenario with a narrower or more localized superiority region.

* **Scenario F**
  Alternative scenario representing a setting where a fixed cutoff is more clearly inadequate.

## Simulation outputs

For each simulated trial, the code records quantities such as:

* whether superiority was declared
* whether non-inferiority was declared
* selected superiority threshold (`c2_hat`)
* selected non-inferiority threshold (`c1_hat`)
* smoothed odds-ratio estimates by biomarker level
* analysis sample sizes for the superiority and non-inferiority comparisons

## Software requirements

The code is written in **R**.

Packages used in the project include:

* `survival`
* `survminer`
* `cmprsk`
* `dplyr`
* `foreach`
* `ggplot2`
* `nleqslv`
* `reshape2`
* `gt`
* `logistf`
* `tidyr`
* `doParallel`
* `doRNG`

You may install missing packages using, for example:

```r
install.packages(c(
  "survival", "survminer", "cmprsk", "dplyr", "foreach",
  "ggplot2", "nleqslv", "reshape2", "gt", "logistf",
  "tidyr", "doParallel", "doRNG"
))
```

## How to run the simulations

Before running any script, make sure all required packages are installed.

### Scenario C: type I error

```r
source("scripts/run_scenarioC_type1.R")
```

### Scenario D: benchmarked at standard fixed-threshold sample size

```r
source("scripts/run_scenarioD_benchmark.R")
```

### Scenario D: targeted-power sample size

```r
source("scripts/run_scenarioD_targeted_power.R")
```

### Scenario E: targeted-power sample size

```r
source("scripts/run_scenarioE_targeted_power.R")
```

### Scenario F: targeted-power sample size

```r
source("scripts/run_scenarioF_targeted_power.R")
```

Simulation results are saved in the `results/` folder.

## How to reproduce summaries and figures

After the simulation outputs have been generated and loaded into the workspace:

```r
source("R/90_scenario_results.R")
source("R/95_figures.R")
```

This creates summary objects and figure objects such as:

* `selection_ha_benchmark`
* `selection_ha_targeted`
* `selection_case1`
* `selection_case2`
* `selection_cum_ha_benchmark`
* `selection_cum_ha_targeted`
* `selection_cum_case1`
* `selection_cum_case2`

## Example operating characteristics

Simple operating characteristics are summarized as:

* **Type I error** under Scenario C:

  * `sum(h0_sim2000$c2_significant == TRUE, na.rm = TRUE)`
  * `sum(h0_sim2000$c1_significant == TRUE, na.rm = TRUE)`

* **Liberal power** under alternative scenarios:

  * `sum(file$c2_significant == TRUE, na.rm = TRUE)`
  * `sum(file$c1_significant == TRUE, na.rm = TRUE)`

Threshold-selection frequencies and cumulative inclusion percentages are prepared in `90_scenario_results.R`.

## Notes on threshold selection outputs

The code uses the following mappings from selected threshold indices to biomarker midpoints:

* **Superiority**

  * `4 -> 0.35`
  * `5 -> 0.45`
  * `6 -> 0.55`
  * `7 -> 0.65`
  * `8 -> 0.75`

* **Non-inferiority**

  * `4 -> 0.35`
  * `5 -> 0.45`
  * `6 -> 0.55`
  * `7 -> 0.65`

Selections coded as `1` for superiority or `0` for non-inferiority represent no selected threshold in the plotting summaries.

## Reproducibility notes

* Parallel simulations are run using `foreach`, `doParallel`, and `doRNG`.
* Random seeds are set explicitly in the scenario runner scripts.
* Candidate thresholds are evaluated on ordinal biomarker bins.
* Smoothing parameters and gatekeeping thresholds follow the values used in the manuscript simulation study.
* Scenario C is used for type I error evaluation, whereas Scenarios D-F are used to evaluate operating characteristics under alternative truths.

## Intended use

This repository is intended to reproduce the simulation results for the FTB design manuscript. The code is organized to separate:

* scenario definitions
* simulation logic
* result summaries
* figure generation

This structure is intended to make the workflow transparent and reproducible.

## Suggested workflow

A typical workflow is:

1. Source one of the scenario runner scripts in `scripts/`
2. Save the simulation result object in `results/`
3. Source `R/90_scenario_results.R`
4. Source `R/95_figures.R`
5. Print or save the figure objects

## Contact / manuscript note

This code repository accompanies the simulation study for the flexible threshold biomarker-based design manuscript. If this repository is shared publicly, this section can be updated to include the manuscript citation, authors, and DOI.







