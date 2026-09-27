
# Reproducibility Code

This repository contains the R code used to reproduce the simulation studies and data applications reported in the main manuscript and Supplementary Material.

The repository is organized into three analysis components:

Semi-supervised inference
Causal inference
Missing-covariate analysis

Supporting functions are separated from simulation and application scripts.

Some R scripts reproduce more than one manuscript table or figure. This is intentional. When several manuscript outputs are based on the same Monte Carlo experiment or data analysis, they are generated from a single script rather than rerunning the same analysis separately for each table.

Generated tables and figures are not included in the repository because the final numerical results are reported in the manuscript and Supplementary Material. Running the corresponding scripts reproduces those results.
## Repository Structure

```text
Semi-supervised/
├── functions/
│   ├── nhanes_functions.R
│   └── semi_supervised_functions.R
│
├── simulation/
│   ├── figure1_semi_supervised.R
│   └── supplementary/
│       └── tableS3_mcar_or2.R
│
└── application/
    └── nhanes/
        └── table6_tableS8_tableS9_nhanes_analysis.R


Causal Inference/
├── functions/
│   ├── ate_functions.R
│   ├── stress_test_functions.R
│   └── lalonde_functions.R
│
├── simulation/
│   ├── table3_causal.R
│   └── supplementary/
│       ├── tableS1_causal_n2000.R
│       ├── tableS2_weight_diagnostics_n1000.R
│       ├── tableS5_tableS6_overlap_stress_test.R
│       └── tableS7_overlap_weight_growth.R
│
└── application/
    └── table5_lalonde.R


Missing Covariate/
├── functions/
│   └── missing_covariate_functions.R
│
└── simulation/
    └── table4_tableS4_missing_covariate_simulation.R
```

## Semi-Supervised NHANES Analysis

The NHANES application is reproduced using:

```text
Semi-supervised/application/nhanes/nhanes_semi_supervised_analysis.R
```

The supporting functions are stored in:

```text
Semi-supervised/functions/nhanes_functions.R
```

This analysis produces:

* Main-text Table 6: compact repeated-split NHANES summary.
* Supplementary Table S8: detailed results from the selected single label-deletion split.
* Supplementary Table S9: coefficient-level results across 1,000 label-deletion splits.

Generated results are saved under:

```text
Semi-supervised/results/nhanes/
```

The NHANES data are downloaded directly in the R script using the `nhanesA` package.

## Severe-Overlap Causal Stress Test

The severe-overlap simulation is reproduced using:

```text
Causal Inference/application/overlap_stress_test/severe_overlap_stress_test.R
```

Supporting functions are stored in:

```text
Causal Inference/functions/stress_test_functions.R
```

This analysis produces:

* Supplementary Table S5: true propensity-score diagnostics.
* Supplementary Table S6: complete estimator-performance results.
* Supplementary figure: mean maximum treated-arm weight across sample sizes and target treatment rates.

The additional maximum-weight growth analysis is reproduced using:

```text
Causal Inference/application/overlap_stress_test/tableS7_overlap_weight_growth.R
```

This script produces:

* Supplementary Table S7: growth of the mean maximum treated-arm weight and the normalized quantity MaxW1 / sqrt(N).

Generated results are saved under:

```text
Causal Inference/results/overlap_stress_test/
```

## Software

The analyses were implemented in R.

Important packages include:

```text
dplyr
tidyr
ggplot2
sandwich
MASS
mgcv
CVXR
caret
nhanesA
numDeriv
```

Additional package dependencies may be loaded by individual analysis scripts.

## Running the Code

Clone or download the repository and set the repository root as the working directory.

The scripts use relative paths. Therefore, scripts should be run from the repository root rather than from individual subdirectories.

For example:

```r
source("Semi-supervised/application/nhanes/nhanes_semi_supervised_analysis.R")
```

and

```r
source("Causal Inference/application/overlap_stress_test/severe_overlap_stress_test.R")
```

For Supplementary Table S7:

```r
source("Causal Inference/application/overlap_stress_test/tableS7_overlap_weight_growth.R")
```

The required output directories are created automatically by the scripts.

## Monte Carlo Settings

The severe-overlap simulations use 1,000 Monte Carlo replications.

The repeated-split NHANES evaluation uses 1,000 label-deletion splits.

Random seeds used for reproducibility are specified directly in the corresponding analysis scripts.

## Output Files

The repository separates analysis code from generated results.

Final manuscript tables are stored in table-specific folders, while raw Monte Carlo results and intermediate diagnostic results are stored in `diagnostics` folders.

This structure is intended to make it clear which files directly reproduce manuscript tables and which files are retained for additional verification.

## Reproducibility

For the most reliable reproduction of the results:

1. Use the repository root as the working directory.
2. Install all required R packages.
3. Run the relevant application script.
4. Check the corresponding folder under `results/` for the generated tables and figures.

Package and R version information used for the final analyses is provided in `sessionInfo.txt`.
