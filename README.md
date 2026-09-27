
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

Semi-Supervised Analyses
Simulation Study

The main semi-supervised simulation figure is reproduced using:

Semi-supervised/simulation/figure1_semi_supervised.R

The supplementary MCAR simulation is reproduced using:

Semi-supervised/simulation/supplementary/tableS3_mcar_or2.R

These scripts use supporting functions from:

Semi-supervised/functions/semi_supervised_functions.R
NHANES Application

The NHANES application is reproduced using:

Semi-supervised/application/nhanes/table6_tableS8_tableS9_nhanes_analysis.R

Supporting functions are stored in:

Semi-supervised/functions/nhanes_functions.R

A single NHANES script is used because the three reported tables are derived from the same underlying application.

The script reproduces:

Main-text Table 6: compact summary of the repeated label-deletion analysis across 1,000 splits.
Supplementary Table S8: detailed results from the selected single label-deletion split.
Supplementary Table S9: coefficient-level repeated-split results, including Monte Carlo variability and analytic standard-error diagnostics.

Thus, Table S8 is the single-split analysis, whereas Tables 6 and S9 summarize the repeated-split analysis.

The NHANES 2017--2018 data are downloaded directly within the R script using the nhanesA package.
