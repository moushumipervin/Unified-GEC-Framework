
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
...

Manuscript Output Map
Output	R script
Figure 1	Semi-supervised/simulation/figure1_semi_supervised.R
Table 3	Causal Inference/simulation/table3_causal.R
Table 4	Missing Covariate/simulation/table4_tableS4_missing_covariate_simulation.R
Table 5	Causal Inference/application/table5_lalonde.R
Table 6	Semi-supervised/application/nhanes/table6_tableS8_tableS9_nhanes_analysis.R
Table S1	Causal Inference/simulation/supplementary/tableS1_causal_n2000.R
Table S2	Causal Inference/simulation/supplementary/tableS2_weight_diagnostics_n1000.R
Table S3	Semi-supervised/simulation/supplementary/tableS3_mcar_or2.R
Table S4	Missing Covariate/simulation/table4_tableS4_missing_covariate_simulation.R
Table S5	Causal Inference/simulation/supplementary/tableS5_tableS6_overlap_stress_test.R
Table S6	Causal Inference/simulation/supplementary/tableS5_tableS6_overlap_stress_test.R
Table S7	Causal Inference/simulation/supplementary/tableS7_overlap_weight_growth.R
Table S8	Semi-supervised/application/nhanes/table6_tableS8_tableS9_nhanes_analysis.R
Table S9	Semi-supervised/application/nhanes/table6_tableS8_tableS9_nhanes_analysis.R
Supplementary overlap figure	Causal Inference/simulation/supplementary/tableS5_tableS6_overlap_stress_test.R
