###############################################################################
# SUPPLEMENTARY TABLE S7
# GROWTH OF MAXIMUM TREATED-ARM WEIGHTS UNDER OVERLAP STRESS
#
# This script reproduces Supplementary Table S7.
#
# The table reports the mean maximum treated-arm weight and its normalization
# by sqrt(N) under increasingly weak treatment overlap.
#
# Propensity-score scenario:
#   OR2PS1 overlap stress-test
#
# Target treatment rates:
#   0.50, 0.25, 0.10
#
# Sample sizes:
#   N = 500, 2000, 8000
#
# Methods reported:
#   IPW
#   EL
#   CE
#
# Monte Carlo replications:
#   M = 1000
#
# Run this script from the ROOT of the GitHub repository.
###############################################################################


###############################################################################
# 0. CLEAN WORKSPACE
###############################################################################

rm(list = ls())


###############################################################################
# 1. SOURCE STRESS-TEST FUNCTIONS
###############################################################################

source(
  "Causal Inference/functions/stress_test_functions.R"
)


###############################################################################
# 2. SETTINGS
###############################################################################

# Target treatment rates
target_rates <- c(
  0.50,
  0.25,
  0.10
)


# Propensity-score strength for overlap stress-test
ps_strength <- 2


# Monte Carlo replications
M <- 1000


# Cross-fitting folds
K <- 4


# Methods shown in Supplementary Table S7
tableS7_methods <- c(
  "IPW",
  "EL",
  "CE"
)


###############################################################################
# 3. OUTPUT PATHS
###############################################################################

results_root <-
  "Causal Inference/results/overlap_stress_test"


tableS7_dir <-
  file.path(
    results_root,
    "tableS7"
  )


dir.create(
  tableS7_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


###############################################################################
# 4. CALIBRATE PROPENSITY-SCORE INTERCEPTS
###############################################################################

ps_intercepts <-
  sapply(
    target_rates,
    find_ps1_intercept,
    ps_strength = ps_strength
  )


names(ps_intercepts) <-
  as.character(
    target_rates
  )


cat(
  "\nCalibrated propensity-score intercepts:\n"
)


print(
  ps_intercepts
)


###############################################################################
# 5. DEFINE SUPPLEMENTARY TABLE S7 SIMULATION GRID
###############################################################################

tableS7_grid <-
  expand.grid(

    N = c(
      500,
      2000,
      8000
    ),

    target_rate = c(
      0.50,
      0.25,
      0.10
    ),

    stringsAsFactors = FALSE
  )


cat(
  "\nSupplementary Table S7 simulation grid:\n"
)


print(
  tableS7_grid
)


###############################################################################
# 6. RUN FULL OVERLAP STRESS-TEST GRID
###############################################################################

stress_results_list <-
  vector(
    "list",
    nrow(tableS7_grid)
  )


for (j in seq_len(nrow(tableS7_grid))) {

  cat(
    "\n============================================================\n",
    "SUPPLEMENTARY TABLE S7: OVERLAP STRESS-TEST CONDITION\n",
    "N = ",
    tableS7_grid$N[j],
    "\nTarget treatment rate = ",
    tableS7_grid$target_rate[j],
    "\nMonte Carlo replications = ",
    M,
    "\n============================================================\n",
    sep = ""
  )


  stress_results_list[[j]] <-
    run_A2_stress_condition(

      N =
        tableS7_grid$N[j],

      target_rate =
        tableS7_grid$target_rate[j],

      M =
        M,

      ps_strength =
        ps_strength,

      K =
        K,

      numerical_jacobian =
        FALSE
    )
}


###############################################################################
# 7. COMBINE RAW STRESS-TEST RESULTS
###############################################################################

stress_results <-
  dplyr::bind_rows(
    stress_results_list
  )


###############################################################################
# 8. CONVERT SELECTED METHODS TO LONG FORMAT
###############################################################################

tableS7_long <-
  dplyr::bind_rows(

    lapply(

      tableS7_methods,

      function(method_name) {

        make_A2_method_rows(
          stress_results,
          method_name
        )
      }
    )
  )


###############################################################################
# 9. IDENTIFY VALID TREATED-ARM WEIGHT DIAGNOSTICS
#
# Mean MaxW1 is calculated only from replications with finite treated-arm
# effective sample size and finite maximum treated-arm weight.
###############################################################################

tableS7_long <-
  tableS7_long |>

  dplyr::mutate(

    weight1_ok =
      is.finite(
        ESS1
      ) &
      is.finite(
        MAXW1
      )
  )


###############################################################################
# 10. SUMMARIZE MEAN MAXIMUM TREATED-ARM WEIGHT
###############################################################################

tableS7_summary <-
  tableS7_long |>

  dplyr::group_by(
    N,
    target_rate,
    method
  ) |>

  dplyr::summarise(

    Mean_MAXW1 =
      mean(
        MAXW1[
          weight1_ok
        ],
        na.rm = TRUE
      ),

    .groups =
      "drop"
  )


###############################################################################
# 11. COMPUTE N^(1/2)-NORMALIZED MAXIMUM WEIGHT
#
# Supplementary Table S7 reports:
#
#   Mean MaxW1
#
# and
#
#   Mean MaxW1 / sqrt(N)
#
###############################################################################

tableS7_summary <-
  tableS7_summary |>

  dplyr::mutate(

    MaxW1_sqrtN =
      Mean_MAXW1 /
      sqrt(N)
  )


###############################################################################
# 12. ORDER TREATMENT RATES, SAMPLE SIZES, AND METHODS
###############################################################################

tableS7_summary <-
  tableS7_summary |>

  dplyr::mutate(

    target_rate =
      factor(
        target_rate,
        levels = c(
          0.50,
          0.25,
          0.10
        )
      ),

    method =
      factor(
        method,
        levels = tableS7_methods
      )
  ) |>

  dplyr::arrange(
    target_rate,
    N,
    method
  )


###############################################################################
# 13. CREATE MANUSCRIPT-FORMATTED SUPPLEMENTARY TABLE S7
#
# Rounding:
#
#   MaxW1          : 1 decimal
#   MaxW1 / sqrt(N): 2 decimals
###############################################################################

tableS7 <-
  tableS7_summary |>

  dplyr::mutate(

    Mean_MAXW1 =
      round(
        Mean_MAXW1,
        1
      ),

    MaxW1_sqrtN =
      round(
        MaxW1_sqrtN,
        2
      )
  ) |>

  dplyr::select(
    target_rate,
    N,
    method,
    Mean_MAXW1,
    MaxW1_sqrtN
  ) |>

  tidyr::pivot_wider(

    names_from =
      method,

    values_from =
      c(
        Mean_MAXW1,
        MaxW1_sqrtN
      ),

    names_glue =
      "{method}_{.value}"
  ) |>

  dplyr::rename(
    Target_Rate =
      target_rate
  )


###############################################################################
# 14. SAVE SUPPLEMENTARY TABLE S7
###############################################################################

write.csv(

  tableS7,

  file =
    file.path(
      tableS7_dir,
      "tableS7_overlap_weight_growth.csv"
    ),

  row.names = FALSE
)


###############################################################################
# 15. SAVE RAW MONTE CARLO RESULTS
###############################################################################

saveRDS(

  stress_results,

  file =
    file.path(
      tableS7_dir,
      "tableS7_overlap_weight_growth_raw.rds"
    )
)


###############################################################################
# 16. OPTIONAL: SAVE LONG-FORM DIAGNOSTIC RESULTS
###############################################################################

write.csv(

  tableS7_long,

  file =
    file.path(
      tableS7_dir,
      "tableS7_overlap_weight_growth_diagnostics.csv"
    ),

  row.names = FALSE
)


###############################################################################
# 17. PRINT SUPPLEMENTARY TABLE S7
###############################################################################

cat(
  "\n============================================================\n"
)

cat(
  "SUPPLEMENTARY TABLE S7\n"
)

cat(
  "GROWTH OF MAXIMUM TREATED-ARM WEIGHTS\n"
)

cat(
  "============================================================\n\n"
)


print(
  tableS7,
  row.names = FALSE
)


###############################################################################
# 18. PRINT OUTPUT LOCATIONS
###############################################################################

cat(
  "\n============================================================\n"
)

cat(
  "SUPPLEMENTARY TABLE S7 ANALYSIS COMPLETE\n"
)

cat(
  "============================================================\n"
)


cat(
  "\nManuscript table:\n",
  file.path(
    tableS7_dir,
    "tableS7_overlap_weight_growth.csv"
  ),
  "\n",
  sep = ""
)


cat(
  "\nRaw Monte Carlo results:\n",
  file.path(
    tableS7_dir,
    "tableS7_overlap_weight_growth_raw.rds"
  ),
  "\n",
  sep = ""
)


cat(
  "\nDiagnostics:\n",
  file.path(
    tableS7_dir,
    "tableS7_overlap_weight_growth_diagnostics.csv"
  ),
  "\n",
  sep = ""
)


cat(
  "\n============================================================\n"
)
