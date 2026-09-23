###############################################################################
# MAIN TABLE 8
# Growth of mean maximum treated-arm weight across sample sizes
# OR2PS1 overlap stress-test
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
A2_rates <- c(
  0.50,
  0.25,
  0.10
)

# Stress-test propensity-score strength
A2_ps_strength <- 2

# Monte Carlo replications
M <- 1000

# Cross-fitting folds
K <- 4


###############################################################################
# 3. CALIBRATE PROPENSITY-SCORE INTERCEPTS
###############################################################################

A2_intercepts <-
  sapply(
    A2_rates,
    find_ps1_intercept,
    ps_strength = A2_ps_strength
  )

names(A2_intercepts) <-
  as.character(
    A2_rates
  )


cat(
  "\nCalibrated propensity-score intercepts:\n"
)

print(
  A2_intercepts
)


###############################################################################
# 4. FULL TABLE 8 GRID
###############################################################################
#
# Table 8 requires:
#
#   target rate = 0.50, 0.25, 0.10
#   N           = 500, 2000, 8000
#
###############################################################################

A2_grid <-
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


print(
  A2_grid
)


###############################################################################
# 5. RUN FULL STRESS-TEST GRID
###############################################################################

A2_all_list <-
  vector(
    "list",
    nrow(A2_grid)
  )


for (j in seq_len(nrow(A2_grid))) {

  cat(
    "\n============================================================\n",
    "TABLE 8 STRESS-TEST CONDITION\n",
    "N = ",
    A2_grid$N[j],
    "\nTarget treatment rate = ",
    A2_grid$target_rate[j],
    "\nMonte Carlo replications = ",
    M,
    "\n============================================================\n",
    sep = ""
  )


  A2_all_list[[j]] <-
    run_A2_stress_condition(

      N =
        A2_grid$N[j],

      target_rate =
        A2_grid$target_rate[j],

      M =
        M,

      ps_strength =
        A2_ps_strength,

      K =
        K,

      numerical_jacobian =
        FALSE
    )
}


A2_all_results <-
  dplyr::bind_rows(
    A2_all_list
  )


###############################################################################
# 6. METHODS USED IN MAIN TABLE 8
###############################################################################

table8_methods <- c(
  "IPW",
  "EL",
  "CE"
)


###############################################################################
# 7. CONVERT RESULTS TO LONG FORMAT
###############################################################################

A2_long_table8 <-
  dplyr::bind_rows(
    lapply(
      table8_methods,
      function(m) {

        make_A2_method_rows(
          A2_all_results,
          m
        )
      }
    )
  )


###############################################################################
# 8. IDENTIFY VALID TREATED-ARM WEIGHT DIAGNOSTICS
###############################################################################
#
# Same logic as the original stress-test implementation:
# Mean MaxW1 should use replications with finite treated-arm weight diagnostics.
###############################################################################

A2_long_table8 <-
  A2_long_table8 |>
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
# 9. COMPUTE MEAN MAXIMUM TREATED-ARM WEIGHT
###############################################################################

A2_table8_summary <-
  A2_long_table8 |>

  dplyr::group_by(
    N,
    target_rate,
    method
  ) |>

  dplyr::summarise(

    Mean_MAXW1 =
      mean(
        MAXW1[weight1_ok],
        na.rm = TRUE
      ),

    .groups = "drop"
  )


###############################################################################
# 10. COMPUTE NORMALIZED MAXIMUM WEIGHT
###############################################################################
#
# Table 8 reports:
#
#   Mean MaxW1
#
# and
#
#   Mean MaxW1 / sqrt(N)
#
###############################################################################

A2_table8_summary <-
  A2_table8_summary |>
  dplyr::mutate(

    MaxW1_sqrtN =
      Mean_MAXW1 /
      sqrt(N)
  )


###############################################################################
# 11. ORDER ROWS AND METHODS
###############################################################################

A2_table8_summary <-
  A2_table8_summary |>
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
        levels = c(
          "IPW",
          "EL",
          "CE"
        )
      )
  ) |>

  dplyr::arrange(
    target_rate,
    N,
    method
  )


###############################################################################
# 12. CREATE WIDE MANUSCRIPT TABLE
###############################################################################

###############################################################################
# 12. CREATE WIDE MANUSCRIPT TABLE
###############################################################################

table8 <-
  A2_table8_summary |>

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
    names_from = method,
    values_from = c(
      Mean_MAXW1,
      MaxW1_sqrtN
    ),
    names_glue = "{method}_{.value}"
  ) |>

  dplyr::rename(
    Target_Rate = target_rate
  )


###############################################################################
# 13. SAVE
###############################################################################

dir.create(
  "Causal Inference/results/table8",
  recursive = TRUE,
  showWarnings = FALSE
)


write.csv(
  table8,
  "Causal Inference/results/table8/table8_overlap_damping.csv",
  row.names = FALSE
)


saveRDS(
  A2_all_results,
  "Causal Inference/results/table8/table8_overlap_damping_raw.rds"
)


###############################################################################
# 14. PRINT
###############################################################################

print(
  table8,
  row.names = FALSE
)
