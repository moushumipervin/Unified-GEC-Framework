###############################################################################
# MAIN TABLE 7
# ATE estimator performance under increasingly limited overlap
# OR2PS1 stress-test scenario
###############################################################################

rm(list = ls())


###############################################################################
# 1. SOURCE FUNCTIONS
###############################################################################

source(
  "Causal Inference/functions/stress_test_functions.R"
)


###############################################################################
# 2. TABLE 7 SETTINGS
###############################################################################

# Target treatment rates shown in Main Table 7
A2_rates <- c(
  0.50,
  0.25,
  0.10
)

# Strength of propensity-score slopes
# Original PS1 = 1
# Stress-test setting used for Table 7 = 2
A2_ps_strength <- 2


# Number of Monte Carlo replications
M <- 1000


# Sample size used in Main Table 7
N_table7 <- 2000


# Number of cross-fitting folds
K <- 4


# True ATE under OR2
A2_true_ATE <- 10


###############################################################################
# 3. CALIBRATE PROPENSITY-SCORE INTERCEPTS
###############################################################################
#
# For each desired treatment rate, find the intercept producing approximately
# that marginal treatment probability.
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
# 4. TABLE 7 DESIGN GRID
###############################################################################
#
# Main Table 7 uses ONLY:
#
#   N = 500
#
# with target treatment rates:
#
#   0.50
#   0.25
#   0.10
#
###############################################################################

A2_grid <-
  data.frame(
    N = rep(
      N_table7,
      length(A2_rates)
    ),
    target_rate = A2_rates,
    stringsAsFactors = FALSE
  )


print(
  A2_grid
)


###############################################################################
# 5. RUN MONTE CARLO STRESS TEST
###############################################################################

A2_all_list <-
  vector(
    "list",
    nrow(A2_grid)
  )


for (j in seq_len(nrow(A2_grid))) {

  cat(
    "\n============================================================\n",
    "TABLE 7 STRESS-TEST CONDITION\n",
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
# 6. METHODS USED IN MAIN TABLE 7
###############################################################################
#
# Table 7 contains:
#
#   IPW
#   AIPW-GAM
#   ET
#   HD
#   CE
#
###############################################################################

table7_methods <- c(
  "IPW",
  "AIPW_GAM",
  "ET",
  "HD",
  "CE"
)


###############################################################################
# 7. CONVERT MONTE CARLO RESULTS TO LONG FORMAT
###############################################################################

A2_long <-
  dplyr::bind_rows(
    lapply(
      table7_methods,
      function(m) {

        make_A2_method_rows(
          A2_all_results,
          m
        )
      }
    )
  )


###############################################################################
# 8. DEFINE VALID REPLICATIONS
###############################################################################
#
# Keep the same definitions used in the original stress-test implementation.
###############################################################################

A2_long <-
  A2_long |>
  dplyr::mutate(

    # ------------------------------------------------------------
    # Finite point estimate
    # Used for Bias and RMSE
    # ------------------------------------------------------------

    est_ok =
      is.finite(
        Estimate
      ),


    # ------------------------------------------------------------
    # Successful analytic inference
    # Used for Coverage
    # ------------------------------------------------------------

    analytic_ok =
      is.finite(
        Estimate
      ) &
      is.finite(
        SE
      ) &
      Success == 1,


    # ------------------------------------------------------------
    # Valid treated-arm weight diagnostics
    # Used for ESS_1
    # ------------------------------------------------------------

    weight1_ok =
      is.finite(
        ESS1
      ) &
      is.finite(
        MAXW1
      )
  )


###############################################################################
# 9. SUMMARIZE TABLE 7 QUANTITIES
###############################################################################
#
# Bias:
#   mean(estimator - true ATE)
#   using all finite point estimates
#
# RMSE:
#   sqrt(mean((estimator - true ATE)^2))
#   using all finite point estimates
#
# Coverage:
#   empirical analytic 95% CI coverage
#   among successful analytic-inference replications
#
# ESS1:
#   mean treated-arm effective sample size
#   among replications with valid treated-arm weight diagnostics
#
###############################################################################

A2_summary <-
  A2_long |>

  dplyr::group_by(
    N,
    target_rate,
    method
  ) |>

  dplyr::summarise(

    # Number of attempted Monte Carlo replications
    N_total =
      M,


    # Number of finite point estimates
    N_est =
      sum(
        est_ok,
        na.rm = TRUE
      ),


    # Number with successful analytic inference
    N_valid =
      sum(
        analytic_ok,
        na.rm = TRUE
      ),


    # Bias
    Bias =
      mean(
        Estimate[est_ok] -
          A2_true_ATE,
        na.rm = TRUE
      ),


    # Monte Carlo SD
    MC_SD =
      sd(
        Estimate[est_ok],
        na.rm = TRUE
      ),


    # RMSE
    RMSE =
      sqrt(
        mean(
          (
            Estimate[est_ok] -
              A2_true_ATE
          )^2,
          na.rm = TRUE
        )
      ),


    # Mean analytic SE
    Avg_SE =
      mean(
        SE[analytic_ok],
        na.rm = TRUE
      ),


    # Empirical analytic coverage
    Coverage =
      mean(
        Coverage[analytic_ok],
        na.rm = TRUE
      ),


    # Mean treated-arm ESS
    Mean_ESS1 =
      mean(
        ESS1[weight1_ok],
        na.rm = TRUE
      ),


    .groups = "drop"
  )


###############################################################################
# 10. CONSTRUCT MAIN MANUSCRIPT TABLE 7
###############################################################################

table7 <-
  A2_summary |>

  dplyr::filter(
    N == N_table7,
    method %in% table7_methods
  ) |>

  dplyr::mutate(

    # Preserve manuscript ordering of target treatment rates
    target_rate =
      factor(
        target_rate,
        levels = c(
          0.50,
          0.25,
          0.10
        )
      ),


    # Preserve manuscript method ordering
    method =
      factor(
        method,
        levels = table7_methods
      )
  ) |>

  dplyr::arrange(
    target_rate,
    method
  ) |>

  dplyr::mutate(

    # Method names exactly as shown in manuscript
    Method =
      dplyr::recode(
        as.character(method),
        "IPW" = "IPW",
        "AIPW_GAM" = "AIPW-GAM",
        "ET" = "ET",
        "HD" = "HD",
        "CE" = "CE"
      ),


    # ------------------------------------------------------------
    # Manuscript precision
    # ------------------------------------------------------------

    Bias =
      round(
        Bias,
        3
      ),

    RMSE =
      round(
        RMSE,
        3
      ),

    Coverage =
      round(
        Coverage,
        3
      ),

    ESS1 =
      round(
        Mean_ESS1,
        1
      )
  ) |>

  dplyr::select(
    Target_Rate = target_rate,
    Method,
    Bias,
    RMSE,
    Coverage,
    ESS1
  )


###############################################################################
# 11. SAVE RESULTS
###############################################################################

dir.create(
  "Causal Inference/results/table7",
  recursive = TRUE,
  showWarnings = FALSE
)


# Final manuscript table
write.csv(
  table7,
  "Causal Inference/results/table7/table7_overlap_stress.csv",
  row.names = FALSE
)


# Raw Monte Carlo output
saveRDS(
  A2_all_results,
  "Causal Inference/results/table7/table7_overlap_stress_raw.rds"
)


###############################################################################
# 12. PRINT FINAL TABLE 7
###############################################################################

cat(
  "\n============================================================\n"
)

cat(
  "MAIN TABLE 7: OVERLAP STRESS TEST\n"
)

cat(
  "OR2PS1, N = 500, M = 1000\n"
)

cat(
  "============================================================\n\n"
)


print(
  table7,
  row.names = FALSE
)
