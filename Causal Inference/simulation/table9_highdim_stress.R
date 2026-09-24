###############################################################################
# MAIN TABLE 9
# High-dimensional calibration stress test
###############################################################################

rm(list = ls())


###############################################################################
# 1. LOAD EXISTING FUNCTIONS
###############################################################################

source(
  "Causal Inference/functions/highdim_stress_functions.R"
)


###############################################################################
# 3. A3 DESIGN
#
# Fixed sample size
# Increasing number of balancing functions q
#
# Primary scenario:
#   OR2PS1
###############################################################################

A3_N <- 2000

A3_M <- 1000

A3_K <- 4


A3_q_grid <- c(
  4,
  10,
  20,
  30,
  40,
  50
)


A3_entropies <- c(
  "SL",
  "EL",
  "ET",
  "HD",
  "CE"
)


###############################################################################
# 4. RUN A3 HIGH-DIMENSIONAL STRESS TEST
#
# outcome_model = 2 : OR2
# ps_model      = 1 : PS1
#
# Keep K = 4.
# K is the number of cross-fitting folds and is unrelated to q.
###############################################################################

A3_results <-
  run_A3_stress_test(

    n =
      A3_N,

    p =
      4,

    M =
      A3_M,

    K =
      A3_K,

    q_grid =
      A3_q_grid,

    outcome_model =
      2,

    ps_model =
      1,

    gec_from =
      "gam",

    entropies =
      A3_entropies,

    maxit =
      1000,

    progress =
      TRUE
  )


###############################################################################
# 5. SAVE RAW MONTE CARLO RESULTS
###############################################################################

dir.create(
  "Causal Inference/results/table9",
  recursive = TRUE,
  showWarnings = FALSE
)


saveRDS(
  A3_results,
  file =
    "Causal Inference/results/table9/table9_highdim_stress_raw.rds"
)


###############################################################################
# 6. QUICK CHECK
###############################################################################

cat(
  "\nDimensions of raw A3 results:\n"
)

print(
  dim(
    A3_results
  )
)


cat(
  "\nNumber of replications by q and method:\n"
)

print(
  table(
    A3_results$q,
    A3_results$method
  )
)


###############################################################################
# 7. A3 SUMMARY
###############################################################################
#
# IMPORTANT:
#
# These calculations are kept from the original implementation.
#
# RMSE is calculated among successful replications.
#
# Median log10 condition number uses all finite recorded
# Max_Log10_Condition values.
###############################################################################

A3_summary <-
  A3_results %>%

  dplyr::group_by(
    q,
    dual_dimension,
    method
  ) %>%

  dplyr::summarise(

    ###########################################################################
    # Number attempted
    ###########################################################################

    N_total =
      dplyr::n(),


    ###########################################################################
    # Successful exact-calibration solves
    ###########################################################################

    N_success =
      sum(
        Success == 1,
        na.rm = TRUE
      ),


    ###########################################################################
    # Solver failure rate
    ###########################################################################

    Failure_Rate =
      mean(
        Failure == 1,
        na.rm = TRUE
      ),


    ###########################################################################
    # Bias among successfully solved replications
    ###########################################################################

    Bias =
      if (
        sum(
          Success == 1 &
            is.finite(
              Estimate
            )
        ) >
          0
      ) {

        mean(
          Error[
            Success == 1 &
              is.finite(
                Estimate
              )
          ],
          na.rm = TRUE
        )

      } else {

        NA_real_
      },


    ###########################################################################
    # Monte Carlo SD among successful solves
    ###########################################################################

    MC_SD =
      if (
        sum(
          Success == 1 &
            is.finite(
              Estimate
            )
        ) >
          1
      ) {

        sd(
          Estimate[
            Success == 1 &
              is.finite(
                Estimate
              )
          ],
          na.rm = TRUE
        )

      } else {

        NA_real_
      },


    ###########################################################################
    # RMSE among successfully solved replications
    ###########################################################################

    RMSE =
      if (
        sum(
          Success == 1 &
            is.finite(
              Squared_Error
            )
        ) >
          0
      ) {

        sqrt(
          mean(
            Squared_Error[
              Success == 1 &
                is.finite(
                  Squared_Error
                )
            ],
            na.rm = TRUE
          )
        )

      } else {

        NA_real_
      },


    ###########################################################################
    # Median finite condition number
    ###########################################################################

    Median_Condition =
      if (
        any(
          is.finite(
            Max_Condition
          )
        )
      ) {

        median(
          Max_Condition[
            is.finite(
              Max_Condition
            )
          ],
          na.rm = TRUE
        )

      } else {

        NA_real_
      },


    ###########################################################################
    # 95th percentile finite condition number
    ###########################################################################

    P95_Condition =
      if (
        any(
          is.finite(
            Max_Condition
          )
        )
      ) {

        as.numeric(
          quantile(
            Max_Condition[
              is.finite(
                Max_Condition
              )
            ],
            probs = 0.95,
            na.rm = TRUE
          )
        )

      } else {

        NA_real_
      },


    ###########################################################################
    # Median log10 condition number
    ###########################################################################

    Median_Log10_Condition =
      if (
        any(
          is.finite(
            Max_Log10_Condition
          )
        )
      ) {

        median(
          Max_Log10_Condition[
            is.finite(
              Max_Log10_Condition
            )
          ],
          na.rm = TRUE
        )

      } else {

        NA_real_
      },


    ###########################################################################
    # 95th percentile log10 condition number
    ###########################################################################

    P95_Log10_Condition =
      if (
        any(
          is.finite(
            Max_Log10_Condition
          )
        )
      ) {

        as.numeric(
          quantile(
            Max_Log10_Condition[
              is.finite(
                Max_Log10_Condition
              )
            ],
            probs = 0.95,
            na.rm = TRUE
          )
        )

      } else {

        NA_real_
      },


    ###########################################################################
    # Frequency of effectively infinite/singular Hessian
    ###########################################################################

    Infinite_Condition_Rate =
      mean(
        is.infinite(
          Max_Condition
        ),
        na.rm = TRUE
      ),


    ###########################################################################
    # Exact-calibration residual
    ###########################################################################

    Median_Balance_Residual =
      if (
        any(
          is.finite(
            Max_Balance_Residual
          )
        )
      ) {

        median(
          Max_Balance_Residual[
            is.finite(
              Max_Balance_Residual
            )
          ],
          na.rm = TRUE
        )

      } else {

        NA_real_
      },


    ###########################################################################
    # 95th percentile calibration residual
    ###########################################################################

    P95_Balance_Residual =
      if (
        any(
          is.finite(
            Max_Balance_Residual
          )
        )
      ) {

        as.numeric(
          quantile(
            Max_Balance_Residual[
              is.finite(
                Max_Balance_Residual
              )
            ],
            probs = 0.95,
            na.rm = TRUE
          )
        )

      } else {

        NA_real_
      },


    ###########################################################################
    # Rank-deficiency rate
    ###########################################################################

    Rank_Deficiency_Rate =
      mean(
        (
          Rank_Deficient1 == 1 |
            Rank_Deficient0 == 1
        ),
        na.rm = TRUE
      ),


    ###########################################################################
    # Newton iterations
    ###########################################################################

    Mean_Iterations =
      mean(
        pmax(
          Iterations1,
          Iterations0,
          na.rm = TRUE
        ),
        na.rm = TRUE
      ),


    .groups =
      "drop"
  )


###############################################################################
# 8. FULL A3 DIAGNOSTIC TABLE
#
# This preserves the larger table produced by the original implementation.
###############################################################################

A3_table <-
  A3_summary %>%

  dplyr::select(

    q,

    dual_dimension,

    method,

    N_total,

    N_success,

    Failure_Rate,

    Bias,

    MC_SD,

    RMSE,

    Median_Log10_Condition,

    P95_Log10_Condition,

    Infinite_Condition_Rate,

    Median_Balance_Residual,

    P95_Balance_Residual,

    Rank_Deficiency_Rate,

    Mean_Iterations
  ) %>%

  dplyr::mutate(

    Failure_Rate =
      round(
        Failure_Rate,
        3
      ),

    Bias =
      round(
        Bias,
        4
      ),

    MC_SD =
      round(
        MC_SD,
        4
      ),

    RMSE =
      round(
        RMSE,
        4
      ),

    Median_Log10_Condition =
      round(
        Median_Log10_Condition,
        3
      ),

    P95_Log10_Condition =
      round(
        P95_Log10_Condition,
        3
      ),

    Infinite_Condition_Rate =
      round(
        Infinite_Condition_Rate,
        3
      ),

    Median_Balance_Residual =
      signif(
        Median_Balance_Residual,
        3
      ),

    P95_Balance_Residual =
      signif(
        P95_Balance_Residual,
        3
      ),

    Rank_Deficiency_Rate =
      round(
        Rank_Deficiency_Rate,
        3
      ),

    Mean_Iterations =
      round(
        Mean_Iterations,
        1
      )
  ) %>%

  dplyr::arrange(
    q,
    method
  )


###############################################################################
# 9. CREATE MAIN MANUSCRIPT TABLE 9
###############################################################################
#
# The simulation was run for:
#
#   q = 4, 10, 20, 30, 40, 50
#
# Main Table 9 displays:
#
#   q = 4, 20, 30, 50
#
# Therefore q = 10 and q = 40 are retained in the simulation but are not
# displayed in the manuscript table.
###############################################################################

table9_q <- c(
  4,
  20,
  30,
  50
)


table9_method_order <- c(
  "CE",
  "EL",
  "ET",
  "HD",
  "SL"
)


table9 <-
  A3_summary %>%

  dplyr::filter(
    q %in% table9_q
  ) %>%

  dplyr::mutate(

    q =
      factor(
        q,
        levels = table9_q
      ),

    method =
      factor(
        method,
        levels = table9_method_order
      ),

    Failure_Rate =
      round(
        Failure_Rate,
        3
      ),

    RMSE =
      round(
        RMSE,
        4
      ),

    Median_Log10_Condition =
      round(
        Median_Log10_Condition,
        3
      )
  ) %>%

  dplyr::arrange(
    q,
    method
  ) %>%

  dplyr::select(

    q,

    Entropy =
      method,

    Failure_Rate,

    RMSE,

    Median_Log10_Condition
  )


###############################################################################
# 10. SAVE MAIN TABLE 9
###############################################################################

write.csv(
  table9,
  "Causal Inference/results/table9/table9_highdim_stress.csv",
  row.names = FALSE
)


###############################################################################
# 11. SAVE FULL A3 SUMMARY
###############################################################################

write.csv(
  A3_table,
  "Causal Inference/results/table9/table9_highdim_stress_diagnostics.csv",
  row.names = FALSE
)


###############################################################################
# 12. PRINT MAIN TABLE 9
###############################################################################

cat(
  "\n============================================================\n"
)

cat(
  "MAIN TABLE 9: HIGH-DIMENSIONAL CALIBRATION STRESS TEST\n"
)

cat(
  "N = 2000, M = 1000, OR2PS1\n"
)

cat(
  "============================================================\n\n"
)


print(
  table9,
  row.names = FALSE
)
