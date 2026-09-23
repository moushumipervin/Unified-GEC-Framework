###############################################################################
# TABLE 4: SCORE-AUGMENTED CALIBRATION UNDER OR2PS1
###############################################################################

rm(list = ls())

source("Causal Inference/functions/ate_functions.R")


###############################################################################
# SIMULATION SETTINGS
###############################################################################

B <- 1000
n <- 1000
p <- 4
K <- 4


###############################################################################
# RUN SIMULATION
###############################################################################

table4_results <- vector(
  "list",
  B
)

for (b in seq_len(B)) {

  if (
    b == 1 ||
    b %% 25 == 0
  ) {

    message(
      "Table 4 OR2PS1: replication ",
      b,
      " / ",
      B
    )
  }


  table4_results[[b]] <- tryCatch(

    run_one_A6_ET(
      rep_id = b,
      n = n,
      p = p,
      K = K
    ),

    error = function(e) {

      message(
        "Table 4 replication ",
        b,
        " failed: ",
        e$message
      )

      data.frame(
        AIPW = NA_real_,
        AIPW_SE = NA_real_,
        AIPW_COV = NA_real_,
        ET = NA_real_,
        ET_SE = NA_real_,
        ET_COV = NA_real_,
        ET_SCORE = NA_real_,
        ET_SCORE_SE = NA_real_,
        ET_SCORE_COV = NA_real_
      )
    }
  )
}


table4_results <-
  dplyr::bind_rows(
    table4_results
  )


###############################################################################
# SUMMARY FUNCTION
###############################################################################

summarize_table4 <- function(
    estimate,
    se,
    coverage,
    truth = 10
) {

  point_ok <-
    is.finite(
      estimate
    )

  inference_ok <-
    point_ok &
    is.finite(se)

  est_point <-
    estimate[
      point_ok
    ]

  se_inf <-
    se[
      inference_ok
    ]

  cov_inf <-
    coverage[
      inference_ok
    ]

  MC_SD <-
    if (length(est_point) > 1) {
      sd(est_point)
    } else {
      NA_real_
    }


  data.frame(

    Mean_Estimate =
      if (length(est_point) > 0) {
        mean(est_point)
      } else {
        NA_real_
      },

    Bias =
      if (length(est_point) > 0) {
        mean(
          est_point - truth
        )
      } else {
        NA_real_
      },

    MC_SD =
      MC_SD,

    RMSE =
      if (length(est_point) > 0) {
        sqrt(
          mean(
            (
              est_point -
                truth
            )^2
          )
        )
      } else {
        NA_real_
      },

    Mean_SE =
      if (length(se_inf) > 0) {
        mean(se_inf)
      } else {
        NA_real_
      },

    SE_MC_Ratio =
      if (
        is.finite(MC_SD) &&
        MC_SD > 0 &&
        length(se_inf) > 0
      ) {
        mean(se_inf) / MC_SD
      } else {
        NA_real_
      },

    Coverage =
      if (length(cov_inf) > 0) {
        mean(
          cov_inf,
          na.rm = TRUE
        )
      } else {
        NA_real_
      },

    Point_Success =
      sum(point_ok),

    Inference_Success =
      sum(inference_ok),

    M_fail =
      sum(!point_ok)
  )
}


###############################################################################
# FULL SUMMARY
###############################################################################

table4_summary <- dplyr::bind_rows(

  `AIPW (estimated PS)` =
    summarize_table4(
      estimate =
        table4_results$AIPW,

      se =
        table4_results$AIPW_SE,

      coverage =
        table4_results$AIPW_COV
    ),


  `Standard ET` =
    summarize_table4(
      estimate =
        table4_results$ET,

      se =
        table4_results$ET_SE,

      coverage =
        table4_results$ET_COV
    ),


  `ET + PS score` =
    summarize_table4(
      estimate =
        table4_results$ET_SCORE,

      se =
        table4_results$ET_SCORE_SE,

      coverage =
        table4_results$ET_SCORE_COV
    ),

  .id = "Method"
)


###############################################################################
# MANUSCRIPT TABLE 4
###############################################################################

table4 <-
  table4_summary |>
  dplyr::select(
    Method,
    Bias,
    MC_SD,
    RMSE,
    M_fail
  ) |>
  dplyr::mutate(
    Bias = round(Bias, 4),
    MC_SD = round(MC_SD, 4),
    RMSE = round(RMSE, 4),
    M_fail = as.integer(M_fail)
  )


###############################################################################
# SAVE RESULTS
###############################################################################

dir.create(
  "Causal Inference/results/table4",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  table4,
  "Causal Inference/results/table4/table4_score_augmented.csv",
  row.names = FALSE
)

saveRDS(
  table4_results,
  "Causal Inference/results/table4/table4_raw_simulations.rds"
)


###############################################################################
# PRINT TABLE 4
###############################################################################

print(table4)
