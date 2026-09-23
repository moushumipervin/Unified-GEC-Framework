# =============================================================================
# Table 4: Score-augmented calibration under OR2PS1
# =============================================================================
#
# This script reproduces Table 4 of the manuscript.
#
# Methods:
#   1. AIPW (estimated PS)
#   2. Standard ET
#   3. ET + PS score
#
# Simulation setting:
#   Scenario = OR2PS1
#   N = 1000
#   Monte Carlo replications = 1000
#
# =============================================================================

rm(list = ls())

# -----------------------------------------------------------------------------
# Source functions
# -----------------------------------------------------------------------------

source("Causal Inference/functions/ate_functions.R")


# -----------------------------------------------------------------------------
# Simulation settings
# -----------------------------------------------------------------------------

B <- 1000
n <- 1000
p <- 4
K <- 4


# -----------------------------------------------------------------------------
# Run Monte Carlo simulation
# -----------------------------------------------------------------------------

table4_results <- vector(
  "list",
  B
)

for (b in seq_len(B)) {

  if (b == 1 || b %% 25 == 0) {
    message(
      "Table 4, OR2PS1: replication ",
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

      NULL
    }
  )
}

table4_results <- dplyr::bind_rows(
  table4_results
)


# -----------------------------------------------------------------------------
# Monte Carlo summary helper
# -----------------------------------------------------------------------------

summarize_table4 <- function(
    estimate,
    se,
    coverage,
    truth = 10,
    B = 1000
) {

  point_ok <- is.finite(estimate)

  inference_ok <-
    point_ok &
    is.finite(se)

  est_point <- estimate[point_ok]

  est_inf <- estimate[inference_ok]

  se_inf <- se[inference_ok]

  cov_inf <- coverage[inference_ok]

  MC_SD <-
    if (length(est_point) > 1) {
      sd(est_point)
    } else {
      NA_real_
    }

  data.frame(

    Bias =
      mean(
        est_point - truth
      ),

    MC_SD =
      MC_SD,

    RMSE =
      sqrt(
        mean(
          (est_point - truth)^2
        )
      ),

    M_fail =
      B - sum(point_ok)
  )
}


# =============================================================================
# Table 4
# =============================================================================

table4 <- dplyr::bind_rows(

  `AIPW (estimated PS)` =
    summarize_table4(
      estimate = table4_results$AIPW,
      se = table4_results$AIPW_SE,
      coverage = table4_results$AIPW_COV,
      B = B
    ),

  `Standard ET` =
    summarize_table4(
      estimate = table4_results$ET,
      se = table4_results$ET_SE,
      coverage = table4_results$ET_COV,
      B = B
    ),

  `ET + PS score` =
    summarize_table4(
      estimate = table4_results$ET_SCORE,
      se = table4_results$ET_SCORE_SE,
      coverage = table4_results$ET_SCORE_COV,
      B = B
    ),

  .id = "Method"
) |>
  dplyr::mutate(
    Bias = round(Bias, 4),
    MC_SD = round(MC_SD, 4),
    RMSE = round(RMSE, 4),
    M_fail = as.integer(M_fail)
  )

# -----------------------------------------------------------------------------
# Save outputs
# -----------------------------------------------------------------------------

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

print(table4)
