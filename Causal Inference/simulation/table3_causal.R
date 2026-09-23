# =============================================================================
# Table 3: Causal inference simulation
# =============================================================================
#
# This script reproduces Table 3 of the manuscript.
#
# Scenarios:
#   OR1PS1
#   OR1PS2
#   OR2PS1
#   OR2PS2
#
# Methods:
#   IPW, EBPS, oCBPS, CBPS, EBCW,
#   AIPW_LM, AIPW_GAM, ET, HD, CE
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

n <- 1000
p <- 4
m <- 1000
K <- 4

methods <- c(
  "IPW",
  "EBPS",
  "oCBPS",
  "CBPS",
  "EBCW",
  "AIPW_LM",
  "AIPW_GAM",
  "ET",
  "HD",
  "CE"
)


# -----------------------------------------------------------------------------
# Run all four scenarios
# -----------------------------------------------------------------------------

results <- run_all_scenarios(
  n = n,
  p = p,
  m = m,
  K = K,
  run_lm = TRUE,
  run_gam = TRUE,
  gec_from = "gam",
  progress = TRUE,
  run_gec = TRUE
)


# -----------------------------------------------------------------------------
# Monte Carlo summary
# -----------------------------------------------------------------------------

summary_table <- make_table(
  results,
  methods = methods
)


# =============================================================================
# Table 3
# =============================================================================

table3_full <- summary_table |>
  dplyr::select(
    Scenario,
    Method,
    Bias,
    SD = MC_SD,
    RMSE,
    Coverage = Coverage_Analytic,
    Failures = N_Failure,
    ESS1 = Mean_ESS_Treated,
    ESS0 = Mean_ESS_Control,
    MaxW1 = Mean_MaxW_Treated,
    MaxW0 = Mean_MaxW_Control
  ) |>
  dplyr::mutate(

    Bias = round(Bias, 4),
    SD = round(SD, 4),
    RMSE = round(RMSE, 4),
    Coverage = round(Coverage, 3),

    Failures = as.integer(Failures),

    ESS1 = round(ESS1, 1),
    ESS0 = round(ESS0, 1),

    MaxW1 = round(MaxW1, 2),
    MaxW0 = round(MaxW0, 2),

    Method = factor(
      Method,
      levels = methods
    ),

    Scenario = factor(
      Scenario,
      levels = c(
        "OR1PS1",
        "OR1PS2",
        "OR2PS1",
        "OR2PS2"
      )
    )
  ) |>
  dplyr::arrange(
    Scenario,
    Method
  )


# -----------------------------------------------------------------------------
# Full diagnostic table
# -----------------------------------------------------------------------------

table3_diagnostics <- summary_table |>
  dplyr::select(
    Scenario,
    Method,
    Bias,
    MC_SD,
    Avg_SE = Avg_SE_Analytic,
    SE_Ratio = SEratio_Analytic_MC,
    RMSE,
    Coverage = Coverage_Analytic,
    Failures = N_Failure,
    ESS1 = Mean_ESS_Treated,
    ESS0 = Mean_ESS_Control,
    MaxW1 = Mean_MaxW_Treated,
    MaxW0 = Mean_MaxW_Control
  ) |>
  dplyr::mutate(

    Bias = round(Bias, 4),
    MC_SD = round(MC_SD, 4),
    Avg_SE = round(Avg_SE, 4),
    SE_Ratio = round(SE_Ratio, 3),
    RMSE = round(RMSE, 4),
    Coverage = round(Coverage, 3),

    Failures = as.integer(Failures),

    ESS1 = round(ESS1, 1),
    ESS0 = round(ESS0, 1),

    MaxW1 = round(MaxW1, 2),
    MaxW0 = round(MaxW0, 2),

    Method = factor(
      Method,
      levels = methods
    ),

    Scenario = factor(
      Scenario,
      levels = c(
        "OR1PS1",
        "OR1PS2",
        "OR2PS1",
        "OR2PS2"
      )
    )
  ) |>
  dplyr::arrange(
    Scenario,
    Method
  )

table3 <- summary_table |>
  dplyr::select(
    Scenario,
    Method,
    Bias,
    MC_SD,
    RMSE,
    Coverage = Coverage_Analytic,
    M_fail = N_Failure
  )
# -----------------------------------------------------------------------------
# Save outputs
# -----------------------------------------------------------------------------
write.csv(
  table3,
  "Causal Inference/results/table3/table3_causal.csv",
  row.names = FALSE
)

write.csv(
  table3_diagnostics,
  "Causal Inference/results/table3/table3_causal_diagnostics.csv",
  row.names = FALSE
)

saveRDS(
  results,
  "Causal Inference/results/table3/table3_raw_simulations.rds"
)


# -----------------------------------------------------------------------------
# Display table
# -----------------------------------------------------------------------------

print(table3)
