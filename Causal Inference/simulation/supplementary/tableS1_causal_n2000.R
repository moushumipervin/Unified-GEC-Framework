# =============================================================================
# Supplementary Table S1: Causal inference simulation, N = 2000
# =============================================================================
#
# This script reproduces Supplementary Table S1.
#
# Sample size:
#   N = 2000
#
# Monte Carlo replications:
#   M = 1000
#
# Scenarios:
#   OR1PS1
#   OR1PS2
#   OR2PS1
#   OR2PS2
#
# Methods reported in Table S1:
#   ET, HD, CE
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

n <- 2000
p <- 4
m <- 1000
K <- 4


# All methods are computed by the simulation framework
methods_all <- c(
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


# Methods reported in Supplementary Table S1
methods_S1 <- c(
  "ET",
  "HD",
  "CE"
)


# -----------------------------------------------------------------------------
# Run all four scenarios
# -----------------------------------------------------------------------------

results_S1 <- run_all_scenarios(
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

summary_S1 <- make_table(
  results_S1,
  methods = methods_all
)


# =============================================================================
# Supplementary Table S1
# =============================================================================

table_S1 <- summary_S1 |>
  dplyr::filter(
    Method %in% methods_S1
  ) |>
  dplyr::select(
    Scenario,
    Method,
    Bias,
    RMSE,
    Coverage = Coverage_Analytic
  ) |>
  dplyr::mutate(

    Bias = round(Bias, 4),

    RMSE = round(RMSE, 4),

    Coverage = round(Coverage, 3),

    Method = factor(
      Method,
      levels = methods_S1
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
# Save outputs
# -----------------------------------------------------------------------------

dir.create(
  "Causal Inference/results/supplementary",
  recursive = TRUE,
  showWarnings = FALSE
)


write.csv(
  table_S1,
  "Causal Inference/results/supplementary/tableS1_causal_n2000.csv",
  row.names = FALSE
)


saveRDS(
  results_S1,
  "Causal Inference/results/supplementary/tableS1_raw_simulations_n2000.rds"
)


# -----------------------------------------------------------------------------
# Display Supplementary Table S1
# -----------------------------------------------------------------------------

print(table_S1)
