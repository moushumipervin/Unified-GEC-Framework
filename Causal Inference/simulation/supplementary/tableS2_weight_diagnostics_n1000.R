rm(list = ls())

source("Causal Inference/functions/ate_functions.R")

results <- readRDS(
  "Causal Inference/results/table3/table3_raw_simulations.rds"
)

summary_table <- make_table(
  results,
  methods = c(
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
)

weight_diagnostics_n1000 <- summary_table |>
  dplyr::filter(
    Method %in% c("IPW", "ET", "HD", "CE")
  ) |>
  dplyr::select(
    Scenario,
    Method,
    ESS1 = Mean_ESS_Treated,
    ESS0 = Mean_ESS_Control,
    MaxW1 = Mean_MaxW_Treated,
    MaxW0 = Mean_MaxW_Control
  ) |>
  dplyr::mutate(
    ESS1 = round(ESS1, 1),
    ESS0 = round(ESS0, 1),
    MaxW1 = round(MaxW1, 2),
    MaxW0 = round(MaxW0, 2)
  )

dir.create(
  "Causal Inference/results/supplementary",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  weight_diagnostics_n1000,
  "Causal Inference/results/supplementary/weight_diagnostics_n1000.csv",
  row.names = FALSE
)

print(weight_diagnostics_n1000)
