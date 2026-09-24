###############################################################################
# SUPPLEMENTARY SEVERE-OVERLAP STRESS TEST
#
# Produces:
#   Table S4: True propensity-score diagnostics
#   Table S5: Complete estimator performance
#   Supplementary figure: Mean maximum treated-arm weight
###############################################################################

rm(list = ls())


###############################################################################
# 1. SOURCE FUNCTIONS
###############################################################################

source(
  "Causal Inference/functions/stress_test_functions.R"
)


###############################################################################
# 2. SETTINGS
###############################################################################

# Desired target treatment rates
A2_rates <- c(
  0.50,
  0.25,
  0.10
)

# Propensity-score slope strength used for severe-overlap stress test
A2_ps_strength <- 2

# Monte Carlo replications
M <- 1000

# Cross-fitting folds
K <- 4

# True ATE under OR2
A2_true_ATE <- 10


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


###############################################################################
# 4. FULL STRESS-TEST GRID
###############################################################################

A2_grid <-
  expand.grid(
    N = c(
      500,
      2000,
      8000
    ),
    target_rate = c(
      0.10,
      0.25,
      0.50
    ),
    stringsAsFactors = FALSE
  )


###############################################################################
# 5. RUN MONTE CARLO SIMULATION
###############################################################################

A2_all_list <-
  vector(
    "list",
    nrow(A2_grid)
  )


for (j in seq_len(nrow(A2_grid))) {

  cat(
    "\n========================================\n",
    "A2 CONDITION\n",
    "N = ",
    A2_grid$N[j],
    ", target rate = ",
    A2_grid$target_rate[j],
    "\n========================================\n",
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
# 6. CONVERT RESULTS TO LONG FORMAT
###############################################################################

A2_long <-
  dplyr::bind_rows(
    lapply(
      A2_methods,
      function(m) {

        make_A2_method_rows(
          A2_all_results,
          m
        )
      }
    )
  )


###############################################################################
# 7. DEFINE VALID REPLICATIONS
###############################################################################

A2_long <-
  A2_long %>%
  dplyr::mutate(

    # Valid point estimate
    est_ok =
      is.finite(
        Estimate
      ),

    # Valid analytic inference
    analytic_ok =
      is.finite(
        Estimate
      ) &
      is.finite(
        SE
      ) &
      Success == 1,

    # Valid treated-arm weight diagnostics
    weight1_ok =
      is.finite(
        ESS1
      ) &
      is.finite(
        MAXW1
      ),

    # Valid control-arm weight diagnostics
    weight0_ok =
      is.finite(
        ESS0
      ) &
      is.finite(
        MAXW0
      )
  )


###############################################################################
# 8. COMPLETE A2 SUMMARY
###############################################################################

A2_summary <-
  A2_long %>%

  dplyr::group_by(
    N,
    target_rate,
    method
  ) %>%

  dplyr::summarise(

    ###########################################################################
    # Number attempted
    ###########################################################################

    N_total =
      M,


    ###########################################################################
    # Number with finite point estimates
    ###########################################################################

    N_est =
      sum(
        est_ok,
        na.rm = TRUE
      ),


    ###########################################################################
    # Number with valid analytic inference
    ###########################################################################

    N_valid =
      sum(
        analytic_ok,
        na.rm = TRUE
      ),


    ###########################################################################
    # Point-estimation failure rate
    ###########################################################################

    Est_Failure_Rate =
      1 -
      N_est / M,


    ###########################################################################
    # Analytic-inference failure rate
    ###########################################################################

    Failure_Rate =
      1 -
      N_valid / M,


    ###########################################################################
    # Average realized treatment rate
    ###########################################################################

    Mean_Treat_Rate =
      mean(
        TREAT_RATE,
        na.rm = TRUE
      ),


    ###########################################################################
    # Mean true propensity
    ###########################################################################

    Mean_True_PS =
      mean(
        TRUE_PS_MEAN,
        na.rm = TRUE
      ),


    ###########################################################################
    # Mean estimate
    ###########################################################################

    Mean_Estimate =
      mean(
        Estimate[est_ok],
        na.rm = TRUE
      ),


    ###########################################################################
    # Bias
    ###########################################################################

    Bias =
      mean(
        Estimate[est_ok] -
          A2_true_ATE,
        na.rm = TRUE
      ),


    ###########################################################################
    # Monte Carlo SD
    ###########################################################################

    MC_SD =
      sd(
        Estimate[est_ok],
        na.rm = TRUE
      ),


    ###########################################################################
    # RMSE
    ###########################################################################

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


    ###########################################################################
    # Mean analytic standard error
    ###########################################################################

    Avg_SE =
      mean(
        SE[analytic_ok],
        na.rm = TRUE
      ),


    ###########################################################################
    # MC SD among analytic-valid replications
    ###########################################################################

    MC_SD_Analytic =
      sd(
        Estimate[analytic_ok],
        na.rm = TRUE
      ),


    ###########################################################################
    # SE / MC SD ratio
    ###########################################################################

    SE_MC_Ratio =
      Avg_SE /
      MC_SD_Analytic,


    ###########################################################################
    # Coverage
    ###########################################################################

    Coverage =
      mean(
        Coverage[analytic_ok],
        na.rm = TRUE
      ),


    ###########################################################################
    # Treated-arm ESS
    ###########################################################################

    Mean_ESS1 =
      mean(
        ESS1[weight1_ok],
        na.rm = TRUE
      ),


    ###########################################################################
    # Control-arm ESS
    ###########################################################################

    Mean_ESS0 =
      mean(
        ESS0[weight0_ok],
        na.rm = TRUE
      ),


    ###########################################################################
    # Mean maximum treated-arm weight
    ###########################################################################

    Mean_MAXW1 =
      mean(
        MAXW1[weight1_ok],
        na.rm = TRUE
      ),


    ###########################################################################
    # Mean maximum control-arm weight
    ###########################################################################

    Mean_MAXW0 =
      mean(
        MAXW0[weight0_ok],
        na.rm = TRUE
      ),


    ###########################################################################
    # Median maximum treated-arm weight
    ###########################################################################

    Median_MAXW1 =
      median(
        MAXW1[weight1_ok],
        na.rm = TRUE
      ),


    ###########################################################################
    # Median maximum control-arm weight
    ###########################################################################

    Median_MAXW0 =
      median(
        MAXW0[weight0_ok],
        na.rm = TRUE
      ),


    ###########################################################################
    # Worst maximum treated-arm weight
    ###########################################################################

    Worst_MAXW1 =
      max(
        MAXW1[weight1_ok],
        na.rm = TRUE
      ),


    ###########################################################################
    # Worst maximum control-arm weight
    ###########################################################################

    Worst_MAXW0 =
      max(
        MAXW0[weight0_ok],
        na.rm = TRUE
      ),


    .groups =
      "drop"
  )


###############################################################################
# 9. CREATE OUTPUT DIRECTORIES
###############################################################################

dir.create(
  "Causal Inference/results/supplementary/overlap_stress",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "Causal Inference/results/supplementary/overlap_stress/figures",
  recursive = TRUE,
  showWarnings = FALSE
)


###############################################################################
# 10. SAVE RAW MONTE CARLO RESULTS
###############################################################################

saveRDS(
  A2_all_results,
  "Causal Inference/results/supplementary/overlap_stress/overlap_stress_raw.rds"
)


###############################################################################
# 11. SUPPLEMENTARY TABLE S4
# TRUE PROPENSITY-SCORE DIAGNOSTICS
###############################################################################

tableS4 <-
  A2_all_results %>%

  dplyr::group_by(
    N_DESIGN,
    RATE_DESIGN
  ) %>%

  dplyr::summarise(

    Mean_PS =
      mean(
        TRUE_PS_MEAN,
        na.rm = TRUE
      ),

    Mean_Min_PS =
      mean(
        TRUE_PS_MIN,
        na.rm = TRUE
      ),

    Mean_P001 =
      mean(
        TRUE_PS_P001,
        na.rm = TRUE
      ),

    Mean_P01 =
      mean(
        TRUE_PS_P01,
        na.rm = TRUE
      ),

    Mean_P05 =
      mean(
        TRUE_PS_P05,
        na.rm = TRUE
      ),

    .groups =
      "drop"
  ) %>%

  dplyr::mutate(

    RATE_DESIGN =
      factor(
        RATE_DESIGN,
        levels = c(
          0.50,
          0.25,
          0.10
        )
      )
  ) %>%

  dplyr::arrange(
    RATE_DESIGN,
    N_DESIGN
  ) %>%

  dplyr::rename(
    N =
      N_DESIGN,
    Rate =
      RATE_DESIGN
  )


write.csv(
  tableS4,
  "Causal Inference/results/supplementary/overlap_stress/tableS4_propensity_diagnostics.csv",
  row.names = FALSE
)


###############################################################################
# 12. SUPPLEMENTARY TABLE S5
# COMPLETE ESTIMATOR PERFORMANCE
###############################################################################

tableS5 <-
  A2_summary %>%

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
          "AIPW_GAM",
          "AIPW_LM",
          "SL",
          "EL",
          "ET",
          "HD",
          "CE"
        )
      )
  ) %>%

  dplyr::arrange(
    target_rate,
    N,
    method
  ) %>%

  dplyr::mutate(

    Method =
      dplyr::recode(
        as.character(method),

        "IPW" =
          "IPW",

        "AIPW_GAM" =
          "AIPW-GAM",

        "AIPW_LM" =
          "AIPW-LM",

        "SL" =
          "SL",

        "EL" =
          "EL",

        "ET" =
          "ET",

        "HD" =
          "HD",

        "CE" =
          "CE"
      ),

    Bias =
      round(
        Bias,
        4
      ),

    RMSE =
      round(
        RMSE,
        4
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
      ),

    MaxW1 =
      round(
        Mean_MAXW1,
        1
      ),

    Failure =
      round(
        Failure_Rate,
        3
      )
  ) %>%

  dplyr::select(
    N,
    Rate =
      target_rate,
    Method,
    Bias,
    RMSE,
    Coverage,
    ESS1,
    MaxW1,
    Failure
  )


write.csv(
  tableS5,
  "Causal Inference/results/supplementary/overlap_stress/tableS5_complete_performance.csv",
  row.names = FALSE
)


###############################################################################
# 13. SUPPLEMENTARY FIGURE
# MEAN MAXIMUM TREATED-ARM WEIGHT
###############################################################################

A2_plot_data <-
  A2_summary %>%

  dplyr::filter(
    !method %in% c(
      "AIPW_LM",
      "AIPW_GAM"
    )
  )


figureS_overlap <-
  ggplot(
    A2_plot_data,
    aes(
      x = N,
      y = Mean_MAXW1,
      group = method,
      color = method,
      linetype = method,
      shape = method
    )
  ) +

  geom_line(
    linewidth = 1
  ) +

  geom_point(
    size = 3
  ) +

  facet_wrap(
    ~ target_rate,
    scales = "fixed",
    labeller = label_both
  ) +

  scale_color_manual(
    values = c(
      "IPW" = "#000000",
      "SL"  = "#6A3D9A",
      "EL"  = "#1B9E77",
      "ET"  = "#D95F02",
      "HD"  = "#377EB8",
      "CE"  = "#A6761D"
    )
  ) +

  scale_shape_manual(
    values = c(
      "IPW" = 16,
      "SL"  = 18,
      "EL"  = 17,
      "ET"  = 15,
      "HD"  = 3,
      "CE"  = 8
    )
  ) +

  labs(
    x =
      "Sample size",

    y =
      "Mean maximum treated-arm weight",

    color =
      "Method",

    linetype =
      "Method",

    shape =
      "Method"
  ) +

  theme_bw() +

  theme(
    legend.position =
      "right",

    panel.grid.minor =
      element_blank(),

    strip.text =
      element_text(
        face = "bold"
      )
  )


print(
  figureS_overlap
)


ggsave(
  filename =
    "Causal Inference/results/supplementary/overlap_stress/figures/figureS_overlap_max_weight.pdf",

  plot =
    figureS_overlap,

  width =
    8,

  height =
    5
)


###############################################################################
# 14. PRINT TABLES
###############################################################################

cat(
  "\n============================================================\n"
)

cat(
  "SUPPLEMENTARY TABLE S4\n"
)

cat(
  "TRUE PROPENSITY-SCORE DIAGNOSTICS\n"
)

cat(
  "============================================================\n\n"
)

print(
  tableS4,
  n = Inf
)


cat(
  "\n============================================================\n"
)

cat(
  "SUPPLEMENTARY TABLE S5\n"
)

cat(
  "COMPLETE ESTIMATOR PERFORMANCE\n"
)

cat(
  "============================================================\n\n"
)

print(
  tableS5,
  n = Inf
)
