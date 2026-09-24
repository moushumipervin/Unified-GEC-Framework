################################################################################
# SUPPLEMENTARY TABLE S3
# MCAR simulation results for the semi-supervised regression setting under OR2
#
# This script uses the Monte Carlo results generated for Figure 1.
#
# It reports:
#   Bias
#   Monte Carlo SD
#   Average analytic SE
#   Analytic SE / Monte Carlo SD
#   RMSE
#   Empirical coverage
#   ARE relative to the supervised estimator
#
# Methods shown in Supplementary Table S3:
#   Supervised
#   DRESS
#   PSSE
#   ET
#   HD
#   CE
#
# Run from the root directory of:
#   Unified-GEC-Framework/
################################################################################


rm(list = ls())


################################################################################
# 1. REQUIRED PACKAGES
################################################################################

required_packages <- c(
  "dplyr"
)


missing_packages <-
  required_packages[
    !vapply(
      required_packages,
      requireNamespace,
      quietly = TRUE,
      FUN.VALUE = logical(1)
    )
  ]


if (length(missing_packages) > 0) {

  install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org"
  )
}


suppressPackageStartupMessages({

  library(dplyr)

})


################################################################################
# 2. READ RAW MONTE CARLO RESULTS FROM FIGURE 1
################################################################################

simulation_results_df <-
  readRDS(
    "Semi Supervised/results/figure1/figure1_mcar_or2_raw.rds"
  )


################################################################################
# 3. CHECK REQUIRED VARIABLES
################################################################################

required_columns <- c(
  "method",
  "parameter",
  "estimate",
  "truth",
  "analytic_se",
  "covered"
)


missing_columns <-
  setdiff(
    required_columns,
    names(simulation_results_df)
  )


if (length(missing_columns) > 0) {

  stop(
    paste0(
      "The raw Monte Carlo file is missing required columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
}


################################################################################
# 4. METHODS INCLUDED IN SUPPLEMENTARY TABLE S3
################################################################################

tableS3_methods <- c(
  "Supervised",
  "DRESS",
  "PSSE",
  "ET",
  "HD",
  "CE"
)


################################################################################
# 5. PARAMETER ORDER
################################################################################

parameter_order <- c(
  "(Intercept)",
  "X1",
  "X2",
  "X3",
  "X4"
)


################################################################################
# 6. FILTER TO METHODS USED IN TABLE S3
################################################################################

tableS3_raw <-
  simulation_results_df %>%

  dplyr::filter(
    method %in%
      tableS3_methods
  )


################################################################################
# 7. MONTE CARLO SUMMARY
################################################################################

tableS3_summary <-
  tableS3_raw %>%

  dplyr::group_by(
    method,
    parameter
  ) %>%

  dplyr::summarise(

    N_success =
      sum(
        is.finite(
          estimate
        )
      ),

    Mean_Estimate =
      mean(
        estimate,
        na.rm = TRUE
      ),

    Bias =
      mean(
        estimate - truth,
        na.rm = TRUE
      ),

    MC_SD =
      sd(
        estimate,
        na.rm = TRUE
      ),

    Ana_SE =
      mean(
        analytic_se,
        na.rm = TRUE
      ),

    SE_MC =
      Ana_SE /
      MC_SD,

    RMSE =
      sqrt(
        mean(
          (estimate - truth)^2,
          na.rm = TRUE
        )
      ),

    Coverage =
      mean(
        covered,
        na.rm = TRUE
      ),

    .groups =
      "drop"
  )


################################################################################
# 8. RELATIVE EFFICIENCY
#
# ARE = Var(Supervised) / Var(Method)
#     = MC_SD(Supervised)^2 / MC_SD(Method)^2
#
# Therefore:
#   ARE > 1 means the method is more efficient than Supervised.
################################################################################

supervised_reference <-
  tableS3_summary %>%

  dplyr::filter(
    method ==
      "Supervised"
  ) %>%

  dplyr::select(
    parameter,
    Supervised_MC_SD =
      MC_SD
  )


tableS3_summary <-
  tableS3_summary %>%

  dplyr::left_join(
    supervised_reference,
    by = "parameter"
  ) %>%

  dplyr::mutate(

    ARE =
      (
        Supervised_MC_SD^2
      ) /
      (
        MC_SD^2
      )
  )


################################################################################
# 9. ORDER METHODS AND PARAMETERS
################################################################################

tableS3_summary <-
  tableS3_summary %>%

  dplyr::mutate(

    method =
      factor(
        method,
        levels =
          tableS3_methods
      ),

    parameter =
      factor(
        parameter,
        levels =
          parameter_order
      )
  ) %>%

  dplyr::arrange(
    method,
    parameter
  )


################################################################################
# 10. FULL-PRECISION TABLE S3
#
# This is the preferred numerical reproducibility file.
# Do NOT round this file.
################################################################################

tableS3_mcar_or2 <-
  tableS3_summary %>%

  dplyr::select(
    Method = method,
    Parameter = parameter,
    Bias,
    MC_SD,
    Ana_SE,
    SE_MC,
    RMSE,
    Coverage,
    ARE
  )


################################################################################
# 11. PRINT FULL-PRECISION TABLE
################################################################################

print(
  tableS3_mcar_or2,
  n = Inf
)


################################################################################
# 12. DISPLAY VERSION
#
# Manuscript-style rounding.
# The numerical CSV above remains full precision.
################################################################################

tableS3_display <-
  tableS3_mcar_or2 %>%

  dplyr::mutate(

    Bias =
      round(
        Bias,
        3
      ),

    MC_SD =
      round(
        MC_SD,
        3
      ),

    Ana_SE =
      round(
        Ana_SE,
        3
      ),

    SE_MC =
      round(
        SE_MC,
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

    ARE =
      round(
        ARE,
        3
      )
  )


################################################################################
# 13. PRINT MANUSCRIPT DISPLAY VERSION
################################################################################

print(
  tableS3_display,
  n = Inf
)


################################################################################
# 14. CREATE OUTPUT DIRECTORY
################################################################################

dir.create(

  "Semi Supervised/results/supplementary/tableS3",

  recursive =
    TRUE,

  showWarnings =
    FALSE
)


################################################################################
# 15. SAVE FULL-PRECISION TABLE
################################################################################

write.csv(

  tableS3_mcar_or2,

  file =
    "Semi Supervised/results/supplementary/tableS3/tableS3_mcar_or2.csv",

  row.names =
    FALSE
)


################################################################################
# 16. SAVE MANUSCRIPT DISPLAY VERSION
################################################################################

write.csv(

  tableS3_display,

  file =
    "Semi Supervised/results/supplementary/tableS3/tableS3_mcar_or2_display.csv",

  row.names =
    FALSE
)


################################################################################
# 17. SAVE SUMMARY OBJECT
################################################################################

saveRDS(

  tableS3_summary,

  file =
    "Semi Supervised/results/supplementary/tableS3/tableS3_mcar_or2_summary.rds"
)


################################################################################
# 18. FINAL MESSAGE
################################################################################

cat(
  "\n============================================================\n"
)

cat(
  "SUPPLEMENTARY TABLE S3 COMPLETE: OR2 + MCAR\n"
)

cat(
  "============================================================\n"
)

cat(
  "\nFull-precision numerical table:\n"
)

cat(
  "Semi Supervised/results/supplementary/tableS3/tableS3_mcar_or2.csv\n"
)

cat(
  "\nManuscript display table:\n"
)

cat(
  "Semi Supervised/results/supplementary/tableS3/tableS3_mcar_or2_display.csv\n"
)

cat(
  "\nR summary object:\n"
)

cat(
  "Semi Supervised/results/supplementary/tableS3/tableS3_mcar_or2_summary.rds\n"
)

cat(
  "\n============================================================\n"
)
