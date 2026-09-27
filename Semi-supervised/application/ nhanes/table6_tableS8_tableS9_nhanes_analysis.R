###############################################################################
# TABLE 6 + TABLE S8 + TABLE S9
# NHANES SEMI-SUPERVISED APPLICATION
#
# This script reproduces:
#
#   Table 6  : Compact repeated-split NHANES summary
#   Table S8 : Single NHANES label-deletion split
#   Table S9 : Full repeated-split coefficient-level NHANES results
#
# Repeated-split design:
#   B = 1000 label-deletion splits
#   50% of the originally observed valid fasting-glucose labels are retained
#   in each split.
#
# Repository structure:
#
#   Semi-supervised/
#   ├── functions/
#   │   └── nhanes_functions.R
#   ├── application/
#   │   └── nhanes/
#   │       └── table6_tableS8_tableS9_nhanes_analysis.R
#   └── results/
#       └── nhanes/
#
# Run this script from the ROOT of the GitHub repository.
###############################################################################


###############################################################################
# 0. CLEAN WORKSPACE
###############################################################################

rm(list = ls())


###############################################################################
# 1. REQUIRED PACKAGES
###############################################################################

required_packages <- c(
  "dplyr",
  "tidyr",
  "sandwich",
  "MASS",
  "mgcv",
  "CVXR",
  "caret",
  "nhanesA"
)


missing_packages <- required_packages[
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
  library(tidyr)
  library(sandwich)
  library(MASS)
  library(mgcv)
  library(CVXR)
  library(caret)
  library(nhanesA)

})


###############################################################################
# 2. SOURCE NHANES FUNCTIONS
###############################################################################

source(
  "Semi-supervised/functions/nhanes_functions.R"
)


###############################################################################
# 3. SETTINGS
###############################################################################

B <- 1000

label_fraction <- 0.50

base_seed <- 123

single_split_seed <- 124


method_order <- c(
  "Supervised",
  "DRESS",
  "PSSE",
  "ET",
  "HD",
  "CE"
)


###############################################################################
# 4. OUTPUT PATHS
###############################################################################

results_root <-
  "Semi-supervised/results/nhanes"


table6_dir <-
  file.path(
    results_root,
    "table6"
  )


tableS8_dir <-
  file.path(
    results_root,
    "tableS8"
  )


tableS9_dir <-
  file.path(
    results_root,
    "tableS9"
  )


diagnostics_dir <-
  file.path(
    results_root,
    "diagnostics"
  )


###############################################################################
# 5. CREATE OUTPUT DIRECTORIES
###############################################################################

dir.create(
  table6_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  tableS8_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  tableS9_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  diagnostics_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


###############################################################################
# 6. DOWNLOAD NHANES 2017-2018 DATA
###############################################################################

demo <-
  nhanesA::nhanes(
    "DEMO_J"
  )


bmx <-
  nhanesA::nhanes(
    "BMX_J"
  )


bpx <-
  nhanesA::nhanes(
    "BPX_J"
  )


glu <-
  nhanesA::nhanes(
    "GLU_J"
  ) %>%

  dplyr::select(
    SEQN,
    Y = LBXGLU,
    WTSAF2YR
  )


###############################################################################
# 7. MERGE NHANES COMPONENTS
###############################################################################

dat <-
  demo %>%

  dplyr::transmute(

    SEQN,

    age =
      RIDAGEYR,

    sex =
      factor(
        RIAGENDR,
        labels = c(
          "Male",
          "Female"
        )
      ),

    race =
      factor(
        RIDRETH3
      )
  ) %>%

  dplyr::left_join(

    bmx %>%
      dplyr::select(
        SEQN,
        BMI = BMXBMI
      ),

    by = "SEQN"
  ) %>%

  dplyr::left_join(

    bpx %>%
      dplyr::select(
        SEQN,
        SBP = BPXSY1,
        DBP = BPXDI1
      ),

    by = "SEQN"
  ) %>%

  dplyr::left_join(

    glu,

    by = "SEQN"
  )


###############################################################################
# 8. DEFINE ANALYTIC COHORT
#
# 1. Restrict to age >= 12.
#
# 2. Retain participants with complete covariates.
#
# 3. A glucose value is considered a valid fasting-glucose label only when
#    WTSAF2YR > 0.
#
# 4. Participants with a released glucose value who do not satisfy the
#    fasting criterion remain in the analytic cohort but are treated as
#    unlabeled.
###############################################################################

datX <-
  dat %>%

  dplyr::filter(

    age >= 12,

    !dplyr::if_any(

      c(
        age,
        sex,
        race,
        BMI,
        SBP,
        DBP
      ),

      is.na
    )
  ) %>%

  dplyr::mutate(

    Y =
      dplyr::if_else(

        !is.na(Y) &
          !is.na(WTSAF2YR) &
          WTSAF2YR > 0,

        Y,

        NA_real_
      )
  ) %>%

  dplyr::select(

    SEQN,

    age,
    sex,
    race,
    BMI,
    SBP,
    DBP,

    Y
  )


###############################################################################
# 9. ANALYTIC SAMPLE CHECKS
###############################################################################

cat(
  "\n============================================================\n"
)

cat(
  "NHANES ANALYTIC SAMPLE\n"
)

cat(
  "============================================================\n"
)


cat(
  "Total N: ",
  nrow(datX),
  "\n",
  sep = ""
)


cat(
  "Originally labeled: ",
  sum(
    !is.na(
      datX$Y
    )
  ),
  "\n",
  sep = ""
)


cat(
  "Originally unlabeled: ",
  sum(
    is.na(
      datX$Y
    )
  ),
  "\n",
  sep = ""
)


cat(
  "Unlabeled fraction: ",
  round(
    mean(
      is.na(
        datX$Y
      )
    ),
    4
  ),
  "\n\n",
  sep = ""
)


###############################################################################
# 10. RACE / ETHNICITY REFERENCE CATEGORY CHECK
#
# Mexican American must be the reference category.
###############################################################################

if (!is.factor(datX$race)) {

  stop(
    "NHANES race variable is not a factor."
  )
}


if (
  levels(datX$race)[1] !=
  "Mexican American"
) {

  stop(
    paste0(
      "Race reference category is ",
      levels(datX$race)[1],
      "; expected Mexican American."
    )
  )
}


cat(
  "NHANES race reference category confirmed: ",
  levels(datX$race)[1],
  "\n\n",
  sep = ""
)


###############################################################################
# 11. ENCODE COVARIATES
###############################################################################

encode_data <- function(df) {

  Xmm <-
    model.matrix(

      ~ age +
        sex +
        race +
        BMI +
        SBP +
        DBP,

      data = df

    )[, -1, drop = FALSE]


  Xmm <-
    as.data.frame(
      Xmm
    )


  colnames(Xmm) <-
    gsub(
      "[^A-Za-z0-9_]",
      "_",
      colnames(Xmm)
    )


  colnames(Xmm) <-
    gsub(
      "_+",
      "_",
      colnames(Xmm)
    )


  colnames(Xmm) <-
    gsub(
      "_$",
      "",
      colnames(Xmm)
    )


  Xmm
}


###############################################################################
# 12. FULL-DATA SUPERVISED BENCHMARK
#
# Uses every originally available valid fasting-glucose observation.
###############################################################################

full_labeled <-
  datX %>%

  dplyr::filter(
    !is.na(Y)
  )


full_labeled_final <-
  cbind(

    Y =
      full_labeled$Y,

    encode_data(
      full_labeled
    )
  )


full_formula <-
  as.formula(

    paste(

      "Y ~",

      paste(

        sprintf(
          "`%s`",
          colnames(
            full_labeled_final
          )[-1]
        ),

        collapse = " + "
      )
    )
  )


full_fit <-
  lm(

    full_formula,

    data =
      as.data.frame(
        full_labeled_final
      )
  )


theta_full <-
  coef(
    full_fit
  )


parameter_names <-
  names(
    theta_full
  )


cat(
  "\nFULL-DATA SUPERVISED BENCHMARK:\n"
)


print(
  theta_full
)


###############################################################################
# 13. RUN B = 1000 LABEL-DELETION SPLITS
#
# run_one_split() must be defined in:
#
#   Semi-supervised/functions/nhanes_functions.R
#
###############################################################################

all_results <-
  vector(
    "list",
    B
  )


for (b in seq_len(B)) {

  cat(
    "\rRunning split ",
    b,
    " of ",
    B,
    sep = ""
  )


  all_results[[b]] <-

    tryCatch(

      run_one_split(

        split_seed =
          base_seed + b,

        datX =
          datX
      ),

      error = function(e) {

        message(
          "\nSplit ",
          base_seed + b,
          " failed: ",
          conditionMessage(e)
        )

        NULL
      }
    )
}


cat(
  "\n"
)


###############################################################################
# 14. COMBINE SUCCESSFUL SPLITS
###############################################################################

results_1000 <-
  dplyr::bind_rows(
    all_results
  )


if (
  nrow(
    results_1000
  ) == 0
) {

  stop(
    "No label-deletion split completed successfully."
  )
}


cat(
  "\nSuccessful split results available for ",
  length(
    unique(
      results_1000$split
    )
  ),
  " split(s).\n",
  sep = ""
)


###############################################################################
# 15. SAVE RAW REPEATED-SPLIT RESULTS
###############################################################################

saveRDS(

  results_1000,

  file =
    file.path(
      diagnostics_dir,
      "nhanes_repeated_split_raw.rds"
    )
)


###############################################################################
# 16. PARAMETER-LEVEL REPEATED-SPLIT SUMMARY
#
# Common source for:
#
#   Table 6
#   Table S9
#
###############################################################################

parameter_summary <-
  results_1000 %>%

  dplyr::group_by(
    method,
    parameter
  ) %>%

  dplyr::summarise(

    Full_Data_Benchmark =
      dplyr::first(
        benchmark
      ),


    MC_Mean =
      mean(
        estimate,
        na.rm = TRUE
      ),


    MC_SD =
      sd(
        estimate,
        na.rm = TRUE
      ),


    Mean_CI_Width =
      mean(
        ci_width,
        na.rm = TRUE
      ),


    SD_CI_Width =
      sd(
        ci_width,
        na.rm = TRUE
      ),


    Inclusion =
      mean(
        covered,
        na.rm = TRUE
      ),


    Mean_Analytic_SE =
      mean(
        analytic_se,
        na.rm = TRUE
      ),


    .groups =
      "drop"
  ) %>%

  dplyr::mutate(

    SE_to_MCSD =
      Mean_Analytic_SE /
      MC_SD,


    method =
      factor(
        method,
        levels = method_order
      ),


    parameter =
      factor(
        parameter,
        levels = parameter_names
      )
  ) %>%

  dplyr::arrange(
    method,
    parameter
  )


###############################################################################
# 17. SAVE PARAMETER-LEVEL DIAGNOSTIC SUMMARY
###############################################################################

write.csv(

  parameter_summary,

  file =
    file.path(
      diagnostics_dir,
      "nhanes_parameter_summary.csv"
    ),

  row.names = FALSE
)


###############################################################################
# 18. TABLE S9
# FULL REPEATED-SPLIT COEFFICIENT-LEVEL RESULTS
###############################################################################


###############################################################################
# 18A. FULL-DATA BENCHMARK ROW
###############################################################################

tableS9_benchmark_row <-
  parameter_summary %>%

  dplyr::select(
    parameter,
    Full_Data_Benchmark
  ) %>%

  dplyr::distinct() %>%

  dplyr::arrange(
    parameter
  ) %>%

  tidyr::pivot_wider(

    names_from =
      parameter,

    values_from =
      Full_Data_Benchmark
  ) %>%

  dplyr::mutate(

    Method =
      "",

    Statistic =
      "Full-data benchmark",

    .before = 1
  )


###############################################################################
# 18B. METHOD-SPECIFIC ROWS
###############################################################################

tableS9_method_rows <-
  parameter_summary %>%

  dplyr::select(

    method,
    parameter,

    MC_Mean,
    MC_SD,
    Mean_CI_Width,
    SD_CI_Width,
    Inclusion,
    Mean_Analytic_SE,
    SE_to_MCSD
  ) %>%

  tidyr::pivot_longer(

    cols =
      c(
        MC_Mean,
        MC_SD,
        Mean_CI_Width,
        SD_CI_Width,
        Inclusion,
        Mean_Analytic_SE,
        SE_to_MCSD
      ),

    names_to =
      "Statistic",

    values_to =
      "Value"
  ) %>%

  dplyr::mutate(

    Statistic =
      dplyr::recode(

        Statistic,

        MC_Mean =
          "MC mean",

        MC_SD =
          "MC SD",

        Mean_CI_Width =
          "Mean CI width",

        SD_CI_Width =
          "SD CI width",

        Inclusion =
          "Empirical coverage",

        Mean_Analytic_SE =
          "Mean analytic SE",

        SE_to_MCSD =
          "SE/MCSD"
      ),


    Method =
      factor(
        method,
        levels = method_order
      ),


    parameter =
      factor(
        parameter,
        levels = parameter_names
      ),


    Statistic =
      factor(

        Statistic,

        levels = c(
          "MC mean",
          "MC SD",
          "Mean CI width",
          "SD CI width",
          "Empirical coverage",
          "Mean analytic SE",
          "SE/MCSD"
        )
      )
  ) %>%

  dplyr::arrange(
    Method,
    Statistic,
    parameter
  ) %>%

  tidyr::pivot_wider(

    names_from =
      parameter,

    values_from =
      Value
  ) %>%

  dplyr::select(
    -method
  )


###############################################################################
# 18C. COMBINE TABLE S9
###############################################################################

tableS9_nhanes <-
  dplyr::bind_rows(

    tableS9_benchmark_row,
    tableS9_method_rows
  )


###############################################################################
# 18D. ROUND TABLE S9
#
# Empirical coverage: 3 decimals
# Everything else:     2 decimals
###############################################################################

tableS9_nhanes_rounded <-
  tableS9_nhanes %>%

  dplyr::mutate(

    dplyr::across(

      -c(
        Method,
        Statistic
      ),

      ~ dplyr::if_else(

        as.character(
          Statistic
        ) ==
          "Empirical coverage",

        round(
          .x,
          3
        ),

        round(
          .x,
          2
        )
      )
    )
  )


###############################################################################
# 18E. SAVE TABLE S9
###############################################################################

write.csv(

  tableS9_nhanes_rounded,

  file =
    file.path(
      tableS9_dir,
      "tableS9_nhanes_repeated_split.csv"
    ),

  row.names = FALSE
)


###############################################################################
# 19. TABLE 6
# COMPACT REPEATED-SPLIT SUMMARY
#
# Entries are medians across the 11 regression coefficients.
#
# CIW ratio:
#
#   Mean CI width for method /
#   Mean CI width for supervised estimator
#
# CIW dispersion:
#
#   SD CI width /
#   Mean CI width
#
# Inclusion:
#
#   Frequency that the split interval contains the
#   full-data supervised estimate.
###############################################################################


###############################################################################
# 19A. SUPERVISED CI-WIDTH REFERENCE
###############################################################################

supervised_ciw <-
  parameter_summary %>%

  dplyr::filter(
    method == "Supervised"
  ) %>%

  dplyr::select(

    parameter,

    Supervised_CIW =
      Mean_CI_Width
  )


###############################################################################
# 19B. CREATE TABLE 6
###############################################################################

table6_nhanes <-
  parameter_summary %>%

  dplyr::left_join(

    supervised_ciw,

    by = "parameter"
  ) %>%

  dplyr::mutate(

    CIW_ratio =
      Mean_CI_Width /
      Supervised_CIW,


    CIW_dispersion =
      SD_CI_Width /
      Mean_CI_Width
  ) %>%

  dplyr::group_by(
    method
  ) %>%

  dplyr::summarise(

    CIW_ratio =
      median(
        CIW_ratio,
        na.rm = TRUE
      ),


    CIW_dispersion =
      median(
        CIW_dispersion,
        na.rm = TRUE
      ),


    Inclusion =
      median(
        Inclusion,
        na.rm = TRUE
      ),


    .groups =
      "drop"
  ) %>%

  dplyr::mutate(

    Method =
      factor(
        method,
        levels = method_order
      )
  ) %>%

  dplyr::arrange(
    Method
  ) %>%

  dplyr::select(

    Method,

    CIW_ratio,

    CIW_dispersion,

    Inclusion
  )


###############################################################################
# 19C. ROUND TABLE 6
###############################################################################

table6_nhanes_rounded <-
  table6_nhanes %>%

  dplyr::mutate(

    CIW_ratio =
      round(
        CIW_ratio,
        2
      ),


    CIW_dispersion =
      round(
        CIW_dispersion,
        2
      ),


    Inclusion =
      round(
        Inclusion,
        3
      )
  )


###############################################################################
# 19D. SAVE TABLE 6
###############################################################################

write.csv(

  table6_nhanes_rounded,

  file =
    file.path(
      table6_dir,
      "table6_nhanes_compact_summary.csv"
    ),

  row.names = FALSE
)


###############################################################################
# 20. TABLE S8
# SINGLE LABEL-DELETION SPLIT
#
# Seed = 124
#
# Statistics:
#
#   Est
#   SE
#   CIW
#   VR
#
# VR =
#
#   supervised analytic variance /
#   method analytic variance
#
###############################################################################


###############################################################################
# 20A. RUN SINGLE SPLIT
###############################################################################

single_result <-
  run_one_split(

    split_seed =
      single_split_seed,

    datX =
      datX
  )


single_result <-
  single_result %>%

  dplyr::filter(
    method %in%
      method_order
  )


###############################################################################
# 20B. SUPERVISED VARIANCE REFERENCE
###############################################################################

sup_reference <-
  single_result %>%

  dplyr::filter(
    method ==
      "Supervised"
  ) %>%

  dplyr::transmute(

    parameter,

    Supervised_Variance =
      analytic_se^2
  )


###############################################################################
# 20C. SINGLE-SPLIT PARAMETER SUMMARY
###############################################################################

single_run_summary <-
  single_result %>%

  dplyr::left_join(

    sup_reference,

    by = "parameter"
  ) %>%

  dplyr::mutate(

    VR =
      Supervised_Variance /
      (analytic_se^2),


    method =
      factor(
        method,
        levels = method_order
      ),


    parameter =
      factor(
        parameter,
        levels = parameter_names
      )
  ) %>%

  dplyr::arrange(
    method,
    parameter
  ) %>%

  dplyr::select(

    method,

    parameter,

    benchmark,

    estimate,

    analytic_se,

    ci_width,

    VR
  )


###############################################################################
# 20D. FULL-DATA BENCHMARK ROW
###############################################################################

tableS8_benchmark_row <-
  single_run_summary %>%

  dplyr::select(
    parameter,
    benchmark
  ) %>%

  dplyr::distinct() %>%

  dplyr::arrange(
    parameter
  ) %>%

  tidyr::pivot_wider(

    names_from =
      parameter,

    values_from =
      benchmark
  ) %>%

  dplyr::mutate(

    Method =
      "",

    Statistic =
      "Full-data benchmark",

    .before = 1
  )


###############################################################################
# 20E. METHOD-SPECIFIC ROWS
###############################################################################

tableS8_method_rows <-
  single_run_summary %>%

  dplyr::select(

    method,

    parameter,

    estimate,

    analytic_se,

    ci_width,

    VR
  ) %>%

  tidyr::pivot_longer(

    cols =
      c(
        estimate,
        analytic_se,
        ci_width,
        VR
      ),

    names_to =
      "Statistic",

    values_to =
      "Value"
  ) %>%

  dplyr::mutate(

    Statistic =
      dplyr::recode(

        Statistic,

        estimate =
          "Est",

        analytic_se =
          "SE",

        ci_width =
          "CIW",

        VR =
          "VR"
      ),


    Statistic =
      factor(

        Statistic,

        levels = c(
          "Est",
          "SE",
          "CIW",
          "VR"
        )
      ),


    Method =
      factor(
        method,
        levels = method_order
      ),


    parameter =
      factor(
        parameter,
        levels = parameter_names
      )
  ) %>%

  dplyr::arrange(
    Method,
    Statistic,
    parameter
  ) %>%

  tidyr::pivot_wider(

    names_from =
      parameter,

    values_from =
      Value
  ) %>%

  dplyr::select(
    -method
  )


###############################################################################
# 20F. COMBINE TABLE S8
###############################################################################

tableS8_nhanes <-
  dplyr::bind_rows(

    tableS8_benchmark_row,

    tableS8_method_rows
  )


###############################################################################
# 20G. ROUND TABLE S8
###############################################################################

tableS8_nhanes_rounded <-
  tableS8_nhanes %>%

  dplyr::mutate(

    dplyr::across(

      -c(
        Method,
        Statistic
      ),

      ~ round(
        .x,
        2
      )
    )
  )


###############################################################################
# 20H. SAVE TABLE S8
###############################################################################

write.csv(

  tableS8_nhanes_rounded,

  file =
    file.path(
      tableS8_dir,
      "tableS8_nhanes_single_split.csv"
    ),

  row.names = FALSE
)


###############################################################################
# 21. SAVE SINGLE-SPLIT PARAMETER-LEVEL DIAGNOSTICS
###############################################################################

write.csv(

  single_run_summary,

  file =
    file.path(
      diagnostics_dir,
      "nhanes_single_split_parameter_results.csv"
    ),

  row.names = FALSE
)


###############################################################################
# 22. PRINT TABLE 6
###############################################################################

cat(
  "\n============================================================\n"
)

cat(
  "TABLE 6: COMPACT REPEATED-SPLIT SUMMARY\n"
)

cat(
  "============================================================\n\n"
)


print(
  table6_nhanes_rounded,
  n = Inf,
  width = Inf
)


###############################################################################
# 23. PRINT TABLE S8
###############################################################################

cat(
  "\n============================================================\n"
)

cat(
  "TABLE S8: SINGLE-SPLIT NHANES RESULTS\n"
)

cat(
  "Seed = ",
  single_split_seed,
  "\n",
  sep = ""
)

cat(
  "============================================================\n\n"
)


print(
  tableS8_nhanes_rounded,
  n = Inf,
  width = Inf
)


###############################################################################
# 24. PRINT TABLE S9
###############################################################################

cat(
  "\n============================================================\n"
)

cat(
  "TABLE S9: REPEATED-SPLIT COEFFICIENT-LEVEL RESULTS\n"
)

cat(
  "============================================================\n\n"
)


print(
  tableS9_nhanes_rounded,
  n = Inf,
  width = Inf
)


###############################################################################
# 25. FINAL OUTPUT LOCATIONS
###############################################################################

cat(
  "\n============================================================\n"
)

cat(
  "NHANES ANALYSIS COMPLETE\n"
)

cat(
  "============================================================\n"
)


cat(
  "\nTable 6:\n",
  file.path(
    table6_dir,
    "table6_nhanes_compact_summary.csv"
  ),
  "\n",
  sep = ""
)


cat(
  "\nTable S8:\n",
  file.path(
    tableS8_dir,
    "tableS8_nhanes_single_split.csv"
  ),
  "\n",
  sep = ""
)


cat(
  "\nTable S9:\n",
  file.path(
    tableS9_dir,
    "tableS9_nhanes_repeated_split.csv"
  ),
  "\n",
  sep = ""
)


cat(
  "\nRaw repeated-split results:\n",
  file.path(
    diagnostics_dir,
    "nhanes_repeated_split_raw.rds"
  ),
  "\n",
  sep = ""
)


cat(
  "\nParameter-level repeated-split diagnostics:\n",
  file.path(
    diagnostics_dir,
    "nhanes_parameter_summary.csv"
  ),
  "\n",
  sep = ""
)


cat(
  "\nSingle-split diagnostics:\n",
  file.path(
    diagnostics_dir,
    "nhanes_single_split_parameter_results.csv"
  ),
  "\n",
  sep = ""
)


cat(
  "\n============================================================\n"
)
