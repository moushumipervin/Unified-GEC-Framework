###############################################################################
# TABLE 5: MISSING-COVARIATE REGRESSION SIMULATION
#
# Reproduces Main Table 5:
# Simulation results for the missing-covariate regression setting
# based on 1,000 Monte Carlo replications.
#
# Run this script from the ROOT of:
#   Unified-GEC-Framework/
#
# Functions:
#   Missing Covariates/functions/missing_covariate_functions.R
#
# Outputs:
#   Missing Covariates/results/table5/
###############################################################################

rm(list = ls())


###############################################################################
# 0. REQUIRED PACKAGES
###############################################################################

required_packages <- c(
  "dplyr",
  "tidyr"
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
  library(tidyr)

})


###############################################################################
# 1. SOURCE MISSING-COVARIATE FUNCTIONS
###############################################################################

source(
  "Missing Covariates/functions/missing_covariate_functions.R"
)


###############################################################################
# 2. SETTINGS
###############################################################################

seed <- 1234

m <- 1000

k <- 3

n <- 1000


###############################################################################
# TRUE TARGET PARAMETERS
###############################################################################

beta_true_OR1 <- c(
  1,
  1,
  2
)

beta_true_OR2 <- c(
  0.500,
  0.795,
  1.000
)

parameter_names <- c(
  "beta0",
  "beta1",
  "beta2"
)


###############################################################################
# MANUSCRIPT TABLE CONVENTION
#
# Bias, MC SD, Mean analytic SE, RMSE, and Mean CI width
# are multiplied by 10.
###############################################################################

scale_factor <- 10


###############################################################################
# 3. SIMULATION SCENARIOS
###############################################################################

scenario_grid <- data.frame(

  OR = c(
    1,
    1,
    2,
    2
  ),

  PS = c(
    1,
    2,
    1,
    2
  ),

  Scenario = c(
    "OR1PS1",
    "OR1PS2",
    "OR2PS1",
    "OR2PS2"
  ),

  stringsAsFactors = FALSE
)


###############################################################################
# 4. HELPER FOR FULL-DATA AND COMPLETE-CASE METHODS
###############################################################################

make_lm_rows <- function(
    fit,
    method,
    truth,
    replication,
    scenario) {

  if (is.null(fit)) {
    return(NULL)
  }


  tab <-
    summary(
      fit
    )$coefficients


  estimate <-
    as.numeric(
      tab[, "Estimate"]
    )


  se <-
    as.numeric(
      tab[, "Std. Error"]
    )


  lower <-
    estimate -
    qnorm(0.975) *
    se


  upper <-
    estimate +
    qnorm(0.975) *
    se


  data.frame(

    replication =
      replication,

    Scenario =
      scenario,

    method =
      method,

    parameter =
      parameter_names,

    estimate =
      estimate,

    truth =
      truth,

    error =
      estimate -
      truth,

    squared_error =
      (
        estimate -
        truth
      )^2,

    analytic_se =
      se,

    ci_lower =
      lower,

    ci_upper =
      upper,

    ci_width =
      upper -
      lower,

    covered =
      as.numeric(
        lower <= truth &
          truth <= upper
      ),

    ESS =
      NA_real_,

    max_weight =
      NA_real_,

    estimator_success =
      1L,

    lambda_success =
      NA_integer_,

    calibration_residual =
      NA_real_,

    variance_success =
      1L,

    stringsAsFactors =
      FALSE
  )
}


###############################################################################
# 5. HELPER FOR IPW / AIPW / ET / HD / CE
###############################################################################

process_method_safe <- function(
    fit,
    method,
    truth,
    replication,
    scenario,
    entropy = NULL) {

  if (is.null(fit)) {
    return(NULL)
  }


  out <-
    tryCatch(

      process_missingcov_fit(

        fit =
          fit,

        method =
          method,

        truth =
          truth,

        entropy =
          entropy,

        parameter_names =
          parameter_names,

        replication =
          replication
      ),

      error = function(e) {

        message(
          method,
          " variance/output failed in replication ",
          replication,
          ": ",
          conditionMessage(e)
        )

        NULL
      }
    )


  if (is.null(out)) {
    return(NULL)
  }


  out$Scenario <-
    scenario


  # process_missingcov_fit() returns analytic SE/CI
  # when sandwich variance estimation succeeds.

  out$variance_success <-
    as.integer(
      is.finite(
        out$analytic_se
      )
    )


  out
}


###############################################################################
# 6. STORAGE
###############################################################################

all_results <- list()

row_counter <- 1L


###############################################################################
# 7. MAIN MONTE CARLO SIMULATION
###############################################################################

for (
  s in
  seq_len(
    nrow(
      scenario_grid
    )
  )
) {


  current_OR <-
    scenario_grid$OR[s]


  current_PS <-
    scenario_grid$PS[s]


  current_scenario <-
    scenario_grid$Scenario[s]


  beta_true <-
    if (
      current_OR == 1
    ) {

      beta_true_OR1

    } else {

      beta_true_OR2
    }


  cat(
    "\n============================================================\n"
  )

  cat(
    "Running scenario: ",
    current_scenario,
    "\n",
    sep = ""
  )

  cat(
    "Target: ",
    paste(
      beta_true,
      collapse = ", "
    ),
    "\n",
    sep = ""
  )

  cat(
    "============================================================\n"
  )


  ###########################################################################
  # MONTE CARLO REPLICATIONS
  ###########################################################################

  for (
    i in
    seq_len(
      m
    )
  ) {


    cat(
      "\r",
      current_scenario,
      ": replication ",
      i,
      " of ",
      m,
      sep = ""
    )


    set.seed(
      seed + i
    )


    ###########################################################################
    # A. GENERATE DATA
    ###########################################################################

    dat <-
      generate_data(

        n =
          n,

        OR =
          current_OR,

        PS =
          current_PS
      )


    ###########################################################################
    # B. FULL-DATA BENCHMARK
    ###########################################################################

    fit_full <-
      tryCatch(

        lm(
          y ~ x + z,
          data = dat
        ),

        error = function(e) {
          NULL
        }
      )


    out_full <-
      make_lm_rows(

        fit =
          fit_full,

        method =
          "Full",

        truth =
          beta_true,

        replication =
          i,

        scenario =
          current_scenario
      )


    if (
      !is.null(
        out_full
      )
    ) {

      all_results[[row_counter]] <-
        out_full

      row_counter <-
        row_counter + 1L
    }


    ###########################################################################
    # C. COMPLETE-CASE ESTIMATOR
    ###########################################################################

    fit_cc <-
      tryCatch(

        lm(

          y ~ x + z,

          data =
            dat[
              dat$D == 1,
              ,
              drop = FALSE
            ]
        ),

        error = function(e) {
          NULL
        }
      )


    out_cc <-
      make_lm_rows(

        fit =
          fit_cc,

        method =
          "CC",

        truth =
          beta_true,

        replication =
          i,

        scenario =
          current_scenario
      )


    if (
      !is.null(
        out_cc
      )
    ) {

      all_results[[row_counter]] <-
        out_cc

      row_counter <-
        row_counter + 1L
    }


    ###########################################################################
    # D. CROSS-FIT z.hat ONCE FOR IPW AND AIPW
    ###########################################################################

    data_all <-
      tryCatch(

        prepare_missingcov_data(

          data_full =
            dat,

          K =
            k,

          seed =
            seed + i
        ),

        error = function(e) {

          message(
            "\nCross-fitting failed in replication ",
            i,
            ": ",
            conditionMessage(e)
          )

          NULL
        }
      )


    if (
      !is.null(
        data_all
      )
    ) {


      #########################################################################
      # E. IPW
      #########################################################################

      fit_IPW <-
        tryCatch(

          estimate_ipw_missingcov(
            data_all =
              data_all
          ),

          error = function(e) {

            message(
              "\nIPW failed in replication ",
              i,
              ": ",
              conditionMessage(e)
            )

            NULL
          }
        )


      out_IPW <-
        process_method_safe(

          fit =
            fit_IPW,

          method =
            "IPW",

          truth =
            beta_true,

          replication =
            i,

          scenario =
            current_scenario
        )


      if (
        !is.null(
          out_IPW
        )
      ) {

        all_results[[row_counter]] <-
          out_IPW

        row_counter <-
          row_counter + 1L
      }


      #########################################################################
      # F. AIPW
      #
      # Point estimator uses the same estimating equation as
      # aipw_sandwich_missingcov().
      #########################################################################

      fit_AIPW <-
        tryCatch(

          estimate_aipw_missingcov(
            data_all =
              data_all
          ),

          error = function(e) {

            message(
              "\nAIPW failed in replication ",
              i,
              ": ",
              conditionMessage(e)
            )

            NULL
          }
        )


      out_AIPW <-
        process_method_safe(

          fit =
            fit_AIPW,

          method =
            "AIPW",

          truth =
            beta_true,

          replication =
            i,

          scenario =
            current_scenario
        )


      if (
        !is.null(
          out_AIPW
        )
      ) {

        all_results[[row_counter]] <-
          out_AIPW

        row_counter <-
          row_counter + 1L
      }


      #########################################################################
      # INITIAL VALUE FOR ET / HD / CE
      #########################################################################

      theta_start <-
        if (
          !is.null(
            fit_AIPW
          ) &&
          all(
            is.finite(
              fit_AIPW$theta
            )
          )
        ) {

          as.numeric(
            fit_AIPW$theta
          )

        } else if (
          !is.null(
            fit_IPW
          ) &&
          all(
            is.finite(
              fit_IPW$theta
            )
          )
        ) {

          as.numeric(
            fit_IPW$theta
          )

        } else {

          c(
            0,
            0,
            0
          )
        }


      #########################################################################
      # G. ET
      #########################################################################

      fit_ET <-
        tryCatch(

          estimate_theta_EM_kfold_dual_ET(

            th =
              theta_start,

            data_full =
              dat,

            K =
              k,

            seed =
              seed + i,

            max.iter =
              50,

            eps =
              1e-6,

            lambda_maxit =
              1000,

            lambda_tol =
              1e-10,

            damping =
              1
          ),

          error = function(e) {

            message(
              "\nET failed in replication ",
              i,
              ": ",
              conditionMessage(e)
            )

            NULL
          }
        )


      out_ET <-
        process_method_safe(

          fit =
            fit_ET,

          method =
            "ET",

          truth =
            beta_true,

          replication =
            i,

          scenario =
            current_scenario,

          entropy =
            "ET"
        )


      if (
        !is.null(
          out_ET
        )
      ) {

        all_results[[row_counter]] <-
          out_ET

        row_counter <-
          row_counter + 1L
      }


      #########################################################################
      # H. HD
      #########################################################################

      fit_HD <-
        tryCatch(

          estimate_theta_EM_kfold_dual_HD(

            th =
              theta_start,

            data_full =
              dat,

            K =
              k,

            seed =
              seed + i,

            max.iter =
              50,

            eps =
              1e-6,

            lambda_maxit =
              1000,

            lambda_tol =
              1e-10,

            damping =
              1
          ),

          error = function(e) {

            message(
              "\nHD failed in replication ",
              i,
              ": ",
              conditionMessage(e)
            )

            NULL
          }
        )


      out_HD <-
        process_method_safe(

          fit =
            fit_HD,

          method =
            "HD",

          truth =
            beta_true,

          replication =
            i,

          scenario =
            current_scenario,

          entropy =
            "HD"
        )


      if (
        !is.null(
          out_HD
        )
      ) {

        all_results[[row_counter]] <-
          out_HD

        row_counter <-
          row_counter + 1L
      }


      #########################################################################
      # I. CE
      #########################################################################

      fit_CE <-
        tryCatch(

          estimate_theta_EM_kfold_dual_CE(

            th =
              theta_start,

            data_full =
              dat,

            K =
              k,

            seed =
              seed + i,

            max.iter =
              50,

            eps =
              1e-6,

            lambda_maxit =
              1000,

            lambda_tol =
              1e-8,

            damping =
              0.1
          ),

          error = function(e) {

            message(
              "\nCE failed in replication ",
              i,
              ": ",
              conditionMessage(e)
            )

            NULL
          }
        )


      out_CE <-
        process_method_safe(

          fit =
            fit_CE,

          method =
            "CE",

          truth =
            beta_true,

          replication =
            i,

          scenario =
            current_scenario,

          entropy =
            "CE"
        )


      if (
        !is.null(
          out_CE
        )
      ) {

        all_results[[row_counter]] <-
          out_CE

        row_counter <-
          row_counter + 1L
      }
    }
  }


  cat(
    "\n"
  )
}


###############################################################################
# 8. COMBINE RAW REPLICATION-LEVEL RESULTS
###############################################################################

raw_results <-
  dplyr::bind_rows(
    all_results
  )


###############################################################################
# 9. MONTE CARLO SUMMARY
###############################################################################

final_table <-
  raw_results %>%

  dplyr::group_by(
    Scenario,
    method,
    parameter
  ) %>%

  dplyr::summarise(


    ###########################################################################
    # NUMBER OF SUCCESSFUL POINT ESTIMATES
    ###########################################################################

    N_Estimate =
      sum(
        is.finite(
          estimate
        )
      ),


    ###########################################################################
    # NUMBER WITH FINITE ANALYTIC SE
    ###########################################################################

    N_Variance =
      sum(
        is.finite(
          analytic_se
        )
      ),


    ###########################################################################
    # MEAN ESTIMATE
    ###########################################################################

    Mean_Estimate =
      mean(
        estimate,
        na.rm = TRUE
      ),


    ###########################################################################
    # BIAS x 10
    ###########################################################################

    Bias =
      mean(
        estimate -
          truth,
        na.rm = TRUE
      ) *
      scale_factor,


    ###########################################################################
    # MONTE CARLO SD x 10
    ###########################################################################

    MC_SD =
      sd(
        estimate,
        na.rm = TRUE
      ) *
      scale_factor,


    ###########################################################################
    # MEAN ANALYTIC SE x 10
    ###########################################################################

    Mean_Analytic_SE =
      mean(
        analytic_se,
        na.rm = TRUE
      ) *
      scale_factor,


    ###########################################################################
    # ANALYTIC SE / MC SD
    ###########################################################################

    SE_to_MCSD =
      (
        mean(
          analytic_se,
          na.rm = TRUE
        ) /
        sd(
          estimate,
          na.rm = TRUE
        )
      ),


    ###########################################################################
    # RMSE x 10
    ###########################################################################

    RMSE =
      sqrt(
        mean(
          (
            estimate -
              truth
          )^2,
          na.rm = TRUE
        )
      ) *
      scale_factor,


    ###########################################################################
    # EMPIRICAL COVERAGE
    ###########################################################################

    Coverage =
      mean(
        covered,
        na.rm = TRUE
      ),


    ###########################################################################
    # MEAN CI WIDTH x 10
    ###########################################################################

    Mean_CI_Width =
      mean(
        ci_width,
        na.rm = TRUE
      ) *
      scale_factor,


    ###########################################################################
    # MEAN ESS
    ###########################################################################

    Mean_ESS =
      if (
        all(
          is.na(
            ESS
          )
        )
      ) {

        NA_real_

      } else {

        mean(
          ESS,
          na.rm = TRUE
        )
      },


    ###########################################################################
    # MEAN MAXIMUM WEIGHT
    ###########################################################################

    Mean_MaxWeight =
      if (
        all(
          is.na(
            max_weight
          )
        )
      ) {

        NA_real_

      } else {

        mean(
          max_weight,
          na.rm = TRUE
        )
      },


    ###########################################################################
    # MEDIAN MAXIMUM WEIGHT
    ###########################################################################

    Median_MaxWeight =
      if (
        all(
          is.na(
            max_weight
          )
        )
      ) {

        NA_real_

      } else {

        median(
          max_weight,
          na.rm = TRUE
        )
      },


    ###########################################################################
    # 95TH PERCENTILE MAXIMUM WEIGHT
    ###########################################################################

    P95_MaxWeight =
      if (
        all(
          is.na(
            max_weight
          )
        )
      ) {

        NA_real_

      } else {

        as.numeric(
          quantile(
            max_weight,
            probs = 0.95,
            na.rm = TRUE
          )
        )
      },


    ###########################################################################
    # ESTIMATOR FAILURE RATE
    ###########################################################################

    Estimator_Failure_Rate =
      1 -
      mean(
        estimator_success,
        na.rm = TRUE
      ),


    ###########################################################################
    # VARIANCE FAILURE RATE
    ###########################################################################

    Variance_Failure_Rate =
      1 -
      mean(
        variance_success,
        na.rm = TRUE
      ),


    ###########################################################################
    # LAMBDA FAILURE RATE
    ###########################################################################

    Lambda_Failure_Rate =
      if (
        all(
          is.na(
            lambda_success
          )
        )
      ) {

        NA_real_

      } else {

        1 -
        mean(
          lambda_success,
          na.rm = TRUE
        )
      },


    ###########################################################################
    # MEAN CALIBRATION RESIDUAL
    ###########################################################################

    Mean_Calibration_Residual =
      if (
        all(
          is.na(
            calibration_residual
          )
        )
      ) {

        NA_real_

      } else {

        mean(
          calibration_residual,
          na.rm = TRUE
        )
      },


    .groups =
      "drop"
  )


###############################################################################
# 10. CLEAN METHOD/PARAMETER LABELS
###############################################################################

final_table <-
  final_table %>%

  dplyr::rename(

    Method =
      method,

    Parameter =
      parameter
  )


###############################################################################
# 11. ARRANGE COMPLETE SUMMARY FOR INSPECTION
###############################################################################

final_table_print <-
  final_table %>%

  dplyr::mutate(

    Scenario =
      factor(
        Scenario,
        levels = c(
          "OR1PS1",
          "OR1PS2",
          "OR2PS1",
          "OR2PS2"
        )
      ),

    Method =
      factor(
        Method,
        levels = c(
          "Full",
          "CC",
          "IPW",
          "AIPW",
          "ET",
          "HD",
          "CE"
        )
      ),

    Parameter =
      factor(
        Parameter,
        levels = c(
          "beta0",
          "beta1",
          "beta2"
        )
      )
  ) %>%

  dplyr::arrange(
    Scenario,
    Method,
    Parameter
  ) %>%

  dplyr::mutate(

    dplyr::across(
      where(is.numeric),
      ~ round(
        .x,
        4
      )
    )
  )


print(
  final_table_print,
  n = Inf
)


###############################################################################
# 12. CREATE TABLE 5 OUTPUT DIRECTORY
###############################################################################

dir.create(
  "Missing Covariates/results/table5",
  recursive = TRUE,
  showWarnings = FALSE
)


###############################################################################
# 13. SAVE RAW MONTE CARLO RESULTS
###############################################################################

saveRDS(
  raw_results,
  "Missing Covariates/results/table5/table5_missing_covariates_raw.rds"
)


###############################################################################
# 14. SAVE COMPLETE MONTE CARLO DIAGNOSTIC SUMMARY
###############################################################################

write.csv(
  final_table,
  "Missing Covariates/results/table5/table5_missing_covariates_diagnostics.csv",
  row.names = FALSE
)


saveRDS(
  final_table,
  "Missing Covariates/results/table5/table5_missing_covariates_diagnostics.rds"
)


###############################################################################
# 15. CREATE THE EXACT MAIN TABLE 5 DATA
#
# Paper columns:
#
# Model | Method |
# beta0: Bias, MC SD, RMSE, Cov.
# beta1: Bias, MC SD, RMSE, Cov.
# beta2: Bias, MC SD, RMSE, Cov.
###############################################################################

table5_missing_covariates <-
  final_table %>%

  dplyr::select(
    Scenario,
    Method,
    Parameter,
    Bias,
    MC_SD,
    RMSE,
    Coverage
  ) %>%


  ###########################################################################
  # Scenario, method, and parameter ordering
  ###########################################################################

  dplyr::mutate(

    Scenario =
      factor(
        Scenario,
        levels = c(
          "OR1PS1",
          "OR1PS2",
          "OR2PS1",
          "OR2PS2"
        )
      ),

    Method =
      factor(
        Method,
        levels = c(
          "Full",
          "CC",
          "IPW",
          "AIPW",
          "ET",
          "HD",
          "CE"
        )
      ),

    Parameter =
      factor(
        Parameter,
        levels = c(
          "beta0",
          "beta1",
          "beta2"
        )
      )
  ) %>%


  dplyr::arrange(
    Scenario,
    Method,
    Parameter
  ) %>%


  ###########################################################################
  # Put beta0, beta1, beta2 next to one another
  ###########################################################################

  tidyr::pivot_wider(

    names_from =
      Parameter,

    values_from =
      c(
        Bias,
        MC_SD,
        RMSE,
        Coverage
      ),

    names_glue =
      "{Parameter}_{.value}"
  ) %>%


  ###########################################################################
  # Exact manuscript column order
  ###########################################################################

  dplyr::select(

    Scenario,
    Method,

    beta0_Bias,
    beta0_MC_SD,
    beta0_RMSE,
    beta0_Coverage,

    beta1_Bias,
    beta1_MC_SD,
    beta1_RMSE,
    beta1_Coverage,

    beta2_Bias,
    beta2_MC_SD,
    beta2_RMSE,
    beta2_Coverage
  ) %>%


  dplyr::arrange(
    Scenario,
    Method
  )


###############################################################################
# 16. PRINT TABLE 5
###############################################################################

print(
  table5_missing_covariates,
  n = Inf
)


###############################################################################
# 17. SAVE TABLE 5
###############################################################################

write.csv(
  table5_missing_covariates,
  "Missing Covariates/results/table5/table5_missing_covariates.csv",
  row.names = FALSE
)


###############################################################################
# 18. OPTIONAL: COMPLETE MANUSCRIPT DIAGNOSTIC TABLE
#
# This is not Main Table 5.
# It keeps additional quantities useful for reproducibility/reviewer checks.
###############################################################################

table5_detailed_diagnostics <-
  final_table %>%

  dplyr::select(

    Scenario,
    Method,
    Parameter,

    Bias,
    MC_SD,

    Mean_Analytic_SE,
    SE_to_MCSD,

    RMSE,
    Coverage,

    Mean_CI_Width,

    Mean_ESS,
    Mean_MaxWeight,
    Median_MaxWeight,
    P95_MaxWeight,

    Estimator_Failure_Rate,
    Variance_Failure_Rate,
    Lambda_Failure_Rate,
    Mean_Calibration_Residual
  ) %>%

  dplyr::mutate(

    Scenario =
      factor(
        Scenario,
        levels = c(
          "OR1PS1",
          "OR1PS2",
          "OR2PS1",
          "OR2PS2"
        )
      ),

    Method =
      factor(
        Method,
        levels = c(
          "Full",
          "CC",
          "IPW",
          "AIPW",
          "ET",
          "HD",
          "CE"
        )
      ),

    Parameter =
      factor(
        Parameter,
        levels = c(
          "beta0",
          "beta1",
          "beta2"
        )
      )
  ) %>%

  dplyr::arrange(
    Scenario,
    Method,
    Parameter
  )


write.csv(
  table5_detailed_diagnostics,
  "Missing Covariates/results/table5_missing_covariates_detailed_diagnostics.csv",
  row.names = FALSE
)


###############################################################################
# 19. FINAL MESSAGE
###############################################################################

cat(
  "\n============================================================\n"
)

cat(
  "TABLE 5 SIMULATION COMPLETE\n"
)

cat(
  "============================================================\n"
)

cat(
  "\nMain manuscript table:\n"
)

cat(
  "Missing Covariates/results/table5_missing_covariates.csv\n"
)

cat(
  "\nRaw Monte Carlo results:\n"
)

cat(
  "Missing Covariates/results/table5_missing_covariates_raw.rds\n"
)

cat(
  "\nComplete diagnostics:\n"
)

cat(
  "Missing Covariates/results/table5_missing_covariates_diagnostics.csv\n"
)

cat(
  "\nDetailed diagnostics:\n"
)

cat(
  "Missing Covariates/results/table5_missing_covariates_detailed_diagnostics.csv\n"
)

cat(
  "\n============================================================\n"
)
