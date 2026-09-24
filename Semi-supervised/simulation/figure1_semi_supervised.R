################################################################################
# FIGURE 1
# Semi-supervised regression simulation under OR2 + MCAR
#
# This script:
#   1. Loads required packages
#   2. Sources the semi-supervised functions
#   3. Sets the OR2 + MCAR simulation
#   4. Computes the target regression parameter
#   5. Runs 1000 Monte Carlo replications
#   6. Fits:
#        Supervised
#        PI
#        EASE
#        DRESS
#        PSSE
#        ET
#        HD
#        CE
#   7. Saves the raw Monte Carlo results
#   8. Recreates Figure 1 as boxplots
#
# Run this script from the root directory of:
# Unified-GEC-Framework/
################################################################################


rm(list = ls())


################################################################################
# 1. REQUIRED PACKAGES
################################################################################

required_packages <- c(
  "MASS",
  "mgcv",
  "CVXR",
  "caret",
  "dplyr",
  "tidyr",
  "purrr",
  "ggplot2",
  "ggh4x",
  "sandwich"
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

  library(MASS)
  library(mgcv)
  library(CVXR)
  library(caret)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(ggh4x)
  library(sandwich)

})


################################################################################
# 2. SOURCE SEMI-SUPERVISED FUNCTIONS
################################################################################

source(
  "Semi-supervised/functions/semi_supervised_functions.R"
)


################################################################################
# 3. GLOBAL SIMULATION PARAMETERS
################################################################################

n <- 1000
# labeled sample size

N <- 1500
# unlabeled sample size

p <- 4
# number of predictors

rep <- 1000
# number of Monte Carlo replications

polyOrder <- 4
# polynomial order used by PSSE and DRESS

OR <- 2
# OR = 1: linear outcome model
# OR = 2: nonlinear outcome model

MAR <- 0
# MAR = 1: MAR
# MAR = 0: MCAR

q <- 5
# number of folds used by DRESS


################################################################################
# 4. TRUE TARGET PARAMETER
################################################################################

set.seed(1234)


n_target <- 10^7


alpha0 <- 1

alpha1 <-
  rep(
    1,
    p
  )

alpha2 <-
  rep(
    1,
    p
  )


X_target <-
  MASS::mvrnorm(

    n =
      n_target,

    mu =
      rep(
        0,
        p
      ),

    Sigma =
      diag(
        p
      )
  )


if (OR == 1) {

  Y_target <-
    alpha0 +
    X_target %*%
    alpha1 +
    rnorm(
      n_target,
      0,
      1
    )

} else {

  Y_target <-
    alpha0 +
    X_target %*%
    alpha1 +
    (
      X_target^3 -
      X_target^2 +
      exp(
        X_target
      )
    ) %*%
    alpha2 +
    rnorm(
      n_target,
      0,
      2
    )
}


target_parameter <-
  as.numeric(
    coef(
      lm(
        Y_target ~ X_target
      )
    )
  )


parameter_names <-
  c(
    "(Intercept)",
    paste0(
      "X",
      seq_len(
        p
      )
    )
  )


cat(
  "\nTarget parameter:\n"
)

print(
  target_parameter
)


################################################################################
# 5. STORAGE
################################################################################

simulation_results <-
  vector(
    "list",
    rep
  )


################################################################################
# 6. MONTE CARLO SIMULATION
################################################################################

for (
  k in
  seq_len(
    rep
  )
) {


  message(
    "Replication ",
    k,
    "/",
    rep
  )


  set.seed(
    k + 20220122
  )


  one_rep <-
    list()


  ##############################################################################
  # 6.1 GENERATE DATA
  ##############################################################################

  DesiredData <-
    GenerateData(

      n =
        n,

      N =
        N,

      p =
        p,

      OR =
        OR,

      MAR =
        MAR
    )


  data_labelled <-
    DesiredData$Data.labelled


  data_unlabelled <-
    DesiredData$Data.unlabelled


  data_full <-
    DesiredData$data_full


  ##############################################################################
  # 6.2 INITIAL OLS
  ##############################################################################

  formula_lm <-
    as.formula(
      paste0(

        "Y ~ ",

        paste0(
          colnames(
            data_labelled[
              ,
              -1
            ]
          ),
          collapse = " + "
        )
      )
    )


  fit_initial <-
    lm(

      formula_lm,

      data =
        as.data.frame(
          data_labelled
        )
    )


  theta_init <-
    as.numeric(
      coef(
        fit_initial
      )
    )


  ##############################################################################
  # 6.3 SUPERVISED
  ##############################################################################

  fit_sup <-
    tryCatch(

      SupervisedEst(
        data_labelled
      ),

      error =
        function(e) {

          message(
            "Supervised failed in replication ",
            k,
            ": ",
            conditionMessage(e)
          )

          NULL
        }
    )


  if (
    !is.null(
      fit_sup
    )
  ) {


    sup_est <-
      as.numeric(
        fit_sup$Est.coef
      )


    sup_se <-
      sqrt(
        diag(
          sandwich::vcovHC(
            fit_initial,
            type = "HC0"
          )
        )
      )


    one_rep[["Supervised"]] <-
      make_simulation_rows(

        rep_id =
          k,

        method =
          "Supervised",

        estimate =
          sup_est,

        truth =
          target_parameter,

        se =
          sup_se,

        parameter_names =
          parameter_names
      )
  }


  ##############################################################################
  # 6.4 PSSE
  ##############################################################################

  fit_PSSE <-
    tryCatch(

      PSSE1(

        data_labelled,

        data_unlabelled,

        type =
          "linear",

        sd =
          TRUE,

        alpha =
          polyOrder
      ),

      error =
        function(e) {

          message(
            "PSSE failed in replication ",
            k,
            ": ",
            conditionMessage(e)
          )

          NULL
        }
    )


  if (
    !is.null(
      fit_PSSE
    )
  ) {


    PSSE_est <-
      as.numeric(
        fit_PSSE$Hattheta
      )


    PSSE_se <-
      if (
        !is.null(
          fit_PSSE$sd.of.hattheta
        )
      ) {

        as.numeric(
          fit_PSSE$sd.of.hattheta
        )

      } else {

        rep(
          NA_real_,
          length(
            target_parameter
          )
        )
      }


    one_rep[["PSSE"]] <-
      make_simulation_rows(

        rep_id =
          k,

        method =
          "PSSE",

        estimate =
          PSSE_est,

        truth =
          target_parameter,

        se =
          PSSE_se,

        parameter_names =
          parameter_names
      )
  }


  ##############################################################################
  # 6.5 PI
  ##############################################################################

  fit_PI <-
    tryCatch(

      PI(

        data_labelled,

        data_unlabelled
      ),

      error =
        function(e) {

          message(
            "PI failed in replication ",
            k,
            ": ",
            conditionMessage(e)
          )

          NULL
        }
    )


  if (
    !is.null(
      fit_PI
    )
  ) {


    PI_est <-
      as.numeric(
        fit_PI$Hattheta
      )


    PI_se <-
      if (
        !is.null(
          fit_PI$sd.of.hattheta
        )
      ) {

        as.numeric(
          fit_PI$sd.of.hattheta
        )

      } else {

        rep(
          NA_real_,
          length(
            target_parameter
          )
        )
      }


    one_rep[["PI"]] <-
      make_simulation_rows(

        rep_id =
          k,

        method =
          "PI",

        estimate =
          PI_est,

        truth =
          target_parameter,

        se =
          PI_se,

        parameter_names =
          parameter_names
      )
  }


  ##############################################################################
  # 6.6 EASE
  ##############################################################################

  fit_EASE <-
    tryCatch(

      EASE(

        data_labelled,

        data_unlabelled,

        K =
          2,

        H =
          2,

        r =
          2
      ),

      error =
        function(e) {

          message(
            "EASE failed in replication ",
            k,
            ": ",
            conditionMessage(e)
          )

          NULL
        }
    )


  if (
    !is.null(
      fit_EASE
    )
  ) {


    EASE_est <-
      as.numeric(
        fit_EASE$Hattheta
      )


    EASE_se <-
      if (
        !is.null(
          fit_EASE$sd.of.hattheta
        )
      ) {

        as.numeric(
          fit_EASE$sd.of.hattheta
        )

      } else {

        rep(
          NA_real_,
          length(
            target_parameter
          )
        )
      }


    one_rep[["EASE"]] <-
      make_simulation_rows(

        rep_id =
          k,

        method =
          "EASE",

        estimate =
          EASE_est,

        truth =
          target_parameter,

        se =
          EASE_se,

        parameter_names =
          parameter_names
      )
  }


  ##############################################################################
  # 6.7 DRESS
  ##############################################################################

  fit_DRESS <-
    tryCatch(

      DRESS1(

        data_labelled,

        data_unlabelled,

        L =
          polyOrder,

        Kfolds =
          q
      ),

      error =
        function(e) {

          message(
            "DRESS failed in replication ",
            k,
            ": ",
            conditionMessage(e)
          )

          NULL
        }
    )


  if (
    !is.null(
      fit_DRESS
    )
  ) {


    DRESS_est <-
      as.numeric(
        fit_DRESS$Hattheta
      )


    DRESS_se <-
      if (
        !is.null(
          fit_DRESS$sd.of.hattheta
        )
      ) {

        as.numeric(
          fit_DRESS$sd.of.hattheta
        )

      } else {

        rep(
          NA_real_,
          length(
            target_parameter
          )
        )
      }


    one_rep[["DRESS"]] <-
      make_simulation_rows(

        rep_id =
          k,

        method =
          "DRESS",

        estimate =
          DRESS_est,

        truth =
          target_parameter,

        se =
          DRESS_se,

        parameter_names =
          parameter_names
      )
  }


  ##############################################################################
  # 6.8 ET
  ##############################################################################

  ET_fit <-
    tryCatch(

      estimate_theta_EM_kfold_dual_ET(

        th =
          theta_init,

        data_full =
          data_full,

        K =
          5,

        seed =
          k + 20220122,

        max.iter =
          500,

        eps =
          1e-4,

        damping =
          0.1
      ),

      error =
        function(e) {

          message(
            "ET estimator failed in replication ",
            k,
            ": ",
            conditionMessage(e)
          )

          NULL
        }
    )


  one_rep[["ET"]] <-
    process_gec_simulation(

      fit_gec =
        ET_fit,

      entropy =
        "ET",

      rep_id =
        k,

      target_parameter =
        target_parameter,

      parameter_names =
        parameter_names
    )


  ##############################################################################
  # 6.9 HD
  ##############################################################################

  HD_fit <-
    tryCatch(

      estimate_theta_EM_kfold_dual_HD(

        th =
          theta_init,

        data_full =
          data_full,

        K =
          5,

        seed =
          k + 20220122,

        max.iter =
          500,

        eps =
          1e-4,

        damping =
          0.1
      ),

      error =
        function(e) {

          message(
            "HD estimator failed in replication ",
            k,
            ": ",
            conditionMessage(e)
          )

          NULL
        }
    )


  one_rep[["HD"]] <-
    process_gec_simulation(

      fit_gec =
        HD_fit,

      entropy =
        "HD",

      rep_id =
        k,

      target_parameter =
        target_parameter,

      parameter_names =
        parameter_names
    )


  ##############################################################################
  # 6.10 CE
  ##############################################################################

  CE_fit <-
    tryCatch(

      estimate_theta_EM_kfold_CE(

        th =
          theta_init,

        data_full =
          data_full,

        K =
          5,

        seed =
          k + 20220122,

        max.iter =
          500,

        eps =
          1e-4,

        damping =
          0.1
      ),

      error =
        function(e) {

          message(
            "CE estimator failed in replication ",
            k,
            ": ",
            conditionMessage(e)
          )

          NULL
        }
    )


  one_rep[["CE"]] <-
    process_gec_simulation(

      fit_gec =
        CE_fit,

      entropy =
        "CE",

      rep_id =
        k,

      target_parameter =
        target_parameter,

      parameter_names =
        parameter_names
    )


  ##############################################################################
  # 6.11 SAVE THIS REPLICATION
  ##############################################################################

  simulation_results[[k]] <-
    dplyr::bind_rows(
      one_rep
    )

}


################################################################################
# 7. COMBINE MONTE CARLO RESULTS
################################################################################

simulation_results_df <-
  dplyr::bind_rows(
    simulation_results
  )


################################################################################
# 8. ORDER METHODS AND CREATE PARAMETER LABELS FOR FIGURE 1
################################################################################

simulation_results_df <-
  simulation_results_df %>%

  dplyr::mutate(

    beta_idx =
      match(
        parameter,
        parameter_names
      ) - 1,

    beta_lab =
      factor(

        paste0(
          "beta[",
          beta_idx,
          "]"
        ),

        levels =
          paste0(
            "beta[",
            0:(
              length(
                parameter_names
              ) - 1
            ),
            "]"
          )
      ),

    method =
      factor(

        method,

        levels =
          c(
            "Supervised",
            "PI",
            "EASE",
            "DRESS",
            "PSSE",
            "ET",
            "HD",
            "CE"
          )
      )
  )


################################################################################
# 9. TARGET VALUES FOR REFERENCE LINES
################################################################################

target_df <-
  data.frame(

    beta_lab =
      factor(

        paste0(
          "beta[",
          0:(
            length(
              target_parameter
            ) - 1
          ),
          "]"
        ),

        levels =
          paste0(
            "beta[",
            0:(
              length(
                target_parameter
              ) - 1
            ),
            "]"
          )
      ),

    target =
      target_parameter
  )


################################################################################
# 10. CREATE FIGURE 1
################################################################################

figure1 <-
  ggplot(

    simulation_results_df,

    aes(
      x =
        method,
      y =
        estimate
    )
  ) +

  geom_boxplot(

    outlier.size =
      0.5,

    fill =
      "grey80",

    na.rm =
      TRUE
  ) +

  geom_hline(

    data =
      target_df,

    aes(
      yintercept =
        target
    ),

    color =
      "red",

    linewidth =
      0.6
  ) +

  facet_wrap(

    ~ beta_lab,

    nrow =
      1,

    labeller =
      label_parsed,

    scales =
      "free_y"
  ) +

  labs(

    x =
      "Method",

    y =
      "Point estimate"
  ) +

  theme_bw() +

  theme(

    text =
      element_text(
        size = 12
      ),

    panel.grid =
      element_blank(),

    axis.text.x =
      element_text(

        angle =
          45,

        hjust =
          1
      ),

    strip.text =
      element_text(
        size = 11
      )
  )


################################################################################
# 11. DISPLAY FIGURE
################################################################################

print(
  figure1
)


################################################################################
# 12. CREATE RESULTS DIRECTORY
################################################################################

dir.create(

  "Semi Supervised/results/figure1",

  recursive =
    TRUE,

  showWarnings =
    FALSE
)


################################################################################
# 13. SAVE FIGURE 1
################################################################################

ggsave(

  filename =
    "Semi Supervised/results/figure1/figure1_boxplots_mcar_or2.pdf",

  plot =
    figure1,

  width =
    15,

  height =
    4
)


################################################################################
# 14. SAVE RAW MONTE CARLO RESULTS
################################################################################

saveRDS(

  simulation_results_df,

  file =
    "Semi Supervised/results/figure1/figure1_mcar_or2_raw.rds"
)


write.csv(

  simulation_results_df,

  file =
    "Semi Supervised/results/figure1/figure1_mcar_or2_raw.csv",

  row.names =
    FALSE
)


################################################################################
# 15. SAVE TARGET PARAMETER
################################################################################

target_parameter_df <-
  data.frame(

    Parameter =
      parameter_names,

    Target =
      target_parameter
  )


write.csv(

  target_parameter_df,

  file =
    "Semi Supervised/results/figure1/figure1_mcar_or2_target.csv",

  row.names =
    FALSE
)


################################################################################
# 16. FINAL MESSAGE
################################################################################

cat(
  "\n============================================================\n"
)

cat(
  "FIGURE 1 SIMULATION COMPLETE: OR2 + MCAR\n"
)

cat(
  "============================================================\n"
)

cat(
  "\nFigure:\n"
)

cat(
  "Semi Supervised/results/figure1/figure1_boxplots_mcar_or2.pdf\n"
)

cat(
  "\nRaw Monte Carlo results:\n"
)

cat(
  "Semi Supervised/results/figure1/figure1_mcar_or2_raw.rds\n"
)

cat(
  "Semi Supervised/results/figure1/figure1_mcar_or2_raw.csv\n"
)

cat(
  "\nTarget parameter:\n"
)

cat(
  "Semi Supervised/results/figure1/figure1_mcar_or2_target.csv\n"
)

cat(
  "\n============================================================\n"
)
