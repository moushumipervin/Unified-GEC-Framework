###############################################################################
# TABLE 10: NSW / PSID / CPS REAL-DATA APPLICATION
###############################################################################
#
# This script reproduces Table 10 of the manuscript.
#
# Data:
#   NSW experimental sample
#   PSID observational controls
#   CPS observational controls
#
# Methods:
#   Unweighted
#   IPW
#   EBPS
#   CBPS
#   EBCW
#   AIPW-LM
#   AIPW-GAM
#   ET
#   HD
#   EL
#
###############################################################################


rm(list = ls())


###############################################################################
# 1. REQUIRED PACKAGES
###############################################################################

required_packages <- c(
  "dplyr",
  "nnet",
  "mgcv",
  "DAAG",
  "haven",
  "numDeriv",
  "MASS",
  "CVXR",
  "CBPS",
  "WeightIt"
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
  library(nnet)
  library(mgcv)
  library(DAAG)
  library(haven)
  library(numDeriv)
})


###############################################################################
# 2. SOURCE LALONDE FUNCTIONS
###############################################################################

source(
  "lalonde_functions.R"
)


###############################################################################
# 3. SETTINGS
###############################################################################

k <- 4
K <- 4

seed <- 2024202

xvars <- c(
  "age",
  "education",
  "black",
  "hispanic",
  "married",
  "nodegree",
  "re74",
  "re75"
)


###############################################################################
# 4. LOAD DATA
###############################################################################

# NSW experimental data
exp_dat <- haven::read_dta(
  "http://www.nber.org/~rdehejia/data/nsw_dw.dta"
)


# PSID and CPS data from DAAG
data(
  "psid1",
  package = "DAAG"
)

data(
  "cps1",
  package = "DAAG"
)


###############################################################################
# 5. PREPARE CPS
###############################################################################

cps1_renamed <- cps1


names(cps1_renamed)[
  names(cps1_renamed) == "trt"
] <- "treat"

names(cps1_renamed)[
  names(cps1_renamed) == "educ"
] <- "education"

names(cps1_renamed)[
  names(cps1_renamed) == "hisp"
] <- "hispanic"

names(cps1_renamed)[
  names(cps1_renamed) == "marr"
] <- "married"

names(cps1_renamed)[
  names(cps1_renamed) == "nodeg"
] <- "nodegree"


# Preserve original construction
cps1_renamed$data_id <-
  seq_len(
    nrow(cps1_renamed)
  )


cps1_renamed <-
  cps1_renamed[
    ,
    c(
      "data_id",
      "treat",
      "age",
      "education",
      "black",
      "hispanic",
      "married",
      "nodegree",
      "re74",
      "re75",
      "re78"
    )
  ]


###############################################################################
# 6. PREPARE PSID
###############################################################################

psid1_renamed <- psid1


names(psid1_renamed)[
  names(psid1_renamed) == "trt"
] <- "treat"

names(psid1_renamed)[
  names(psid1_renamed) == "educ"
] <- "education"

names(psid1_renamed)[
  names(psid1_renamed) == "hisp"
] <- "hispanic"

names(psid1_renamed)[
  names(psid1_renamed) == "marr"
] <- "married"

names(psid1_renamed)[
  names(psid1_renamed) == "nodeg"
] <- "nodegree"


psid1_renamed$data_id <-
  seq_len(
    nrow(psid1_renamed)
  )


psid1_renamed <-
  psid1_renamed[
    ,
    c(
      "data_id",
      "treat",
      "age",
      "education",
      "black",
      "hispanic",
      "married",
      "nodegree",
      "re74",
      "re75",
      "re78"
    )
  ]


###############################################################################
# 7. PREPARE NSW
###############################################################################

nsw_treat <-
  exp_dat |>
  dplyr::filter(
    treat == 1
  ) |>
  dplyr::mutate(
    G = 1L
  )


nsw_ctrl <-
  exp_dat |>
  dplyr::filter(
    treat == 0
  ) |>
  dplyr::mutate(
    G = 2L
  )


# Preserve original group coding
psid1_renamed$G <- 3

cps1_renamed$G <- 4


###############################################################################
# 8. COMBINE DATA
###############################################################################
#
# IMPORTANT:
# Preserve the original base rbind() construction and row order.
###############################################################################

dat_all <- rbind(
  nsw_treat,
  nsw_ctrl,
  psid1_renamed,
  cps1_renamed
)


dat_all2 <- dat_all


d4 <-
  dat_all2 |>
  dplyr::filter(
    G %in% c(
      1,
      2,
      3,
      4
    )
  ) |>
  dplyr::mutate(
    G = factor(G)
  )


d4$ID <-
  seq_len(
    nrow(d4)
  )


###############################################################################
# 9. CHECK GROUP COUNTS
###############################################################################

cat(
  "\nGroup counts:\n"
)

print(
  d4 |>
    dplyr::count(G)
)


###############################################################################
# 10. NSW EXPERIMENTAL BENCHMARK
###############################################################################

benchmark <-
  mean(
    d4$re78[
      d4$G == "1"
    ]
  ) -
  mean(
    d4$re78[
      d4$G == "2"
    ]
  )


cat(
  "\nNSW experimental benchmark = ",
  benchmark,
  "\n\n",
  sep = ""
)


###############################################################################
# 11. MULTINOMIAL SOURCE-MEMBERSHIP MODEL
###############################################################################

fml_multi <-
  as.formula(
    paste(
      "G ~",
      paste(
        xvars,
        collapse = "+"
      )
    )
  )


fit_multi <-
  nnet::multinom(
    fml_multi,
    data = dat_all2,
    trace = FALSE
  )


P <-
  predict(
    fit_multi,
    type = "probs"
  )


P <- as.matrix(P)


p1 <- P[, "1"]
p2 <- P[, "2"]
p3 <- P[, "3"]
p4 <- P[, "4"]


pNSW <- p1 + p2


# Original transport propensity ratios
w1 <- pNSW / p1
w2 <- pNSW / p2
w3 <- pNSW / p3
w4 <- pNSW / p4


###############################################################################
# 12. MULTINOMIAL OBJECT FOR JOINT SANDWICH VARIANCE
###############################################################################

X_multi <-
  model.matrix(
    fml_multi,
    data = dat_all2
  )


B_multi <-
  as.matrix(
    coef(fit_multi)
  )


B_multi <-
  B_multi[
    c(
      "2",
      "3",
      "4"
    ),
    ,
    drop = FALSE
  ]


ps_obj <- list(
  
  fit = fit_multi,
  
  formula = fml_multi,
  
  X = X_multi,
  
  B = B_multi,
  
  phi_hat =
    as.vector(
      t(B_multi)
    ),
  
  P = P,
  
  xvars = xvars
)


###############################################################################
# 13. UNWEIGHTED ESTIMATORS
###############################################################################

uw_nsw <-
  unweighted_function(
    d4,
    cat1 = 1,
    cat0 = 2
  )


uw_psid <-
  unweighted_function(
    d4,
    cat1 = 1,
    cat0 = 3
  )


uw_cps <-
  unweighted_function(
    d4,
    cat1 = 1,
    cat0 = 4
  )


###############################################################################
# 14. IPW
###############################################################################

ipw_nsw <-
  IPW_multinom_joint_final(
    d4 = d4,
    donor_g = 2,
    ps_obj = ps_obj,
    benchmark = benchmark
  )


ipw_psid <-
  IPW_multinom_joint_final(
    d4 = d4,
    donor_g = 3,
    ps_obj = ps_obj,
    benchmark = benchmark
  )


ipw_cps <-
  IPW_multinom_joint_final(
    d4 = d4,
    donor_g = 4,
    ps_obj = ps_obj,
    benchmark = benchmark
  )


###############################################################################
# 15. AIPW ORIGINAL POINT ESTIMATES
###############################################################################
#
# The original code first obtains the historical point estimates using
# run_aipw_cv(). The revised joint-sandwich functions are then used for SEs.
###############################################################################

res_nsw_lm_old <-
  run_aipw_cv(
    cat0 = 2,
    w_use = w2,
    t0_model = "lm"
  )


res_nsw_gam_old <-
  run_aipw_cv(
    cat0 = 2,
    w_use = w2,
    t0_model = "gam"
  )


res_psid_lm_old <-
  run_aipw_cv(
    cat0 = 3,
    w_use = w3,
    t0_model = "lm"
  )


res_psid_gam_old <-
  run_aipw_cv(
    cat0 = 3,
    w_use = w3,
    t0_model = "gam"
  )


res_cps_lm_old <-
  run_aipw_cv(
    cat0 = 4,
    w_use = w4,
    t0_model = "lm"
  )


res_cps_gam_old <-
  run_aipw_cv(
    cat0 = 4,
    w_use = w4,
    t0_model = "gam"
  )


###############################################################################
# 16. AIPW JOINT-SANDWICH STANDARD ERRORS
###############################################################################

aipw_nsw_lm <-
  AIPW_multinom_joint_final(
    d4 = d4,
    donor_g = 2,
    ps_obj = ps_obj,
    model = "lm",
    K = K,
    seed = seed,
    benchmark = benchmark
  )


aipw_nsw_gam <-
  AIPW_multinom_joint_final(
    d4 = d4,
    donor_g = 2,
    ps_obj = ps_obj,
    model = "gam",
    K = K,
    seed = seed,
    benchmark = benchmark
  )


aipw_psid_lm <-
  AIPW_multinom_joint_final(
    d4 = d4,
    donor_g = 3,
    ps_obj = ps_obj,
    model = "lm",
    K = K,
    seed = seed,
    benchmark = benchmark
  )


aipw_psid_gam <-
  AIPW_multinom_joint_final(
    d4 = d4,
    donor_g = 3,
    ps_obj = ps_obj,
    model = "gam",
    K = K,
    seed = seed,
    benchmark = benchmark
  )


aipw_cps_lm <-
  AIPW_multinom_joint_final(
    d4 = d4,
    donor_g = 4,
    ps_obj = ps_obj,
    model = "lm",
    K = K,
    seed = seed,
    benchmark = benchmark
  )


aipw_cps_gam <-
  AIPW_multinom_joint_final(
    d4 = d4,
    donor_g = 4,
    ps_obj = ps_obj,
    model = "gam",
    K = K,
    seed = seed,
    benchmark = benchmark
  )


###############################################################################
# 17. PRESERVE ORIGINAL AIPW POINT ESTIMATES
###############################################################################
#
# Point estimate = original run_aipw_cv()
# SE             = revised joint multinomial sandwich
###############################################################################

aipw_nsw_lm$ATE <-
  res_nsw_lm_old$ATE

aipw_nsw_gam$ATE <-
  res_nsw_gam_old$ATE


aipw_psid_lm$ATE <-
  res_psid_lm_old$ATE

aipw_psid_gam$ATE <-
  res_psid_gam_old$ATE


aipw_cps_lm$ATE <-
  res_cps_lm_old$ATE

aipw_cps_gam$ATE <-
  res_cps_gam_old$ATE


###############################################################################
# 18. PACKAGE COMPARATOR SOURCE OBJECTS
###############################################################################

source_nsw <-
  make_transport_source_data(
    d4 = d4,
    donor_g = 2,
    xvars = xvars
  )


source_psid <-
  make_transport_source_data(
    d4 = d4,
    donor_g = 3,
    xvars = xvars
  )


source_cps <-
  make_transport_source_data(
    d4 = d4,
    donor_g = 4,
    xvars = xvars
  )


###############################################################################
# 19. EBPS / CBPS / EBCW
###############################################################################

# NSW
ebps_nsw <-
  transport_ebps(
    source_nsw,
    benchmark
  )

cbps_nsw <-
  transport_cbps(
    source_nsw,
    benchmark
  )

ebcw_nsw <-
  transport_ebcw(
    source_nsw,
    benchmark
  )


# PSID
ebps_psid <-
  transport_ebps(
    source_psid,
    benchmark
  )

cbps_psid <-
  transport_cbps(
    source_psid,
    benchmark
  )

ebcw_psid <-
  transport_ebcw(
    source_psid,
    benchmark
  )


# CPS
ebps_cps <-
  transport_ebps(
    source_cps,
    benchmark
  )

cbps_cps <-
  transport_cbps(
    source_cps,
    benchmark
  )

ebcw_cps <-
  transport_ebcw(
    source_cps,
    benchmark
  )


###############################################################################
# 20. ET / HD / EL
###############################################################################

run_three_gec <- function(
    donor_g
) {
  
  ans <-
    lapply(
      c(
        "ET",
        "HD",
        "EL"
      ),
      function(ent) {
        
        tryCatch(
          
          run_GEC_transport_simulation_style(
            d4 = d4,
            donor_g = donor_g,
            ps_obj = ps_obj,
            entropy = ent,
            K = K,
            seed = seed
          ),
          
          error = function(e) {
            
            list(
              ATE = NA_real_,
              SE = NA_real_,
              variance = NA_real_,
              CI = c(
                NA_real_,
                NA_real_
              ),
              success = FALSE,
              ESS = NA_real_,
              MaxW = NA_real_,
              CalibrationMax = NA_real_,
              variance_type =
                paste0(
                  "failed: ",
                  conditionMessage(e)
                )
            )
          }
        )
      }
    )
  
  names(ans) <-
    c(
      "ET",
      "HD",
      "EL"
    )
  
  ans
}


gec_nsw <-
  run_three_gec(
    donor_g = 2
  )


gec_psid <-
  run_three_gec(
    donor_g = 3
  )


gec_cps <-
  run_three_gec(
    donor_g = 4
  )


###############################################################################
# 21. BUILD LONG-FORM RESULTS
###############################################################################

rows <- list()


add_uw <- function(
    donor,
    x
) {
  
  make_final_real_row(
    donor = donor,
    method = "Unweighted",
    est = as.numeric(
      x[1]
    ),
    se = as.numeric(
      x[2]
    ),
    benchmark = benchmark,
    variance_type =
      "ordinary unweighted variance"
  )
}


add_fit <- function(
    donor,
    method,
    x
) {
  
  make_final_real_row(
    
    donor = donor,
    
    method = method,
    
    est = x$ATE,
    
    se = x$SE,
    
    benchmark = benchmark,
    
    ess =
      if (!is.null(x$ESS)) {
        x$ESS
      } else if (!is.null(x$ESS0)) {
        x$ESS0
      } else {
        NA_real_
      },
    
    maxw =
      if (!is.null(x$MaxW)) {
        x$MaxW
      } else if (!is.null(x$MaxW0)) {
        x$MaxW0
      } else {
        NA_real_
      },
    
    cal =
      if (!is.null(x$CalibrationMax)) {
        x$CalibrationMax
      } else {
        NA_real_
      },
    
    success =
      if (!is.null(x$success)) {
        isTRUE(x$success)
      } else {
        is.finite(x$ATE) &&
          is.finite(x$SE)
      },
    
    variance_type =
      if (!is.null(x$variance_type)) {
        x$variance_type
      } else {
        ""
      }
  )
}


###############################################################################
# 22. NSW RESULTS
###############################################################################

rows[[length(rows) + 1L]] <-
  add_uw(
    "NSW",
    uw_nsw
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "NSW",
    "IPW",
    ipw_nsw
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "NSW",
    "EBPS",
    ebps_nsw
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "NSW",
    "CBPS",
    cbps_nsw
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "NSW",
    "EBCW",
    ebcw_nsw
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "NSW",
    "AIPW (LM)",
    aipw_nsw_lm
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "NSW",
    "AIPW (GAM)",
    aipw_nsw_gam
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "NSW",
    "ET",
    gec_nsw[["ET"]]
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "NSW",
    "HD",
    gec_nsw[["HD"]]
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "NSW",
    "EL",
    gec_nsw[["EL"]]
  )


###############################################################################
# 23. PSID RESULTS
###############################################################################

rows[[length(rows) + 1L]] <-
  add_uw(
    "PSID",
    uw_psid
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "PSID",
    "IPW",
    ipw_psid
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "PSID",
    "EBPS",
    ebps_psid
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "PSID",
    "CBPS",
    cbps_psid
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "PSID",
    "EBCW",
    ebcw_psid
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "PSID",
    "AIPW (LM)",
    aipw_psid_lm
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "PSID",
    "AIPW (GAM)",
    aipw_psid_gam
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "PSID",
    "ET",
    gec_psid[["ET"]]
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "PSID",
    "HD",
    gec_psid[["HD"]]
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "PSID",
    "EL",
    gec_psid[["EL"]]
  )


###############################################################################
# 24. CPS RESULTS
###############################################################################

rows[[length(rows) + 1L]] <-
  add_uw(
    "CPS",
    uw_cps
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "CPS",
    "IPW",
    ipw_cps
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "CPS",
    "EBPS",
    ebps_cps
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "CPS",
    "CBPS",
    cbps_cps
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "CPS",
    "EBCW",
    ebcw_cps
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "CPS",
    "AIPW (LM)",
    aipw_cps_lm
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "CPS",
    "AIPW (GAM)",
    aipw_cps_gam
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "CPS",
    "ET",
    gec_cps[["ET"]]
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "CPS",
    "HD",
    gec_cps[["HD"]]
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "CPS",
    "EL",
    gec_cps[["EL"]]
  )


###############################################################################
# 25. COMBINE LONG-FORM RESULTS
###############################################################################

final_long <-
  dplyr::bind_rows(
    rows
  )


method_order <- c(
  "Unweighted",
  "IPW",
  "EBPS",
  "CBPS",
  "EBCW",
  "AIPW (LM)",
  "AIPW (GAM)",
  "ET",
  "HD",
  "EL"
)


final_long$Method <-
  factor(
    final_long$Method,
    levels = method_order
  )


final_long <-
  final_long |>
  dplyr::mutate(
    Donor =
      factor(
        Donor,
        levels = c(
          "NSW",
          "PSID",
          "CPS"
        )
      )
  ) |>
  dplyr::arrange(
    Method,
    Donor
  )


###############################################################################
# 26. BUILD MANUSCRIPT TABLE 10
###############################################################################
#
# The manuscript defines evaluation bias relative to the rounded NSW
# experimental benchmark of 1794.
###############################################################################

paper_benchmark <- 1794


get_result <- function(
    donor_name,
    method_name
) {
  
  z <-
    final_long |>
    dplyr::filter(
      as.character(Donor) ==
        donor_name,
      as.character(Method) ==
        method_name
    )
  
  
  if (nrow(z) != 1L) {
    
    return(
      c(
        Estimate = NA_real_,
        SE = NA_real_
      )
    )
  }
  
  
  c(
    Estimate =
      z$Estimate[1],
    
    SE =
      z$SE[1]
  )
}


table10_rows <-
  lapply(
    method_order,
    function(m) {
      
      nsw <-
        get_result(
          "NSW",
          m
        )
      
      psid <-
        get_result(
          "PSID",
          m
        )
      
      cps <-
        get_result(
          "CPS",
          m
        )
      
      
      data.frame(
        
        Estimator = m,
        
        NSW_Est =
          round(
            nsw["Estimate"],
            0
          ),
        
        NSW_SE =
          round(
            nsw["SE"],
            0
          ),
        
        PSID_Est =
          round(
            psid["Estimate"],
            0
          ),
        
        PSID_EB =
          round(
            psid["Estimate"] -
              paper_benchmark,
            0
          ),
        
        PSID_SE =
          round(
            psid["SE"],
            0
          ),
        
        CPS_Est =
          round(
            cps["Estimate"],
            0
          ),
        
        CPS_EB =
          round(
            cps["Estimate"] -
              paper_benchmark,
            0
          ),
        
        CPS_SE =
          round(
            cps["SE"],
            0
          ),
        
        stringsAsFactors = FALSE
      )
    }
  )


table10 <-
  dplyr::bind_rows(
    table10_rows
  )


###############################################################################
# 27. LABEL METHODS EXACTLY AS IN MANUSCRIPT
###############################################################################

table10 <-
  table10 |>
  dplyr::mutate(
    Estimator =
      dplyr::recode(
        Estimator,
        "AIPW (LM)" = "AIPW-LM",
        "AIPW (GAM)" = "AIPW-GAM"
      )
  )


###############################################################################
# 28. SAVE TABLE 10
###############################################################################

dir.create(
  "Causal Inference/results/table10",
  recursive = TRUE,
  showWarnings = FALSE
)


write.csv(
  table10,
  "Causal Inference/results/table10/table10_lalonde.csv",
  row.names = FALSE
)


###############################################################################
# 29. OPTIONAL: SAVE FULL UNROUNDED RESULTS FOR CHECKING
###############################################################################

write.csv(
  final_long,
  "Causal Inference/results/table10/table10_lalonde_diagnostics.csv",
  row.names = FALSE
)


###############################################################################
# 30. PRINT TABLE 10
###############################################################################

cat(
  "\n============================================================\n"
)

cat(
  "TABLE 10: NSW / PSID / CPS REAL-DATA APPLICATION\n"
)

cat(
  "============================================================\n\n"
)


print(
  table10,
  row.names = FALSE
)
