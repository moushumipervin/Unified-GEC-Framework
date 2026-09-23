# =============================================================================
# ATE simulation functions: fully commented version
# =============================================================================
# Purpose of this script
# ----------------------
# This script collects all helper functions used in the ATE simulation study.
# The main goal here is not to change the code, but to explain clearly what
# each function is doing so that another reader can follow the workflow by
# reading this file alone.
#
# Broad structure of the script
# -----------------------------
# 1. Generate data under four simulation scenarios:
#      OR1PS1, OR1PS2, OR2PS1, OR2PS2
#    where OR = outcome regression model and PS = propensity score model.
#
# 2. Build cross-fitting folds and fit nuisance models:
#      - linear model (LM)
#      - generalized additive model (GAM)
#
# 3. Compute ATE estimators:
#      - IPW
#      - AIPW with cross-fitting
#      - oCBPS / CBPS
#      - entropy balancing (EBPS in the stored object name)
#      - EBCW
#      - proposed HD and ET estimators
#
# 4. Run one replication, one scenario, or all scenarios.
#
# 5. Build the final 2 x 2 boxplot figure for method comparison.
#
# Required packages
# -----------------
# install.packages(c("mgcv",  "CBPS", "ATE", "WeightIt",
#                    "dplyr", "tidyr", "ggplot2", "patchwork"))
# =============================================================================

suppressPackageStartupMessages({
  library(mgcv)
  library(CBPS)
  library(ATE)
  library(WeightIt)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# ---------------------------------------------------
# 1. Data-generating mechanism
# ---------------------------------------------------

# Function: generate_ate_data
# ---------------------------
# Generates one simulated dataset for a chosen outcome model and propensity
# score model. The function returns:
#   - the simulated observed dataset
#   - the true ATE used in that scenario
#   - a scenario label such as OR1PS1
#
# Scenario meaning:
#   outcome_model = 1  -> linear outcome model
#   outcome_model = 2  -> nonlinear outcome model
#   ps_model      = 1  -> linear logistic propensity score
#   ps_model      = 2  -> nonlinear logistic propensity score
#
# The observed outcome is y = D*y1 + (1-D)*y0.
# A fitted propensity score pi.hat is also stored in the returned data frame.
generate_ate_data <- function(n , p , outcome_model , ps_model ,
                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(p == 4)
  
  Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(Z) <- paste0("x", seq_len(p))
  
  eta <- rnorm(n)
  alpha1 <- rep(1, p)
  alpha2 <- rep(1, p)
  
  if (outcome_model == 1) {
    beta1 <- rep(0.5, p)
    y1 <- 1 + as.numeric(Z %*% beta1) + rnorm(n)
    y0 <-     as.numeric(Z %*% beta1) + rnorm(n)
    true_ate <- 1
  } else {
    Zexp <- exp(pmax(pmin(Z, 3), -3))
    Zt <- (Z - 1)^3 - Z^2 + Z / (1 + Zexp) + 10
    y1 <- 10 + as.numeric(Z %*% alpha1) + 0.5 * as.numeric(Zt %*% alpha2) + eta
    y0 <-      as.numeric(Z %*% alpha1) + 0.5 * as.numeric(Zt %*% alpha2) + eta
    true_ate <- 10
  }
  
  if (ps_model == 1) {
    px <- plogis(
      -(0.25 +
          Z[,1] +
          0.5 * Z[,2] -
          0.5 * Z[,3] -
          0.1 * Z[,4])
    )
  } else {
    px <- plogis(
      -(Z[,1] -
          0.5 * Z[,2] * Z[,1] -
          Z[,3]^2 +
          0.5 * Z[,4]^3)
    )
  }
  
  D <- rbinom(n, size = 1, prob = px)
  y <- D * y1 + (1 - D) * y0
  
  dat <- data.frame(y = y, Z, D = D)
  dat$pi.hat <- fitted(glm(D ~ ., data = dat[, c(paste0("x", 1:p), "D")], family = binomial()))
  
  list(
    data = dat,
    true_ate = true_ate,
    scenario = paste0("OR", outcome_model, "PS", ps_model)
  )
}


# ==============================================================
# TRUE conditional outcome regression for OR2
# ==============================================================

true_m_or2 <- function(dat) {
  
  xvars <- paste0("x", 1:4)
  
  Z <- as.matrix(
    dat[, xvars, drop = FALSE]
  )
  
  # EXACTLY reproduce the transformation in generate_ate_data()
  Zexp <- exp(
    pmax(
      pmin(Z, 3),
      -3
    )
  )
  
  Zt <- (Z - 1)^3 -
    Z^2 +
    Z / (1 + Zexp) +
    10
  
  m0 <- rowSums(Z) +
    0.5 * rowSums(Zt)
  
  m1 <- 10 +
    rowSums(Z) +
    0.5 * rowSums(Zt)
  
  list(
    m1 = as.numeric(m1),
    m0 = as.numeric(m0)
  )
}
# ---------------------------------------------------
# 2. Cross-fitting helpers
# ---------------------------------------------------

# Function: build_SU_folds
# ------------------------
# Builds K roughly equal folds for two index sets:
#   idx_S  = the observations treated as the labeled/training side
#   idx_U0 = the observations treated as the complementary side
#
# The function keeps track of both row indices and IDs for each fold.
# This is later used by the cross-fitting routines.
build_SU_folds <- function(
    df, K, id_col = "ID",
    idx_S, idx_U0,
    seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  df$ID <- seq_len(nrow(df))
  
  id <- df[[id_col]]
  
  nS  <- length(idx_S)
  nU0 <- length(idx_U0)
  if (nS < K || nU0 < K) stop("Need at least K labeled and K unlabeled observations.")
  
  split_K <- function(idx, K) {
    if (length(idx) == 0) return(rep(list(integer(0)), K))
    idx <- sample(idx)
    split(idx, rep(1:K, length.out = length(idx)))
  }
  
  S_parts  <- split_K(idx_S,  K)
  U0_parts <- split_K(idx_U0, K)
  
  folds <- vector("list", K)
  fold_id_S <- integer(nrow(df))
  fold_id_U <- integer(nrow(df))
  
  for (k in 1:K) {
    S_k_idx  <- S_parts[[k]]
    U0_k_idx <- U0_parts[[k]]
    U_k_idx  <- c(S_k_idx, U0_k_idx)
    
    folds[[k]] <- list(
      S_k_ids = id[S_k_idx],
      U_k_ids = id[U_k_idx],
      S_k_idx = S_k_idx,
      U_k_idx = U_k_idx
    )
    
    fold_id_S[S_k_idx] <- k
    fold_id_U[U_k_idx] <- k
  }
  
  list(
    folds = folds,
    fold_id_S = fold_id_S,
    fold_id_U = fold_id_U
  )
}


# Function: k_fold_function_lm
# ----------------------------
# Performs cross-fitting using a linear model for the outcome regression.
#
# For each fold:
#   1. Fit the outcome model on labeled training observations only
#   2. Predict outcomes on the validation fold
#   3. Store the fold-level data together with y.hat
#
# The returned object is a list of K validation-fold data frames.
k_fold_function_lm <- function(df, K, idx_S, idx_U0, seed) {
  res <- build_SU_folds(df, K, id_col = "ID", idx_S, idx_U0, seed = seed)
  df$ID <- seq_len(nrow(df))
  lab_idx_all <- idx_S
  
  data_unlabeled <- list()
  
  for (k in 1:K) {
    train_idx_lab <- setdiff(lab_idx_all, res$folds[[k]]$S_k_idx)
    validation_data_labeled <- df[res$folds[[k]]$S_k_idx, ]
    validation_data_unlabeled <- df[res$folds[[k]]$U_k_idx, ]
    
    train_data_labeled <- df[train_idx_lab, ]
    
    formula <- as.formula(
      paste0(
        "y~",
        paste0(
          colnames(subset(train_data_labeled, select = -c(D, y, ID, pi.hat))),
          collapse = "+"
        )
      )
    )
    
    lm_model <- lm(
      formula,
      data = as.data.frame(subset(train_data_labeled, select = -c(D, ID, pi.hat)))
    )
    
    y_hat <- predict(lm_model, newdata = as.data.frame(validation_data_unlabeled))
    validation_data_unlabeled$y.hat <- y_hat
    
    data_unlabeled[[k]] <- validation_data_unlabeled
  }
  
  return(data_unlabeled)
}


# Function: k_fold_function_gam
# -----------------------------
# Same idea as k_fold_function_lm, but uses a GAM with smooth terms for all
# covariates instead of an ordinary linear model.
k_fold_function_gam <- function(df, K, idx_S, idx_U0, seed) {
  res <- build_SU_folds(df, K, id_col = "ID", idx_S, idx_U0, seed = seed)
  df$ID <- seq_len(nrow(df))
  lab_idx_all <- idx_S
  
  data_unlabeled <- list()
  
  for (k in 1:K) {
    train_idx_lab <- setdiff(lab_idx_all, res$folds[[k]]$S_k_idx)
    validation_data_labeled <- df[res$folds[[k]]$S_k_idx, ]
    validation_data_unlabeled <- df[res$folds[[k]]$U_k_idx, ]
    
    train_data_labeled <- df[train_idx_lab, ]
    
    formula <- as.formula(
      paste0(
        "y~",
        paste0(
          "s(",
          colnames(subset(train_data_labeled, select = -c(D, y, ID, pi.hat))),
          ")",
          collapse = "+"
        )
      )
    )
    
    gam_model <- mgcv::gam(
      formula,
      data = as.data.frame(subset(train_data_labeled, select = -c(D, ID, pi.hat)))
    )
    
    y_hat <- predict(gam_model, newdata = as.data.frame(validation_data_unlabeled))
    validation_data_unlabeled$y.hat <- y_hat
    
    data_unlabeled[[k]] <- validation_data_unlabeled
  }
  
  return(data_unlabeled)
}



# ==============================================================
# Cross-fitted AIPW inference
# Joint sandwich for (phi, tau)
#
# Accounts for estimation of propensity-score parameter phi.
# Cross-fitted m1(X), m0(X) are treated as fixed.
# ==============================================================

estimate_aipw_inference <- function(
    fold_T1,
    fold_T0,
    true_ATE = NA_real_
) {
  
  # ============================================================
  # 1. Reconstruct cross-fitted datasets
  # ============================================================
  
  dat1 <- do.call(
    rbind,
    fold_T1
  )
  
  dat0 <- do.call(
    rbind,
    fold_T0
  )
  
  
  # ------------------------------------------------------------
  # Align both prediction datasets by subject ID
  # ------------------------------------------------------------
  
  dat1 <- dat1[
    order(dat1$ID),
  ]
  
  dat0 <- dat0[
    order(dat0$ID),
  ]
  
  
  stopifnot(
    identical(
      as.integer(dat1$ID),
      as.integer(dat0$ID)
    )
  )
  
  
  # ============================================================
  # 2. Basic quantities
  # ============================================================
  
  y <- as.numeric(dat1$y)
  
  D <- as.numeric(dat1$D)
  
  m1 <- as.numeric(dat1$y.hat)
  
  m0 <- as.numeric(dat0$y.hat)
  
  N <- length(y)
  
  
  # ============================================================
  # 3. Fit propensity-score model
  #
  # Same working logistic model used in the GEC analysis
  # ============================================================
  
  xvars <- grep(
    "^x\\d+$",
    names(dat1),
    value = TRUE
  )
  
  
  ps_formula <- as.formula(
    paste(
      "D ~",
      paste(
        xvars,
        collapse = " + "
      )
    )
  )
  
  
  ps_fit <- glm(
    ps_formula,
    data = dat1,
    family = binomial()
  )
  
  
  X <- model.matrix(
    ps_fit
  )
  
  
  phi_hat <- as.numeric(
    coef(ps_fit)
  )
  
  
  pi_hat <- plogis(
    as.vector(
      X %*% phi_hat
    )
  )
  
  
  # Numerical protection
  pi_hat <- pmin(
    pmax(
      pi_hat,
      1e-8
    ),
    1 - 1e-8
  )
  
  
  # ============================================================
  # 4. AIPW point estimator
  # ============================================================
  
  psi1 <-
    m1 +
    D / pi_hat *
    (y - m1)
  
  
  psi0 <-
    m0 +
    (1 - D) /
    (1 - pi_hat) *
    (y - m0)
  
  
  pseudo_ate <-
    psi1 - psi0
  
  
  tau_hat <-
    mean(
      pseudo_ate
    )
  
  
  # ============================================================
  # 5. Joint estimating equations
  #
  # beta = (phi, tau)
  #
  # Psi_phi,i = X_i (D_i - pi_i)
  #
  # Psi_tau,i =
  #    m1_i - m0_i
  #  + D_i/pi_i (Y_i-m1_i)
  #  - (1-D_i)/(1-pi_i) (Y_i-m0_i)
  #  - tau
  # ============================================================
  
  
  # ------------------------------------------------------------
  # Propensity-score estimating equation
  # ------------------------------------------------------------
  
  Psi_phi <-
    X *
    as.numeric(
      D - pi_hat
    )
  
  
  # ------------------------------------------------------------
  # AIPW estimating equation
  # ------------------------------------------------------------
  
  Psi_tau <-
    pseudo_ate -
    tau_hat
  
  
  # ------------------------------------------------------------
  # Full observation-level estimating-function matrix
  # ------------------------------------------------------------
  
  Psi <- cbind(
    Psi_phi,
    tau = Psi_tau
  )
  
  
  # ============================================================
  # 6. Bread matrix
  #
  # A = - 1/N sum d Psi_i / d beta'
  # ============================================================
  
  r <- ncol(X)
  
  p_total <- r + 1
  
  
  A <- matrix(
    0,
    nrow = p_total,
    ncol = p_total
  )
  
  
  idx_phi <- seq_len(r)
  
  idx_tau <- r + 1
  
  
  # ------------------------------------------------------------
  # A_phi,phi
  #
  # -d[X(D-pi)]/dphi'
  # = X X' pi(1-pi)
  # ------------------------------------------------------------
  
  A_phi_phi <-
    crossprod(
      X,
      X *
        as.numeric(
          pi_hat *
            (1 - pi_hat)
        )
    ) / N
  
  
  A[
    idx_phi,
    idx_phi
  ] <- A_phi_phi
  
  
  # ------------------------------------------------------------
  # A_phi,tau = 0
  # ------------------------------------------------------------
  
  A[
    idx_phi,
    idx_tau
  ] <- 0
  
  
  # ------------------------------------------------------------
  # A_tau,phi
  #
  # Derivative of AIPW estimating equation with respect to phi.
  #
  # d/dphi [
  # D/pi (Y-m1)
  # -
  # (1-D)/(1-pi) (Y-m0)
  # ]
  #
  # =
  # -D (1-pi)/pi (Y-m1) X
  # -(1-D) pi/(1-pi) (Y-m0) X
  #
  # Bread uses minus derivative.
  # ------------------------------------------------------------
  
  term_treated <-
    D *
    ((1 - pi_hat) / pi_hat) *
    (y - m1)
  
  
  term_control <-
    (1 - D) *
    (pi_hat / (1 - pi_hat)) *
    (y - m0)
  
  
  A_tau_phi <-
    colMeans(
      X *
        as.numeric(
          term_treated +
            term_control
        )
    )
  
  
  A[
    idx_tau,
    idx_phi
  ] <- A_tau_phi
  
  
  # ------------------------------------------------------------
  # A_tau,tau
  #
  # Psi_tau contains -tau,
  # therefore -d Psi_tau/dtau = 1
  # ------------------------------------------------------------
  
  A[
    idx_tau,
    idx_tau
  ] <- 1
  
  
  # ============================================================
  # 7. Meat
  #
  # B = 1/N sum Psi_i Psi_i'
  # ============================================================
  
  B <-
    crossprod(
      Psi
    ) / N
  
  
  # ============================================================
  # 8. Sandwich covariance
  #
  # Var(beta_hat)
  # =
  # 1/N A^{-1} B A^{-T}
  # ============================================================
  
  A_inv <- tryCatch(
    solve(A),
    error = function(e) NULL
  )
  
  
  if (is.null(A_inv) ||
      any(!is.finite(A_inv))) {
    
    return(
      list(
        
        ATE = tau_hat,
        
        SE = NA_real_,
        
        variance = NA_real_,
        
        coverage = NA_real_,
        
        CI_lower = NA_real_,
        
        CI_upper = NA_real_,
        
        ESS1 = NA_real_,
        ESS0 = NA_real_,
        
        max_weight1 = NA_real_,
        max_weight0 = NA_real_,
        
        A = A,
        B = B,
        Psi = Psi,
        
        success = FALSE
      )
    )
  }
  
  
  V_beta <-
    (
      A_inv %*%
        B %*%
        t(A_inv)
    ) / N
  
  
  # ============================================================
  # 9. Extract variance of tau
  # ============================================================
  
  var_tau <-
    as.numeric(
      V_beta[
        idx_tau,
        idx_tau
      ]
    )
  
  
  if (!is.finite(var_tau) ||
      var_tau < 0) {
    
    se_tau <- NA_real_
    
  } else {
    
    se_tau <-
      sqrt(
        var_tau
      )
  }
  
  
  # ============================================================
  # 10. Confidence interval
  # ============================================================
  
  if (is.finite(se_tau)) {
    
    CI_lower <-
      tau_hat -
      qnorm(0.975) *
      se_tau
    
    
    CI_upper <-
      tau_hat +
      qnorm(0.975) *
      se_tau
    
  } else {
    
    CI_lower <- NA_real_
    
    CI_upper <- NA_real_
  }
  
  
  # ============================================================
  # 11. Coverage
  # ============================================================
  
  if (is.finite(true_ATE) &&
      is.finite(CI_lower) &&
      is.finite(CI_upper)) {
    
    coverage <-
      as.numeric(
        CI_lower <= true_ATE &&
          true_ATE <= CI_upper
      )
    
  } else {
    
    coverage <- NA_real_
  }
  
  
  # ============================================================
  # 12. Propensity-weight diagnostics
  #
  # These are PS-weight diagnostics, not AIPW
  # calibration weights.
  # ============================================================
  
  w1 <-
    1 /
    pi_hat[D == 1]
  
  
  w0 <-
    1 /
    (1 - pi_hat[D == 0])
  
  
  ESS1 <-
    sum(w1)^2 /
    sum(w1^2)
  
  
  ESS0 <-
    sum(w0)^2 /
    sum(w0^2)
  
  
  max_weight1 <-
    max(
      w1
    )
  
  
  max_weight0 <-
    max(
      w0
    )
  
  
  # ============================================================
  # 13. Return
  # ============================================================
  
  list(
    
    ATE = tau_hat,
    
    SE = se_tau,
    
    variance = var_tau,
    
    CI_lower = CI_lower,
    
    CI_upper = CI_upper,
    
    coverage = coverage,
    
    ESS1 = ESS1,
    
    ESS0 = ESS0,
    
    max_weight1 =
      max_weight1,
    
    max_weight0 =
      max_weight0,
    
    phi_hat =
      phi_hat,
    
    pi_hat =
      pi_hat,
    
    A = A,
    
    B = B,
    
    V_beta =
      V_beta,
    
    Psi = Psi,
    
    success =
      is.finite(se_tau)
  )
}




# ==============================================================
# Common helper:
# convert estimate + SE + weights into the same quantities
# stored for ET / HD / CE
# ==============================================================

make_benchmark_output <- function(
    ate_hat,
    se_hat,
    dat,
    weights = NULL,
    true_ATE = NA_real_
) {
  
  # ------------------------------------------------------------
  # Coverage
  # ------------------------------------------------------------
  
  if (is.finite(ate_hat) &&
      is.finite(se_hat) &&
      is.finite(true_ATE)) {
    
    lower <-
      ate_hat -
      qnorm(0.975) * se_hat
    
    upper <-
      ate_hat +
      qnorm(0.975) * se_hat
    
    coverage <-
      as.numeric(
        lower <= true_ATE &&
          true_ATE <= upper
      )
    
  } else {
    
    coverage <- NA_real_
  }
  
  
  # ------------------------------------------------------------
  # Weight diagnostics
  # ------------------------------------------------------------
  
  if (!is.null(weights) &&
      length(weights) == nrow(dat) &&
      all(is.finite(weights))) {
    
    w1 <- weights[
      dat$D == 1
    ]
    
    w0 <- weights[
      dat$D == 0
    ]
    
    
    ESS1 <-
      sum(w1)^2 /
      sum(w1^2)
    
    
    ESS0 <-
      sum(w0)^2 /
      sum(w0^2)
    
    
    max_weight1 <-
      max(w1)
    
    
    max_weight0 <-
      max(w0)
    
  } else {
    
    ESS1 <- NA_real_
    ESS0 <- NA_real_
    
    max_weight1 <- NA_real_
    max_weight0 <- NA_real_
  }
  
  
  list(
    
    ATE = ate_hat,
    
    SE = se_hat,
    
    coverage = coverage,
    
    ESS1 = ESS1,
    ESS0 = ESS0,
    
    max_weight1 = max_weight1,
    max_weight0 = max_weight0,
    
    success =
      is.finite(ate_hat) &&
      is.finite(se_hat)
  )
}


estimate_ebps_inference <- function(
    dat,
    true_ATE = NA_real_
) {
  
  xvars <- grep(
    "^x\\d+$",
    names(dat),
    value = TRUE
  )
  
  form <- as.formula(
    paste(
      "D ~",
      paste(
        xvars,
        collapse = " + "
      )
    )
  )
  
  
  # ------------------------------------------------------------
  # Entropy-balancing weights
  # ------------------------------------------------------------
  
  fit_w <- tryCatch(
    
    WeightIt::weightit(
      form,
      data = dat,
      method = "ebal",
      estimand = "ATE"
    ),
    
    error = function(e) NULL
  )
  
  
  if (is.null(fit_w)) {
    
    return(
      make_benchmark_output(
        ate_hat = NA_real_,
        se_hat = NA_real_,
        dat = dat,
        true_ATE = true_ATE
      )
    )
  }
  
  
  w <- as.numeric(
    fit_w$weights
  )
  
  
  # ------------------------------------------------------------
  # Weighted outcome model
  # ------------------------------------------------------------
  fit_y <- WeightIt::lm_weightit(
    y ~ D,
    data = dat,
    weightit = fit_w
  )
  
  ate_hat <- as.numeric(
    coef(fit_y)["D"]
  )
  
  se_hat <- sqrt(
    vcov(fit_y)["D", "D"]
  )
  
  
  
  make_benchmark_output(
    ate_hat = ate_hat,
    se_hat = se_hat,
    dat = dat,
    weights = w,
    true_ATE = true_ATE
  )
}



estimate_cbps_inference <- function(
    dat,
    true_ATE = NA_real_
) {
  
  xvars <- grep(
    "^x\\d+$",
    names(dat),
    value = TRUE
  )
  
  form <- reformulate(
    xvars,
    response = "D"
  )
  
  fit_w <- tryCatch(
    WeightIt::weightit(
      form,
      data = dat,
      method = "cbps",
      estimand = "ATE",
      over = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit_w)) {
    return(
      make_benchmark_output(
        ate_hat = NA_real_,
        se_hat = NA_real_,
        dat = dat,
        true_ATE = true_ATE
      )
    )
  }
  
  fit_y <- tryCatch(
    WeightIt::lm_weightit(
      y ~ D,
      data = dat,
      weightit = fit_w,
      vcov = "asympt"
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit_y)) {
    return(
      make_benchmark_output(
        ate_hat = NA_real_,
        se_hat = NA_real_,
        dat = dat,
        true_ATE = true_ATE
      )
    )
  }
  
  ate_hat <- as.numeric(
    coef(fit_y)["D"]
  )
  
  se_hat <- sqrt(
    vcov(fit_y)["D", "D"]
  )
  
  w <- as.numeric(
    fit_w$weights
  )
  
  make_benchmark_output(
    ate_hat = ate_hat,
    se_hat = se_hat,
    dat = dat,
    weights = w,
    true_ATE = true_ATE
  )
}

estimate_ocbps_inference <- function(
    dat,
    true_ATE = NA_real_
) {
  
  xvars <- grep(
    "^x\\d+$",
    names(dat),
    value = TRUE
  )
  
  X2 <- as.matrix(
    cbind(
      1,
      dat[, xvars, drop = FALSE]
    )
  )
  
  
  # ------------------------------------------------------------
  # 1. Fit oCBPS
  # ------------------------------------------------------------
  
  ocbps_fit <- tryCatch(
    CBPS::CBPS(
      dat$D ~ X2,
      ATT = 0,
      method = "exact",
      baseline.formula = ~ X2,
      diff.formula = ~ X2
    ),
    error = function(e) NULL
  )
  
  
  if (is.null(ocbps_fit)) {
    
    return(
      make_benchmark_output(
        ate_hat = NA_real_,
        se_hat = NA_real_,
        dat = dat,
        true_ATE = true_ATE
      )
    )
  }
  
  
  # ------------------------------------------------------------
  # 2. Extract oCBPS weights
  # ------------------------------------------------------------
  
  w <- as.numeric(
    ocbps_fit$weights
  )
  
  
  # ------------------------------------------------------------
  # 3. Natural weighted ATE
  #
  # Use normalized weighted means within each treatment arm
  # ------------------------------------------------------------
  
  mu1_hat <- sum(
    w[dat$D == 1] *
      dat$y[dat$D == 1]
  ) /
    sum(
      w[dat$D == 1]
    )
  
  
  mu0_hat <- sum(
    w[dat$D == 0] *
      dat$y[dat$D == 0]
  ) /
    sum(
      w[dat$D == 0]
    )
  
  
  ate_hat <-
    mu1_hat -
    mu0_hat
  
  
  # ------------------------------------------------------------
  # 4. oCBPS package variance
  # ------------------------------------------------------------
  
  ocbps_inf <- tryCatch(
    CBPS::AsyVar(
      Y = dat$y,
      CBPS_obj = ocbps_fit,
      method = "oCBPS",
      CI = 0.95
    ),
    error = function(e) NULL
  )
  
  
  if (is.null(ocbps_inf)) {
    
    return(
      make_benchmark_output(
        ate_hat = ate_hat,
        se_hat = NA_real_,
        dat = dat,
        weights = w,
        true_ATE = true_ATE
      )
    )
  }
  
  
  # ------------------------------------------------------------
  # 5. Extract SE only
  #
  # IMPORTANT:
  # keep the ATE from the oCBPS weights above.
  # Do NOT replace it with ocbps_inf$mu.hat.
  # ------------------------------------------------------------
  
  se_hat <- as.numeric(
    ocbps_inf$std.err
  )
  
  
  # ------------------------------------------------------------
  # 6. Return
  # ------------------------------------------------------------
  
  make_benchmark_output(
    ate_hat = ate_hat,
    se_hat = se_hat,
    dat = dat,
    weights = w,
    true_ATE = true_ATE
  )
}
estimate_ocbps_inference1<- function(
    dat,
    true_ATE = NA_real_
) {
  
  xvars <- grep(
    "^x\\d+$",
    names(dat),
    value = TRUE
  )
  
  Xcov <- as.matrix(
    dat[, xvars, drop = FALSE]
  )
  
  # ------------------------------------------------------------
  # 1. Fit oCBPS
  # ------------------------------------------------------------
  
  X2 <- cbind(
    1,
    Xcov
  )
  
  ocbps_fit <- tryCatch(
    CBPS::CBPS(
      dat$D ~ X2,
      ATT = 0,
      method = "exact",
      baseline.formula = ~ X2,
      diff.formula = ~ X2,
      standardize = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(ocbps_fit)) {
    
    return(
      make_benchmark_output(
        ate_hat = NA_real_,
        se_hat = NA_real_,
        dat = dat,
        true_ATE = true_ATE
      )
    )
  }
  
  
  # ------------------------------------------------------------
  # 2. Basic quantities
  # ------------------------------------------------------------
  
  y <- as.numeric(
    dat$y
  )
  
  D <- as.numeric(
    dat$D
  )
  
  N <- nrow(
    dat
  )
  
  beta_hat <- as.numeric(
    ocbps_fit$coefficients
  )
  
  
  # Design matrix actually corresponding to beta
  X <- model.matrix(
    ~ x1 + x2 + x3 + x4,
    data = dat
  )
  
  if (
    length(beta_hat) != ncol(X)
  ) {
    
    stop(
      "Dimension of oCBPS coefficient vector does not match X."
    )
  }
  
  
  # ------------------------------------------------------------
  # 3. Fitted propensity scores
  # ------------------------------------------------------------
  
  pi_hat <- plogis(
    as.vector(
      X %*% beta_hat
    )
  )
  
  pi_hat <- pmin(
    pmax(
      pi_hat,
      1e-6
    ),
    1 - 1e-6
  )
  
  
  # ------------------------------------------------------------
  # 4. ATE point estimate
  #
  # Horvitz-Thompson style because standardize = FALSE
  # ------------------------------------------------------------
  
  mu1_hat <- mean(
    D * y / pi_hat
  )
  
  mu0_hat <- mean(
    (1 - D) * y /
      (1 - pi_hat)
  )
  
  ate_hat <-
    mu1_hat -
    mu0_hat
  
  
  # ------------------------------------------------------------
  # 5. Joint parameter
  #
  # theta = (beta, mu1, mu0)
  # ------------------------------------------------------------
  
  theta_hat <- c(
    beta_hat,
    mu1_hat,
    mu0_hat
  )
  
  p_beta <- length(
    beta_hat
  )
  
  
  # ------------------------------------------------------------
  # 6. Observation-level estimating equations
  # ------------------------------------------------------------
  
  psi_fun <- function(theta) {
    
    beta <- theta[
      seq_len(
        p_beta
      )
    ]
    
    mu1 <- theta[
      p_beta + 1
    ]
    
    mu0 <- theta[
      p_beta + 2
    ]
    
    
    pi <- plogis(
      as.vector(
        X %*% beta
      )
    )
    
    pi <- pmin(
      pmax(
        pi,
        1e-6
      ),
      1 - 1e-6
    )
    
    
    # ----------------------------------------------------------
    # oCBPS balance moments
    #
    # For your baseline and diff formula both using X2,
    # construct the two sets of moments used by optimal CBPS.
    # ----------------------------------------------------------
    
    # baseline balance component
    h1 <- (
      D / pi -
        (1 - D) /
        (1 - pi)
    )
    
    Psi_base <-
      X *
      as.numeric(
        h1
      )
    
    
    # treatment-effect / difference component
    h2 <-
      D / pi -
      1
    
    Psi_diff <-
      X *
      as.numeric(
        h2
      )
    
    
    # ----------------------------------------------------------
    # Outcome mean equations
    # ----------------------------------------------------------
    
    Psi_mu1 <-
      D / pi *
      (
        y -
          mu1
      )
    
    Psi_mu0 <-
      (1 - D) /
      (1 - pi) *
      (
        y -
          mu0
      )
    
    
    cbind(
      Psi_base,
      Psi_diff,
      mu1 = Psi_mu1,
      mu0 = Psi_mu0
    )
  }
  
  
  # ------------------------------------------------------------
  # 7. Evaluate estimating functions
  # ------------------------------------------------------------
  
  Psi_hat <- psi_fun(
    theta_hat
  )
  
  
  # ------------------------------------------------------------
  # 8. Meat
  #
  # B = 1/N sum psi_i psi_i'
  # ------------------------------------------------------------
  
  B <- crossprod(
    Psi_hat
  ) / N
  
  
  # ------------------------------------------------------------
  # 9. Numerical Jacobian
  #
  # A = - d mean(psi) / d theta'
  # ------------------------------------------------------------
  
  J <- numDeriv::jacobian(
    func = function(theta) {
      
      colMeans(
        psi_fun(
          theta
        )
      )
    },
    
    x = theta_hat
  )
  
  
  A <- -J
  
  
  # ------------------------------------------------------------
  # 10. Generalized inverse
  #
  # oCBPS can be overidentified because there are more
  # moment equations than parameters.
  # ------------------------------------------------------------
  
  A_inv <- tryCatch(
    MASS::ginv(
      A
    ),
    error = function(e) NULL
  )
  
  if (
    is.null(A_inv) ||
    any(
      !is.finite(
        A_inv
      )
    )
  ) {
    
    return(
      make_benchmark_output(
        ate_hat = ate_hat,
        se_hat = NA_real_,
        dat = dat,
        weights = as.numeric(
          ocbps_fit$weights
        ),
        true_ATE = true_ATE
      )
    )
  }
  
  
  # ------------------------------------------------------------
  # 11. Sandwich covariance
  # ------------------------------------------------------------
  
  V <- (
    A_inv %*%
      B %*%
      t(
        A_inv
      )
  ) / N
  
  
  # ------------------------------------------------------------
  # 12. ATE contrast
  #
  # tau = mu1 - mu0
  # ------------------------------------------------------------
  
  contrast <- rep(
    0,
    length(
      theta_hat
    )
  )
  
  contrast[
    p_beta + 1
  ] <- 1
  
  contrast[
    p_beta + 2
  ] <- -1
  
  
  var_ate <- as.numeric(
    t(
      contrast
    ) %*%
      V %*%
      contrast
  )
  
  
  se_hat <- if (
    is.finite(
      var_ate
    ) &&
    var_ate >= 0
  ) {
    
    sqrt(
      var_ate
    )
    
  } else {
    
    NA_real_
  }
  
  
  # ------------------------------------------------------------
  # 13. Original oCBPS weights for diagnostics
  # ------------------------------------------------------------
  
  w <- as.numeric(
    ocbps_fit$weights
  )
  
  
  # ------------------------------------------------------------
  # 14. Return
  # ------------------------------------------------------------
  
  make_benchmark_output(
    ate_hat = ate_hat,
    se_hat = se_hat,
    dat = dat,
    weights = w,
    true_ATE = true_ATE
  )
}
estimate_ocbps_inference_old <- function(
    dat,
    true_ATE = NA_real_
) {
  
  xvars <- grep(
    "^x\\d+$",
    names(dat),
    value = TRUE
  )
  
  X2 <- as.matrix(
    cbind(
      1,
      dat[, xvars, drop = FALSE]
    )
  )
  
  
  # ------------------------------------------------------------
  # 1. Estimate oCBPS weights
  # ------------------------------------------------------------
  
  ocbps_fit <- tryCatch(
    CBPS::CBPS(
      dat$D ~ X2,
      ATT = 0,
      method = "exact",
      baseline.formula = ~ X2,
      diff.formula = ~ X2
    ),
    error = function(e) NULL
  )
  
  
  if (is.null(ocbps_fit)) {
    
    return(
      make_benchmark_output(
        ate_hat = NA_real_,
        se_hat = NA_real_,
        dat = dat,
        true_ATE = true_ATE
      )
    )
  }
  
  
  w <- as.numeric(
    ocbps_fit$weights
  )
  
  
  # ------------------------------------------------------------
  # 2. Weighted ATE model
  # ------------------------------------------------------------
  
  fit_y <- tryCatch(
    lm(
      y ~ D,
      data = dat,
      weights = w
    ),
    error = function(e) NULL
  )
  
  
  if (is.null(fit_y)) {
    
    return(
      make_benchmark_output(
        ate_hat = NA_real_,
        se_hat = NA_real_,
        dat = dat,
        true_ATE = true_ATE
      )
    )
  }
  
  
  # ------------------------------------------------------------
  # 3. ATE
  # ------------------------------------------------------------
  
  ate_hat <- as.numeric(
    coef(fit_y)["D"]
  )
  
  
  # ------------------------------------------------------------
  # 4. Robust sandwich variance
  # ------------------------------------------------------------
  
  V <- tryCatch(
    sandwich::vcovHC(
      fit_y,
      type = "HC0"
    ),
    error = function(e) NULL
  )
  
  
  if (is.null(V)) {
    
    se_hat <- NA_real_
    
  } else {
    
    var_hat <- as.numeric(
      V["D", "D"]
    )
    
    se_hat <- if (
      is.finite(var_hat) &&
      var_hat >= 0
    ) {
      sqrt(var_hat)
    } else {
      NA_real_
    }
  }
  
  
  # ------------------------------------------------------------
  # 5. Return
  # ------------------------------------------------------------
  
  make_benchmark_output(
    ate_hat = ate_hat,
    se_hat = se_hat,
    dat = dat,
    weights = w,
    true_ATE = true_ATE
  )
}



estimate_ocbps_inference_old <- function(
    dat,
    true_ATE = NA_real_
) {
  
  xvars <- grep(
    "^x\\d+$",
    names(dat),
    value = TRUE
  )
  
  
  X2 <- as.matrix(
    cbind(
      1,
      dat[
        ,
        xvars,
        drop = FALSE
      ]
    )
  )
  
  
  ocbps_fit <- tryCatch(
    
    CBPS::CBPS(
      dat$D ~ X2,
      ATT = 0,
      method = "exact",
      baseline.formula = ~ X2,
      diff.formula = ~ X2
    ),
    
    error = function(e) NULL
  )
  
  
  if (is.null(ocbps_fit)) {
    
    return(
      make_benchmark_output(
        ate_hat = NA_real_,
        se_hat = NA_real_,
        dat = dat,
        true_ATE = true_ATE
      )
    )
  }
  
  
  # ------------------------------------------------------------
  # Package-based oCBPS ATE and variance
  # ------------------------------------------------------------
  
  ocbps_inf <- tryCatch(
    
    CBPS::AsyVar(
      Y = dat$y,
      CBPS_obj = ocbps_fit,
      method = "oCBPS",
      CI = 0.95
    ),
    
    error = function(e) NULL
  )
  
  
  if (is.null(ocbps_inf)) {
    
    return(
      make_benchmark_output(
        ate_hat = NA_real_,
        se_hat = NA_real_,
        dat = dat,
        true_ATE = true_ATE
      )
    )
  }
  
  
  # oCBPS ATE from CBPS package
  ate_hat <- as.numeric(
    ocbps_inf$mu.hat
  )
  
  # oCBPS SE from CBPS package
  se_hat <- as.numeric(
    ocbps_inf$std.err
  )
  
  # Optional variance
  var_hat <- as.numeric(
    ocbps_inf$var
  )
  
  
  # Keep weights only for ESS / max-weight diagnostics
  w <- as.numeric(
    ocbps_fit$weights
  ) 
  make_benchmark_output(
    ate_hat = ate_hat,
    se_hat = se_hat,
    dat = dat,
    weights = w,
    true_ATE = true_ATE
  )
}


estimate_ebcw_inference <- function(
    dat,
    true_ATE = NA_real_
) {
  
  xvars <- grep(
    "^x\\d+$",
    names(dat),
    value = TRUE
  )
  
  X_ebcw <- data.frame(
    dat[, xvars, drop = FALSE]
  )
  
  fit <- tryCatch(
    ATE::ATE(
      as.numeric(dat$y),
      as.integer(dat$D),
      X_ebcw,
      ATT = FALSE
    ),
    error = function(e) {
      message("EBCW error: ", e$message)
      NULL
    }
  )
  
  if (is.null(fit)) {
    return(
      make_benchmark_output(
        ate_hat = NA_real_,
        se_hat = NA_real_,
        dat = dat,
        weights = NULL,
        true_ATE = true_ATE
      )
    )
  }
  
  sm <- summary(fit)
  
  ate_hat <- as.numeric(
    sm$Estimate["ATE", "Estimate"]
  )
  
  se_hat <- as.numeric(
    sm$Estimate["ATE", "Std. Error"]
  )
  w1 <- as.numeric(sm$weights.treat)
  w0 <- as.numeric(sm$weights.placebo)
  
  w <- w1 + w0
  make_benchmark_output(
    ate_hat = ate_hat,
    se_hat = se_hat,
    dat = dat,
    weights = w,
    true_ATE = true_ATE
  )
}
estimate_ebcw_inference2 <- function(
    dat,
    true_ATE = NA_real_
) {
  
  xvars <- grep(
    "^x\\d+$",
    names(dat),
    value = TRUE
  )
  
  # ------------------------------------------------------------
  # Covariates supplied directly to EBCW
  # No propensity-score term
  # ------------------------------------------------------------
  
  X_ebcw <- as.matrix(
    dat[, xvars, drop = FALSE]
  )
  
  
  # ------------------------------------------------------------
  # Fit EBCW
  # ------------------------------------------------------------
  
  fit <- tryCatch(
    ATE::ATE(
      Y = dat$y,
      Ti = dat$D,
      X = X_ebcw,
      ATT = FALSE
    ),
    error = function(e) NULL
  )
  
  
  if (is.null(fit)) {
    
    return(
      make_benchmark_output(
        ate_hat = NA_real_,
        se_hat = NA_real_,
        dat = dat,
        true_ATE = true_ATE
      )
    )
  }
  
  
  # ------------------------------------------------------------
  # Package summary
  # ------------------------------------------------------------
  
  sm <- summary(fit)
  
  
  # ATE
  ate_hat <- as.numeric(
    sm$Estimate["ATE", "Estimate"]
  )
  
  
  # Package variance-based SE
  se_hat <- as.numeric(
    sm$Estimate["ATE", "Std. Error"]
  )
  
  
  # ------------------------------------------------------------
  # Return
  # ------------------------------------------------------------
  
  make_benchmark_output(
    ate_hat = ate_hat,
    se_hat = se_hat,
    dat = dat,
    weights = NULL,
    true_ATE = true_ATE
  )
}


estimate_ebcw_inference_old <- function(
    dat,
    true_ATE = NA_real_
) {
  
  xvars <- grep(
    "^x\\d+$",
    names(dat),
    value = TRUE
  )
  
  
  X_chan <- data.frame(
    
    dat[
      ,
      xvars,
      drop = FALSE
    ],
    
    log_pi_hat =
      log(
        pmax(
          dat$pi.hat,
          1e-8
        )
      )
  )
  
  
  fit <- tryCatch(
    
    ATE::ATE(
      dat$y,
      dat$D,
      X_chan,
      ATT = FALSE
    ),
    
    error = function(e) NULL
  )
  
  
  if (is.null(fit)) {
    
    return(
      make_benchmark_output(
        ate_hat = NA_real_,
        se_hat = NA_real_,
        dat = dat,
        true_ATE = true_ATE
      )
    )
  }
  
  
  sm <- summary(fit)
  
  
  ate_hat <-
    sm$Estimate[
      3,
      1
    ]
  
  
  # Check your summary(fit) object once.
  # If column 2 is the reported SE:
  se_hat <-
    sm$Estimate[
      3,
      2
    ]
  
  
  make_benchmark_output(
    ate_hat = ate_hat,
    se_hat = se_hat,
    dat = dat,
    weights = NULL,
    true_ATE = true_ATE
  )
}

# ==============================================================
# IPW + joint sandwich inference
# ==============================================================

estimate_ipw_inference <- function(
    dat,
    true_ATE = NA_real_
) {
  
  # ------------------------------------------------------------
  # 1. Propensity-score model
  # ------------------------------------------------------------
  
  xvars <- grep("^x\\d+$", names(dat), value = TRUE)
  
  ps_formula <- as.formula(
    paste(
      "D ~",
      paste(xvars, collapse = " + ")
    )
  )
  
  ps_fit <- glm(
    ps_formula,
    data = dat,
    family = binomial()
  )
  
  X <- model.matrix(ps_fit)
  
  phi_hat <- as.numeric(coef(ps_fit))
  
  pi_hat <- pmin(
    pmax(fitted(ps_fit), 1e-8),
    1 - 1e-8
  )
  
  y <- as.numeric(dat$y)
  D <- as.numeric(dat$D)
  N <- nrow(dat)
  
  
  # ------------------------------------------------------------
  # 2. Standard HT-style IPW means
  # ------------------------------------------------------------
  
  mu1_hat <- mean(
    D * y / pi_hat
  )
  
  mu0_hat <- mean(
    (1 - D) * y / (1 - pi_hat)
  )
  
  ate_hat <- mu1_hat - mu0_hat
  
  
  # ------------------------------------------------------------
  # 3. Joint parameter
  # theta = (phi, mu1, mu0)
  # ------------------------------------------------------------
  
  theta_hat <- c(
    phi_hat,
    mu1_hat,
    mu0_hat
  )
  
  r <- ncol(X)
  
  
  # ------------------------------------------------------------
  # 4. Observation-level estimating equations
  # ------------------------------------------------------------
  
  psi_fun <- function(theta) {
    
    phi <- theta[seq_len(r)]
    
    mu1 <- theta[r + 1]
    mu0 <- theta[r + 2]
    
    pi <- plogis(
      as.vector(X %*% phi)
    )
    
    pi <- pmin(
      pmax(pi, 1e-8),
      1 - 1e-8
    )
    
    # propensity-score score
    Psi_phi <-
      X * as.numeric(D - pi)
    
    # HT IPW mean equations
    Psi_mu1 <-
      D * y / pi - mu1
    
    Psi_mu0 <-
      (1 - D) * y / (1 - pi) - mu0
    
    cbind(
      Psi_phi,
      mu1 = Psi_mu1,
      mu0 = Psi_mu0
    )
  }
  
  
  # ------------------------------------------------------------
  # 5. Sandwich variance
  # ------------------------------------------------------------
  
  Psi_hat <- psi_fun(theta_hat)
  
  B <- crossprod(Psi_hat) / N
  
  J <- numDeriv::jacobian(
    func = function(theta) {
      colMeans(
        psi_fun(theta)
      )
    },
    x = theta_hat
  )
  
  A <- -J
  
  A_inv <- tryCatch(
    solve(A),
    error = function(e) NULL
  )
  
  if (is.null(A_inv) ||
      any(!is.finite(A_inv))) {
    
    return(
      list(
        ATE = ate_hat,
        SE = NA_real_,
        coverage = NA_real_,
        ESS1 = NA_real_,
        ESS0 = NA_real_,
        max_weight1 = NA_real_,
        max_weight0 = NA_real_,
        success = FALSE
      )
    )
  }
  
  V <- (
    A_inv %*%
      B %*%
      t(A_inv)
  ) / N
  
  
  # ------------------------------------------------------------
  # 6. ATE = mu1 - mu0
  # ------------------------------------------------------------
  
  contrast <- rep(
    0,
    length(theta_hat)
  )
  
  contrast[r + 1] <- 1
  contrast[r + 2] <- -1
  
  var_ate <- as.numeric(
    t(contrast) %*%
      V %*%
      contrast
  )
  
  se_ate <- if (
    is.finite(var_ate) &&
    var_ate >= 0
  ) {
    sqrt(var_ate)
  } else {
    NA_real_
  }
  
  
  # ------------------------------------------------------------
  # 7. Coverage
  # ------------------------------------------------------------
  
  coverage <- if (
    is.finite(true_ATE) &&
    is.finite(se_ate)
  ) {
    
    lower <-
      ate_hat -
      qnorm(0.975) * se_ate
    
    upper <-
      ate_hat +
      qnorm(0.975) * se_ate
    
    as.numeric(
      lower <= true_ATE &&
        true_ATE <= upper
    )
    
  } else {
    
    NA_real_
  }
  
  
  # ------------------------------------------------------------
  # 8. Weight diagnostics
  # ------------------------------------------------------------
  
  w1 <- 1 / pi_hat[D == 1]
  
  w0 <- 1 / (1 - pi_hat[D == 0])
  
  ESS1 <-
    sum(w1)^2 /
    sum(w1^2)
  
  ESS0 <-
    sum(w0)^2 /
    sum(w0^2)
  
  
  # ------------------------------------------------------------
  # 9. Return
  # ------------------------------------------------------------
  
  list(
    ATE = ate_hat,
    SE = se_ate,
    coverage = coverage,
    
    ESS1 = ESS1,
    ESS0 = ESS0,
    
    max_weight1 = max(w1),
    max_weight0 = max(w0),
    
    pi_hat = pi_hat,
    phi_hat = phi_hat,
    
    success = is.finite(se_ate)
  )
}


###############################################################################
# ANALYTIC SANDWICH VARIANCE FOR CAUSAL GEC ATE
###############################################################################

gec_ate_sandwich <- function(
    y,
    T,
    X,                 # PS design matrix INCLUDING intercept
    yhat1,
    yhat0,
    lambda1,
    lambda0,
    entropy = c("SL", "EL", "ET", "HD", "CE"),
    phi_hat = NULL,
    pi_hat = NULL,
    tol = 1e-10
) {
  
  entropy <- match.arg(entropy)
  
  y <- as.numeric(y)
  T <- as.numeric(T)
  X <- as.matrix(X)
  
  yhat1 <- as.numeric(yhat1)
  yhat0 <- as.numeric(yhat0)
  
  lambda1 <- as.numeric(lambda1)
  lambda0 <- as.numeric(lambda0)
  
  N <- length(y)
  
  stopifnot(
    length(T) == N,
    nrow(X) == N,
    length(yhat1) == N,
    length(yhat0) == N
  )
  
  
  ###########################################################################
  # 1. Propensity score
  ###########################################################################
  
  # If phi_hat is supplied, reconstruct pi_hat from the logistic model.
  if (!is.null(phi_hat)) {
    
    phi_hat <- as.numeric(phi_hat)
    
    if (length(phi_hat) != ncol(X)) {
      stop("length(phi_hat) must equal ncol(X).")
    }
    
    pi_hat <- plogis(
      as.vector(X %*% phi_hat)
    )
    
  } else {
    
    if (is.null(pi_hat)) {
      stop("Supply either phi_hat or pi_hat.")
    }
    
    pi_hat <- as.numeric(pi_hat)
  }
  
  # numerical protection
  pi_hat <- pmin(
    pmax(pi_hat, 1e-8),
    1 - 1e-8
  )
  
  
  ###########################################################################
  # 2. Entropy-specific debiasing covariate q_G(p)
  ###########################################################################
  
  q_fun <- function(p) {
    
    if (entropy == "SL") {
      return(1 / p - 1)
    }
    
    if (entropy == "EL") {
      return(-p)
    }
    
    if (entropy == "ET") {
      return(log(1 / p))
    }
    
    if (entropy == "HD") {
      return(-sqrt(p) / 2)
    }
    
    if (entropy == "CE") {
      return(log1p(-p))
    }
  }
  
  
  ###########################################################################
  # 3. Derivative q_G'(p)
  ###########################################################################
  
  q_prime <- function(p) {
    
    if (entropy == "SL") {
      return(-1 / p^2)
    }
    
    if (entropy == "EL") {
      return(rep(-1, length(p)))
    }
    
    if (entropy == "ET") {
      return(-1 / p)
    }
    
    if (entropy == "HD") {
      return(-1 / (4 * sqrt(p)))
    }
    
    if (entropy == "CE") {
      return(-1 / (1 - p))
    }
  }
  
  
  ###########################################################################
  # 4. Entropy-specific weight omega(eta)
  ###########################################################################
  
  weight_fun <- function(eta) {
    
    if (entropy == "SL") {
      return(1 + eta)
    }
    
    if (entropy == "EL") {
      return(-1 / eta)
    }
    
    if (entropy == "ET") {
      return(exp(eta))
    }
    
    if (entropy == "HD") {
      return(1 / (4 * eta^2))
    }
    
    if (entropy == "CE") {
      return(-1 / expm1(eta))
    }
  }
  
  
  ###########################################################################
  # 5. omega_dot = d omega / d eta
  ###########################################################################
  
  weight_deriv <- function(eta) {
    
    if (entropy == "SL") {
      return(rep(1, length(eta)))
    }
    
    if (entropy == "EL") {
      return(1 / eta^2)
    }
    
    if (entropy == "ET") {
      return(exp(eta))
    }
    
    if (entropy == "HD") {
      return(-1 / (2 * eta^3))
    }
    
    if (entropy == "CE") {
      return(
        exp(eta) / expm1(eta)^2
      )
    }
  }
  
  
  ###########################################################################
  # 6. Arm-specific probabilities
  ###########################################################################
  
  p1 <- pi_hat
  p0 <- 1 - pi_hat
  
  D1 <- T
  D0 <- 1 - T
  
  
  ###########################################################################
  # 7. Calibration vectors S_ti
  #
  #    S_ti = (1, mhat_t(X_i), q_G(p_ti))'
  ###########################################################################
  
  S1 <- cbind(
    intercept = 1,
    yhat = yhat1,
    g_pi = q_fun(p1)
  )
  
  S0 <- cbind(
    intercept = 1,
    yhat = yhat0,
    g_pi = q_fun(p0)
  )
  
  q1 <- ncol(S1)
  q0 <- ncol(S0)
  r  <- ncol(X)
  
  if (length(lambda1) != q1) {
    stop("lambda1 dimension does not match S1.")
  }
  
  if (length(lambda0) != q0) {
    stop("lambda0 dimension does not match S0.")
  }
  
  
  ###########################################################################
  # 8. eta_ti and weights
  ###########################################################################
  
  eta1 <- as.vector(
    S1 %*% lambda1
  )
  
  eta0 <- as.vector(
    S0 %*% lambda0
  )
  
  w1 <- weight_fun(eta1)
  w0 <- weight_fun(eta0)
  
  wdot1 <- weight_deriv(eta1)
  wdot0 <- weight_deriv(eta0)
  
  
  ###########################################################################
  # 9. Closed-form theta estimates
  ###########################################################################
  
  theta1_hat <- sum(
    D1 * w1 * y
  ) / N
  
  theta0_hat <- sum(
    D0 * w0 * y
  ) / N
  
  ate_hat <- theta1_hat - theta0_hat
  
  
  ###########################################################################
  # 10. Construct J_ti = d S_ti / d phi'
  #
  # Only the final component of S_ti depends on phi.
  ###########################################################################
  
  # logistic derivative:
  #
  # d pi / d phi' = pi(1-pi) X'
  #
  
  dpi <- pi_hat * (1 - pi_hat)
  
  # dimensions:
  # J1_array[i,,] = q1 x r
  # J0_array[i,,] = q0 x r
  
  J1_array <- array(
    0,
    dim = c(N, q1, r)
  )
  
  J0_array <- array(
    0,
    dim = c(N, q0, r)
  )
  
  # d q(p1) / d phi
  dq1_dphi <- (
    q_prime(p1) * dpi
  ) * X
  
  # p0 = 1-pi, hence
  # d p0 / d phi = - d pi / d phi
  dq0_dphi <- (
    q_prime(p0) * (-dpi)
  ) * X
  
  for (i in seq_len(N)) {
    
    J1_array[i, q1, ] <- dq1_dphi[i, ]
    
    J0_array[i, q0, ] <- dq0_dphi[i, ]
  }
  
  
  ###########################################################################
  # 11. ANALYTIC BREAD
  ###########################################################################
  
  p_total <- r + q1 + 1 + q0 + 1
  
  A <- matrix(
    0,
    nrow = p_total,
    ncol = p_total
  )
  
  
  # Parameter positions
  idx_phi <- seq_len(r)
  
  idx_lam1 <- (
    max(idx_phi) + 1
  ):(
    max(idx_phi) + q1
  )
  
  idx_theta1 <- max(idx_lam1) + 1
  
  idx_lam0 <- (
    idx_theta1 + 1
  ):(
    idx_theta1 + q0
  )
  
  idx_theta0 <- max(idx_lam0) + 1
  
  
  ###########################################################################
  # A_phi,phi
  ###########################################################################
  
  A_phi_phi <- crossprod(
    X,
    X * as.numeric(pi_hat * (1 - pi_hat))
  ) / N
  
  A[
    idx_phi,
    idx_phi
  ] <- A_phi_phi
  
  
  ###########################################################################
  # Initialize arm-specific bread blocks
  ###########################################################################
  
  A_lam1_lam1 <- matrix(0, q1, q1)
  A_lam0_lam0 <- matrix(0, q0, q0)
  
  A_lam1_phi <- matrix(0, q1, r)
  A_lam0_phi <- matrix(0, q0, r)
  
  A_theta1_lam1 <- matrix(0, 1, q1)
  A_theta0_lam0 <- matrix(0, 1, q0)
  
  A_theta1_phi <- matrix(0, 1, r)
  A_theta0_phi <- matrix(0, 1, r)
  
  
  ###########################################################################
  # Compute observation-wise analytic derivatives
  ###########################################################################
  
  for (i in seq_len(N)) {
    
    s1 <- matrix(
      S1[i, ],
      ncol = 1
    )
    
    s0 <- matrix(
      S0[i, ],
      ncol = 1
    )
    
    J1 <- matrix(
      J1_array[i, , ],
      nrow = q1,
      ncol = r
    )
    
    J0 <- matrix(
      J0_array[i, , ],
      nrow = q0,
      ncol = r
    )
    
    
    #########################################################################
    # A_lambda1,lambda1
    #########################################################################
    
    A_lam1_lam1 <-
      A_lam1_lam1 -
      D1[i] *
      wdot1[i] *
      (s1 %*% t(s1)) / N
    
    
    #########################################################################
    # A_lambda0,lambda0
    #########################################################################
    
    A_lam0_lam0 <-
      A_lam0_lam0 -
      D0[i] *
      wdot0[i] *
      (s0 %*% t(s0)) / N
    
    
    #########################################################################
    # A_theta1,lambda1
    #########################################################################
    
    A_theta1_lam1 <-
      A_theta1_lam1 -
      D1[i] *
      (y[i] - theta1_hat) *
      wdot1[i] *
      t(s1) / N
    
    
    #########################################################################
    # A_theta0,lambda0
    #########################################################################
    
    A_theta0_lam0 <-
      A_theta0_lam0 -
      D0[i] *
      (y[i] - theta0_hat) *
      wdot0[i] *
      t(s0) / N
    
    
    #########################################################################
    # A_lambda1,phi
    #########################################################################
    
    term1_1 <-
      (D1[i] * w1[i] - 1) * J1
    
    term2_1 <-
      D1[i] *
      wdot1[i] *
      s1 %*%
      (t(lambda1) %*% J1)
    
    A_lam1_phi <-
      A_lam1_phi -
      (term1_1 + term2_1) / N
    
    
    #########################################################################
    # A_lambda0,phi
    #########################################################################
    
    term1_0 <-
      (D0[i] * w0[i] - 1) * J0
    
    term2_0 <-
      D0[i] *
      wdot0[i] *
      s0 %*%
      (t(lambda0) %*% J0)
    
    A_lam0_phi <-
      A_lam0_phi -
      (term1_0 + term2_0) / N
    
    
    #########################################################################
    # A_theta1,phi
    #########################################################################
    
    A_theta1_phi <-
      A_theta1_phi -
      D1[i] *
      (y[i] - theta1_hat) *
      wdot1[i] *
      (t(lambda1) %*% J1) / N
    
    
    #########################################################################
    # A_theta0,phi
    #########################################################################
    
    A_theta0_phi <-
      A_theta0_phi -
      D0[i] *
      (y[i] - theta0_hat) *
      wdot0[i] *
      (t(lambda0) %*% J0) / N
  }
  
  
  ###########################################################################
  # A_theta,theta empirical weighted Jacobian
  ###########################################################################
  
  A_theta1_theta1 <- sum(
    D1 * w1
  ) / N
  
  A_theta0_theta0 <- sum(
    D0 * w0
  ) / N
  
  
  ###########################################################################
  # Fill full A matrix
  ###########################################################################
  
  A[
    idx_lam1,
    idx_phi
  ] <- A_lam1_phi
  
  A[
    idx_lam1,
    idx_lam1
  ] <- A_lam1_lam1
  
  A[
    idx_theta1,
    idx_phi
  ] <- A_theta1_phi
  
  A[
    idx_theta1,
    idx_lam1
  ] <- A_theta1_lam1
  
  A[
    idx_theta1,
    idx_theta1
  ] <- A_theta1_theta1
  
  
  A[
    idx_lam0,
    idx_phi
  ] <- A_lam0_phi
  
  A[
    idx_lam0,
    idx_lam0
  ] <- A_lam0_lam0
  
  A[
    idx_theta0,
    idx_phi
  ] <- A_theta0_phi
  
  A[
    idx_theta0,
    idx_lam0
  ] <- A_theta0_lam0
  
  A[
    idx_theta0,
    idx_theta0
  ] <- A_theta0_theta0
  
  
  ###########################################################################
  # 12. Observation-level estimating functions Psi_i
  ###########################################################################
  
  Psi_phi <- X *
    as.numeric(T - pi_hat)
  
  Psi_lam1 <- S1 *
    as.numeric(D1 * w1 - 1)
  
  Psi_theta1 <- D1 *
    w1 *
    (y - theta1_hat)
  
  Psi_lam0 <- S0 *
    as.numeric(D0 * w0 - 1)
  
  Psi_theta0 <- D0 *
    w0 *
    (y - theta0_hat)
  
  
  Psi <- cbind(
    Psi_phi,
    Psi_lam1,
    theta1 = Psi_theta1,
    Psi_lam0,
    theta0 = Psi_theta0
  )
  
  
  ###########################################################################
  # 13. Meat
  ###########################################################################
  
  B <- crossprod(Psi) / N
  
  
  ###########################################################################
  # 14. Sandwich covariance
  ###########################################################################
  
  A_inv <- tryCatch(
    solve(A),
    error = function(e) NULL
  )
  
  if (is.null(A_inv)) {
    
    warning("Analytic bread matrix is singular.")
    
    return(
      list(
        ATE = ate_hat,
        theta1 = theta1_hat,
        theta0 = theta0_hat,
        SE = NA_real_,
        variance = NA_real_,
        A = A,
        B = B,
        Psi = Psi
      )
    )
  }
  
  V_beta <- (
    A_inv %*%
      B %*%
      t(A_inv)
  ) / N
  
  
  ###########################################################################
  # 15. Delta-method contrast for ATE = theta1 - theta0
  ###########################################################################
  
  contrast <- rep(
    0,
    p_total
  )
  
  contrast[idx_theta1] <- 1
  contrast[idx_theta0] <- -1
  
  var_ate <- as.numeric(
    t(contrast) %*%
      V_beta %*%
      contrast
  )
  
  se_ate <- sqrt(
    max(var_ate, 0)
  )
  
  ci_lower <- ate_hat -
    qnorm(0.975) * se_ate
  
  ci_upper <- ate_hat +
    qnorm(0.975) * se_ate
  
  
  ###########################################################################
  # Diagnostics
  ###########################################################################
  
  calibration1 <- max(
    abs(
      colSums(
        S1 *
          as.numeric(D1 * w1 - 1)
      )
    )
  )
  
  calibration0 <- max(
    abs(
      colSums(
        S0 *
          as.numeric(D0 * w0 - 1)
      )
    )
  )
  
  ###########################################################################
  # Weight diagnostics: max weight and effective sample size
  ###########################################################################
  
  w1_obs <- w1[D1 == 1]
  w0_obs <- w0[D0 == 1]
  
  max_weight1 <- max(w1_obs)
  max_weight0 <- max(w0_obs)
  
  ESS1 <- sum(w1_obs)^2 / sum(w1_obs^2)
  ESS0 <- sum(w0_obs)^2 / sum(w0_obs^2)
  ###########################################################################
  # Return
  ###########################################################################
  
  list(
    
    ATE = ate_hat,
    
    theta1 = theta1_hat,
    theta0 = theta0_hat,
    
    variance = var_ate,
    SE = se_ate,
    
    CI_lower = ci_lower,
    CI_upper = ci_upper,
    
    A = A,
    B = B,
    V_beta = V_beta,
    Psi = Psi,
    
    S1 = S1,
    S0 = S0,
    
    eta1 = eta1,
    eta0 = eta0,
    
    weight1 = w1,
    weight0 = w0,
    max_weight1 = max_weight1,
    max_weight0 = max_weight0,
    
    ESS1 = ESS1,
    ESS0 = ESS0,
    weight_deriv1 = wdot1,
    weight_deriv0 = wdot0,
    
    J1 = J1_array,
    J0 = J0_array,
    
    calibration_raw1 = calibration1,
    calibration_raw0 = calibration0,
    
    theta1_bread = A_theta1_theta1,
    theta0_bread = A_theta0_theta0,
    
    entropy = entropy
  )
}
# Function: compute_crossfit_aipw
# --------------------------------
# Combines the foldwise predictions from the treated and control sides to form
# the cross-fitted AIPW estimator of the ATE.
compute_crossfit_aipw <- function(fold_T1, fold_T0) {
  data_T1 <- do.call(rbind, fold_T1)
  yhat_T1 <- data_T1$y.hat
  y_T1 <- data_T1$y
  pi_hat_T1 <- data_T1$pi.hat
  T1 <- data_T1$D
  
  data_T0 <- do.call(rbind, fold_T0)
  yhat_T0 <- data_T0$y.hat
  y_T0 <- data_T0$y
  pi_hat_T0 <- data_T0$pi.hat
  T0 <- data_T0$D
  
  mean(yhat_T1 + (T1 / pi_hat_T1) * (y_T1 - yhat_T1)) -
    mean(yhat_T0 + ((1 - T0) / (1 - pi_hat_T0)) * (y_T0 - yhat_T0))
}











newton_backtracking_dual <- function(
    S,
    D,
    weight_fun,
    weight_deriv,
    F_fun,
    valid_eta,
    lambda_start,
    tol = 1e-8,
    maxit ,
    min_alpha = 1e-10,
    armijo = 1e-4,
    ridge = 1e-10) {
  
  N <- nrow(S)
  q <- ncol(S)
  
  lambda <- as.numeric(lambda_start)
  
  # ------------------------------------------------------------
  # Objective
  # ------------------------------------------------------------
  
  objective <- function(lambda) {
    
    eta <- as.numeric(S %*% lambda)
    
    # domain only needed where D = 1 for the nonlinear weight map
    if (!valid_eta(eta[D == 1])) {
      return(Inf)
    }
    
    (
      sum(F_fun(eta[D == 1])) -
        sum(eta)
    ) / N
  }
  
  
  # ------------------------------------------------------------
  # Iterate Newton + backtracking
  # ------------------------------------------------------------
  
  for (iter in seq_len(maxit)) {
    
    eta <- as.numeric(
      S %*% lambda
    )
    
    if (!valid_eta(eta[D == 1])) {
      return(list(
        converged = FALSE,
        reason = "current iterate outside dual domain"
      ))
    }
    
    
    # ----------------------------------------------------------
    # weights
    # ----------------------------------------------------------
    
    w <- numeric(N)
    
    w[D == 1] <-
      weight_fun(
        eta[D == 1]
      )
    
    
    # ----------------------------------------------------------
    # Calibration score
    #
    # 1/N sum S_i (D_i w_i - 1)
    # ----------------------------------------------------------
    
    score <- colMeans(
      S *
        as.numeric(
          D * w - 1
        )
    )
    
    max_score <- max(
      abs(score)
    )
    # actual calibration-equation residuals
    raw_score <- colSums(
      S * as.numeric(D * w - 1)
    )
    
    max_raw_score <- max(abs(raw_score))
    
    # ----------------------------------------------------------
    # Stop based on the ACTUAL calibration equation
    # ----------------------------------------------------------
    
    if (is.finite(max_raw_score) &&
        max_raw_score < tol) {
      
      return(list(
        converged = TRUE,
        reason = "converged",
        lambda = lambda,
        weights_all = w,
        weights = w[D == 1],
        eta = eta,
        
        score = score,
        max_score = max_score,
        
        raw_score = raw_score,
        max_raw_score = max_raw_score,
        
        iterations = iter
      ))
    }
    # ----------------------------------------------------------
    # Hessian / Newton Jacobian
    # ----------------------------------------------------------
    
    wp <- numeric(N)
    
    wp[D == 1] <-
      weight_deriv(
        eta[D == 1]
      )
    
    
    H <- crossprod(
      S,
      S * as.numeric(D * wp)
    ) / N
    
    
    # tiny ridge only for numerical stability
    H <- H +
      diag(ridge, q)
    
    
    step <- tryCatch(
      solve(H, score),
      error = function(e) NULL
    )
    
    if (is.null(step) ||
        any(!is.finite(step))) {
      
      return(list(
        converged = FALSE,
        reason = "singular Newton Hessian",
        lambda = lambda,
        score = score,
        max_score = max_score
      ))
    }
    
    
    # ----------------------------------------------------------
    # Backtracking
    # ----------------------------------------------------------
    
    Q_old <- objective(lambda)
    
    alpha <- 1
    
    accepted <- FALSE
    
    
    while (alpha >= min_alpha) {
      
      lambda_new <-
        lambda -
        alpha * step
      
      eta_new <-
        as.numeric(
          S %*% lambda_new
        )
      
      
      # 1. Must stay inside entropy-specific dual domain
      domain_ok <-
        valid_eta(
          eta_new[D == 1]
        )
      
      
      if (domain_ok) {
        
        Q_new <-
          objective(
            lambda_new
          )
        
        
        # Armijo decrease
        decrease_ok <-
          is.finite(Q_new) &&
          Q_new <=
          Q_old -
          armijo *
          alpha *
          sum(score * step)
        
        
        if (decrease_ok) {
          
          lambda <- lambda_new
          
          accepted <- TRUE
          
          break
        }
      }
      
      
      # Shorten Newton step
      alpha <- alpha / 2
    }
    
    
    if (!accepted) {
      
      return(list(
        converged = FALSE,
        reason = "backtracking failed",
        lambda = lambda,
        score = score,
        max_score = max_score,
        iterations = iter
      ))
    }
  }
  
  
  # ------------------------------------------------------------
  # Maximum iterations reached
  # ------------------------------------------------------------
  
  eta <- as.numeric(
    S %*% lambda
  )
  
  w <- numeric(N)
  
  if (valid_eta(eta[D == 1])) {
    
    w[D == 1] <-
      weight_fun(
        eta[D == 1]
      )
    
    score <- colMeans(
      S *
        as.numeric(
          D * w - 1
        )
    )
    
    max_score <- max(abs(score))
    raw_score <- colSums(
      S * as.numeric(D * w - 1)
    )
    
    max_raw_score <- max(abs(raw_score))
    
  } else {
    
    score <- rep(NA_real_, q)
    max_score <- Inf
    raw_score <- rep(NA_real_, q)
    max_raw_score <- Inf
  }
  
  
  list(
    converged = FALSE,
    reason = "maximum iterations reached",
    lambda = lambda,
    weights_all = w,
    weights = w[D == 1],
    eta = eta,
    score = score,
    max_score = max_score,
    raw_score = raw_score,
    max_raw_score = max_raw_score,
    iterations = maxit
  )
}



solve_lambda_dual <- function(
    b_mat,
    pi_hat,
    D,
    entropy = c("SL", "EL", "ET", "HD", "CE"),
    lambda_start = NULL,
    maxit = 1000,
    
    eps_domain = 1e-8) {
  
  entropy <- match.arg(entropy)
  
  b_mat <- as.matrix(b_mat)
  pi_hat <- as.numeric(pi_hat)
  D <- as.numeric(D)
  
  N <- length(D)
  
  stopifnot(
    nrow(b_mat) == N,
    length(pi_hat) == N,
    length(D) == N
  )
  
  # Avoid exact 0/1 propensity values
  pi_hat <- pmin(
    pmax(pi_hat, 1e-8),
    1 - 1e-8
  )
  
  
  # ============================================================
  # 1. Entropy-specific debiasing covariate g(pi^{-1})
  # ============================================================
  
  g_pi <- switch(
    entropy,
    
    # Squared loss:
    # g(omega) = omega - 1
    # g(pi^{-1}) = pi^{-1} - 1
    SL = 1 / pi_hat - 1,
    
    # Empirical likelihood:
    # g(omega) = -1/omega
    # g(pi^{-1}) = -pi
    EL = -pi_hat,
    
    # Exponential tilting:
    # g(omega) = log(omega)
    # g(pi^{-1}) = log(pi^{-1})
    ET = log(1 / pi_hat),
    
    # Hellinger:
    # g(omega) = -1/(2 sqrt(omega))
    # g(pi^{-1}) = -sqrt(pi)/2
    HD = -sqrt(pi_hat) / 2,
    
    # Contrast entropy:
    # g(pi^{-1}) = log(1 - pi)
    CE = log1p(-pi_hat)
  )
  
  weight_deriv <- switch(
    entropy,
    
    SL = function(eta) {
      rep(1, length(eta))
    },
    
    EL = function(eta) {
      1 / eta^2
    },
    
    ET = function(eta) {
      exp(eta)
    },
    
    HD = function(eta) {
      -1 / (2 * eta^3)
    },
    
    CE = function(eta) {
      exp(eta) / expm1(eta)^2
    }
  )
  # ============================================================
  # 2. Calibration vector
  #
  # Include 1 explicitly so normalization is part of the
  # calibration equations:
  #
  # sum_i D_i omega_i = N
  # ============================================================
  
  S <- cbind(
    intercept = 1,
    b_mat,
    g_pi = g_pi
  )
  
  q <- ncol(S)
  
  
  # ============================================================
  # 3. Entropy-specific inverse dual map omega = g^{-1}(eta)
  # ============================================================
  
  weight_fun <- switch(
    entropy,
    
    SL = function(eta) {
      1 + eta
    },
    
    EL = function(eta) {
      -1 / eta
    },
    
    ET = function(eta) {
      exp(eta)
    },
    
    HD = function(eta) {
      1 / (4 * eta^2)
    },
    
    CE = function(eta) {
      -1 / expm1(eta)
    }
  )
  
  
  # ============================================================
  # 4. Valid dual domain
  # ============================================================
  
  valid_eta <- switch(
    entropy,
    
    SL = function(eta) {
      all(eta > -1 + eps_domain)
    },
    
    EL = function(eta) {
      all(eta < -eps_domain)
    },
    
    ET = function(eta) {
      all(is.finite(eta))
    },
    
    HD = function(eta) {
      all(eta < -eps_domain)
    },
    
    CE = function(eta) {
      all(eta < -eps_domain)
    }
  )
  
  
  # ============================================================
  # 5. Primitive F satisfying F'(eta) = omega(eta)
  # ============================================================
  
  F_fun <- switch(
    entropy,
    
    SL = function(eta) {
      eta + 0.5 * eta^2
    },
    
    EL = function(eta) {
      -log(-eta)
    },
    
    ET = function(eta) {
      exp(eta)
    },
    
    HD = function(eta) {
      -1 / (4 * eta)
    },
    
    CE = function(eta) {
      eta - log(-expm1(eta))
    }
  )
  
  
  # ============================================================
  # 6. Starting value
  # ============================================================
  
  if (is.null(lambda_start)) {
    
    lambda_start <- rep(0, q)
    
    if (entropy == "ET") {
      
      # Correct-PS motivated ET start
      lambda_start[q] <- 1
      
    } else if (entropy %in% c("EL", "HD", "CE")) {
      
      # Need eta < 0
      # Intercept term makes this convenient
      lambda_start[1] <- -1
      
    } else if (entropy == "SL") {
      
      lambda_start[1] <- 0
    }
  }
  
  # ============================================================
  # 9. Optimization
  # ============================================================
  
  fit <- newton_backtracking_dual(
    S = S,
    D = D,
    weight_fun = weight_fun,
    weight_deriv = weight_deriv,
    F_fun = F_fun,
    valid_eta = valid_eta,
    lambda_start = lambda_start,
    tol = 1e-8,
    maxit = maxit
  )
  
  # ============================================================
  # 10. Final weights
  # ============================================================
  
  # ============================================================
  # Attach information needed later
  # ============================================================
  
  fit$S <- S
  fit$g_pi <- g_pi
  fit$entropy <- entropy
  
  return(fit)
}



estimate_ATE_dual <- function(
    fold_t1,
    fold_t0,
    entropy = c("SL", "EL", "ET", "HD", "CE"),
    maxit = 1000) {
  
  entropy <- match.arg(entropy)
  
  
  # ============================================================
  # Reconstruct full cross-fitted datasets
  # ============================================================
  
  dat1 <- do.call(
    rbind,
    fold_t1
  )
  
  dat0 <- do.call(
    rbind,
    fold_t0
  )
  
  
  # ============================================================
  # Treated arm
  #
  # D1 = I(T = 1)
  # pi1 = P(T = 1 | X)
  # b1 = cross-fitted prediction of Y(1)
  # ============================================================
  
  D1 <- as.numeric(
    dat1$D == 1
  )
  
  fit1 <- solve_lambda_dual(
    b_mat = matrix(
      dat1$y.hat,
      ncol = 1
    ),
    pi_hat = dat1$pi.hat,
    D = D1,
    entropy = entropy,
    maxit = maxit
  )
  
  
  # ============================================================
  # Control arm
  #
  # D0 = I(T = 0)
  # pi0 = P(T = 0 | X) = 1 - pi1
  # b0 = cross-fitted prediction of Y(0)
  # ============================================================
  
  D0 <- as.numeric(
    dat0$D == 0
  )
  
  fit0 <- solve_lambda_dual(
    b_mat = matrix(
      dat0$y.hat,
      ncol = 1
    ),
    pi_hat = 1 - dat0$pi.hat,
    D = D0,
    entropy = entropy,
    maxit = maxit
  )
  
  
  
  # HARD FAILURE
  if (is.null(fit1$weights_all) ||
      is.null(fit0$weights_all) ||
      any(!is.finite(fit1$weights_all)) ||
      any(!is.finite(fit0$weights_all))) {
    
    return(
      list(
        ATE = NA_real_,
        theta1 = NA_real_,
        theta0 = NA_real_,
        fit1 = fit1,
        fit0 = fit0,
        hard_failure = TRUE
      )
    )
  }
  
  # SOFT FAILURE / tolerance warning
  if (!isTRUE(fit1$converged) ||
      !isTRUE(fit0$converged)) {
    
    warning(
      paste0(
        entropy,
        " calibration warning: ",
        "raw_score1 = ", fit1$max_raw_score,
        ", raw_score0 = ", fit0$max_raw_score,
        ", reason1 = ", fit1$reason,
        ", reason0 = ", fit0$reason
      )
    )
  }
  # ============================================================
  # Closed-form ATE estimator
  #
  # theta_t =
  # 1/N sum_i D_ti omega_ti Y_i
  #
  # No theta iteration.
  # ============================================================
  
  N1 <- nrow(dat1)
  N0 <- nrow(dat0)
  
  theta1_hat <- sum(
    D1 *
      fit1$weights_all *
      dat1$y
  ) / N1
  
  
  theta0_hat <- sum(
    D0 *
      fit0$weights_all *
      dat0$y
  ) / N0
  
  
  ATE_hat <-
    theta1_hat -
    theta0_hat
  
  
  # ============================================================
  # Diagnostics
  # ============================================================
  
  list(
    ATE = ATE_hat,
    
    theta1 = theta1_hat,
    theta0 = theta0_hat,
    
    fit1 = fit1,
    fit0 = fit0,
    
    normalization1 =
      sum(D1 * fit1$weights_all) / N1,
    
    normalization0 =
      sum(D0 * fit0$weights_all) / N0,
    
    balance_yhat1 =
      sum(D1 *
            fit1$weights_all *
            dat1$y.hat) -
      sum(dat1$y.hat),
    
    balance_yhat0 =
      sum(D0 *
            fit0$weights_all *
            dat0$y.hat) -
      sum(dat0$y.hat),
    
    max_score1 =
      fit1$max_score,
    
    max_score0 =
      fit0$max_score
  )
}


# ---------------------------------------------------
# 4. Benchmark estimators
# ---------------------------------------------------

# Function: estimate_ipw
# ----------------------
# Standard normalized inverse probability weighted estimator of the ATE.
estimate_ipw <- function(dat) {
  sum(dat$D * dat$y / dat$pi.hat) / sum(dat$D / dat$pi.hat) -
    sum((1 - dat$D) * dat$y / (1 - dat$pi.hat)) / sum((1 - dat$D) / (1 - dat$pi.hat))
}

# Function: estimate_cbps_pair
# ----------------------------
# Computes two CBPS-based benchmark estimators:
#   - oCBPS
#   - CBPS
# and returns both in a named vector.
estimate_cbps_pair <- function(dat) {
  xvars <- grep("^x\\d+$", names(dat), value = TRUE)
  X2 <- as.matrix(cbind(1, dat[, xvars, drop = FALSE]))
  
  ocbps_model <- CBPS::CBPS(dat$D ~ X2, ATT = 0, method = "exact",
                            baseline.formula = ~ X2, diff.formula = ~ X2)
  ocbps <- coef(lm(dat$y ~ dat$D, weights = ocbps_model$weights))["dat$D"]
  
  cbps_model <- CBPS::CBPS(dat$D ~ X2, ATT = 0, method = "exact")
  cbps <- coef(lm(dat$y ~ dat$D, weights = cbps_model$weights))["dat$D"]
  
  c(oCBPS = unname(ocbps), CBPS = unname(cbps))
}

# Function: estimate_ebal
# -----------------------
# Computes the entropy balancing benchmark using WeightIt and then forms the
# weighted mean difference in outcomes between treated and control groups.
estimate_ebal <- function(dat) {
  xvars <- grep("^x\\d+$", names(dat), value = TRUE)
  form <- as.formula(paste("D ~", paste(xvars, collapse = " + ")))
  ebal_fit <- WeightIt::weightit(form, data = dat, method = "ebal", estimand = "ATE")
  w <- ebal_fit$weights
  sum(w * dat$D * dat$y) / sum(w * dat$D) -
    sum(w * (1 - dat$D) * dat$y) / sum(w * (1 - dat$D))
}

# Function: estimate_ebcw
# -----------------------
# Computes the EBCW benchmark through the ATE package.
estimate_ebcw <- function(dat) {
  xvars <- grep("^x\\d+$", names(dat), value = TRUE)
  X_chan <- data.frame(dat[, xvars, drop = FALSE], log_pi_hat = log(dat$pi.hat))
  fit <- ATE::ATE(dat$y, dat$D, X_chan, ATT = FALSE)
  summary(fit)$Estimate[3, 1]
}

# ---------------------------------------------------
# 5. One-replication wrapper
# ---------------------------------------------------


psi_bar_ate <- function(
    beta,
    y,
    T,
    X,
    yhat1,
    yhat0,
    entropy
) {
  
  N <- length(y)
  
  X <- as.matrix(X)
  
  r <- ncol(X)
  
  # Our S has:
  # intercept + yhat + debiasing variable
  q1 <- 3
  q0 <- 3
  
  
  ##########################################################################
  # 1. Unpack beta
  #
  # beta = (phi, lambda1, theta1, lambda0, theta0)
  ##########################################################################
  
  idx_phi <- seq_len(r)
  
  idx_lam1 <- (r + 1):(r + q1)
  
  idx_theta1 <- r + q1 + 1
  
  idx_lam0 <- (idx_theta1 + 1):(idx_theta1 + q0)
  
  idx_theta0 <- max(idx_lam0) + 1
  
  
  phi <- beta[idx_phi]
  
  lambda1 <- beta[idx_lam1]
  
  theta1 <- beta[idx_theta1]
  
  lambda0 <- beta[idx_lam0]
  
  theta0 <- beta[idx_theta0]
  
  
  ##########################################################################
  # 2. Recompute propensity scores from candidate phi
  ##########################################################################
  
  pi <- plogis(
    as.vector(X %*% phi)
  )
  
  pi <- pmin(
    pmax(pi, 1e-8),
    1 - 1e-8
  )
  
  p1 <- pi
  p0 <- 1 - pi
  
  
  ##########################################################################
  # 3. Entropy-specific debiasing function
  ##########################################################################
  
  q_fun <- function(p) {
    
    if (entropy == "SL")
      return(1 / p - 1)
    
    if (entropy == "EL")
      return(-p)
    
    if (entropy == "ET")
      return(log(1 / p))
    
    if (entropy == "HD")
      return(-sqrt(p) / 2)
    
    if (entropy == "CE")
      return(log1p(-p))
  }
  
  
  ##########################################################################
  # 4. Entropy-specific weight function
  ##########################################################################
  
  weight_fun <- function(eta) {
    
    if (entropy == "SL")
      return(1 + eta)
    
    if (entropy == "EL")
      return(-1 / eta)
    
    if (entropy == "ET")
      return(exp(eta))
    
    if (entropy == "HD")
      return(1 / (4 * eta^2))
    
    if (entropy == "CE")
      return(-1 / expm1(eta))
  }
  
  
  ##########################################################################
  # 5. Build S1 and S0 at candidate phi
  ##########################################################################
  
  S1 <- cbind(
    intercept = 1,
    yhat = yhat1,
    g_pi = q_fun(p1)
  )
  
  S0 <- cbind(
    intercept = 1,
    yhat = yhat0,
    g_pi = q_fun(p0)
  )
  
  
  ##########################################################################
  # 6. Candidate weights
  ##########################################################################
  
  eta1 <- as.vector(
    S1 %*% lambda1
  )
  
  eta0 <- as.vector(
    S0 %*% lambda0
  )
  
  
  # Domain protection for numerical differentiation
  if (entropy == "SL") {
    
    if (any(eta1[D1 == 1] <= -1) ||
        any(eta0[D0 == 1] <= -1)) {
      
      return(
        rep(NA_real_, length(beta))
      )
    }
  }
  
  if (entropy %in% c("EL", "HD", "CE")) {
    
    if (any(eta1[D1 == 1] >= 0) ||
        any(eta0[D0 == 1] >= 0)) {
      
      return(
        rep(NA_real_, length(beta))
      )
    }
  }
  
  
  w1 <- weight_fun(eta1)
  w0 <- weight_fun(eta0)
  
  D1 <- T
  D0 <- 1 - T
  
  
  ##########################################################################
  # 7. Eq. (27) propensity-score estimating function
  #
  # For logistic PS:
  #
  # ((T/pi)-1) h(phi)
  # with h(phi)=pi X
  #
  # exactly equals X(T-pi).
  ##########################################################################
  
  Psi_phi <- X *
    as.numeric(T - pi)
  
  
  ##########################################################################
  # 8. Calibration estimating equations
  ##########################################################################
  
  Psi_lam1 <- S1 *
    as.numeric(D1 * w1 - 1)
  
  Psi_lam0 <- S0 *
    as.numeric(D0 * w0 - 1)
  
  
  ##########################################################################
  # 9. Outcome estimating equations
  ##########################################################################
  
  Psi_theta1 <-
    D1 *
    w1 *
    (y - theta1)
  
  Psi_theta0 <-
    D0 *
    w0 *
    (y - theta0)
  
  
  ##########################################################################
  # 10. Return mean joint estimating equation
  ##########################################################################
  
  Psi <- cbind(
    Psi_phi,
    Psi_lam1,
    theta1 = Psi_theta1,
    Psi_lam0,
    theta0 = Psi_theta0
  )
  
  colMeans(Psi)
}





gec_ate_inference <- function(
    dat,
    yhat1,
    yhat0,
    fit1,
    fit0,
    ps_fit,
    entropy,
    true_ATE = NA_real_,
    numerical_jacobian = TRUE
) {
  
  N <- nrow(dat)
  
  #####################################################################
  # Align objects
  #####################################################################
  
  X_ps <- model.matrix(ps_fit)
  
  phi_hat <- coef(ps_fit)
  
  lambda1_hat <- fit1$lambda
  lambda0_hat <- fit0$lambda
  
  
  #####################################################################
  # Require successful calibration
  #####################################################################
  
  if (!isTRUE(fit1$converged) ||
      !isTRUE(fit0$converged)) {
    
    return(
      list(
        success = FALSE,
        entropy = entropy
      )
    )
  }
  
  
  #####################################################################
  # 1. Analytic sandwich
  #####################################################################
  
  sw <- tryCatch(
    
    gec_ate_sandwich(
      y = dat$y,
      T = dat$D,
      X = X_ps,
      yhat1 = yhat1,
      yhat0 = yhat0,
      lambda1 = lambda1_hat,
      lambda0 = lambda0_hat,
      entropy = entropy,
      phi_hat = phi_hat
    ),
    
    error = function(e) NULL
  )
  
  
  if (is.null(sw) ||
      !is.finite(sw$ATE) ||
      !is.finite(sw$SE)) {
    
    return(
      list(
        success = FALSE,
        entropy = entropy
      )
    )
  }
  
  
  #####################################################################
  # Point estimate
  #####################################################################
  
  ate_hat <- sw$ATE
  
  se_analytic <- sw$SE
  var_analytic <- sw$variance
  
  
  #####################################################################
  # Analytic CI and coverage
  #####################################################################
  
  lower_analytic <- ate_hat -
    qnorm(0.975) * se_analytic
  
  upper_analytic <- ate_hat +
    qnorm(0.975) * se_analytic
  
  
  if (is.finite(true_ATE)) {
    
    cover_analytic <- as.integer(
      lower_analytic <= true_ATE &&
        true_ATE <= upper_analytic
    )
    
  } else {
    
    cover_analytic <- NA_integer_
  }
  
  
  #####################################################################
  # Defaults for numerical Jacobian
  #####################################################################
  
  se_numeric <- NA_real_
  var_numeric <- NA_real_
  
  lower_numeric <- NA_real_
  upper_numeric <- NA_real_
  
  cover_numeric <- NA_integer_
  
  jacobian_max_diff <- NA_real_
  
  
  #####################################################################
  # 2. Numerical Jacobian sandwich
  #####################################################################
  
  if (numerical_jacobian) {
    
    beta_hat <- c(
      phi_hat,
      lambda1_hat,
      sw$theta1,
      lambda0_hat,
      sw$theta0
    )
    
    
    J_numeric <- tryCatch(
      
      numDeriv::jacobian(
        func = function(beta) {
          
          psi_bar_ate(
            beta = beta,
            y = dat$y,
            T = dat$D,
            X = X_ps,
            yhat1 = yhat1,
            yhat0 = yhat0,
            entropy = entropy
          )
          
        },
        
        x = beta_hat,
        
        method = "Richardson"
      ),
      
      error = function(e) NULL
    )
    
    
    if (!is.null(J_numeric) &&
        all(is.finite(J_numeric))) {
      
      # Bread = minus Jacobian
      A_numeric <- -J_numeric
      
      # Same meat for analytic and numerical versions
      B_hat <- sw$B
      
      
      Ainv_numeric <- tryCatch(
        solve(A_numeric),
        error = function(e) NULL
      )
      
      
      if (!is.null(Ainv_numeric)) {
        
        V_numeric <-
          (
            Ainv_numeric %*%
              B_hat %*%
              t(Ainv_numeric)
          ) / N
        
        
        #################################################################
        # ATE contrast theta1 - theta0
        #################################################################
        
        par_names <- c(
          paste0("phi_", colnames(X_ps)),
          paste0("lambda1_", colnames(sw$S1)),
          "theta1",
          paste0("lambda0_", colnames(sw$S0)),
          "theta0"
        )
        
        contrast <- rep(
          0,
          length(beta_hat)
        )
        
        contrast[
          which(par_names == "theta1")
        ] <- 1
        
        contrast[
          which(par_names == "theta0")
        ] <- -1
        
        
        var_numeric <- as.numeric(
          t(contrast) %*%
            V_numeric %*%
            contrast
        )
        
        
        if (is.finite(var_numeric) &&
            var_numeric >= 0) {
          
          se_numeric <- sqrt(
            var_numeric
          )
          
          
          lower_numeric <-
            ate_hat -
            qnorm(0.975) * se_numeric
          
          upper_numeric <-
            ate_hat +
            qnorm(0.975) * se_numeric
          
          
          if (is.finite(true_ATE)) {
            
            cover_numeric <- as.integer(
              lower_numeric <= true_ATE &&
                true_ATE <= upper_numeric
            )
          }
        }
        
        
        #################################################################
        # Analytic-vs-numerical bread check
        #################################################################
        
        if (all(dim(sw$A) ==
                dim(A_numeric))) {
          
          jacobian_max_diff <-
            max(
              abs(
                sw$A -
                  A_numeric
              )
            )
        }
      }
    }
  }
  
  
  #####################################################################
  # Return one row's worth of information
  #####################################################################
  
  list(
    
    success = TRUE,
    
    entropy = entropy,
    
    ATE = ate_hat,
    
    SE_analytic = se_analytic,
    SE_numeric = se_numeric,
    
    Var_analytic = var_analytic,
    Var_numeric = var_numeric,
    
    lower_analytic = lower_analytic,
    upper_analytic = upper_analytic,
    
    lower_numeric = lower_numeric,
    upper_numeric = upper_numeric,
    
    cover_analytic = cover_analytic,
    cover_numeric = cover_numeric,
    
    max_weight1 = sw$max_weight1,
    max_weight0 = sw$max_weight0,
    
    ESS1 = sw$ESS1,
    ESS0 = sw$ESS0,
    
    jacobian_max_diff = jacobian_max_diff,
    
    calibration1 = sw$calibration_raw1,
    calibration0 = sw$calibration_raw0
  )
}



inference_to_row <- function(x, replication) {
  
  if (!isTRUE(x$success)) {
    
    return(
      data.frame(
        replication = replication,
        entropy = x$entropy,
        success = FALSE,
        
        ATE = NA_real_,
        SE_analytic = NA_real_,
        SE_numeric = NA_real_,
        
        Var_analytic = NA_real_,
        Var_numeric = NA_real_,
        
        cover_analytic = NA_real_,
        cover_numeric = NA_real_,
        
        max_weight1 = NA_real_,
        max_weight0 = NA_real_,
        
        ESS1 = NA_real_,
        ESS0 = NA_real_,
        
        jacobian_max_diff = NA_real_,
        
        calibration1 = NA_real_,
        calibration0 = NA_real_
      )
    )
  }
  
  
  data.frame(
    
    replication = replication,
    
    entropy = x$entropy,
    
    success = TRUE,
    
    ATE = x$ATE,
    
    SE_analytic = x$SE_analytic,
    
    SE_numeric = x$SE_numeric,
    
    Var_analytic = x$Var_analytic,
    
    Var_numeric = x$Var_numeric,
    
    cover_analytic = x$cover_analytic,
    
    cover_numeric = x$cover_numeric,
    
    max_weight1 = x$max_weight1,
    
    max_weight0 = x$max_weight0,
    
    ESS1 = x$ESS1,
    
    ESS0 = x$ESS0,
    
    jacobian_max_diff =
      x$jacobian_max_diff,
    
    calibration1 =
      x$calibration1,
    
    calibration0 =
      x$calibration0
  )
}







# ---------------------------------------------------
# 5. One-replication wrapper
#    USING YOUR ORIGINAL FOLD STRUCTURE
# ---------------------------------------------------

# Function: run_one_replication
# -----------------------------
# Runs a single Monte Carlo replication:
#   1. generate one dataset
#   2. compute benchmark estimators
#   3. compute cross-fitted AIPW estimators
#   4. compute the proposed estimators
#
# Returns one row of estimates.
run_one_replication_old <- function(
    n = 1000,
    p = 4,
    K = 4,
    rep_id = 1,
    outcome_model = 1,
    ps_model = 1,
    run_lm = TRUE,
    run_gam = TRUE,
    gec_from = c("gam", "lm"),
    numerical_jacobian = TRUE,oracle_m=FALSE
) {
  
  gec_from <- match.arg(gec_from)
  
  
  # ==============================================================
  # 1. Generate data
  # ==============================================================
  
  sim <- generate_ate_data(
    n = n,
    p = p,
    outcome_model = outcome_model,
    ps_model = ps_model,
    seed = rep_id + 20242022
  )
  
  dat <- sim$data
  
  # Make sure ID exists
  if (!"ID" %in% names(dat)) {
    dat$ID <- seq_len(nrow(dat))
  }
  
  D <- dat$D
  
  
  # ==============================================================
  # 2. True ATE
  #
  # OR1: true ATE = 1
  # OR2: true ATE = 10
  # ==============================================================
  
  true_ATE <- if (outcome_model == 1) {
    1
  } else {
    10
  }
  
  
  # ==============================================================
  # 3. Existing estimators
  # ==============================================================
  
  estimates <- c()
  
  estimates["IPW"] <- estimate_ipw(dat)
  
  estimates["EBPS"] <- estimate_ebal(dat)
  
  estimates[c("oCBPS", "CBPS")] <-
    estimate_cbps_pair(dat)
  
  estimates["EBCW"] <- estimate_ebcw(dat)
  
  
  # ==============================================================
  # 4. Cross-fitting objects
  # ==============================================================
  
  fold_T1_lm <- NULL
  fold_T0_lm <- NULL
  
  fold_T1_gam <- NULL
  fold_T0_gam <- NULL
  
  
  # --------------------------------------------------------------
  # LM
  # --------------------------------------------------------------
  
  if (run_lm) {
    
    fold_T1_lm <- k_fold_function_lm(
      df = dat,
      K = K,
      idx_S = which(D == 1),
      idx_U0 = which(D == 0),
      seed = rep_id + 20242022
    )
    
    fold_T0_lm <- k_fold_function_lm(
      df = dat,
      K = K,
      idx_S = which(D == 0),
      idx_U0 = which(D == 1),
      seed = rep_id + 20242022
    )
    
    estimates["AIPW_LM"] <-
      compute_crossfit_aipw(
        fold_T1_lm,
        fold_T0_lm
      )
  }
  
  
  # --------------------------------------------------------------
  # GAM
  # --------------------------------------------------------------
  
  if (run_gam) {
    
    fold_T1_gam <- k_fold_function_gam(
      df = dat,
      K = K,
      idx_S = which(D == 1),
      idx_U0 = which(D == 0),
      seed = rep_id + 20242022
    )
    
    fold_T0_gam <- k_fold_function_gam(
      df = dat,
      K = K,
      idx_S = which(D == 0),
      idx_U0 = which(D == 1),
      seed = rep_id + 20242022
    )
    
    estimates["AIPW_GAM"] <-
      compute_crossfit_aipw(
        fold_T1_gam,
        fold_T0_gam
      )
  }
  
  
  
  
  # ==============================================================
  # 5. Choose calibration functions for GEC
  #
  # oracle_m = FALSE:
  #     use the existing cross-fitted LM/GAM predictions
  #
  # oracle_m = TRUE:
  #     use the TRUE m1(X), m0(X)
  #
  # This oracle option is only a variance diagnostic.
  # ==============================================================
  
  if (oracle_m) {
    
    # ------------------------------------------------------------
    # TRUE conditional outcome functions
    # ------------------------------------------------------------
    
    if (outcome_model != 2) {
      stop("oracle_m diagnostic currently set up for OR2 only.")
    }
    
    oracle <- true_m_or2(dat)
    
    # Make two complete copies because estimate_ATE_dual()
    # expects separate treatment-arm calibration objects.
    
    dat1_oracle <- dat
    dat0_oracle <- dat
    
    dat1_oracle$y.hat <- oracle$m1
    dat0_oracle$y.hat <- oracle$m0
    
    # estimate_ATE_dual() expects lists of data frames.
    # No cross-fitting is necessary because m_t(X) is known exactly.
    fold_t1_gec <- list(dat1_oracle)
    fold_t0_gec <- list(dat0_oracle)
    
  } else {
    
    # ============================================================
    # EXISTING estimated calibration-function pathway
    # ============================================================
    
    if (gec_from == "gam") {
      
      if (is.null(fold_T1_gam) ||
          is.null(fold_T0_gam)) {
        
        fold_T1_gam <- k_fold_function_gam(
          df = dat,
          K = K,
          idx_S = which(D == 1),
          idx_U0 = which(D == 0),
          seed = rep_id + 20242022
        )
        
        fold_T0_gam <- k_fold_function_gam(
          df = dat,
          K = K,
          idx_S = which(D == 0),
          idx_U0 = which(D == 1),
          seed = rep_id + 20242022
        )
      }
      
      fold_t1_gec <- fold_T1_gam
      fold_t0_gec <- fold_T0_gam
      
    } else {
      
      if (is.null(fold_T1_lm) ||
          is.null(fold_T0_lm)) {
        
        fold_T1_lm <- k_fold_function_lm(
          df = dat,
          K = K,
          idx_S = which(D == 1),
          idx_U0 = which(D == 0),
          seed = rep_id + 20242022
        )
        
        fold_T0_lm <- k_fold_function_lm(
          df = dat,
          K = K,
          idx_S = which(D == 0),
          idx_U0 = which(D == 1),
          seed = rep_id + 20242022
        )
      }
      
      fold_t1_gec <- fold_T1_lm
      fold_t0_gec <- fold_T0_lm
    }
  }
  
  # ==============================================================
  # 6. Reconstruct full cross-fitted datasets
  #
  # IMPORTANT for sandwich:
  # align treated/control predictions by ID.
  # ==============================================================
  
  dat1_gec <- do.call(
    rbind,
    fold_t1_gec
  )
  
  dat0_gec <- do.call(
    rbind,
    fold_t0_gec
  )
  
  dat1_gec <- dat1_gec[
    order(dat1_gec$ID),
  ]
  
  dat0_gec <- dat0_gec[
    order(dat0_gec$ID),
  ]
  
  dat <- dat[
    order(dat$ID),
  ]
  
  
  stopifnot(
    identical(
      as.integer(dat1_gec$ID),
      as.integer(dat$ID)
    ),
    identical(
      as.integer(dat0_gec$ID),
      as.integer(dat$ID)
    )
  )
  
  
  yhat1 <- as.numeric(
    dat1_gec$y.hat
  )
  
  yhat0 <- as.numeric(
    dat0_gec$y.hat
  )
  
  
  # ==============================================================
  # 7. Fit the SAME working propensity model used for pi.hat
  #
  # Needed because sandwich requires phi_hat and model matrix X.
  # ==============================================================
  
  xvars <- setdiff(
    names(dat),
    c(
      "D",
      "y",
      "ID",
      "pi.hat",
      "y.hat"
    )
  )
  
  ps_formula <- as.formula(
    paste(
      "D ~",
      paste(
        xvars,
        collapse = " + "
      )
    )
  )
  
  ps_fit <- glm(
    ps_formula,
    data = dat,
    family = binomial()
  )
  
  phi_hat <- coef(ps_fit)
  
  X_ps <- model.matrix(
    ps_fit
  )
  
  
  # Optional consistency check:
  # this should be tiny if generate_ate_data() used the same PS fit.
  if ("pi.hat" %in% names(dat)) {
    
    pi_difference <- max(
      abs(
        fitted(ps_fit) -
          dat$pi.hat
      )
    )
    
  } else {
    
    pi_difference <- NA_real_
  }
  
  
  # ==============================================================
  # 8. Container for GEC inference results
  # ==============================================================
  
  gec_output <- c()
  
  
  # ==============================================================
  # 9. Proposed GEC estimators
  #
  # Kim: SL / EL / ET / HD / CE
  # ==============================================================
  
  for (ent in c(
    "ET",
    "HD",
    "CE"
  )) {
    
    
    # ------------------------------------------------------------
    # 9a. Point estimation
    # ------------------------------------------------------------
    
    fit_gec <- tryCatch(
      
      estimate_ATE_dual(
        fold_t1 = fold_t1_gec,
        fold_t0 = fold_t0_gec,
        entropy = ent
      ),
      
      error = function(e) NULL
    )
    
    
    # ------------------------------------------------------------
    # Hard failure
    # ------------------------------------------------------------
    
    if (is.null(fit_gec) ||
        !is.finite(fit_gec$ATE) ||
        is.null(fit_gec$fit1) ||
        is.null(fit_gec$fit0)) {
      
      estimates[ent] <- NA_real_
      
      gec_output[
        paste0(ent, "_SE_ANA")
      ] <- NA_real_
      
      gec_output[
        paste0(ent, "_SE_NUM")
      ] <- NA_real_
      
      gec_output[
        paste0(ent, "_COV_ANA")
      ] <- NA_real_
      
      gec_output[
        paste0(ent, "_COV_NUM")
      ] <- NA_real_
      
      gec_output[
        paste0(ent, "_ESS1")
      ] <- NA_real_
      
      gec_output[
        paste0(ent, "_ESS0")
      ] <- NA_real_
      
      gec_output[
        paste0(ent, "_MAXW1")
      ] <- NA_real_
      
      gec_output[
        paste0(ent, "_MAXW0")
      ] <- NA_real_
      
      gec_output[
        paste0(ent, "_JACDIFF")
      ] <- NA_real_
      
      gec_output[
        paste0(ent, "_SUCCESS")
      ] <- 0
      
      next
    }
    
    
    estimates[ent] <-
      fit_gec$ATE
    
    
    # ------------------------------------------------------------
    # 9b. Sandwich inference
    # ------------------------------------------------------------
    
    inf_gec <- tryCatch(
      
      gec_ate_inference(
        dat = dat,
        yhat1 = yhat1,
        yhat0 = yhat0,
        fit1 = fit_gec$fit1,
        fit0 = fit_gec$fit0,
        ps_fit = ps_fit,
        entropy = ent,
        true_ATE = true_ATE,
        numerical_jacobian =
          numerical_jacobian
      ),
      
      error = function(e) NULL
    )
    
    
    # ------------------------------------------------------------
    # Sandwich failure
    # ------------------------------------------------------------
    
    if (is.null(inf_gec) ||
        !isTRUE(inf_gec$success)) {
      
      gec_output[
        paste0(ent, "_SE_ANA")
      ] <- NA_real_
      
      gec_output[
        paste0(ent, "_SE_NUM")
      ] <- NA_real_
      
      gec_output[
        paste0(ent, "_COV_ANA")
      ] <- NA_real_
      
      gec_output[
        paste0(ent, "_COV_NUM")
      ] <- NA_real_
      
      gec_output[
        paste0(ent, "_ESS1")
      ] <- NA_real_
      
      gec_output[
        paste0(ent, "_ESS0")
      ] <- NA_real_
      
      gec_output[
        paste0(ent, "_MAXW1")
      ] <- NA_real_
      
      gec_output[
        paste0(ent, "_MAXW0")
      ] <- NA_real_
      
      gec_output[
        paste0(ent, "_JACDIFF")
      ] <- NA_real_
      
      gec_output[
        paste0(ent, "_SUCCESS")
      ] <- 0
      
      next
    }
    
    
    # ------------------------------------------------------------
    # 9c. Store analytic SE
    # ------------------------------------------------------------
    
    gec_output[
      paste0(ent, "_SE_ANA")
    ] <- inf_gec$SE_analytic
    
    
    # ------------------------------------------------------------
    # 9d. Store numerical-Jacobian SE
    # ------------------------------------------------------------
    
    gec_output[
      paste0(ent, "_SE_NUM")
    ] <- inf_gec$SE_numeric
    
    
    # ------------------------------------------------------------
    # 9e. Coverage indicators
    #
    # Each replication gives 0 or 1.
    # Monte Carlo mean later = coverage probability.
    # ------------------------------------------------------------
    
    gec_output[
      paste0(ent, "_COV_ANA")
    ] <- inf_gec$cover_analytic
    
    gec_output[
      paste0(ent, "_COV_NUM")
    ] <- inf_gec$cover_numeric
    
    
    # ------------------------------------------------------------
    # 9f. Effective sample sizes
    # ------------------------------------------------------------
    
    gec_output[
      paste0(ent, "_ESS1")
    ] <- inf_gec$ESS1
    
    gec_output[
      paste0(ent, "_ESS0")
    ] <- inf_gec$ESS0
    
    
    # ------------------------------------------------------------
    # 9g. Maximum calibration weights
    # ------------------------------------------------------------
    
    gec_output[
      paste0(ent, "_MAXW1")
    ] <- inf_gec$max_weight1
    
    gec_output[
      paste0(ent, "_MAXW0")
    ] <- inf_gec$max_weight0
    
    
    # ------------------------------------------------------------
    # 9h. Analytic-vs-numerical Jacobian diagnostic
    # ------------------------------------------------------------
    
    gec_output[
      paste0(ent, "_JACDIFF")
    ] <- inf_gec$jacobian_max_diff
    
    
    # ------------------------------------------------------------
    # 9i. Success indicator
    # ------------------------------------------------------------
    
    gec_output[
      paste0(ent, "_SUCCESS")
    ] <- 1
  }
  
  
  # ==============================================================
  # 10. Additional replication information
  # ==============================================================
  
  misc_output <- c(
    TRUE_ATE = true_ATE,
    N_TREATED = sum(dat$D == 1),
    N_CONTROL = sum(dat$D == 0),
    PI_FIT_DIFF = pi_difference
  )
  
  
  # ==============================================================
  # 11. Return ONE wide row
  # ==============================================================
  
  output <- c(
    estimates,
    gec_output,
    misc_output
  )
  
  as.data.frame(
    as.list(output),
    check.names = FALSE
  )
}

run_one_replication <- function(
    n = 1000,
    p = 4,
    K = 4,
    rep_id = 1,
    outcome_model = 1,
    ps_model = 1,
    run_lm = TRUE,
    run_gam = TRUE,
    gec_from = c("gam", "lm"),
    numerical_jacobian = TRUE,
    oracle_m = FALSE,
    run_gec = FALSE          # <---- NEW SWITCH
) {
  
  gec_from <- match.arg(gec_from)
  
  
  # ==============================================================
  # 1. Generate data
  # ==============================================================
  
  sim <- generate_ate_data(
    n = n,
    p = p,
    outcome_model = outcome_model,
    ps_model = ps_model,
    seed = rep_id + 20242022
  )
  
  dat <- sim$data
  
  if (!"ID" %in% names(dat)) {
    dat$ID <- seq_len(nrow(dat))
  }
  
  D <- dat$D
  
  
  # ==============================================================
  # 2. True ATE
  # ==============================================================
  
  true_ATE <- if (outcome_model == 1) {
    1
  } else {
    10
  }
  
  
  # ==============================================================
  # 3. Containers
  # ==============================================================
  
  estimates <- c()
  
  benchmark_output <- c()
  
  gec_output <- c()
  
  
  # ==============================================================
  # 4. IPW inference
  # ==============================================================
  
  ipw_inf <- tryCatch(
    estimate_ipw_inference(
      dat = dat,
      true_ATE = true_ATE
    ),
    error = function(e) NULL
  )
  
  if (!is.null(ipw_inf)) {
    
    estimates["IPW"] <-
      ipw_inf$ATE
    
    benchmark_output["IPW_SE_ANA"] <-
      ipw_inf$SE
    
    benchmark_output["IPW_COV_ANA"] <-
      ipw_inf$coverage
    
    benchmark_output["IPW_ESS1"] <-
      ipw_inf$ESS1
    
    benchmark_output["IPW_ESS0"] <-
      ipw_inf$ESS0
    
    benchmark_output["IPW_MAXW1"] <-
      ipw_inf$max_weight1
    
    benchmark_output["IPW_MAXW0"] <-
      ipw_inf$max_weight0
    
    benchmark_output["IPW_SUCCESS"] <-
      as.integer(ipw_inf$success)
    
  } else {
    
    estimates["IPW"] <- NA_real_
    
    benchmark_output["IPW_SE_ANA"] <- NA_real_
    benchmark_output["IPW_COV_ANA"] <- NA_real_
    benchmark_output["IPW_ESS1"] <- NA_real_
    benchmark_output["IPW_ESS0"] <- NA_real_
    benchmark_output["IPW_MAXW1"] <- NA_real_
    benchmark_output["IPW_MAXW0"] <- NA_real_
    benchmark_output["IPW_SUCCESS"] <- 0
  }
  
  
  # ==============================================================
  # 5. EBPS
  # ==============================================================
  
  ebps_inf <- tryCatch(
    estimate_ebps_inference(
      dat = dat,
      true_ATE = true_ATE
    ),
    error = function(e) NULL
  )
  
  if (!is.null(ebps_inf)) {
    
    estimates["EBPS"] <-
      ebps_inf$ATE
    
    benchmark_output["EBPS_SE_ANA"] <-
      ebps_inf$SE
    
    benchmark_output["EBPS_COV_ANA"] <-
      ebps_inf$coverage
    
    benchmark_output["EBPS_ESS1"] <-
      ebps_inf$ESS1
    
    benchmark_output["EBPS_ESS0"] <-
      ebps_inf$ESS0
    
    benchmark_output["EBPS_MAXW1"] <-
      ebps_inf$max_weight1
    
    benchmark_output["EBPS_MAXW0"] <-
      ebps_inf$max_weight0
    
    benchmark_output["EBPS_SUCCESS"] <-
      as.integer(ebps_inf$success)
    
  } else {
    
    estimates["EBPS"] <- NA_real_
    
    benchmark_output["EBPS_SE_ANA"] <- NA_real_
    benchmark_output["EBPS_COV_ANA"] <- NA_real_
    benchmark_output["EBPS_ESS1"] <- NA_real_
    benchmark_output["EBPS_ESS0"] <- NA_real_
    benchmark_output["EBPS_MAXW1"] <- NA_real_
    benchmark_output["EBPS_MAXW0"] <- NA_real_
    benchmark_output["EBPS_SUCCESS"] <- 0
  }
  
  
  # ==============================================================
  # 6. oCBPS
  # ==============================================================
  
  ocbps_inf <- tryCatch(
    estimate_ocbps_inference(
      dat = dat,
      true_ATE = true_ATE
    ),
    error = function(e) NULL
  )
  
  if (!is.null(ocbps_inf)) {
    
    estimates["oCBPS"] <-
      ocbps_inf$ATE
    
    benchmark_output["oCBPS_SE_ANA"] <-
      ocbps_inf$SE
    
    benchmark_output["oCBPS_COV_ANA"] <-
      ocbps_inf$coverage
    
    benchmark_output["oCBPS_ESS1"] <-
      ocbps_inf$ESS1
    
    benchmark_output["oCBPS_ESS0"] <-
      ocbps_inf$ESS0
    
    benchmark_output["oCBPS_MAXW1"] <-
      ocbps_inf$max_weight1
    
    benchmark_output["oCBPS_MAXW0"] <-
      ocbps_inf$max_weight0
    
    benchmark_output["oCBPS_SUCCESS"] <-
      as.integer(ocbps_inf$success)
    
  } else {
    
    estimates["oCBPS"] <- NA_real_
    
    benchmark_output["oCBPS_SE_ANA"] <- NA_real_
    benchmark_output["oCBPS_COV_ANA"] <- NA_real_
    benchmark_output["oCBPS_ESS1"] <- NA_real_
    benchmark_output["oCBPS_ESS0"] <- NA_real_
    benchmark_output["oCBPS_MAXW1"] <- NA_real_
    benchmark_output["oCBPS_MAXW0"] <- NA_real_
    benchmark_output["oCBPS_SUCCESS"] <- 0
  }
  
  
  # ==============================================================
  # 7. CBPS
  # ==============================================================
  
  cbps_inf <- tryCatch(
    estimate_cbps_inference(
      dat = dat,
      true_ATE = true_ATE
    ),
    error = function(e) NULL
  )
  
  if (!is.null(cbps_inf)) {
    
    estimates["CBPS"] <-
      cbps_inf$ATE
    
    benchmark_output["CBPS_SE_ANA"] <-
      cbps_inf$SE
    
    benchmark_output["CBPS_COV_ANA"] <-
      cbps_inf$coverage
    
    benchmark_output["CBPS_ESS1"] <-
      cbps_inf$ESS1
    
    benchmark_output["CBPS_ESS0"] <-
      cbps_inf$ESS0
    
    benchmark_output["CBPS_MAXW1"] <-
      cbps_inf$max_weight1
    
    benchmark_output["CBPS_MAXW0"] <-
      cbps_inf$max_weight0
    
    benchmark_output["CBPS_SUCCESS"] <-
      as.integer(cbps_inf$success)
    
  } else {
    
    estimates["CBPS"] <- NA_real_
    
    benchmark_output["CBPS_SE_ANA"] <- NA_real_
    benchmark_output["CBPS_COV_ANA"] <- NA_real_
    benchmark_output["CBPS_ESS1"] <- NA_real_
    benchmark_output["CBPS_ESS0"] <- NA_real_
    benchmark_output["CBPS_MAXW1"] <- NA_real_
    benchmark_output["CBPS_MAXW0"] <- NA_real_
    benchmark_output["CBPS_SUCCESS"] <- 0
  }
  
  
  # ==============================================================
  # 8. EBCW
  # ==============================================================
  
  ebcw_inf <- tryCatch(
    estimate_ebcw_inference(
      dat = dat,
      true_ATE = true_ATE
    ),
    error = function(e) NULL
  )
  
  if (!is.null(ebcw_inf)) {
    
    estimates["EBCW"] <-
      ebcw_inf$ATE
    
    benchmark_output["EBCW_SE_ANA"] <-
      ebcw_inf$SE
    
    benchmark_output["EBCW_COV_ANA"] <-
      ebcw_inf$coverage
    
    benchmark_output["EBCW_ESS1"] <-
      ebcw_inf$ESS1
    
    benchmark_output["EBCW_ESS0"] <-
      ebcw_inf$ESS0
    
    benchmark_output["EBCW_MAXW1"] <-
      ebcw_inf$max_weight1
    
    benchmark_output["EBCW_MAXW0"] <-
      ebcw_inf$max_weight0
    
    benchmark_output["EBCW_SUCCESS"] <-
      as.integer(ebcw_inf$success)
    
  } else {
    
    estimates["EBCW"] <- NA_real_
    
    benchmark_output["EBCW_SE_ANA"] <- NA_real_
    benchmark_output["EBCW_COV_ANA"] <- NA_real_
    benchmark_output["EBCW_ESS1"] <- NA_real_
    benchmark_output["EBCW_ESS0"] <- NA_real_
    benchmark_output["EBCW_MAXW1"] <- NA_real_
    benchmark_output["EBCW_MAXW0"] <- NA_real_
    benchmark_output["EBCW_SUCCESS"] <- 0
  }
  
  
  # ==============================================================
  # 9. Cross-fitting objects
  # ==============================================================
  
  fold_T1_lm <- NULL
  fold_T0_lm <- NULL
  
  fold_T1_gam <- NULL
  fold_T0_gam <- NULL
  
  
  # ==============================================================
  # 10. AIPW-LM
  # ==============================================================
  
  if (run_lm) {
    
    fold_T1_lm <- k_fold_function_lm(
      df = dat,
      K = K,
      idx_S = which(D == 1),
      idx_U0 = which(D == 0),
      seed = rep_id + 20242022
    )
    
    fold_T0_lm <- k_fold_function_lm(
      df = dat,
      K = K,
      idx_S = which(D == 0),
      idx_U0 = which(D == 1),
      seed = rep_id + 20242022
    )
    
    
    aipw_lm_inf <- tryCatch(
      estimate_aipw_inference(
        fold_T1 = fold_T1_lm,
        fold_T0 = fold_T0_lm,
        true_ATE = true_ATE
      ),
      error = function(e) NULL
    )
    
    
    if (!is.null(aipw_lm_inf)) {
      
      estimates["AIPW_LM"] <-
        aipw_lm_inf$ATE
      
      benchmark_output["AIPW_LM_SE_ANA"] <-
        aipw_lm_inf$SE
      
      benchmark_output["AIPW_LM_COV_ANA"] <-
        aipw_lm_inf$coverage
      
      benchmark_output["AIPW_LM_ESS1"] <-
        aipw_lm_inf$ESS1
      
      benchmark_output["AIPW_LM_ESS0"] <-
        aipw_lm_inf$ESS0
      
      benchmark_output["AIPW_LM_MAXW1"] <-
        aipw_lm_inf$max_weight1
      
      benchmark_output["AIPW_LM_MAXW0"] <-
        aipw_lm_inf$max_weight0
      
      benchmark_output["AIPW_LM_SUCCESS"] <-
        as.integer(aipw_lm_inf$success)
    }
  }
  
  
  # ==============================================================
  # 11. AIPW-GAM
  # ==============================================================
  
  if (run_gam) {
    
    fold_T1_gam <- k_fold_function_gam(
      df = dat,
      K = K,
      idx_S = which(D == 1),
      idx_U0 = which(D == 0),
      seed = rep_id + 20242022
    )
    
    fold_T0_gam <- k_fold_function_gam(
      df = dat,
      K = K,
      idx_S = which(D == 0),
      idx_U0 = which(D == 1),
      seed = rep_id + 20242022
    )
    
    
    aipw_gam_inf <- tryCatch(
      estimate_aipw_inference(
        fold_T1 = fold_T1_gam,
        fold_T0 = fold_T0_gam,
        true_ATE = true_ATE
      ),
      error = function(e) NULL
    )
    
    
    if (!is.null(aipw_gam_inf)) {
      
      estimates["AIPW_GAM"] <-
        aipw_gam_inf$ATE
      
      benchmark_output["AIPW_GAM_SE_ANA"] <-
        aipw_gam_inf$SE
      
      benchmark_output["AIPW_GAM_COV_ANA"] <-
        aipw_gam_inf$coverage
      
      benchmark_output["AIPW_GAM_ESS1"] <-
        aipw_gam_inf$ESS1
      
      benchmark_output["AIPW_GAM_ESS0"] <-
        aipw_gam_inf$ESS0
      
      benchmark_output["AIPW_GAM_MAXW1"] <-
        aipw_gam_inf$max_weight1
      
      benchmark_output["AIPW_GAM_MAXW0"] <-
        aipw_gam_inf$max_weight0
      
      benchmark_output["AIPW_GAM_SUCCESS"] <-
        as.integer(aipw_gam_inf$success)
    }
  }
  
  
  # ==============================================================
  # 12. GEC: ET / HD / CE
  #
  # TEMPORARILY MUTED WHEN run_gec = FALSE
  # ==============================================================
  
  if (run_gec) {
    
    # ------------------------------------------------------------
    # Choose calibration functions
    # ------------------------------------------------------------
    
    if (oracle_m) {
      
      if (outcome_model != 2) {
        stop(
          "oracle_m diagnostic currently set up for OR2 only."
        )
      }
      
      oracle <- true_m_or2(dat)
      
      dat1_oracle <- dat
      dat0_oracle <- dat
      
      dat1_oracle$y.hat <- oracle$m1
      dat0_oracle$y.hat <- oracle$m0
      
      fold_t1_gec <- list(dat1_oracle)
      fold_t0_gec <- list(dat0_oracle)
      
    } else {
      
      if (gec_from == "gam") {
        
        if (is.null(fold_T1_gam) ||
            is.null(fold_T0_gam)) {
          
          fold_T1_gam <- k_fold_function_gam(
            df = dat,
            K = K,
            idx_S = which(D == 1),
            idx_U0 = which(D == 0),
            seed = rep_id + 20242022
          )
          
          fold_T0_gam <- k_fold_function_gam(
            df = dat,
            K = K,
            idx_S = which(D == 0),
            idx_U0 = which(D == 1),
            seed = rep_id + 20242022
          )
        }
        
        fold_t1_gec <- fold_T1_gam
        fold_t0_gec <- fold_T0_gam
        
      } else {
        
        if (is.null(fold_T1_lm) ||
            is.null(fold_T0_lm)) {
          
          fold_T1_lm <- k_fold_function_lm(
            df = dat,
            K = K,
            idx_S = which(D == 1),
            idx_U0 = which(D == 0),
            seed = rep_id + 20242022
          )
          
          fold_T0_lm <- k_fold_function_lm(
            df = dat,
            K = K,
            idx_S = which(D == 0),
            idx_U0 = which(D == 1),
            seed = rep_id + 20242022
          )
        }
        
        fold_t1_gec <- fold_T1_lm
        fold_t0_gec <- fold_T0_lm
      }
    }
    
    
    # ------------------------------------------------------------
    # Reconstruct GEC prediction datasets
    # ------------------------------------------------------------
    
    dat1_gec <- do.call(
      rbind,
      fold_t1_gec
    )
    
    dat0_gec <- do.call(
      rbind,
      fold_t0_gec
    )
    
    dat1_gec <- dat1_gec[
      order(dat1_gec$ID),
    ]
    
    dat0_gec <- dat0_gec[
      order(dat0_gec$ID),
    ]
    
    dat <- dat[
      order(dat$ID),
    ]
    
    
    stopifnot(
      identical(
        as.integer(dat1_gec$ID),
        as.integer(dat$ID)
      ),
      identical(
        as.integer(dat0_gec$ID),
        as.integer(dat$ID)
      )
    )
    
    
    yhat1 <-
      as.numeric(
        dat1_gec$y.hat
      )
    
    yhat0 <-
      as.numeric(
        dat0_gec$y.hat
      )
    
    
    # ------------------------------------------------------------
    # PS model for GEC sandwich
    # ------------------------------------------------------------
    
    xvars <- setdiff(
      names(dat),
      c(
        "D",
        "y",
        "ID",
        "pi.hat",
        "y.hat"
      )
    )
    
    
    ps_formula <- as.formula(
      paste(
        "D ~",
        paste(
          xvars,
          collapse = " + "
        )
      )
    )
    
    
    ps_fit <- glm(
      ps_formula,
      data = dat,
      family = binomial()
    )
    
    
    phi_hat <- coef(
      ps_fit
    )
    
    
    X_ps <- model.matrix(
      ps_fit
    )
    
    
    pi_difference <-
      max(
        abs(
          fitted(ps_fit) -
            dat$pi.hat
        )
      )
    
    
    # ------------------------------------------------------------
    # ET / HD / CE
    # ------------------------------------------------------------
    
    for (ent in c(
      "ET",
      "HD",
      "CE"
    )) {
      
      fit_gec <- tryCatch(
        
        estimate_ATE_dual(
          fold_t1 = fold_t1_gec,
          fold_t0 = fold_t0_gec,
          entropy = ent
        ),
        
        error = function(e) NULL
      )
      
      
      if (is.null(fit_gec) ||
          !is.finite(fit_gec$ATE) ||
          is.null(fit_gec$fit1) ||
          is.null(fit_gec$fit0)) {
        
        estimates[ent] <- NA_real_
        
        gec_output[paste0(ent, "_SE_ANA")] <- NA_real_
        gec_output[paste0(ent, "_SE_NUM")] <- NA_real_
        gec_output[paste0(ent, "_COV_ANA")] <- NA_real_
        gec_output[paste0(ent, "_COV_NUM")] <- NA_real_
        gec_output[paste0(ent, "_ESS1")] <- NA_real_
        gec_output[paste0(ent, "_ESS0")] <- NA_real_
        gec_output[paste0(ent, "_MAXW1")] <- NA_real_
        gec_output[paste0(ent, "_MAXW0")] <- NA_real_
        gec_output[paste0(ent, "_JACDIFF")] <- NA_real_
        gec_output[paste0(ent, "_SUCCESS")] <- 0
        
        next
      }
      
      
      estimates[ent] <-
        fit_gec$ATE
      
      
      inf_gec <- tryCatch(
        
        gec_ate_inference(
          dat = dat,
          yhat1 = yhat1,
          yhat0 = yhat0,
          fit1 = fit_gec$fit1,
          fit0 = fit_gec$fit0,
          ps_fit = ps_fit,
          entropy = ent,
          true_ATE = true_ATE,
          numerical_jacobian =
            numerical_jacobian
        ),
        
        error = function(e) NULL
      )
      
      
      if (is.null(inf_gec) ||
          !isTRUE(inf_gec$success)) {
        
        gec_output[paste0(ent, "_SE_ANA")] <- NA_real_
        gec_output[paste0(ent, "_SE_NUM")] <- NA_real_
        gec_output[paste0(ent, "_COV_ANA")] <- NA_real_
        gec_output[paste0(ent, "_COV_NUM")] <- NA_real_
        gec_output[paste0(ent, "_ESS1")] <- NA_real_
        gec_output[paste0(ent, "_ESS0")] <- NA_real_
        gec_output[paste0(ent, "_MAXW1")] <- NA_real_
        gec_output[paste0(ent, "_MAXW0")] <- NA_real_
        gec_output[paste0(ent, "_JACDIFF")] <- NA_real_
        gec_output[paste0(ent, "_SUCCESS")] <- 0
        
        next
      }
      
      
      gec_output[paste0(ent, "_SE_ANA")] <-
        inf_gec$SE_analytic
      
      gec_output[paste0(ent, "_SE_NUM")] <-
        inf_gec$SE_numeric
      
      gec_output[paste0(ent, "_COV_ANA")] <-
        inf_gec$cover_analytic
      
      gec_output[paste0(ent, "_COV_NUM")] <-
        inf_gec$cover_numeric
      
      gec_output[paste0(ent, "_ESS1")] <-
        inf_gec$ESS1
      
      gec_output[paste0(ent, "_ESS0")] <-
        inf_gec$ESS0
      
      gec_output[paste0(ent, "_MAXW1")] <-
        inf_gec$max_weight1
      
      gec_output[paste0(ent, "_MAXW0")] <-
        inf_gec$max_weight0
      
      gec_output[paste0(ent, "_JACDIFF")] <-
        inf_gec$jacobian_max_diff
      
      gec_output[paste0(ent, "_SUCCESS")] <- 1
    }
    
  } else {
    
    # ============================================================
    # ET / HD / CE intentionally skipped for now
    # ============================================================
    
    pi_difference <- NA_real_
  }
  
  
  # ==============================================================
  # 13. Miscellaneous
  # ==============================================================
  
  misc_output <- c(
    TRUE_ATE = true_ATE,
    N_TREATED = sum(dat$D == 1),
    N_CONTROL = sum(dat$D == 0),
    PI_FIT_DIFF = pi_difference
  )
  
  
  # ==============================================================
  # 14. Return
  # ==============================================================
  
  output <- c(
    estimates,
    benchmark_output,
    gec_output,
    misc_output
  )
  
  
  as.data.frame(
    as.list(output),
    check.names = FALSE
  )
}

# ---------------------------------------------------
# 6. Scenario runner
# ---------------------------------------------------


# Function: run_simulation_scenario
# ---------------------------------
# Repeats one chosen scenario m times and stacks the results.
run_simulation_scenario_old <- function(n = 1000, p = 4, m = 500, K = 4,
                                        outcome_model = 1, ps_model = 1,
                                        run_lm = TRUE,
                                        run_gam = TRUE,
                                        gec_from = c("gam", "lm"),
                                        numerical_jacobian = TRUE,progress = TRUE, oracle_m = FALSE) {
  gec_from <- match.arg(gec_from)
  
  out <- vector("list", m)
  
  for (rep in seq_len(m)) {
    if (progress && rep %% 25 == 0) {
      message("Scenario OR", outcome_model, "PS", ps_model, ": replication ", rep, "/", m)
    }
    
    out[[rep]] <- tryCatch(
      run_one_replication(
        n = n, p = p, K = K, rep_id = rep,
        outcome_model = outcome_model,
        ps_model = ps_model,
        run_lm = run_lm,
        run_gam = run_gam,
        gec_from = gec_from,numerical_jacobian = numerical_jacobian, oracle_m = oracle_m
      ),
      error = function(e) {
        warning(sprintf("Replication %s failed: %s", rep, e$message))
        NULL
      }
    )
  }
  
  dplyr::bind_rows(out)
}

run_simulation_scenario <- function(
    n = 1000,
    p = 4,
    m = 500,
    K = 4,
    outcome_model = 1,
    ps_model = 1,
    run_lm = TRUE,
    run_gam = TRUE,
    gec_from = c("gam", "lm"),
    numerical_jacobian = TRUE,
    progress = TRUE,
    oracle_m = FALSE,
    run_gec = FALSE
) {
  
  gec_from <- match.arg(gec_from)
  
  out <- vector(
    "list",
    m
  )
  
  for (rep in seq_len(m)) {
    
    if (progress &&
        rep %% 25 == 0) {
      
      message(
        "Scenario OR",
        outcome_model,
        "PS",
        ps_model,
        ": replication ",
        rep,
        "/",
        m
      )
    }
    
    
    out[[rep]] <- tryCatch(
      
      run_one_replication(
        n = n,
        p = p,
        K = K,
        rep_id = rep,
        outcome_model = outcome_model,
        ps_model = ps_model,
        run_lm = run_lm,
        run_gam = run_gam,
        gec_from = gec_from,
        numerical_jacobian =
          numerical_jacobian,
        oracle_m = oracle_m,
        run_gec = run_gec
      ),
      
      error = function(e) {
        
        warning(
          sprintf(
            "Replication %s failed: %s",
            rep,
            e$message
          )
        )
        
        NULL
      }
    )
  }
  
  
  dplyr::bind_rows(
    out
  )
}

# Function: run_all_scenarios
# ---------------------------
# Runs all four combinations of outcome model and propensity score model and
# returns a named list of result tables.
run_all_scenarios_old <- function(n , p , m , K ,
                                  run_lm = TRUE,
                                  run_gam = TRUE,
                                  gec_from = c("gam", "lm"),numerical_jacobian = TRUE,
                                  progress = TRUE, oracle_m = FALSE) {
  gec_from <- match.arg(gec_from)
  
  scenarios <- list(
    OR1PS1 = c(1, 1),
    OR1PS2 = c(1, 2),
    OR2PS1 = c(2, 1),
    OR2PS2 = c(2, 2)
  )
  
  results <- lapply(names(scenarios), function(name) {
    om_ps <- scenarios[[name]]
    run_simulation_scenario(
      n = n, p = p, m = m, K = K,
      outcome_model = om_ps[1],
      ps_model = om_ps[2],
      run_lm = run_lm,
      run_gam = run_gam,
      gec_from = gec_from,numerical_jacobian = numerical_jacobian,
      progress = progress, oracle_m = FALSE
    )
  })
  
  names(results) <- names(scenarios)
  results
}




run_all_scenarios <- function(
    n,
    p,
    m,
    K,
    run_lm = TRUE,
    run_gam = TRUE,
    gec_from = c("gam", "lm"),
    numerical_jacobian = TRUE,
    progress = TRUE,
    oracle_m = FALSE,
    run_gec = FALSE          # NEW
) {
  
  gec_from <- match.arg(gec_from)
  
  scenarios <- list(
    OR1PS1 = c(1, 1),
    OR1PS2 = c(1, 2),
    OR2PS1 = c(2, 1),
    OR2PS2 = c(2, 2)
  )
  
  results <- lapply(
    names(scenarios),
    function(name) {
      
      om_ps <- scenarios[[name]]
      
      run_simulation_scenario(
        n = n,
        p = p,
        m = m,
        K = K,
        
        outcome_model = om_ps[1],
        ps_model = om_ps[2],
        
        run_lm = run_lm,
        run_gam = run_gam,
        
        gec_from = gec_from,
        
        numerical_jacobian =
          numerical_jacobian,
        
        progress = progress,
        
        oracle_m = oracle_m,   # FIXED
        
        run_gec = run_gec      # NEW
      )
    }
  )
  
  names(results) <-
    names(scenarios)
  
  results
}



make_table1 <- function(
    results,
    methods = c("ET", "HD", "CE")
) {
  
  true_ates <- c(
    OR1PS1 = 1,
    OR1PS2 = 1,
    OR2PS1 = 10,
    OR2PS2 = 10
  )
  
  out <- list()
  
  counter <- 1
  
  for (sc in names(results)) {
    
    x <- results[[sc]]
    
    true_ate <- true_ates[sc]
    
    for (meth in methods) {
      
      # ----------------------------------------------------------
      # Point estimates
      # ----------------------------------------------------------
      
      # ============================================================
      # Extract method-specific quantities
      # ============================================================
      
      est <- x[[meth]]
      
      se_ana <- x[[paste0(meth, "_SE_ANA")]]
      se_num <- x[[paste0(meth, "_SE_NUM")]]
      
      ESS1 <- x[[paste0(meth, "_ESS1")]]
      ESS0 <- x[[paste0(meth, "_ESS0")]]
      
      maxw1 <- x[[paste0(meth, "_MAXW1")]]
      maxw0 <- x[[paste0(meth, "_MAXW0")]]
      
      success <- x[[paste0(meth, "_SUCCESS")]]
      
      jacdiff <- x[[paste0(meth, "_JACDIFF")]]
      
      
      # ============================================================
      # COMMON replication set
      #
      # Same observations used for:
      # estimate, bias, MC SD, RMSE, SEs, coverage, ESS and weights
      # ============================================================
      # Main estimator/analytic-inference subset
      analytic_ok <-
        is.finite(est) &
        is.finite(se_ana) &
        success == 1
      
      # Only for analytic-vs-numerical Jacobian comparison
      common_ok <-
        analytic_ok &
        is.finite(se_num)
      
      N_total <- nrow(x)
      N_common <- sum(common_ok)
      
      
      # ============================================================
      # Subset everything to common replications
      # ============================================================
      
      est_c <- est[common_ok]
      
      se_ana_c <- se_ana[common_ok]
      se_num_c <- se_num[common_ok]
      
      ESS1_c <- ESS1[common_ok]
      ESS0_c <- ESS0[common_ok]
      
      maxw1_c <- maxw1[common_ok]
      maxw0_c <- maxw0[common_ok]
      
      jacdiff_c <- jacdiff[common_ok]
      
      N1_c <- x$N_TREATED[common_ok]
      N0_c <- x$N_CONTROL[common_ok]
      
      
      # ============================================================
      # Point-estimator performance
      # ============================================================
      
      Mean_Estimate <- mean(
        est_c,
        na.rm = TRUE
      )
      
      Bias <- mean(
        est_c - true_ate,
        na.rm = TRUE
      )
      
      MC_SD <- sd(
        est_c,
        na.rm = TRUE
      )
      
      RMSE <- sqrt(
        mean(
          (est_c - true_ate)^2,
          na.rm = TRUE
        )
      )
      
      
      # ============================================================
      # Analytic and numerical sandwich SE
      # ============================================================
      
      Avg_SE_Analytic <- mean(
        se_ana_c,
        na.rm = TRUE
      )
      
      Avg_SE_Numeric <- mean(
        se_num_c,
        na.rm = TRUE
      )
      
      
      # ============================================================
      # Coverage -- recomputed on SAME common replications
      # ============================================================
      
      lower_ana <-
        est_c - qnorm(0.975) * se_ana_c
      
      upper_ana <-
        est_c + qnorm(0.975) * se_ana_c
      
      cover_ana_c <-
        as.integer(
          lower_ana <= true_ate &
            true_ate <= upper_ana
        )
      
      Coverage_Analytic <- mean(
        cover_ana_c
      )
      
      
      lower_num <-
        est_c - qnorm(0.975) * se_num_c
      
      upper_num <-
        est_c + qnorm(0.975) * se_num_c
      
      cover_num_c <-
        as.integer(
          lower_num <= true_ate &
            true_ate <= upper_num
        )
      
      Coverage_Numeric <- mean(
        cover_num_c
      )
      
      
      # ============================================================
      # ESS
      # ============================================================
      
      Mean_ESS_Treated <- mean(
        ESS1_c,
        na.rm = TRUE
      )
      
      Mean_ESS_Control <- mean(
        ESS0_c,
        na.rm = TRUE
      )
      
      
      # ============================================================
      # Maximum calibration weights
      # ============================================================
      
      Mean_MaxW_Treated <- mean(
        maxw1_c,
        na.rm = TRUE
      )
      
      Mean_MaxW_Control <- mean(
        maxw0_c,
        na.rm = TRUE
      )
      
      
      # ============================================================
      # Mean arm sizes
      # ============================================================
      
      Mean_N_Treated <- mean(
        N1_c,
        na.rm = TRUE
      )
      
      Mean_N_Control <- mean(
        N0_c,
        na.rm = TRUE
      )
      
      
      # ============================================================
      # Analytic/numeric SE versus Monte Carlo SD
      # ============================================================
      
      SEratio_Analytic_MC <-
        Avg_SE_Analytic / MC_SD
      
      SEratio_Numeric_MC <-
        Avg_SE_Numeric / MC_SD
      
      
      # ============================================================
      # Jacobian check
      # ============================================================
      
      Mean_Jacobian_Diff <- mean(
        jacdiff_c,
        na.rm = TRUE
      )
      
      Max_Jacobian_Diff <- max(
        jacdiff_c,
        na.rm = TRUE
      )
      
      
      # ============================================================
      # Numerical availability / exclusion rate
      # ============================================================
      
      Common_Rate <-
        N_common / N_total
      
      Numerical_Failure_Rate <-
        1 - Common_Rate
      # ----------------------------------------------------------
      # Store row
      # ----------------------------------------------------------
      
      out[[counter]] <- data.frame(
        
        Scenario = sc,
        Method = meth,
        
        N_Total = N_total,
        N_Common = N_common,
        Common_Rate = Common_Rate,
        
        Mean_Estimate = Mean_Estimate,
        
        Bias = Bias,
        MC_SD = MC_SD,
        RMSE = RMSE,
        
        Avg_SE_Analytic =
          Avg_SE_Analytic,
        
        Avg_SE_Numeric =
          Avg_SE_Numeric,
        
        SEratio_Analytic_MC =
          SEratio_Analytic_MC,
        
        SEratio_Numeric_MC =
          SEratio_Numeric_MC,
        
        Coverage_Analytic =
          Coverage_Analytic,
        
        Coverage_Numeric =
          Coverage_Numeric,
        
        Mean_ESS_Treated =
          Mean_ESS_Treated,
        
        Mean_ESS_Control =
          Mean_ESS_Control,
        
        Mean_MaxW_Treated =
          Mean_MaxW_Treated,
        
        Mean_MaxW_Control =
          Mean_MaxW_Control,
        
        Mean_N_Treated =
          Mean_N_Treated,
        
        Mean_N_Control =
          Mean_N_Control,
        
        Numerical_Failure_Rate =
          Numerical_Failure_Rate,
        
        Mean_Jacobian_Diff =
          Mean_Jacobian_Diff,
        
        Max_Jacobian_Diff =
          Max_Jacobian_Diff
      )
      counter <- counter + 1
    }
  }
  
  dplyr::bind_rows(out)
}

make_table <- function(
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
) {
  
  true_ates <- c(
    OR1PS1 = 1,
    OR1PS2 = 1,
    OR2PS1 = 10,
    OR2PS2 = 10
  )
  
  gec_methods <- c(
    "ET",
    "HD",
    "CE"
  )
  
  safe_mean <- function(z) {
    z <- z[is.finite(z)]
    
    if (length(z) > 0L) {
      mean(z)
    } else {
      NA_real_
    }
  }
  
  safe_sd <- function(z) {
    z <- z[is.finite(z)]
    
    if (length(z) > 1L) {
      sd(z)
    } else {
      NA_real_
    }
  }
  
  safe_max <- function(z) {
    z <- z[is.finite(z)]
    
    if (length(z) > 0L) {
      max(z)
    } else {
      NA_real_
    }
  }
  
  out <- list()
  counter <- 1L
  
  for (sc in names(results)) {
    
    x <- results[[sc]]
    
    if (is.null(x) || nrow(x) == 0L) {
      next
    }
    
    true_ate <- true_ates[[sc]]
    
    for (meth in methods) {
      
      est_col <- meth
      se_col <- paste0(meth, "_SE_ANA")
      ess1_col <- paste0(meth, "_ESS1")
      ess0_col <- paste0(meth, "_ESS0")
      maxw1_col <- paste0(meth, "_MAXW1")
      maxw0_col <- paste0(meth, "_MAXW0")
      success_col <- paste0(meth, "_SUCCESS")
      
      needed <- c(
        est_col,
        se_col,
        ess1_col,
        ess0_col,
        maxw1_col,
        maxw0_col,
        success_col
      )
      
      missing_cols <- setdiff(
        needed,
        names(x)
      )
      
      if (length(missing_cols) > 0L) {
        
        warning(
          sprintf(
            "%s / %s skipped; missing columns: %s",
            sc,
            meth,
            paste(
              missing_cols,
              collapse = ", "
            )
          )
        )
        
        next
      }
      
      est <- x[[est_col]]
      se_ana <- x[[se_col]]
      
      ESS1 <- x[[ess1_col]]
      ESS0 <- x[[ess0_col]]
      
      maxw1 <- x[[maxw1_col]]
      maxw0 <- x[[maxw0_col]]
      
      success <- x[[success_col]]
      
      N_total <- nrow(x)
      
      analytic_ok <-
        is.finite(est) &
        is.finite(se_ana) &
        !is.na(success) &
        success == 1
      
      N_analytic <- sum(analytic_ok)
      
      Analytic_Rate <-
        N_analytic /
        N_total
      
      est_a <-
        est[analytic_ok]
      
      se_ana_a <-
        se_ana[analytic_ok]
      
      ESS1_a <-
        ESS1[analytic_ok]
      
      ESS0_a <-
        ESS0[analytic_ok]
      
      maxw1_a <-
        maxw1[analytic_ok]
      
      maxw0_a <-
        maxw0[analytic_ok]
      
      if ("N_TREATED" %in% names(x)) {
        N1_a <-
          x[["N_TREATED"]][analytic_ok]
      } else {
        N1_a <- NA_real_
      }
      
      if ("N_CONTROL" %in% names(x)) {
        N0_a <-
          x[["N_CONTROL"]][analytic_ok]
      } else {
        N0_a <- NA_real_
      }
      
      Mean_Estimate <-
        safe_mean(est_a)
      
      Bias <-
        if (N_analytic > 0L) {
          
          mean(
            est_a -
              true_ate
          )
          
        } else {
          
          NA_real_
        }
      
      MC_SD <-
        safe_sd(est_a)
      
      RMSE <-
        if (N_analytic > 0L) {
          
          sqrt(
            mean(
              (est_a -
                 true_ate)^2
            )
          )
          
        } else {
          
          NA_real_
        }
      
      Avg_SE_Analytic <-
        safe_mean(
          se_ana_a
        )
      
      SEratio_Analytic_MC <-
        if (
          is.finite(MC_SD) &&
          MC_SD > 0 &&
          is.finite(Avg_SE_Analytic)
        ) {
          
          Avg_SE_Analytic /
            MC_SD
          
        } else {
          
          NA_real_
        }
      
      if (N_analytic > 0L) {
        
        lower_ana <-
          est_a -
          qnorm(0.975) *
          se_ana_a
        
        upper_ana <-
          est_a +
          qnorm(0.975) *
          se_ana_a
        
        Coverage_Analytic <-
          mean(
            lower_ana <= true_ate &
              true_ate <= upper_ana
          )
        
      } else {
        
        Coverage_Analytic <-
          NA_real_
      }
      
      Mean_ESS_Treated <-
        safe_mean(
          ESS1_a
        )
      
      Mean_ESS_Control <-
        safe_mean(
          ESS0_a
        )
      
      Mean_MaxW_Treated <-
        safe_mean(
          maxw1_a
        )
      
      Mean_MaxW_Control <-
        safe_mean(
          maxw0_a
        )
      
      Worst_MaxW_Treated <-
        safe_max(
          maxw1_a
        )
      
      Worst_MaxW_Control <-
        safe_max(
          maxw0_a
        )
      
      Mean_N_Treated <-
        safe_mean(
          N1_a
        )
      
      Mean_N_Control <-
        safe_mean(
          N0_a
        )
      
      
      # =====================================================
      # Numerical-Jacobian diagnostics only for ET/HD/CE
      # =====================================================
      
      N_common <- NA_integer_
      Common_Rate <- NA_real_
      
      MC_SD_Common <- NA_real_
      
      Avg_SE_Analytic_Common <- NA_real_
      Avg_SE_Numeric_Common <- NA_real_
      
      SEratio_Num_to_Ana_Common <- NA_real_
      
      SEratio_Analytic_MC_Common <- NA_real_
      SEratio_Numeric_MC_Common <- NA_real_
      
      Coverage_Analytic_Common <- NA_real_
      Coverage_Numeric_Common <- NA_real_
      
      Numerical_Failure_Rate <- NA_real_
      
      Mean_Jacobian_Diff <- NA_real_
      Max_Jacobian_Diff <- NA_real_
      
      
      if (meth %in% gec_methods) {
        
        se_num_col <-
          paste0(
            meth,
            "_SE_NUM"
          )
        
        jacdiff_col <-
          paste0(
            meth,
            "_JACDIFF"
          )
        
        if (
          se_num_col %in% names(x) &&
          jacdiff_col %in% names(x)
        ) {
          
          se_num <-
            x[[se_num_col]]
          
          jacdiff <-
            x[[jacdiff_col]]
          
          common_ok <-
            analytic_ok &
            is.finite(se_num)
          
          N_common <-
            sum(common_ok)
          
          Common_Rate <-
            if (N_analytic > 0L) {
              
              N_common /
                N_analytic
              
            } else {
              
              NA_real_
            }
          
          Numerical_Failure_Rate <-
            if (is.finite(Common_Rate)) {
              
              1 -
                Common_Rate
              
            } else {
              
              NA_real_
            }
          
          est_c <-
            est[common_ok]
          
          se_ana_c <-
            se_ana[common_ok]
          
          se_num_c <-
            se_num[common_ok]
          
          MC_SD_Common <-
            safe_sd(
              est_c
            )
          
          Avg_SE_Analytic_Common <-
            safe_mean(
              se_ana_c
            )
          
          Avg_SE_Numeric_Common <-
            safe_mean(
              se_num_c
            )
          
          SEratio_Num_to_Ana_Common <-
            if (
              is.finite(
                Avg_SE_Analytic_Common
              ) &&
              Avg_SE_Analytic_Common > 0 &&
              is.finite(
                Avg_SE_Numeric_Common
              )
            ) {
              
              Avg_SE_Numeric_Common /
                Avg_SE_Analytic_Common
              
            } else {
              
              NA_real_
            }
          
          SEratio_Analytic_MC_Common <-
            if (
              is.finite(MC_SD_Common) &&
              MC_SD_Common > 0 &&
              is.finite(
                Avg_SE_Analytic_Common
              )
            ) {
              
              Avg_SE_Analytic_Common /
                MC_SD_Common
              
            } else {
              
              NA_real_
            }
          
          SEratio_Numeric_MC_Common <-
            if (
              is.finite(MC_SD_Common) &&
              MC_SD_Common > 0 &&
              is.finite(
                Avg_SE_Numeric_Common
              )
            ) {
              
              Avg_SE_Numeric_Common /
                MC_SD_Common
              
            } else {
              
              NA_real_
            }
          
          if (N_common > 0L) {
            
            lower_ana_c <-
              est_c -
              qnorm(0.975) *
              se_ana_c
            
            upper_ana_c <-
              est_c +
              qnorm(0.975) *
              se_ana_c
            
            Coverage_Analytic_Common <-
              mean(
                lower_ana_c <= true_ate &
                  true_ate <= upper_ana_c
              )
            
            lower_num_c <-
              est_c -
              qnorm(0.975) *
              se_num_c
            
            upper_num_c <-
              est_c +
              qnorm(0.975) *
              se_num_c
            
            Coverage_Numeric_Common <-
              mean(
                lower_num_c <= true_ate &
                  true_ate <= upper_num_c
              )
          }
          
          jacdiff_common <-
            jacdiff[
              common_ok
            ]
          
          Mean_Jacobian_Diff <-
            safe_mean(
              jacdiff_common
            )
          
          Max_Jacobian_Diff <-
            safe_max(
              jacdiff_common
            )
        }
      }
      N_success <- sum(
        !is.na(success) & success == 1
      )
      
      N_failure <- N_total - N_success
      
      Failure_Rate <- N_failure / N_total
      
      out[[counter]] <- data.frame(
        
        Scenario = sc,
        
        Method = meth,
        
        N_Total = N_total,
        N_Success = N_success,
        N_Failure = N_failure,
        Failure_Rate = Failure_Rate,
        
        N_Analytic =
          N_analytic,
        
        Analytic_Rate =
          Analytic_Rate,
        
        N_Common =
          N_common,
        
        Common_Rate =
          Common_Rate,
        
        Mean_Estimate =
          Mean_Estimate,
        
        Bias =
          Bias,
        
        MC_SD =
          MC_SD,
        
        RMSE =
          RMSE,
        
        Avg_SE_Analytic =
          Avg_SE_Analytic,
        
        SEratio_Analytic_MC =
          SEratio_Analytic_MC,
        
        Coverage_Analytic =
          Coverage_Analytic,
        
        Mean_ESS_Treated =
          Mean_ESS_Treated,
        
        Mean_ESS_Control =
          Mean_ESS_Control,
        
        Mean_MaxW_Treated =
          Mean_MaxW_Treated,
        
        Mean_MaxW_Control =
          Mean_MaxW_Control,
        
        Worst_MaxW_Treated =
          Worst_MaxW_Treated,
        
        Worst_MaxW_Control =
          Worst_MaxW_Control,
        
        Mean_N_Treated =
          Mean_N_Treated,
        
        Mean_N_Control =
          Mean_N_Control,
        
        MC_SD_Common =
          MC_SD_Common,
        
        Avg_SE_Analytic_Common =
          Avg_SE_Analytic_Common,
        
        Avg_SE_Numeric_Common =
          Avg_SE_Numeric_Common,
        
        SEratio_Num_to_Ana_Common =
          SEratio_Num_to_Ana_Common,
        
        SEratio_Analytic_MC_Common =
          SEratio_Analytic_MC_Common,
        
        SEratio_Numeric_MC_Common =
          SEratio_Numeric_MC_Common,
        
        Coverage_Analytic_Common =
          Coverage_Analytic_Common,
        
        Coverage_Numeric_Common =
          Coverage_Numeric_Common,
        
        Numerical_Failure_Rate =
          Numerical_Failure_Rate,
        
        Mean_Jacobian_Diff =
          Mean_Jacobian_Diff,
        
        Max_Jacobian_Diff =
          Max_Jacobian_Diff,
        
        stringsAsFactors = FALSE
      )
      
      counter <-
        counter +
        1L
    }
  }
  
  if (length(out) == 0L) {
    return(
      data.frame()
    )
  }
  
  dplyr::bind_rows(
    out
  )
}
make_table_old <- function(
    results,
    methods = c("ET", "HD", "CE")
) {
  
  # ============================================================
  # True ATE in each simulation scenario
  # ============================================================
  
  true_ates <- c(
    OR1PS1 = 1,
    OR1PS2 = 1,
    OR2PS1 = 10,
    OR2PS2 = 10
  )
  
  
  # ============================================================
  # Safe helper functions
  # ============================================================
  
  safe_mean <- function(z) {
    
    z <- z[is.finite(z)]
    
    if (length(z) > 0L) {
      mean(z)
    } else {
      NA_real_
    }
  }
  
  
  safe_sd <- function(z) {
    
    z <- z[is.finite(z)]
    
    if (length(z) > 1L) {
      sd(z)
    } else {
      NA_real_
    }
  }
  
  
  safe_max <- function(z) {
    
    z <- z[is.finite(z)]
    
    if (length(z) > 0L) {
      max(z)
    } else {
      NA_real_
    }
  }
  
  
  # ============================================================
  # Storage
  # ============================================================
  
  out <- list()
  
  counter <- 1L
  
  
  # ============================================================
  # Loop over scenarios
  # ============================================================
  
  for (sc in names(results)) {
    
    x <- results[[sc]]
    
    if (is.null(x) || nrow(x) == 0L) {
      next
    }
    
    
    true_ate <- true_ates[[sc]]
    
    
    # ==========================================================
    # Loop over GEC methods
    # ==========================================================
    
    for (meth in methods) {
      
      
      # ========================================================
      # Check required columns
      # ========================================================
      
      needed <- c(
        meth,
        paste0(meth, "_SE_ANA"),
        paste0(meth, "_SE_NUM"),
        paste0(meth, "_ESS1"),
        paste0(meth, "_ESS0"),
        paste0(meth, "_MAXW1"),
        paste0(meth, "_MAXW0"),
        paste0(meth, "_SUCCESS"),
        paste0(meth, "_JACDIFF"),
        "N_TREATED",
        "N_CONTROL"
      )
      
      
      missing_cols <- setdiff(
        needed,
        names(x)
      )
      
      
      if (length(missing_cols) > 0L) {
        
        warning(
          sprintf(
            "%s / %s skipped; missing columns: %s",
            sc,
            meth,
            paste(
              missing_cols,
              collapse = ", "
            )
          )
        )
        
        next
      }
      
      
      # ========================================================
      # Extract quantities
      # ========================================================
      
      est <-
        x[[meth]]
      
      se_ana <-
        x[[paste0(meth, "_SE_ANA")]]
      
      se_num <-
        x[[paste0(meth, "_SE_NUM")]]
      
      ESS1 <-
        x[[paste0(meth, "_ESS1")]]
      
      ESS0 <-
        x[[paste0(meth, "_ESS0")]]
      
      maxw1 <-
        x[[paste0(meth, "_MAXW1")]]
      
      maxw0 <-
        x[[paste0(meth, "_MAXW0")]]
      
      success <-
        x[[paste0(meth, "_SUCCESS")]]
      
      jacdiff <-
        x[[paste0(meth, "_JACDIFF")]]
      
      
      N_total <- nrow(x)
      
      
      # ========================================================
      # 1. MAIN ESTIMATOR / ANALYTIC-INFERENCE SUBSET
      #
      # Numerical Jacobian availability is NOT required here.
      #
      # This subset is used for:
      #   Mean estimate
      #   Bias
      #   MC SD
      #   RMSE
      #   Analytic sandwich SE
      #   Analytic coverage
      #   ESS
      #   Maximum weights
      # ========================================================
      
      analytic_ok <-
        is.finite(est) &
        is.finite(se_ana) &
        !is.na(success) &
        success == 1
      
      
      N_analytic <-
        sum(analytic_ok)
      
      
      Analytic_Rate <-
        N_analytic / N_total
      
      
      # ========================================================
      # Subset quantities for MAIN estimator comparison
      # ========================================================
      
      est_a <-
        est[analytic_ok]
      
      se_ana_a <-
        se_ana[analytic_ok]
      
      ESS1_a <-
        ESS1[analytic_ok]
      
      ESS0_a <-
        ESS0[analytic_ok]
      
      maxw1_a <-
        maxw1[analytic_ok]
      
      maxw0_a <-
        maxw0[analytic_ok]
      
      N1_a <-
        x[["N_TREATED"]][analytic_ok]
      
      N0_a <-
        x[["N_CONTROL"]][analytic_ok]
      
      
      # ========================================================
      # Main point-estimator performance
      # ========================================================
      
      Mean_Estimate <-
        safe_mean(est_a)
      
      
      Bias <-
        if (N_analytic > 0L) {
          
          mean(
            est_a - true_ate
          )
          
        } else {
          
          NA_real_
        }
      
      
      MC_SD <-
        safe_sd(est_a)
      
      
      RMSE <-
        if (N_analytic > 0L) {
          
          sqrt(
            mean(
              (est_a - true_ate)^2
            )
          )
          
        } else {
          
          NA_real_
        }
      
      
      # ========================================================
      # Main analytic sandwich SE
      # ========================================================
      
      Avg_SE_Analytic <-
        safe_mean(se_ana_a)
      
      
      # ========================================================
      # MAIN ANALYTIC COVERAGE
      #
      # Uses all analytic_ok replications.
      #
      # CI:
      # estimate +/- 1.96 * analytic SE
      # ========================================================
      
      if (N_analytic > 0L) {
        
        lower_ana <-
          est_a -
          qnorm(0.975) * se_ana_a
        
        
        upper_ana <-
          est_a +
          qnorm(0.975) * se_ana_a
        
        
        cover_ana <-
          as.integer(
            lower_ana <= true_ate &
              true_ate <= upper_ana
          )
        
        
        Coverage_Analytic <-
          mean(
            cover_ana
          )
        
      } else {
        
        Coverage_Analytic <-
          NA_real_
        
      }
      
      
      # ========================================================
      # Main ESS diagnostics
      # ========================================================
      
      Mean_ESS_Treated <-
        safe_mean(ESS1_a)
      
      Mean_ESS_Control <-
        safe_mean(ESS0_a)
      
      
      # ========================================================
      # Main maximum-weight diagnostics
      # ========================================================
      
      Mean_MaxW_Treated <-
        safe_mean(maxw1_a)
      
      Mean_MaxW_Control <-
        safe_mean(maxw0_a)
      
      
      Worst_MaxW_Treated <-
        safe_max(maxw1_a)
      
      Worst_MaxW_Control <-
        safe_max(maxw0_a)
      
      
      # ========================================================
      # Mean treated/control sample sizes
      # ========================================================
      
      Mean_N_Treated <-
        safe_mean(N1_a)
      
      Mean_N_Control <-
        safe_mean(N0_a)
      
      
      # ========================================================
      # Analytic SE versus Monte Carlo SD
      #
      # Main scientific check of sandwich performance.
      # ========================================================
      
      SEratio_Analytic_MC <-
        if (
          is.finite(MC_SD) &&
          MC_SD > 0 &&
          is.finite(Avg_SE_Analytic)
        ) {
          
          Avg_SE_Analytic / MC_SD
          
        } else {
          
          NA_real_
        }
      
      
      # ========================================================
      # 2. COMMON SUBSET
      #
      # ONLY for analytic-vs-numerical validation.
      #
      # Do NOT use this subset for main bias/RMSE/ESS/etc.
      # ========================================================
      
      common_ok <-
        analytic_ok &
        is.finite(se_num)
      
      
      N_common <-
        sum(common_ok)
      
      
      # Of analytically valid replications, how many also
      # have a numerical-Jacobian SE?
      Common_Rate <-
        if (N_analytic > 0L) {
          
          N_common / N_analytic
          
        } else {
          
          NA_real_
        }
      
      
      Numerical_Failure_Rate <-
        if (is.finite(Common_Rate)) {
          
          1 - Common_Rate
          
        } else {
          
          NA_real_
        }
      
      
      # ========================================================
      # Common-subset quantities
      # ========================================================
      
      est_c <-
        est[common_ok]
      
      se_ana_c <-
        se_ana[common_ok]
      
      se_num_c <-
        se_num[common_ok]
      
      
      # ========================================================
      # Fair analytic-versus-numerical SE comparison
      #
      # BOTH calculated on exactly the same replications.
      # ========================================================
      
      Avg_SE_Analytic_Common <-
        safe_mean(se_ana_c)
      
      
      Avg_SE_Numeric_Common <-
        safe_mean(se_num_c)
      
      
      # ========================================================
      # Analytic vs numerical SE ratio
      #
      # Should be close to 1 if implementations agree.
      # ========================================================
      
      SEratio_Num_to_Ana_Common <-
        if (
          is.finite(Avg_SE_Analytic_Common) &&
          Avg_SE_Analytic_Common > 0 &&
          is.finite(Avg_SE_Numeric_Common)
        ) {
          
          Avg_SE_Numeric_Common /
            Avg_SE_Analytic_Common
          
        } else {
          
          NA_real_
        }
      
      
      # ========================================================
      # COMMON-SUBSET ANALYTIC COVERAGE
      #
      # Used only for fair comparison with numerical coverage.
      # ========================================================
      
      if (N_common > 0L) {
        
        lower_ana_c <-
          est_c -
          qnorm(0.975) * se_ana_c
        
        
        upper_ana_c <-
          est_c +
          qnorm(0.975) * se_ana_c
        
        
        Coverage_Analytic_Common <-
          mean(
            lower_ana_c <= true_ate &
              true_ate <= upper_ana_c
          )
        
        
        # ======================================================
        # Numerical coverage on SAME common replications
        # ======================================================
        
        lower_num_c <-
          est_c -
          qnorm(0.975) * se_num_c
        
        
        upper_num_c <-
          est_c +
          qnorm(0.975) * se_num_c
        
        
        Coverage_Numeric_Common <-
          mean(
            lower_num_c <= true_ate &
              true_ate <= upper_num_c
          )
        
      } else {
        
        Coverage_Analytic_Common <-
          NA_real_
        
        Coverage_Numeric_Common <-
          NA_real_
        
      }
      
      
      # ========================================================
      # Numerical SE versus MC SD on COMMON subset
      #
      # This is only a secondary diagnostic because MC SD based
      # on the reduced common set may be unstable for CE.
      # ========================================================
      
      MC_SD_Common <-
        safe_sd(est_c)
      
      
      SEratio_Analytic_MC_Common <-
        if (
          is.finite(MC_SD_Common) &&
          MC_SD_Common > 0 &&
          is.finite(Avg_SE_Analytic_Common)
        ) {
          
          Avg_SE_Analytic_Common /
            MC_SD_Common
          
        } else {
          
          NA_real_
        }
      
      
      SEratio_Numeric_MC_Common <-
        if (
          is.finite(MC_SD_Common) &&
          MC_SD_Common > 0 &&
          is.finite(Avg_SE_Numeric_Common)
        ) {
          
          Avg_SE_Numeric_Common /
            MC_SD_Common
          
        } else {
          
          NA_real_
        }
      
      
      # ========================================================
      # Jacobian comparison
      #
      # Only meaningful when numerical Jacobian exists.
      # ========================================================
      
      jacdiff_common <-
        jacdiff[common_ok]
      
      
      Mean_Jacobian_Diff <-
        safe_mean(jacdiff_common)
      
      
      Max_Jacobian_Diff <-
        safe_max(jacdiff_common)
      
      
      # ========================================================
      # Store one scenario-method row
      # ========================================================
      
      out[[counter]] <- data.frame(
        
        Scenario = sc,
        
        Method = meth,
        
        
        # ------------------------------------------------------
        # Number of replications
        # ------------------------------------------------------
        
        N_Total =
          N_total,
        
        N_Analytic =
          N_analytic,
        
        Analytic_Rate =
          Analytic_Rate,
        
        N_Common =
          N_common,
        
        Common_Rate =
          Common_Rate,
        
        
        # ------------------------------------------------------
        # MAIN estimator performance
        # Uses analytic_ok
        # ------------------------------------------------------
        
        Mean_Estimate =
          Mean_Estimate,
        
        Bias =
          Bias,
        
        MC_SD =
          MC_SD,
        
        RMSE =
          RMSE,
        
        Avg_SE_Analytic =
          Avg_SE_Analytic,
        
        SEratio_Analytic_MC =
          SEratio_Analytic_MC,
        
        Coverage_Analytic =
          Coverage_Analytic,
        
        
        # ------------------------------------------------------
        # MAIN weighting diagnostics
        # Uses analytic_ok
        # ------------------------------------------------------
        
        Mean_ESS_Treated =
          Mean_ESS_Treated,
        
        Mean_ESS_Control =
          Mean_ESS_Control,
        
        Mean_MaxW_Treated =
          Mean_MaxW_Treated,
        
        Mean_MaxW_Control =
          Mean_MaxW_Control,
        
        Worst_MaxW_Treated =
          Worst_MaxW_Treated,
        
        Worst_MaxW_Control =
          Worst_MaxW_Control,
        
        Mean_N_Treated =
          Mean_N_Treated,
        
        Mean_N_Control =
          Mean_N_Control,
        
        
        # ------------------------------------------------------
        # ANALYTIC vs NUMERICAL validation
        # Uses common_ok only
        # ------------------------------------------------------
        
        MC_SD_Common =
          MC_SD_Common,
        
        Avg_SE_Analytic_Common =
          Avg_SE_Analytic_Common,
        
        Avg_SE_Numeric_Common =
          Avg_SE_Numeric_Common,
        
        SEratio_Num_to_Ana_Common =
          SEratio_Num_to_Ana_Common,
        
        SEratio_Analytic_MC_Common =
          SEratio_Analytic_MC_Common,
        
        SEratio_Numeric_MC_Common =
          SEratio_Numeric_MC_Common,
        
        Coverage_Analytic_Common =
          Coverage_Analytic_Common,
        
        Coverage_Numeric_Common =
          Coverage_Numeric_Common,
        
        
        # ------------------------------------------------------
        # Numerical-Jacobian availability
        # ------------------------------------------------------
        
        Numerical_Failure_Rate =
          Numerical_Failure_Rate,
        
        
        # ------------------------------------------------------
        # Jacobian check
        # ------------------------------------------------------
        
        Mean_Jacobian_Diff =
          Mean_Jacobian_Diff,
        
        Max_Jacobian_Diff =
          Max_Jacobian_Diff,
        
        
        stringsAsFactors = FALSE
      )
      
      
      counter <-
        counter + 1L
    }
  }
  
  
  # ============================================================
  # Return table
  # ============================================================
  
  if (length(out) == 0L) {
    return(
      data.frame()
    )
  }
  
  
  dplyr::bind_rows(out)
}
# ---------------------------------------------------
# 7. Plotting helpers
# ---------------------------------------------------

# Function: panel_boxplot_gg
# --------------------------
# Makes one boxplot panel for a given scenario.
panel_boxplot_gg <- function(data, true_ate, panel_label, ylim_range = NULL) {
  df_long <- tidyr::pivot_longer(
    as.data.frame(data),
    cols = everything(),
    names_to = "Method",
    values_to = "ATE"
  )
  df_long$Method <- factor(df_long$Method, levels = colnames(data))
  
  p <- ggplot(df_long, aes(x = Method, y = ATE)) +
    geom_boxplot(fill = "grey75", color = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = true_ate, color = "red", linewidth = 0.4) +
    labs(title = panel_label, x = NULL, y = "Point estimate") +
    theme_light() +
    theme(
      text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
  
  if (!is.null(ylim_range)) {
    p <- p + coord_cartesian(ylim = ylim_range)
  }
  p
}

# Function: plot_4panel_boxplots
# ------------------------------
# Combines the four scenario-specific boxplots into the final 2 x 2 figure.




panel_boxplot_gg <- function(
    data,
    true_ate,
    panel_label,
    ylim_range = NULL
) {
  
  df_long <- tidyr::pivot_longer(
    as.data.frame(data),
    cols = everything(),
    names_to = "Method",
    values_to = "ATE"
  )
  
  df_long$Method <- factor(
    df_long$Method,
    levels = c(
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
  
  p <- ggplot(
    df_long,
    aes(
      x = Method,
      y = ATE
    )
  ) +
    geom_boxplot(
      fill = "grey75",
      color = "grey40",
      linewidth = 0.4,
      outlier.size = 0.7
    ) +
    geom_hline(
      yintercept = true_ate,
      color = "red",
      linewidth = 0.5
    ) +
    labs(
      title = panel_label,
      x = NULL,
      y = "Point estimate"
    ) +
    theme_light() +
    theme(
      text = element_text(size = 11),
      plot.title = element_text(
        hjust = 0.5
      ),
      panel.grid = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    )
  
  if (!is.null(ylim_range)) {
    
    p <- p +
      coord_cartesian(
        ylim = ylim_range
      )
  }
  
  p
}





plot_4panel_boxplots <- function(
    result_list,
    true_ates = c(
      OR1PS1 = 1,
      OR1PS2 = 1,
      OR2PS1 = 10,
      OR2PS2 = 10
    ),
    ylim_list = list(
      OR1PS1 = c(0.65, 1.40),
      OR1PS2 = c(0.65, 1.40),
      OR2PS1 = c(7.0, 14.0),
      OR2PS2 = c(4.0, 12.0)
    )
) {
  
  p1 <- panel_boxplot_gg(
    result_list$OR1PS1,
    true_ates[["OR1PS1"]],
    "(a) OR1PS1",
    ylim_list$OR1PS1
  )
  
  p2 <- panel_boxplot_gg(
    result_list$OR1PS2,
    true_ates[["OR1PS2"]],
    "(b) OR1PS2",
    ylim_list$OR1PS2
  )
  
  p3 <- panel_boxplot_gg(
    result_list$OR2PS1,
    true_ates[["OR2PS1"]],
    "(c) OR2PS1",
    ylim_list$OR2PS1
  )
  
  p4 <- panel_boxplot_gg(
    result_list$OR2PS2,
    true_ates[["OR2PS2"]],
    "(d) OR2PS2",
    ylim_list$OR2PS2
  )
  
  (p1 | p2) / (p3 | p4)
}


###############################################################################
# A6: SCORE-AUGMENTED GEC ATE ESTIMATOR
#
# Standard GEC balances:
#   1, yhat_t, g(pi_t^{-1})
#
# A6 additionally balances the estimated propensity-score direction:
#
#   h_t(X) = pi_t(X) X
#
# for the arm-specific logistic response model.
#
# Treated arm:
#   pi_1 = P(T=1|X)
#   h_1  = pi_1 X
#
# Control arm:
#   pi_0 = P(T=0|X) = 1-pi_1
#   h_0  = pi_0 X
#
# Because solve_lambda_dual() already adds the intercept and g(pi^{-1}),
# we append h to b_mat rather than changing solve_lambda_dual().
###############################################################################

estimate_ATE_dual_score_aug <- function(
    fold_t1,
    fold_t0,
    entropy = "CE",
    maxit = 1000
) {
  
  entropy <- match.arg(
    entropy,
    choices = c("SL", "EL", "ET", "HD", "CE")
  )
  
  
  ###########################################################################
  # 1. Reconstruct cross-fitted datasets
  ###########################################################################
  
  dat1 <- do.call(
    rbind,
    fold_t1
  )
  
  dat0 <- do.call(
    rbind,
    fold_t0
  )
  
  
  # Keep subjects aligned
  if ("ID" %in% names(dat1)) {
    dat1 <- dat1[
      order(dat1$ID),
      ,
      drop = FALSE
    ]
  }
  
  if ("ID" %in% names(dat0)) {
    dat0 <- dat0[
      order(dat0$ID),
      ,
      drop = FALSE
    ]
  }
  
  
  ###########################################################################
  # 2. Covariates used in the PS model
  ###########################################################################
  
  xvars <- grep(
    "^x\\d+$",
    names(dat1),
    value = TRUE
  )
  
  if (length(xvars) == 0) {
    stop("No x1, ..., xp variables found.")
  }
  
  
  # Logistic PS design matrix INCLUDING intercept
  X1 <- model.matrix(
    reformulate(xvars),
    data = dat1
  )
  
  X0 <- model.matrix(
    reformulate(xvars),
    data = dat0
  )
  
  
  ###########################################################################
  # 3. Arm-specific response indicators and probabilities
  ###########################################################################
  
  D1 <- as.numeric(
    dat1$D == 1
  )
  
  D0 <- as.numeric(
    dat0$D == 0
  )
  
  
  p1 <- as.numeric(
    dat1$pi.hat
  )
  
  p0 <- as.numeric(
    1 - dat0$pi.hat
  )
  
  
  # Numerical protection
  p1 <- pmin(
    pmax(p1, 1e-8),
    1 - 1e-8
  )
  
  p0 <- pmin(
    pmax(p0, 1e-8),
    1 - 1e-8
  )
  
  
  ###########################################################################
  # 4. Propensity-score directions
  #
  # For logistic response probability p:
  #
  # h(O;phi)
  # = [1/(1-p)] dp/dphi
  # = p X
  ###########################################################################
  
  h1 <- X1 * p1
  
  h0 <- X0 * p0
  
  
  colnames(h1) <- paste0(
    "h1_",
    colnames(X1)
  )
  
  colnames(h0) <- paste0(
    "h0_",
    colnames(X0)
  )
  
  
  ###########################################################################
  # 5. A6 balancing functions
  #
  # IMPORTANT:
  # solve_lambda_dual() automatically adds:
  #
  #   intercept
  #   g(pi^{-1})
  #
  # Therefore b_mat contains:
  #
  #   yhat + propensity-score directions
  ###########################################################################
  
  B1_aug <- cbind(
    yhat = dat1$y.hat,
    h1
  )
  
  B0_aug <- cbind(
    yhat = dat0$y.hat,
    h0
  )
  
  
  ###########################################################################
  # 6. Solve calibration problem
  ###########################################################################
  
  fit1 <- solve_lambda_dual(
    b_mat = B1_aug,
    pi_hat = p1,
    D = D1,
    entropy = entropy,
    maxit = maxit
  )
  
  
  fit0 <- solve_lambda_dual(
    b_mat = B0_aug,
    pi_hat = p0,
    D = D0,
    entropy = entropy,
    maxit = maxit
  )
  
  
  ###########################################################################
  # 7. Check for hard failure
  ###########################################################################
  
  if (
    is.null(fit1$weights_all) ||
    is.null(fit0$weights_all) ||
    any(!is.finite(fit1$weights_all)) ||
    any(!is.finite(fit0$weights_all))
  ) {
    
    return(
      list(
        ATE = NA_real_,
        theta1 = NA_real_,
        theta0 = NA_real_,
        fit1 = fit1,
        fit0 = fit0,
        success = FALSE
      )
    )
  }
  
  
  ###########################################################################
  # 8. Closed-form arm means and ATE
  ###########################################################################
  
  N1 <- nrow(dat1)
  N0 <- nrow(dat0)
  
  
  theta1_hat <- sum(
    D1 *
      fit1$weights_all *
      dat1$y
  ) / N1
  
  
  theta0_hat <- sum(
    D0 *
      fit0$weights_all *
      dat0$y
  ) / N0
  
  
  ATE_hat <-
    theta1_hat -
    theta0_hat
  
  
  ###########################################################################
  # 9. Return
  ###########################################################################
  
  list(
    ATE = ATE_hat,
    
    theta1 = theta1_hat,
    theta0 = theta0_hat,
    
    fit1 = fit1,
    fit0 = fit0,
    
    h1 = h1,
    h0 = h0,
    
    B1_aug = B1_aug,
    B0_aug = B0_aug,
    
    success =
      isTRUE(fit1$converged) &&
      isTRUE(fit0$converged)
  )
}










###############################################################################
# A6: SCORE-AUGMENTED ET ESTIMATOR
###############################################################################

estimate_ATE_ET_score <- function(
    fold_t1,
    fold_t0,
    maxit = 1000
) {
  
  ###########################################################################
  # 1. Reconstruct cross-fitted datasets
  ###########################################################################
  
  dat1 <- do.call(
    rbind,
    fold_t1
  )
  
  dat0 <- do.call(
    rbind,
    fold_t0
  )
  
  dat1 <- dat1[
    order(dat1$ID),
    ,
    drop = FALSE
  ]
  
  dat0 <- dat0[
    order(dat0$ID),
    ,
    drop = FALSE
  ]
  
  
  stopifnot(
    identical(
      as.integer(dat1$ID),
      as.integer(dat0$ID)
    )
  )
  
  
  ###########################################################################
  # 2. Fit the SAME working propensity-score model
  #
  # D ~ x1 + x2 + x3 + x4
  ###########################################################################
  
  xvars <- grep(
    "^x\\d+$",
    names(dat1),
    value = TRUE
  )
  
  ps_formula <- reformulate(
    xvars,
    response = "D"
  )
  
  ps_fit <- glm(
    ps_formula,
    data = dat1,
    family = binomial()
  )
  
  X <- model.matrix(
    ps_fit
  )
  
  pi_hat <- as.numeric(
    fitted(ps_fit)
  )
  
  pi_hat <- pmin(
    pmax(pi_hat, 1e-8),
    1 - 1e-8
  )
  
  
  ###########################################################################
  # 3. Treatment-arm quantities
  ###########################################################################
  
  D1 <- as.numeric(
    dat1$D == 1
  )
  
  D0 <- as.numeric(
    dat0$D == 0
  )
  
  p1 <- pi_hat
  
  p0 <- 1 - pi_hat
  
  
  ###########################################################################
  # 4. Estimated propensity-score direction
  #
  # For logistic PS:
  #
  # h(O; phi)
  # =
  # 1/(1-pi) * d pi / d phi
  # =
  # pi X
  #
  # Therefore:
  # treated arm: h1 = p1 X
  # control arm: h0 = p0 X
  ###########################################################################
  
  h1 <- X * p1
  
  h0 <- X * p0
  
  
  colnames(h1) <- paste0(
    "h1_",
    colnames(X)
  )
  
  colnames(h0) <- paste0(
    "h0_",
    colnames(X)
  )
  
  
  ###########################################################################
  # 5. Augmented calibration functions
  #
  # solve_lambda_dual() automatically adds:
  #
  #   intercept
  #   g(pi^{-1})
  #
  # So here we pass:
  #
  #   yhat + propensity-score directions
  ###########################################################################
  
  B1_aug <- cbind(
    yhat = as.numeric(dat1$y.hat),
    h1
  )
  
  B0_aug <- cbind(
    yhat = as.numeric(dat0$y.hat),
    h0
  )
  
  storage.mode(B1_aug) <- "double"
  
  storage.mode(B0_aug) <- "double"
  
  
  ###########################################################################
  # 6. Solve ET calibration
  ###########################################################################
  
  fit1 <- solve_lambda_dual(
    b_mat = B1_aug,
    pi_hat = p1,
    D = D1,
    entropy = "ET",
    maxit = maxit
  )
  
  
  fit0 <- solve_lambda_dual(
    b_mat = B0_aug,
    pi_hat = p0,
    D = D0,
    entropy = "ET",
    maxit = maxit
  )
  
  
  ###########################################################################
  # 7. Hard-failure check
  ###########################################################################
  
  if (
    is.null(fit1$weights_all) ||
    is.null(fit0$weights_all) ||
    any(!is.finite(fit1$weights_all)) ||
    any(!is.finite(fit0$weights_all))
  ) {
    
    return(
      list(
        success = FALSE,
        
        ATE = NA_real_,
        
        theta1 = NA_real_,
        theta0 = NA_real_,
        
        fit1 = fit1,
        fit0 = fit0,
        
        ps_fit = ps_fit
      )
    )
  }
  
  
  ###########################################################################
  # 8. Closed-form ATE
  ###########################################################################
  
  N <- nrow(dat1)
  
  
  theta1_hat <- sum(
    D1 *
      fit1$weights_all *
      dat1$y
  ) / N
  
  
  theta0_hat <- sum(
    D0 *
      fit0$weights_all *
      dat0$y
  ) / N
  
  
  ATE_hat <-
    theta1_hat -
    theta0_hat
  
  
  ###########################################################################
  # 9. Return
  ###########################################################################
  
  list(
    
    success =
      isTRUE(fit1$converged) &&
      isTRUE(fit0$converged),
    
    ATE = ATE_hat,
    
    theta1 = theta1_hat,
    
    theta0 = theta0_hat,
    
    fit1 = fit1,
    
    fit0 = fit0,
    
    ps_fit = ps_fit,
    
    h1 = h1,
    
    h0 = h0,
    
    B1_aug = B1_aug,
    
    B0_aug = B0_aug
  )
}

###############################################################################
# A6: JOINT SANDWICH INFERENCE FOR SCORE-AUGMENTED ET
###############################################################################

ET_score_inference <- function(
    dat,
    yhat1,
    yhat0,
    fit_score,
    true_ATE = NA_real_
) {
  
  ###########################################################################
  # 1. Align data
  ###########################################################################
  
  dat <- dat[
    order(dat$ID),
    ,
    drop = FALSE
  ]
  
  y <- as.numeric(dat$y)
  T <- as.numeric(dat$D)
  
  yhat1 <- as.numeric(yhat1)
  yhat0 <- as.numeric(yhat0)
  
  N <- nrow(dat)
  
  
  ###########################################################################
  # 2. Propensity-score model
  ###########################################################################
  
  ps_fit <- fit_score$ps_fit
  
  X <- model.matrix(ps_fit)
  
  storage.mode(X) <- "double"
  
  phi_hat <- as.numeric(
    coef(ps_fit)
  )
  
  r <- length(phi_hat)
  
  
  ###########################################################################
  # 3. Calibration parameters
  ###########################################################################
  
  lambda1_hat <- as.numeric(
    fit_score$fit1$lambda
  )
  
  lambda0_hat <- as.numeric(
    fit_score$fit0$lambda
  )
  
  q1 <- length(lambda1_hat)
  q0 <- length(lambda0_hat)
  
  
  ###########################################################################
  # 4. Full parameter vector
  #
  # beta =
  # (phi, lambda1, theta1, lambda0, theta0)
  ###########################################################################
  
  beta_hat <- c(
    phi_hat,
    lambda1_hat,
    fit_score$theta1,
    lambda0_hat,
    fit_score$theta0
  )
  
  
  ###########################################################################
  # 5. Parameter indices
  ###########################################################################
  
  idx_phi <- seq_len(r)
  
  idx_lam1 <-
    (r + 1):(r + q1)
  
  idx_theta1 <-
    r + q1 + 1
  
  idx_lam0 <-
    (idx_theta1 + 1):
    (idx_theta1 + q0)
  
  idx_theta0 <-
    idx_theta1 + q0 + 1
  
  
  ###########################################################################
  # 6. Observation-level estimating functions
  ###########################################################################
  
  psi_fun <- function(beta) {
    
    phi <- beta[idx_phi]
    
    lambda1 <- beta[idx_lam1]
    theta1 <- beta[idx_theta1]
    
    lambda0 <- beta[idx_lam0]
    theta0 <- beta[idx_theta0]
    
    
    #########################################################################
    # Propensity score
    #########################################################################
    
    pi <- plogis(
      as.numeric(
        X %*% phi
      )
    )
    
    pi <- pmin(
      pmax(pi, 1e-8),
      1 - 1e-8
    )
    
    
    p1 <- pi
    p0 <- 1 - pi
    
    
    #########################################################################
    # Estimated PS directions
    #
    # Treated response model:
    # h1 = pi X
    #
    # Control response model:
    # using +p0 X gives the same calibration span as -p0 X.
    # This matches estimate_ATE_ET_score().
    #########################################################################
    
    h1 <- X * p1
    
    h0 <- X * p0
    
    
    #########################################################################
    # ET debiasing coordinates
    #
    # g(p^{-1}) = log(p^{-1})
    #########################################################################
    
    g1 <- log(1 / p1)
    
    g0 <- log(1 / p0)
    
    
    #########################################################################
    # Score-augmented calibration vectors
    #########################################################################
    
    S1 <- cbind(
      intercept = 1,
      yhat = yhat1,
      h1,
      g_pi = g1
    )
    
    S0 <- cbind(
      intercept = 1,
      yhat = yhat0,
      h0,
      g_pi = g0
    )
    
    storage.mode(S1) <- "double"
    storage.mode(S0) <- "double"
    
    
    #########################################################################
    # ET weights
    #
    # omega(eta) = exp(eta)
    #########################################################################
    
    eta1 <- as.numeric(
      S1 %*% lambda1
    )
    
    eta0 <- as.numeric(
      S0 %*% lambda0
    )
    
    
    w1 <- exp(eta1)
    
    w0 <- exp(eta0)
    
    
    D1 <- T
    D0 <- 1 - T
    
    
    #########################################################################
    # PS estimating equation
    #########################################################################
    
    Psi_phi <-
      X *
      as.numeric(
        T - pi
      )
    
    
    #########################################################################
    # Calibration estimating equations
    #########################################################################
    
    Psi_lam1 <-
      S1 *
      as.numeric(
        D1 * w1 - 1
      )
    
    
    Psi_lam0 <-
      S0 *
      as.numeric(
        D0 * w0 - 1
      )
    
    
    #########################################################################
    # Outcome equations
    #########################################################################
    
    Psi_theta1 <-
      D1 *
      w1 *
      (y - theta1)
    
    
    Psi_theta0 <-
      D0 *
      w0 *
      (y - theta0)
    
    
    #########################################################################
    # Stack all estimating equations
    #########################################################################
    
    cbind(
      Psi_phi,
      Psi_lam1,
      theta1 = Psi_theta1,
      Psi_lam0,
      theta0 = Psi_theta0
    )
  }
  
  
  ###########################################################################
  # 7. Evaluate estimating functions
  ###########################################################################
  
  Psi_hat <- psi_fun(
    beta_hat
  )
  
  if (any(!is.finite(Psi_hat))) {
    
    return(
      list(
        success = FALSE,
        ATE = fit_score$ATE,
        SE = NA_real_,
        coverage = NA_real_
      )
    )
  }
  
  
  ###########################################################################
  # 8. Meat
  ###########################################################################
  
  B <- crossprod(
    Psi_hat
  ) / N
  
  
  ###########################################################################
  # 9. Numerical Jacobian
  ###########################################################################
  
  J <- tryCatch(
    
    numDeriv::jacobian(
      
      func = function(beta) {
        
        colMeans(
          psi_fun(beta)
        )
      },
      
      x = beta_hat,
      
      method = "Richardson"
    ),
    
    error = function(e) NULL
  )
  
  
  if (
    is.null(J) ||
    any(!is.finite(J))
  ) {
    
    return(
      list(
        success = FALSE,
        ATE = fit_score$ATE,
        SE = NA_real_,
        coverage = NA_real_
      )
    )
  }
  
  
  ###########################################################################
  # Bread = negative Jacobian
  ###########################################################################
  
  A <- -J
  
  
  A_inv <- tryCatch(
    solve(A),
    error = function(e) NULL
  )
  
  
  if (
    is.null(A_inv) ||
    any(!is.finite(A_inv))
  ) {
    
    return(
      list(
        success = FALSE,
        ATE = fit_score$ATE,
        SE = NA_real_,
        coverage = NA_real_
      )
    )
  }
  
  
  ###########################################################################
  # 10. Sandwich covariance
  ###########################################################################
  
  V <- (
    A_inv %*%
      B %*%
      t(A_inv)
  ) / N
  
  
  ###########################################################################
  # 11. ATE contrast theta1 - theta0
  ###########################################################################
  
  contrast <- rep(
    0,
    length(beta_hat)
  )
  
  contrast[idx_theta1] <- 1
  
  contrast[idx_theta0] <- -1
  
  
  var_ATE <- as.numeric(
    t(contrast) %*%
      V %*%
      contrast
  )
  
  
  if (
    !is.finite(var_ATE) ||
    var_ATE < 0
  ) {
    
    return(
      list(
        success = FALSE,
        ATE = fit_score$ATE,
        SE = NA_real_,
        coverage = NA_real_
      )
    )
  }
  
  
  SE <- sqrt(var_ATE)
  
  ATE <- fit_score$ATE
  
  
  ###########################################################################
  # 12. 95% Wald CI
  ###########################################################################
  
  lower <-
    ATE -
    qnorm(0.975) * SE
  
  upper <-
    ATE +
    qnorm(0.975) * SE
  
  
  coverage <- if (
    is.finite(true_ATE)
  ) {
    
    as.integer(
      lower <= true_ATE &&
        true_ATE <= upper
    )
    
  } else {
    
    NA_integer_
  }
  
  
  ###########################################################################
  # 13. Weight diagnostics
  ###########################################################################
  
  w1_obs <-
    fit_score$fit1$weights_all[
      T == 1
    ]
  
  w0_obs <-
    fit_score$fit0$weights_all[
      T == 0
    ]
  
  
  ESS1 <-
    sum(w1_obs)^2 /
    sum(w1_obs^2)
  
  
  ESS0 <-
    sum(w0_obs)^2 /
    sum(w0_obs^2)
  
  
  ###########################################################################
  # 14. Return
  ###########################################################################
  
  list(
    
    success = TRUE,
    
    ATE = ATE,
    
    SE = SE,
    
    variance = var_ATE,
    
    lower = lower,
    
    upper = upper,
    
    coverage = coverage,
    
    ESS1 = ESS1,
    
    ESS0 = ESS0,
    
    max_weight1 =
      max(w1_obs),
    
    max_weight0 =
      max(w0_obs),
    
    A = A,
    
    B = B,
    
    V = V
  )
}

###############################################################################
# A6: ONE REPLICATION
#
# OR2PS1
#
# Compare:
#   1. Estimated-PS AIPW-GAM
#   2. Standard ET
#   3. Score-augmented ET
###############################################################################

run_one_A6_ET <- function(
    rep_id,
    n = 1000,
    p = 4,
    K = 4
) {
  
  true_ATE <- 10
  
  
  ###########################################################################
  # 1. Generate OR2PS1 data
  ###########################################################################
  
  sim <- generate_ate_data(
    n = n,
    p = p,
    outcome_model = 2,
    ps_model = 1,
    seed = rep_id + 20242022
  )
  
  
  dat <- sim$data
  
  
  if (!"ID" %in% names(dat)) {
    
    dat$ID <-
      seq_len(
        nrow(dat)
      )
  }
  
  
  D <- dat$D
  
  
  ###########################################################################
  # 2. GAM cross-fitting
  ###########################################################################
  
  fold_T1 <- k_fold_function_gam(
    df = dat,
    K = K,
    idx_S = which(D == 1),
    idx_U0 = which(D == 0),
    seed = rep_id + 20242022
  )
  
  
  fold_T0 <- k_fold_function_gam(
    df = dat,
    K = K,
    idx_S = which(D == 0),
    idx_U0 = which(D == 1),
    seed = rep_id + 20242022
  )
  
  
  ###########################################################################
  # Align prediction datasets
  ###########################################################################
  
  dat1 <- do.call(
    rbind,
    fold_T1
  )
  
  dat0 <- do.call(
    rbind,
    fold_T0
  )
  
  
  dat1 <- dat1[
    order(dat1$ID),
    ,
    drop = FALSE
  ]
  
  dat0 <- dat0[
    order(dat0$ID),
    ,
    drop = FALSE
  ]
  
  dat <- dat[
    order(dat$ID),
    ,
    drop = FALSE
  ]
  
  
  yhat1 <-
    as.numeric(
      dat1$y.hat
    )
  
  yhat0 <-
    as.numeric(
      dat0$y.hat
    )
  
  
  ###########################################################################
  # 3. Estimated-propensity AIPW
  ###########################################################################
  
  aipw <- tryCatch(
    
    estimate_aipw_inference(
      fold_T1 = fold_T1,
      fold_T0 = fold_T0,
      true_ATE = true_ATE
    ),
    
    error = function(e) {
      
      message(
        "AIPW error rep ",
        rep_id,
        ": ",
        e$message
      )
      
      NULL
    }
  )
  
  
  ###########################################################################
  # 4. Standard ET point estimate
  ###########################################################################
  
  ET_fit <- tryCatch(
    
    estimate_ATE_dual(
      fold_t1 = fold_T1,
      fold_t0 = fold_T0,
      entropy = "ET"
    ),
    
    error = function(e) {
      
      message(
        "ET error rep ",
        rep_id,
        ": ",
        e$message
      )
      
      NULL
    }
  )
  
  
  ###########################################################################
  # Standard ET inference
  ###########################################################################
  
  ET_inf <- NULL
  
  
  if (
    !is.null(ET_fit) &&
    is.finite(ET_fit$ATE) &&
    !is.null(ET_fit$fit1) &&
    !is.null(ET_fit$fit0)
  ) {
    
    xvars <- grep(
      "^x\\d+$",
      names(dat),
      value = TRUE
    )
    
    
    ps_fit <- glm(
      reformulate(
        xvars,
        response = "D"
      ),
      data = dat,
      family = binomial()
    )
    
    
    ET_inf <- tryCatch(
      
      gec_ate_inference(
        dat = dat,
        yhat1 = yhat1,
        yhat0 = yhat0,
        fit1 = ET_fit$fit1,
        fit0 = ET_fit$fit0,
        ps_fit = ps_fit,
        entropy = "ET",
        true_ATE = true_ATE,
        numerical_jacobian = FALSE
      ),
      
      error = function(e) {
        
        message(
          "ET inference error rep ",
          rep_id,
          ": ",
          e$message
        )
        
        NULL
      }
    )
  }
  
  
  ###########################################################################
  # 5. Score-augmented ET
  ###########################################################################
  
  ET_score_fit <- tryCatch(
    
    estimate_ATE_ET_score(
      fold_t1 = fold_T1,
      fold_t0 = fold_T0
    ),
    
    error = function(e) {
      
      message(
        "ET-score error rep ",
        rep_id,
        ": ",
        e$message
      )
      
      NULL
    }
  )
  
  
  ###########################################################################
  # ET + score inference
  ###########################################################################
  
  ET_score_inf <- NULL
  
  
  if (
    !is.null(ET_score_fit) &&
    isTRUE(ET_score_fit$success)
  ) {
    
    ET_score_inf <- tryCatch(
      
      ET_score_inference(
        dat = dat,
        yhat1 = yhat1,
        yhat0 = yhat0,
        fit_score = ET_score_fit,
        true_ATE = true_ATE
      ),
      
      error = function(e) {
        
        message(
          "ET-score inference error rep ",
          rep_id,
          ": ",
          e$message
        )
        
        NULL
      }
    )
  }
  
  
  ###########################################################################
  # 6. Return one row
  ###########################################################################
  
  data.frame(
    
    replication = rep_id,
    
    
    # AIPW
    AIPW =
      if (!is.null(aipw))
        aipw$ATE
    else NA_real_,
    
    AIPW_SE =
      if (!is.null(aipw))
        aipw$SE
    else NA_real_,
    
    AIPW_COV =
      if (!is.null(aipw))
        aipw$coverage
    else NA_real_,
    
    
    # Standard ET
    ET =
      if (!is.null(ET_inf))
        ET_inf$ATE
    else if (!is.null(ET_fit))
      ET_fit$ATE
    else NA_real_,
    
    ET_SE =
      if (!is.null(ET_inf))
        ET_inf$SE_analytic
    else NA_real_,
    
    ET_COV =
      if (!is.null(ET_inf))
        ET_inf$cover_analytic
    else NA_real_,
    
    
    # ET + PS score
    ET_SCORE =
      if (!is.null(ET_score_inf))
        ET_score_inf$ATE
    else if (!is.null(ET_score_fit))
      ET_score_fit$ATE
    else NA_real_,
    
    ET_SCORE_SE =
      if (!is.null(ET_score_inf))
        ET_score_inf$SE
    else NA_real_,
    
    ET_SCORE_COV =
      if (!is.null(ET_score_inf))
        ET_score_inf$coverage
    else NA_real_,
    
    
    ET_SUCCESS =
      as.integer(
        !is.null(ET_fit) &&
          isTRUE(ET_fit$fit1$converged) &&
          isTRUE(ET_fit$fit0$converged)
      ),
    
    
    ET_SCORE_SUCCESS =
      as.integer(
        !is.null(ET_score_inf) &&
          isTRUE(ET_score_inf$success)
      ),
    
    stringsAsFactors = FALSE
  )
}
