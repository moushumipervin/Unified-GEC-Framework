###############################################################################
# Title: Calibration, Weighting, and Transport Estimation Utilities
# Project: Entropy Calibration / Partially Observed Data Analysis
# Author: Saif Hasan
#
# Description:
# This script collects utility functions used for calibration weighting,
# CBPS, EBCW, AIPW transport estimation, cross-fitting, entropy tilting,
# and bootstrap inference.
#
# Notes:
# - Function bodies are kept unchanged.
# - The script is organized for readability and reproducibility.
# - Required objects such as d4, X, Y, T, exp_dat, idx_nsw, xvars, etc.
#   should be created in the analysis script before calling these functions.
###############################################################################

###############################################################################
# 0. Required Packages
###############################################################################
required_packages <- c(
  "MASS",
  "CVXR",
  "CBPS",
  "nnet",
  "mgcv"
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
  library(MASS)
  library(CVXR)
  library(CBPS)
  library(nnet)
  library(mgcv)
})

###############################################################################
# 1. Calibration Weighting Functions
###############################################################################

cal_tilt_mean_multi <- function(X, target, y, maxit = 2000) {
  X <- as.matrix(X)
  target <- as.numeric(target)
  stopifnot(ncol(X) == length(target))
  stopifnot(length(y) == nrow(X))
  
  obj <- function(lambda) {
    eta <- drop(X %*% lambda)
    m <- max(eta)
    log(sum(exp(eta - m))) + m - sum(target * lambda)
  }
  
  grad <- function(lambda) {
    eta <- drop(X %*% lambda)
    m <- max(eta)
    w <- exp(eta - m)
    w <- w / sum(w)
    drop(crossprod(X, w)) - target
  }
  
  fit <- optim(
    par = rep(0, ncol(X)),
    fn = obj, gr = grad,
    method = "BFGS",
    control = list(maxit = maxit, reltol = 1e-10)
  )
  
  eta <- drop(X %*% fit$par)
  m <- max(eta)
  w <- exp(eta - m)
  w <- w / sum(w)
  
  list(
    EY = sum(w * y),
    w  = w,
    lambda = fit$par,
    converged = (fit$convergence == 0),
    opt = fit
  )
}

cal_tilt_mean_multi_cvxr <- function(X, target, y, solver) {
  X <- as.matrix(X)
  target <- as.numeric(target)
  y <- as.numeric(y)
  n <- nrow(X)
  p <- ncol(X)
  
  w <- CVXR::Variable(n)
  
  obj <- CVXR::Minimize(CVXR::sum_entries(-CVXR::entr(w) - w))
  
  constr <- list(
    w >= 0,
    CVXR::sum_entries(w) == 1,
    t(X) %*% w == target
  )
  
  prob <- CVXR::Problem(obj, constr)
  sol  <- CVXR::solve(prob, solver = solver)
  
  w_hat <- as.numeric(sol$getValue(w))
  mom_err <- drop(t(X) %*% w_hat - target)
  
  list(
    EY = sum(w_hat * y),
    w  = w_hat,
    status = sol$status,
    max_moment_error = max(abs(mom_err)),
    sum_w = sum(w_hat),
    min_w = min(w_hat)
  )
}

cal_tilt_mean_multi_se <- function(X, target, y, fit) {
  w  <- fit$w
  mu <- fit$EY
  
  X <- as.matrix(X)
  target <- as.numeric(target)
  
  Xc <- sweep(X, 2, target, "-")
  A <- t(Xc) %*% (Xc * w)
  b <- drop(t(Xc) %*% ((y - mu) * w))
  
  beta <- tryCatch(
    qr.solve(A, b),
    error = function(e) MASS::ginv(A) %*% b
  )
  beta <- drop(beta)
  
  psi <- w * ((y - mu) - drop(Xc %*% beta))
  
  var_mu <- sum(psi^2)
  se_mu  <- sqrt(var_mu)
  
  list(mu = mu, se = se_mu, var = var_mu, beta = beta, psi = psi)
}

cal_tilt_mean_multi_cvxr <- function(X, target, y, solver ) {
  solver <- match.arg(solver)
  X <- as.matrix(X)
  target <- as.numeric(target)
  y <- as.numeric(y)
  
  n <- nrow(X)
  p <- ncol(X)
  stopifnot(length(target) == p, length(y) == n)
  
  w <- CVXR::Variable(n)
  
  obj <- CVXR::Maximize(CVXR::sum(CVXR::entr(w)))
  
  constr <- list(
    w >= 0,
    CVXR::sum(w) == 1,
    t(X) %*% w == target
  )
  
  prob <- CVXR::Problem(obj, constr)
  sol  <- CVXR::solve(prob, solver = solver)
  
  w_hat <- as.numeric(sol$getValue(w))
  w_hat[w_hat < 0] <- 0
  w_hat <- w_hat / sum(w_hat)
  
  list(
    EY = sum(w_hat * y),
    w  = w_hat,
    converged = (sol$status %in% c("optimal", "optimal_inaccurate")),
    status = sol$status,
    dual_sumw = tryCatch(sol$getDualValue(constr[[2]]), error = function(e) NA),
    dual_X    = tryCatch(sol$getDualValue(constr[[3]]), error = function(e) NA)
  )
}

cal_tilt_mean_multi_cvxr_hellinger <- function(X, target, y, solver = "ECOS") {
  X <- as.matrix(X)
  target <- as.numeric(target)
  y <- as.numeric(y)
  
  n <- nrow(X)
  p <- ncol(X)
  stopifnot(length(target) == p, length(y) == n)
  
  w <- CVXR::Variable(n)
  obj <- CVXR::Minimize(-2 * CVXR::sum_entries(sqrt(w)))
  
  constr <- list(
    w >= 0,
    CVXR::sum_entries(w) == 1,
    t(X) %*% w == target
  )
  
  prob <- CVXR::Problem(obj, constr)
  sol  <- CVXR::solve(prob, solver = solver)
  
  w_hat <- as.numeric(sol$getValue(w))
  list(EY = sum(w_hat * y), w = w_hat, status = sol$status)
}

###############################################################################
# 2. Fixed-Weight IPW Variance
###############################################################################

IPW_var_fixed <- function(d4, ps0, cat1, cat0) {
  y1 <- d4$re78[d4$G == cat1]
  y0 <- d4$re78[d4$G == cat0]
  w0 <- ps0[d4$G == cat0]
  
  n1 <- length(y1)
  
  mu1 <- mean(y1)
  mu0 <- sum(w0 * y0) / sum(w0)
  
  var_mu1 <- var(y1) / n1
  var_mu0 <- sum((w0^2) * (y0 - mu0)^2) / (sum(w0)^2)
  
  var_ate <- var_mu1 + var_mu0
  se_ate  <- sqrt(var_ate)
  
  list(
    ate_hat = mu1 - mu0, se = se_ate, var = var_ate,
    mu1 = mu1, mu0 = mu0, var_mu1 = var_mu1, var_mu0 = var_mu0
  )
}

###############################################################################
# 3. CBPS Functions
###############################################################################

cbps_function <- function(d4, xvars, g_target = c("1","2"), g_donor ) {
  dat <- d4[d4$G %in% c(g_target, g_donor), ]
  dat$S <- as.integer(dat$G %in% g_target)
  
  fml <- as.formula(paste("S ~", paste(xvars, collapse = " + ")))
  fit <- CBPS::CBPS(fml, data = dat, verbose = FALSE)
  ehat <- drop(fit$fitted.values)
  
  w <- rep(NA_real_, nrow(dat))
  w[dat$S == 0] <- ehat[dat$S == 0] / (1 - ehat[dat$S == 0])
  
  n_target <- sum(dat$S == 1)
  w[dat$S == 0] <- w[dat$S == 0] * (n_target / sum(w[dat$S == 0]))
  
  list(dat = dat, w = w, ehat = ehat, fit = fit)
}

CBPS_nsw_with_psid_cps <- function(d4, xvars, g_donor) {
  mu1 <- mean(d4$re78[d4$G == "1"])
  
  cb <- cbps_function(d4, xvars, g_target = c("1","2"), g_donor)
  dat <- cb$dat
  w   <- cb$w
  
  mu0 <- sum(w[dat$S == 0] * dat$re78[dat$S == 0]) / sum(w[dat$S == 0])
  
  c(mu1 = mu1, mu0 = mu0, ATE = mu1 - mu0)
}

CBPS_nsw_with_donor_fixedvar <- function(d4, xvars, g_donor) {
  y1 <- d4$re78[d4$G == "1"]
  mu1 <- mean(y1)
  var_mu1 <- var(y1) / length(y1)
  
  cb <- cbps_function(d4, xvars, g_target = c("1","2"), g_donor = g_donor)
  dat <- cb$dat
  y0  <- dat$re78[dat$S == 0]
  w0  <- cb$w[dat$S == 0]
  
  mu0 <- sum(w0 * y0) / sum(w0)
  
  wtil <- w0 / sum(w0)
  var_mu0 <- sum(wtil^2 * (y0 - mu0)^2)
  
  ate <- mu1 - mu0
  var_ate <- var_mu1 + var_mu0
  se_ate <- sqrt(var_ate)
  
  return(list(mu1 = mu1, mu0 = mu0, ATE = ate, w = w0, SE = se_ate, Var = var_ate))
}

###############################################################################
# 4. EBCW Function
###############################################################################

EBCW_function <- function(cat0, mu1, target, y1) {
  res0 <- cal_tilt_mean_multi(X[T == cat0, , drop = FALSE], target, Y[T == cat0])
  ATE <- mu1 - res0$EY
  var_y1 <- var(y1) / length(y1)
  se0 <- cal_tilt_mean_multi_se(X[T == cat0, , drop = FALSE], target, Y[T == cat0], res0)
  SE_ATE <- sqrt(var_y1 + se0$var)
  return(list(ATE = ATE, w = res0$w, SE = SE_ATE))
}

###############################################################################
# 5. AIPW Transport Functions
###############################################################################

aipw_transport_general <- function(d4, xvars,
                                   cat1 = 1, target_cats = c(1,2),
                                   cat0 = 3,
                                   w0,
                                   yvar = "re78", gvar = "G", tvar = "treat",
                                   outcome_model = c("lm","gam"),
                                   gam_smooth = TRUE,
                                   maxit_gam = NULL) {
  
  outcome_model <- match.arg(outcome_model)
  
  G <- d4[[gvar]]
  Y <- d4[[yvar]]
  
  idx1      <- (G == cat1)
  idxTarget <- (G %in% target_cats)
  idx0      <- (G == cat0)
  
  stopifnot(sum(idx1) > 1, sum(idxTarget) > 1, sum(idx0) > 1)
  stopifnot(length(w0) == nrow(d4))
  
  donor  <- d4[idx0, , drop = FALSE]
  target <- d4[idxTarget, , drop = FALSE]
  
  if (outcome_model == "lm") {
    fml <- stats::reformulate(xvars, response = yvar)
    m0  <- stats::lm(fml, data = donor)
  } else {
    stopifnot(requireNamespace("mgcv", quietly = TRUE))
    df_donor <- donor[, c(yvar, xvars), drop = FALSE]
    
    terms <- vapply(xvars, function(v) {
      if (gam_smooth && is.numeric(df_donor[[v]])) paste0("s(", v, ", bs='cs')")
      else v
    }, character(1))
    
    fml <- as.formula(paste(yvar, "~", paste(terms, collapse = " + ")))
    m0  <- mgcv::gam(fml, data = df_donor, method = "REML")
  }
  
  m0_0      <- as.numeric(stats::predict(m0, newdata = donor))
  m0_target <- as.numeric(stats::predict(m0, newdata = target))
  
  mu1 <- mean(Y[idx1])
  m0bar_target <- mean(m0_target)
  
  r0 <- Y[idx0] - m0_0
  sw <- sum(w0[idx0])
  rbar_w <- sum(w0[idx0] * r0) / sw
  
  mu0_aipw <- m0bar_target + rbar_w
  ate <- mu1 - mu0_aipw
  
  n1   <- sum(idx1)
  nT   <- sum(idxTarget)
  
  IF <- numeric(nrow(d4))
  IF[idx1] <- (Y[idx1] - mu1) / n1
  IF[idxTarget] <- IF[idxTarget] - (m0_target - m0bar_target) / nT
  IF[idx0] <- IF[idx0] - (w0[idx0] * (r0 - rbar_w)) / sw
  
  var_hat <- sum(IF^2)
  se_hat  <- sqrt(var_hat)
  ci <- ate + c(-1, 1) * 1.96 * se_hat
  
  list(
    ATE = ate, se = se_hat, var = var_hat, ci = ci,
    mu1 = mu1, mu0_aipw = mu0_aipw,
    cat1 = cat1, cat0 = cat0, target_cats = target_cats,
    m0 = m0
  )
}

aipw_transport_gam <- function(d4, xvars,
                               cat1 = 1, target_cats = c(1,2),
                               cat0 = 3,
                               w0,
                               yvar = "re78", gvar = "G", tvar = "treat",
                               cont_vars = NULL,
                               cat_vars  = NULL,
                               bs = "cs",
                               k  = -1,
                               method = "REML",
                               donor_controls_only = FALSE) {
  
  stopifnot(requireNamespace("mgcv", quietly = TRUE))
  
  G <- d4[[gvar]]
  Y <- d4[[yvar]]
  T <- d4[[tvar]]
  
  idx1      <- (G == cat1)
  idxTarget <- (G %in% target_cats)
  idx0_all  <- (G == cat0)
  
  stopifnot(sum(idx1) > 1, sum(idxTarget) > 1, sum(idx0_all) > 1)
  stopifnot(length(w0) == nrow(d4))
  
  if (donor_controls_only) {
    idx0 <- idx0_all & (T == 0)
    stopifnot(sum(idx0) > 1)
  } else {
    idx0 <- idx0_all
  }
  
  donor  <- d4[idx0, , drop = FALSE]
  target <- d4[idxTarget, , drop = FALSE]
  
  if (is.null(cont_vars) && is.null(cat_vars)) {
    cont_vars <- xvars[vapply(donor[, xvars, drop = FALSE], is.numeric, logical(1)) |
                         vapply(donor[, xvars, drop = FALSE], is.integer, logical(1))]
    cat_vars  <- setdiff(xvars, cont_vars)
  } else {
    if (is.null(cont_vars)) cont_vars <- setdiff(xvars, cat_vars)
    if (is.null(cat_vars))  cat_vars  <- setdiff(xvars, cont_vars)
  }
  
  for (v in cat_vars) {
    donor[[v]]  <- as.factor(donor[[v]])
    target[[v]] <- as.factor(target[[v]])
    lv <- union(levels(donor[[v]]), levels(target[[v]]))
    donor[[v]]  <- factor(donor[[v]],  levels = lv)
    target[[v]] <- factor(target[[v]], levels = lv)
  }
  
  smooth_terms <- if (length(cont_vars) > 0) {
    if (k == -1) {
      paste0("s(", cont_vars, ", bs='", bs, "')")
    } else {
      paste0("s(", cont_vars, ", bs='", bs, "', k=", k, ")")
    }
  } else character(0)
  
  linear_terms <- cat_vars
  rhs <- paste(c(smooth_terms, linear_terms), collapse = " + ")
  fml <- as.formula(paste(yvar, "~", rhs))
  
  df_donor <- donor[, c(yvar, xvars), drop = FALSE]
  m0 <- mgcv::gam(fml, data = df_donor, method = method)
  
  m0_0      <- as.numeric(stats::predict(m0, newdata = donor))
  m0_target <- as.numeric(stats::predict(m0, newdata = target))
  
  mu1 <- mean(Y[idx1])
  m0bar_target <- mean(m0_target)
  
  r0 <- Y[idx0] - m0_0
  sw <- sum(w0[idx0])
  rbar_w <- sum(w0[idx0] * r0) / sw
  
  mu0_aipw <- m0bar_target + rbar_w
  ate <- mu1 - mu0_aipw
  
  n1 <- sum(idx1)
  nT <- sum(idxTarget)
  
  IF <- numeric(nrow(d4))
  IF[idx1]      <- (Y[idx1] - mu1) / n1
  IF[idxTarget] <- IF[idxTarget] - (m0_target - m0bar_target) / nT
  IF[idx0]      <- IF[idx0] - (w0[idx0] * (r0 - rbar_w)) / sw
  
  var_hat <- sum(IF^2)
  se_hat  <- sqrt(var_hat)
  ci <- ate + c(-1, 1) * 1.96 * se_hat
  
  list(
    ATE = ate, se = se_hat, var = var_hat, ci = ci,
    mu1 = mu1, mu0_aipw = mu0_aipw,
    cat1 = cat1, cat0 = cat0, target_cats = target_cats,
    cont_vars = cont_vars, cat_vars = cat_vars,
    formula = fml,
    m0 = m0,
    y_hat_target = m0_target,
    y_hat_ctrl = m0_0,
    w = w0[idx0]
  )
}

###############################################################################
# 6. Cross-Fitting Functions
###############################################################################

build_SU_folds <- function(df, K, id_col = "ID", idx_S, idx_U0, seed ) {
  set.seed(seed)
  if (!id_col %in% names(df)) df[[id_col]] <- seq_len(nrow(df))
  id <- df[[id_col]]
  
  if (length(idx_S) < K || length(idx_U0) < K)
    stop("Need at least K S and K U0 observations.")
  
  split_K <- function(idx, K) {
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
    
    folds[[k]] <- list(S_k_idx = S_k_idx, U_k_idx = U_k_idx,
                       S_k_ids = id[S_k_idx], U_k_ids = id[U_k_idx])
    
    fold_id_S[S_k_idx] <- k
    fold_id_U[U0_k_idx] <- k
  }
  
  list(folds = folds, fold_id_S = fold_id_S, fold_id_U = fold_id_U)
}

k_fold_predict_yhat <- function(df, K, idx_S, idx_U0, seed,
                                y = "re78",
                                drop = c("data_id","treat","ID","G")) {
  
  df <- df
  df$ID <- seq_len(nrow(df))
  
  res <- build_SU_folds(df, K, id_col = "ID", idx_S = idx_S, idx_U0 = idx_U0, seed = seed)
  xvars <- setdiff(names(df), c(drop, y))
  out_list <- vector("list", K)
  
  for (k in 1:K) {
    S_k_idx <- res$folds[[k]]$S_k_idx
    U_k_idx <- res$folds[[k]]$U_k_idx
    
    train_idx <- setdiff(idx_S, S_k_idx)
    train_df  <- df[train_idx, , drop = FALSE]
    pred_df   <- df[U_k_idx,   , drop = FALSE]
    
    fml <- reformulate(xvars, response = y)
    fit <- lm(fml, data = train_df)
    
    pred_df$y.hat <- predict(fit, newdata = pred_df)
    pred_df$fold  <- k
    out_list[[k]] <- pred_df
  }
  
  out_list
}

k_fold_predict_yhat_gam <- function(df, K, idx_S, idx_U0, seed,
                                    y = "re78",
                                    drop = c("data_id","treat","ID","G"),
                                    bs = "cs",
                                    method = "REML") {
  stopifnot(requireNamespace("mgcv", quietly = TRUE))
  
  df <- df
  df$ID <- seq_len(nrow(df))
  
  res <- build_SU_folds(df, K, id_col = "ID", idx_S = idx_S, idx_U0 = idx_U0, seed = seed)
  
  base_xvars <- setdiff(names(df), c(drop, y))
  out_list <- vector("list", K)
  
  for (k in 1:K) {
    S_k_idx <- res$folds[[k]]$S_k_idx
    U_k_idx <- res$folds[[k]]$U_k_idx
    
    train_idx <- setdiff(idx_S, S_k_idx)
    train_df  <- df[train_idx, , drop = FALSE]
    pred_df   <- df[U_k_idx,   , drop = FALSE]
    
    xvars <- base_xvars
    x_keep <- character(0)
    
    for (v in xvars) {
      if (is.character(train_df[[v]])) train_df[[v]] <- factor(train_df[[v]])
      if (is.character(pred_df[[v]]))  pred_df[[v]]  <- factor(pred_df[[v]])
      
      if (is.numeric(train_df[[v]])) {
        u <- unique(train_df[[v]][!is.na(train_df[[v]])])
        if (length(u) <= 2) {
          train_df[[v]] <- factor(train_df[[v]])
          pred_df[[v]]  <- factor(pred_df[[v]])
        }
      }
      
      if (is.factor(train_df[[v]])) {
        pred_df[[v]] <- factor(pred_df[[v]], levels = levels(train_df[[v]]))
      }
      
      if (is.factor(train_df[[v]])) {
        if (nlevels(droplevels(train_df[[v]])) >= 2) x_keep <- c(x_keep, v)
      } else if (is.numeric(train_df[[v]])) {
        u <- unique(train_df[[v]][!is.na(train_df[[v]])])
        if (length(u) >= 2) x_keep <- c(x_keep, v)
      } else {
        u <- unique(train_df[[v]][!is.na(train_df[[v]])])
        if (length(u) >= 2) x_keep <- c(x_keep, v)
      }
    }
    
    terms <- vapply(x_keep, function(v) {
      if (is.numeric(train_df[[v]])) {
        u <- unique(train_df[[v]][!is.na(train_df[[v]])])
        if (length(u) >= 5) paste0("s(", v, ", bs='", bs, "')") else v
      } else {
        v
      }
    }, character(1))
    
    if (length(terms) == 0) {
      fml <- as.formula(paste(y, "~ 1"))
    } else {
      fml <- as.formula(paste(y, "~", paste(terms, collapse = " + ")))
    }
    
    fit <- mgcv::gam(fml, data = train_df, method = method)
    
    pred_df$y.hat <- as.numeric(predict(fit, newdata = pred_df, type = "response"))
    pred_df$fold  <- k
    out_list[[k]] <- pred_df
  }
  
  out_list
}

###############################################################################
# 7. Variance Utility Functions for AIPW
###############################################################################

aipw_var_correct_ordering <- function(data_T1, data_T0, w1, w3, id = "ID") {
  psid0 <- data_T0[data_T0$G == 3, ]
  nsw0  <- data_T0[data_T0$G %in% c(1,2), ]
  idx3  <- which(data_T0$G == 3)
  
  m0_psid <- psid0$y.hat
  m0_nsw  <- nsw0$y.hat
  w3_3    <- w3[idx3]
  r0      <- psid0$re78 - m0_psid
  
  mu0_aipw <- mean(m0_nsw) + sum(w3_3 * r0) / sum(w3_3)
  
  g1   <- data_T1[data_T1$G == 1, ]
  nsw1 <- data_T1[data_T1$G %in% c(1,2), ]
  idx1 <- which(data_T1$G == 1)
  
  m1_g1  <- g1$y.hat
  m1_nsw <- nsw1$y.hat
  w1_1   <- w1[idx1]
  r1     <- g1$re78 - m1_g1
  
  mu1_aipw <- mean(m1_nsw) + sum(w1_1 * r1) / sum(w1_1)
  ate <- mu1_aipw - mu0_aipw
  
  sw1 <- sum(w1_1)
  r1bar_w <- sum(w1_1 * r1) / sw1
  IF_ratio1 <- (w1_1 * (r1 - r1bar_w)) / sw1
  
  sw3 <- sum(w3_3)
  r0bar_w <- sum(w3_3 * r0) / sw3
  IF_ratio0 <- (w3_3 * (r0 - r0bar_w)) / sw3
  
  nsw1_tbl <- nsw1[, c(id, "y.hat")]
  names(nsw1_tbl) <- c("ID", "m1hat")
  nsw0_tbl <- nsw0[, c(id, "y.hat")]
  names(nsw0_tbl) <- c("ID", "m0hat")
  
  nsw_align <- merge(nsw1_tbl, nsw0_tbl, by = "ID", all = FALSE, sort = FALSE)
  
  nT <- nrow(nsw_align)
  m1bar_T <- mean(nsw_align$m1hat)
  m0bar_T <- mean(nsw_align$m0hat)
  
  IF_target <- (nsw_align$m1hat - m1bar_T) / nT - (nsw_align$m0hat - m0bar_T) / nT
  
  var_hat <- sum(IF_target^2) + sum(IF_ratio1^2) + sum(IF_ratio0^2)
  
  se <- sqrt(var_hat)
  ci <- ate + c(-1, 1) * 1.96 * se
  
  list(ATE = ate, mu1 = mu1_aipw, mu0 = mu0_aipw,
       var = var_hat, se = se, ci = ci,
       n_target = nT)
}

aipw_var_mu1_mean <- function(data_T1, data_T0, w3, cat0, id = "ID") {
  psid0 <- data_T0[data_T0$G == cat0, ]
  nsw0  <- data_T0[data_T0$G %in% c(1,2), ]
  idx3  <- which(data_T0$G == cat0)
  
  m0_psid <- psid0$y.hat
  m0_nsw  <- nsw0$y.hat
  w3_3    <- w3[idx3]
  r0      <- psid0$re78 - m0_psid
  
  mu0_aipw <- mean(m0_nsw) + sum(w3_3 * r0) / sum(w3_3)
  
  g1 <- data_T1[data_T1$G == 1, ]
  mu1 <- mean(g1$re78)
  ate <- mu1 - mu0_aipw
  
  n1 <- nrow(g1)
  IF_mu1 <- (g1$re78 - mu1) / n1
  
  nT0 <- nrow(nsw0)
  m0bar_T0 <- mean(m0_nsw)
  IF_target0 <- -(m0_nsw - m0bar_T0) / nT0
  
  sw3 <- sum(w3_3)
  r0bar_w <- sum(w3_3 * r0) / sw3
  IF_ratio0 <- -(w3_3 * (r0 - r0bar_w)) / sw3
  
  var_hat <- sum(IF_mu1^2) + sum(IF_target0^2) + sum(IF_ratio0^2)
  se <- sqrt(var_hat)
  ci <- ate + c(-1, 1) * 1.96 * se
  
  list(ATE = ate, mu1 = mu1, mu0 = mu0_aipw,
       var = var_hat, se = se, ci = ci,
       n1 = n1, n_target0 = nT0)
}

###############################################################################
# 8. Entropy Tilting (ET) and Hellinger Distance (HD)
###############################################################################

ET_function <- function(cat0, data_T0, y1, w) {
  nsw  <- data_T0[data_T0$G %in% c(1,2), ]
  m0_nsw  <- nsw$y.hat
  idx_nsw <- which(data_T0$G %in% c(1,2))
  idx0 <- which(data_T0$G == cat0)
  
  target0 <- cbind(mean(m0_nsw), mean((log(1 / w)[idx_nsw])))
  y_T0 <- data_T0$re78[idx0]
  m0 <- data_T0$y.hat[idx0]
  X <- cbind(m0, (log(1 / w))[idx0])
  
  res0 <- cal_tilt_mean_multi(X, target0, y_T0)
  ATE <- mean(y1) - res0$EY
  var_y1 <- var(y1) / length(y1)
  se0 <- cal_tilt_mean_multi_se(X, target0, y_T0, res0)
  SE_ATE <- sqrt(var_y1 + se0$var)
  
  return(list(ATE = ATE, w = res0$w, se = SE_ATE))
}

ET_function_no_fold <- function(cat0, data_T0, y1, w, y_hat_target, y_hat_ctrl,
                                gvar = "G", target_cats = c(1,2)) {
  
  stopifnot(length(w) == nrow(data_T0))
  
  idx0 <- which(data_T0[[gvar]] == cat0)
  idx_nsw <- which(data_T0[[gvar]] %in% target_cats)
  
  w_safe <- pmax(w, 1e-12)
  y_T0 <- data_T0$re78[idx0]
  
  X <- cbind(
    y_hat_ctrl,
    log(1 / w_safe)[idx0]
  )
  
  target0 <- c(
    mean(y_hat_target),
    mean(log(1 / w_safe)[idx_nsw])
  )
  
  res0 <- cal_tilt_mean_multi(X, target0, y_T0)
  se0 <- cal_tilt_mean_multi_se(X, target0, y_T0, res0)
  
  list(
    ATE = res0$EY,
    w = res0$w,
    var.ate = se0$var
  )
}

ET_no_fold_nsw <- function(dat, y = "re78", d = "treat",
                           covars = c("age","education","black","hispanic",
                                      "married","nodegree","re74","re75")) {
  
  X <- intersect(covars, names(dat))
  
  ps <- glm(reformulate(X, d), data = dat, family = binomial())
  ehat <- predict(ps, type = "response")
  
  m1 <- lm(reformulate(X, y), data = dat[dat[[d]] == 1, ])
  m0 <- lm(reformulate(X, y), data = dat[dat[[d]] == 0, ])
  m1x <- predict(m1, newdata = dat)
  m0x <- predict(m0, newdata = dat)
  
  idx1 <- which(exp_dat$treat == 1)
  idx0 <- which(exp_dat$treat == 0)
  
  res1 <- ET_function_no_fold(
    cat0 = 1, data_T0 = dat, y1 = exp_dat$re78[idx1],
    w = ehat, y_hat_target = m1x, y_hat_ctrl = m1x[idx1]
  )
  
  res0 <- ET_function_no_fold(
    cat0 = 0, data_T0 = dat, y1 = exp_dat$re78[idx1],
    w = 1 - ehat, y_hat_target = m0x, y_hat_ctrl = m0x[idx0]
  )
  
  ATE <- res1$EY1 - res0$EY1
  se <- sqrt(res1$var.ate + res0$var.ate)
  return(round(c(ATE, se)))
}

HD_function <- function(cat0, data_T0, y1, w) {
  nsw  <- data_T0[data_T0$G %in% c(1,2), ]
  m0_nsw  <- nsw$y.hat
  idx_nsw <- which(data_T0$G %in% c(1,2))
  idx0 <- which(data_T0$G == cat0)
  
  target0 <- cbind(mean(m0_nsw), mean((log(1 / w)[idx_nsw])))
  y_T0 <- data_T0$re78[idx0]
  m0 <- data_T0$y.hat[idx0]
  X <- cbind(m0, (log(1 / w))[idx0])
  
  res0 <- cal_tilt_mean_multi_cvxr_hellinger(X, target0, y_T0, solver = "ECOS")
  ATE <- mean(y1) - res0$EY
  var_y1 <- var(y1) / length(y1)
  se0 <- cal_tilt_mean_multi_se(X, target0, y_T0, res0)
  SE_ATE <- sqrt(var_y1 + se0$var)
  
  return(round(c(ATE, ATE - ate, SE_ATE)))
}

###############################################################################
# 9. Unweighted and HT Estimators
###############################################################################

unweighted_function <- function(d4, cat1, cat0) {
  y1 <- d4$re78[d4$G == cat1]
  y0 <- d4$re78[d4$G == cat0]
  ate_unweighted <- mean(y1) - mean(y0)
  se_unweighted  <- sqrt(var(y1) / length(y1) + var(y0) / length(y0))
  return(round(c(ate_unweighted, se_unweighted)))
}

HT_general_fixedvar <- function(d, yvar, gvar = "G",
                                g1, g0,
                                w1, w0,
                                N_target = NULL) {
  
  y <- d[[yvar]]
  G <- d[[gvar]]
  
  idx1 <- which(G %in% g1)
  idx0 <- which(G %in% g0)
  
  stopifnot(length(idx1) > 1, length(idx0) > 1)
  
  if (is.null(N_target)) N_target <- length(which(G %in% c(g1, g0)))
  
  y1 <- y[idx1]
  y0 <- y[idx0]
  ww1 <- w1[idx1]
  ww0 <- w0[idx0]
  
  mu1_HT <- sum(ww1 * y1) / N_target
  mu0_HT <- sum(ww0 * y0) / N_target
  HT <- mu1_HT - mu0_HT
  
  a1 <- ww1 * y1
  a0 <- ww0 * y0
  
  var_mu1 <- sum((a1 - mean(a1))^2) / ((length(a1) - 1) * N_target^2)
  var_mu0 <- sum((a0 - mean(a0))^2) / ((length(a0) - 1) * N_target^2)
  
  var_HT <- var_mu1 + var_mu0
  se_HT  <- sqrt(var_HT)
  
  list(mu1_HT = mu1_HT, mu0_HT = mu0_HT, HT = HT,
       var_mu1 = var_mu1, var_mu0 = var_mu0,
       var_HT = var_HT, se_HT = se_HT,
       N_target = N_target, n1 = length(idx1), n0 = length(idx0))
}

###############################################################################
# 10. Bootstrap for ET with Multinomial PS and GAM Refit
###############################################################################

bootstrap_ET_refit_multinom_and_gam <- function(d4,
                                                xvars = c("age","education","black","hispanic",
                                                          "married","nodegree","re74","re75"),
                                                donor_cat = 3,
                                                cat1 = 1, target_cats = c(1,2),
                                                yvar = "re78", gvar = "G", tvar = "treat",
                                                cont_vars = c("age","education","re74","re75"),
                                                cat_vars  = c("black","hispanic","married","nodegree"),
                                                bs = "cs",
                                                donor_controls_only = FALSE,
                                                eps = 1e-8,
                                                B = 300, seed = 2026) {
  
  stopifnot(requireNamespace("nnet", quietly = TRUE))
  stopifnot(requireNamespace("mgcv", quietly = TRUE))
  
  set.seed(seed)
  G <- d4[[gvar]]
  
  idx_g1 <- which(G == 1)
  idx_g2 <- which(G == 2)
  idx_g3 <- which(G == 3)
  idx_g4 <- which(G == 4)
  
  stopifnot(length(idx_g1) > 1, length(idx_g2) > 1, length(idx_g3) > 1, length(idx_g4) > 1)
  stopifnot(donor_cat %in% c(3,4))
  
  n1 <- length(idx_g1)
  n2 <- length(idx_g2)
  n3 <- length(idx_g3)
  n4 <- length(idx_g4)
  
  ate_b <- numeric(B)
  
  for (b in seq_len(B)) {
    s1 <- sample(idx_g1, n1, replace = TRUE)
    s2 <- sample(idx_g2, n2, replace = TRUE)
    s3 <- sample(idx_g3, n3, replace = TRUE)
    s4 <- sample(idx_g4, n4, replace = TRUE)
    
    d4_b <- d4[c(s1, s2, s3, s4), , drop = FALSE]
    d4_b[[gvar]] <- factor(d4_b[[gvar]], levels = c(1,2,3,4))
    
    fml_multi_b <- as.formula(paste(gvar, "~", paste(xvars, collapse = "+")))
    fit_multi_b <- nnet::multinom(fml_multi_b, data = d4_b, trace = FALSE)
    
    P <- predict(fit_multi_b, type = "probs")
    if (is.null(colnames(P))) stop("predict(multinom, type='probs') returned no column names.")
    
    p1 <- P[, "1"]
    p2 <- P[, "2"]
    p3 <- P[, "3"]
    p4 <- P[, "4"]
    pNSW <- p1 + p2
    
    if (donor_cat == 3) {
      w_donor <- pNSW / pmax(p3, eps)
    } else {
      w_donor <- pNSW / pmax(p4, eps)
    }
    
    res_gam_b <- aipw_transport_gam(
      d4 = d4_b, xvars = xvars,
      cat1 = cat1, target_cats = target_cats,
      cat0 = donor_cat,
      w0 = w_donor,
      yvar = yvar, gvar = gvar, tvar = tvar,
      cont_vars = cont_vars,
      cat_vars  = cat_vars,
      bs = bs,
      donor_controls_only = donor_controls_only
    )
    
    y_hat_target_b <- res_gam_b$y_hat_target
    y_hat_ctrl_b   <- res_gam_b$y_hat_ctrl
    y1_b <- d4_b[[yvar]][d4_b[[gvar]] == as.character(cat1)]
    
    res_b <- ET_function_no_fold(
      cat0 = donor_cat,
      data_T0 = d4_b,
      y1 = y1_b,
      w  = w_donor,
      y_hat_target = y_hat_target_b,
      y_hat_ctrl   = y_hat_ctrl_b
    )
    
    ate_b[b] <- res_b$ATE
  }
  
  se_boot <- sd(ate_b, na.rm = TRUE)
  
  list(
    ate_boot = ate_b,
    se_boot  = se_boot,
    ci_pct   = quantile(ate_b, probs = c(0.025, 0.975), na.rm = TRUE),
    ci_norm  = mean(ate_b, na.rm = TRUE) + c(-1,1) * 1.96 * se_boot
  )
}




run_et <- function(cat0, idx_S0, idx_U0_0, w, d4, k = 4, seed = 2024202) {
  
  # T0 prediction
  fold_T0 <- k_fold_predict_yhat_gam(
    df = d4,
    K = k,
    idx_S = idx_S0,
    idx_U0 = idx_U0_0,
    seed = seed
  )
  data_T0 <- dplyr::bind_rows(fold_T0)
  
  # T1 prediction
  fold_T1 <- k_fold_predict_yhat_gam(
    df = d4,
    K = k,
    idx_S = which(d4$G == 1),
    idx_U0 = which(d4$G %in% c(2, 3, 4)),
    seed = seed
  )
  data_T1 <- dplyr::bind_rows(fold_T1)
  
  # reorder
  d4$ID <- seq_len(nrow(d4))
  data_T0 <- data_T0[order(data_T0$ID), ]
  data_T1 <- data_T1[order(data_T1$ID), ]
  
  # ET result
  ET_function(cat0 = cat0, data_T0 = data_T0, y1 = y1, w = w)
}


##############################################################################
#10 AIPW functions
##############################################################################
run_aipw_cv <- function(cat0, w_use, t0_model = c("lm", "gam")) {
  
  t0_model <- match.arg(t0_model)
  
  # ----------------------------
  # T0 prediction
  # ----------------------------
  idx_S  <- which(d4$G == cat0)
  idx_U0 <- which(d4$G %in% setdiff(c(1, 2, 3, 4), cat0))
  
  if (t0_model == "lm") {
    fold_T0 <- k_fold_predict_yhat(
      df = d4,
      K = k,
      idx_S = idx_S,
      idx_U0 = idx_U0,
      seed = seed
    )
  } else {
    fold_T0 <- k_fold_predict_yhat_gam(
      df = d4,
      K = k,
      idx_S = idx_S,
      idx_U0 = idx_U0,
      seed = seed
    )
  }
  
  data_T0 <- dplyr::bind_rows(fold_T0)
  
  # ----------------------------
  # T1 prediction: GAM
  # ----------------------------
  idx_S  <- which(d4$G == 1)
  idx_U0 <- which(d4$G %in% c(2, 3, 4))
  
  fold_T1 <- k_fold_predict_yhat_gam(
    df = d4,
    K = k,
    idx_S = idx_S,
    idx_U0 = idx_U0,
    seed = seed
  )
  
  data_T1 <- dplyr::bind_rows(fold_T1)
  
  # ----------------------------
  # Reorder and estimate
  # ----------------------------
  data_T0 <- data_T0[order(data_T0$ID), ]
  data_T1 <- data_T1[order(data_T1$ID), ]
  
  res <- aipw_var_mu1_mean(
    data_T1 = data_T1,
    data_T0 = data_T0,
    w3 = w_use,
    cat0 = cat0,
    id = "ID"
  )
  
  return(res)
}











###############################################################################
#11. Plot function
###############################################################################

w_ecdf_df <- function(x, w = NULL) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  x <- x[ok]
  
  if (is.null(w)) {
    w <- rep(1, length(x))
  } else {
    w <- as.numeric(w)[ok]
  }
  
  # keep nonnegative weights
  #w[w < 0] <- 0
  #if (sum(w) == 0) w <- rep(1, length(w))
  
  o <- order(x)
  x <- x[o]; w <- w[o]
  
  # normalize to sum to 1 for a CDF
  w <- w / sum(w)
  data.frame(x = x, F = cumsum(w))
}

get_cal_weight <- function(cat0, idx_U0_train, log_term) {
  
  # Step 1: fit T0 prediction
  idx_S  <- which(d4$G == cat0)
  idx_U0 <- which(d4$G %in% idx_U0_train)
  
  fold_T0 <- k_fold_predict_yhat_gam(
    df = d4,
    K = k,
    idx_S = idx_S,
    idx_U0 = idx_U0,
    seed = seed
  )
  data_T0 <- dplyr::bind_rows(fold_T0)
  
  # Step 2: fit T1 prediction
  idx_S  <- which(d4$G == 1)
  idx_U0 <- which(d4$G %in% c(2, 3, 4))
  
  fold_T1 <- k_fold_predict_yhat_gam(
    df = d4,
    K = k,
    idx_S = idx_S,
    idx_U0 = idx_U0,
    seed = seed
  )
  data_T1 <- dplyr::bind_rows(fold_T1)
  
  # Step 3: reorder
  d4$ID <- seq_len(nrow(d4))
  data_T0 <- data_T0[order(data_T0$ID), ]
  data_T1 <- data_T1[order(data_T1$ID), ]
  
  # Step 4: build calibration target
  nsw <- data_T0[data_T0$G %in% c(1, 2), ]
  m0_nsw <- nsw$y.hat
  idx_nsw <- which(data_T0$G %in% c(1, 2))
  idx0 <- which(data_T0$G == cat0)
  
  target0 <- cbind(
    mean(m0_nsw),
    mean(log_term[idx_nsw])
  )
  
  y_T0 <- data_T0$re78[idx0]
  m0   <- data_T0$y.hat[idx0]
  X    <- cbind(m0, log_term[idx0])
  
  res0 <- cal_tilt_mean_multi(X, target0, y_T0)
  
  return(res0$w)
}


plot_cdf_panel <- function(x_target, x_donor,
                           w_cal, w_ipw, w_ebcw, w_CBPS,
                           xlab = "", ylab = "CDF",
                           add_legend = TRUE,
                           xlim = NULL) {
  
  # target empirical (NSW control)  -- keep RED
  df_t   <- w_ecdf_df(x_target, NULL)
  
  # donor weighted curves
  df_cal <- w_ecdf_df(x_donor, w_cal)
  df_ipw <- w_ecdf_df(x_donor, w_ipw)
  df_eb  <- w_ecdf_df(x_donor, w_ebcw)
  df_cb  <- w_ecdf_df(x_donor, w_CBPS)
  
  
  xr <- if (is.null(xlim)) range(c(df_t$x, df_cal$x, df_ipw$x, df_eb$x, df_cb$x), finite = TRUE) else xlim
  
  plot(df_t$x, df_t$F, type = "l",
       col = "red", lty = 1, lwd = 2,
       xlab = xlab, ylab = ylab, xlim = xr)
  
  
  # donor lines (different colors)
  lines(df_cal$x, df_cal$F, col = "blue",      lty = 2, lwd = 2)  # calibration
  lines(df_ipw$x, df_ipw$F, col = "black",     lty = 3, lwd = 2)  # IPW
  lines(df_eb$x,  df_eb$F,  col = "darkgreen", lty = 4, lwd = 2)  # EBCW
  lines(df_cb$x,  df_cb$F,  col = "purple",    lty = 5, lwd = 2)  # CBPS
  
  if (add_legend) {
    legend("bottomright",
           legend = c("NSW (empirical)", "ET", "IPW", "EBCW", "CBPS"),
           col    = c("red", "blue", "black", "darkgreen", "purple"),
           lty    = c(1, 2, 3, 4, 5),
           lwd    = c(2, 2, 2, 2, 2),
           bty = "n", cex = 0.6)
  }
}



###############################################################################
# End of Script
###############################################################################


###############################################################################
# ============================================================================
# FINAL ADDITIONS FOR THE REAL-DATA TRANSPORT ANALYSIS
# ============================================================================
#
# These functions DO NOT replace the original cross-fitting code above.
# In particular, k_fold_predict_yhat() and k_fold_predict_yhat_gam() above
# are used exactly as they were in the previous real-data analysis.
#
# Changes made here:
#   1. four-category multinomial source model is retained;
#   2. IPW/AIPW variance includes multinomial-phi uncertainty;
#   3. ET/HD/EL use the same dual/Newton entropy machinery as the simulation;
#   4. ET/HD/EL variance uses a joint sandwich with multinomial phi;
#   5. theta is NOT used in calibration; it is estimated after lambda.
###############################################################################


safe_transport_ratio <- function(num, den, eps = 1e-12) {
  pmax(as.numeric(num), eps) / pmax(as.numeric(den), eps)
}


transport_ess <- function(w) {
  w <- as.numeric(w)
  w <- w[is.finite(w)]
  if (length(w) == 0L || sum(w^2) <= 0) return(NA_real_)
  sum(w)^2 / sum(w^2)
}


transport_ci <- function(est, se) {
  if (!is.finite(est) || !is.finite(se)) return(c(NA_real_, NA_real_))
  est + c(-1, 1) * qnorm(.975) * se
}


###############################################################################
# A. Four-category multinomial model and score
###############################################################################

fit_multinom_transport_final <- function(d4, xvars) {

  dat <- d4
  dat$G <- factor(dat$G, levels = c(1, 2, 3, 4))

  fml <- as.formula(
    paste("G ~", paste(xvars, collapse = " + "))
  )

  fit <- nnet::multinom(
    fml,
    data = dat,
    trace = FALSE
  )

  X <- model.matrix(fml, data = dat)

  B <- as.matrix(coef(fit))
  B <- B[c("2", "3", "4"), , drop = FALSE]

  # phi ordering is blockwise:
  # beta_2, beta_3, beta_4
  phi_hat <- as.vector(t(B))

  P <- predict(fit, type = "probs")
  P <- as.matrix(P)
  P <- P[, c("1", "2", "3", "4"), drop = FALSE]

  list(
    fit = fit,
    formula = fml,
    X = X,
    B = B,
    phi_hat = phi_hat,
    P = P,
    xvars = xvars
  )
}


multinom_prob_final <- function(phi, X) {

  X <- as.matrix(X)
  r <- ncol(X)

  stopifnot(length(phi) == 3L * r)

  B <- matrix(
    phi,
    nrow = 3L,
    ncol = r,
    byrow = TRUE
  )

  eta <- X %*% t(B)

  # Stable baseline-category softmax; G=1 is baseline.
  mx <- pmax(0, apply(eta, 1, max))
  e0 <- exp(-mx)
  ek <- exp(eta - mx)
  den <- e0 + rowSums(ek)

  cbind(
    "1" = e0 / den,
    "2" = ek[, 1] / den,
    "3" = ek[, 2] / den,
    "4" = ek[, 3] / den
  )
}


multinom_score_final <- function(phi, X, G) {

  G <- as.character(G)
  P <- multinom_prob_final(phi, X)

  do.call(
    cbind,
    lapply(
      c("2", "3", "4"),
      function(g) {
        X * as.numeric((G == g) - P[, g])
      }
    )
  )
}


multinom_information_final <- function(phi, X) {

  P <- multinom_prob_final(phi, X)
  Prest <- P[, c("2", "3", "4"), drop = FALSE]

  N <- nrow(X)
  r <- ncol(X)

  A <- matrix(0, 3L * r, 3L * r)

  for (j in 1:3) {
    rr <- ((j - 1L) * r + 1L):(j * r)

    for (k in 1:3) {
      cc <- ((k - 1L) * r + 1L):(k * r)

      s <- Prest[, j] * (
        as.numeric(j == k) - Prest[, k]
      )

      A[rr, cc] <- crossprod(
        X,
        X * as.numeric(s)
      ) / N
    }
  }

  A
}


transport_ratio_final <- function(P, donor_g) {

  donor_g <- as.character(donor_g)

  pNSW <- P[, "1"] + P[, "2"]

  safe_transport_ratio(
    pNSW,
    P[, donor_g]
  )
}


###############################################################################
# B. Exact previous cross-fitted predictions
#
# IMPORTANT:
# These call the ORIGINAL functions above without modifying the GAM formula,
# factor handling, basis dimension, folds, or seed.
###############################################################################

get_original_crossfit_T0 <- function(
    d4,
    donor_g,
    K = 4,
    seed = 2024202,
    model = c("lm", "gam")
) {

  model <- match.arg(model)

  idx_S <- which(
    as.character(d4$G) == as.character(donor_g)
  )

  idx_U0 <- which(
    as.character(d4$G) != as.character(donor_g)
  )

  folds <- if (model == "lm") {
    k_fold_predict_yhat(
      df = d4,
      K = K,
      idx_S = idx_S,
      idx_U0 = idx_U0,
      seed = seed
    )
  } else {
    # EXACT previous GAM routine.
    k_fold_predict_yhat_gam(
      df = d4,
      K = K,
      idx_S = idx_S,
      idx_U0 = idx_U0,
      seed = seed
    )
  }

  dat <- dplyr::bind_rows(folds)
  dat <- dat[order(dat$ID), , drop = FALSE]

  dat
}


###############################################################################
# C. Generic numerical sandwich for unconstrained point estimators
###############################################################################

sandwich_numeric_final <- function(
    beta_hat,
    psi_fun,
    contrast,
    N
) {

  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("Install package 'numDeriv'.")
  }

  Psi <- psi_fun(beta_hat)

  if (any(!is.finite(Psi))) {
    return(list(
      success = FALSE,
      SE = NA_real_,
      variance = NA_real_,
      A = NULL,
      B = NULL,
      V = NULL,
      Psi = Psi
    ))
  }

  B <- crossprod(Psi) / N

  J <- tryCatch(
    numDeriv::jacobian(
      function(beta) colMeans(psi_fun(beta)),
      beta_hat
    ),
    error = function(e) NULL
  )

  if (is.null(J) || any(!is.finite(J))) {
    return(list(
      success = FALSE,
      SE = NA_real_,
      variance = NA_real_,
      A = NULL,
      B = B,
      V = NULL,
      Psi = Psi
    ))
  }

  A <- -J

  Ainv <- tryCatch(
    solve(A),
    error = function(e) NULL
  )

  if (is.null(Ainv) || any(!is.finite(Ainv))) {
    return(list(
      success = FALSE,
      SE = NA_real_,
      variance = NA_real_,
      A = A,
      B = B,
      V = NULL,
      Psi = Psi
    ))
  }

  V <- (Ainv %*% B %*% t(Ainv)) / N

  vv <- as.numeric(
    t(contrast) %*% V %*% contrast
  )

  if (!is.finite(vv) || vv < -1e-8) {
    se <- NA_real_
  } else {
    vv <- max(vv, 0)
    se <- sqrt(vv)
  }

  list(
    success = is.finite(se),
    SE = se,
    variance = vv,
    A = A,
    B = B,
    V = V,
    Psi = Psi
  )
}


###############################################################################
# D. IPW: SAME point estimator as previous code;
#         new variance includes multinomial phi uncertainty.
###############################################################################

IPW_multinom_joint_final <- function(
    d4,
    donor_g,
    ps_obj,
    benchmark = NA_real_
) {

  Y <- as.numeric(d4$re78)
  G <- as.character(d4$G)
  X <- as.matrix(ps_obj$X)

  N <- nrow(d4)

  phi_hat <- as.numeric(ps_obj$phi_hat)
  pphi <- length(phi_hat)

  P_hat <- multinom_prob_final(phi_hat, X)
  rr_hat <- transport_ratio_final(P_hat, donor_g)

  I1 <- as.numeric(G == "1")
  Ig <- as.numeric(G == as.character(donor_g))

  # EXACT old point estimator.
  mu1_hat <- mean(Y[I1 == 1])

  mu0_hat <- sum(
    rr_hat[Ig == 1] * Y[Ig == 1]
  ) / sum(
    rr_hat[Ig == 1]
  )

  ate_hat <- mu1_hat - mu0_hat

  beta_hat <- c(
    phi_hat,
    mu1_hat,
    mu0_hat
  )

  idx_phi <- seq_len(pphi)
  idx_mu1 <- pphi + 1L
  idx_mu0 <- pphi + 2L

  psi_fun <- function(beta) {

    phi <- beta[idx_phi]
    mu1 <- beta[idx_mu1]
    mu0 <- beta[idx_mu0]

    P <- multinom_prob_final(phi, X)
    rr <- transport_ratio_final(P, donor_g)

    cbind(
      multinom_score_final(phi, X, G),
      mu1 = I1 * (Y - mu1),
      mu0 = Ig * rr * (Y - mu0)
    )
  }

  contrast <- rep(0, length(beta_hat))
  contrast[idx_mu1] <- 1
  contrast[idx_mu0] <- -1

  sw <- sandwich_numeric_final(
    beta_hat,
    psi_fun,
    contrast,
    N
  )

  wd <- rr_hat[Ig == 1]

  list(
    ATE = ate_hat,
    SE = sw$SE,
    variance = sw$variance,
    CI = transport_ci(ate_hat, sw$SE),
    ESS = transport_ess(wd),
    MaxW = max(wd),
    success = sw$success,
    A = sw$A,
    B = sw$B,
    V = sw$V,
    variance_type =
      "joint sandwich: multinomial phi + mu1 + donor mean"
  )
}


###############################################################################
# E. AIPW: EXACT previous point-estimation formula and EXACT old cross-fitting.
#
# Cross-fitted m0 is treated as fixed in the sandwich, matching simulation.
#
# beta = (phi_multinom, mu1, mbar_target, residual_bar)
###############################################################################

AIPW_multinom_joint_final <- function(
    d4,
    donor_g,
    ps_obj,
    model = c("lm", "gam"),
    K = 4,
    seed = 2024202,
    benchmark = NA_real_
) {

  model <- match.arg(model)

  data_T0 <- get_original_crossfit_T0(
    d4 = d4,
    donor_g = donor_g,
    K = K,
    seed = seed,
    model = model
  )

  # The exact old run_aipw_cv() uses the same object and the same formula
  # implemented by aipw_var_mu1_mean().
  G0 <- as.character(data_T0$G)

  Y0 <- as.numeric(data_T0$re78)
  m0 <- as.numeric(data_T0$y.hat)

  # Align multinomial probabilities by ID because data_T0 is reordered by ID.
  ord <- match(data_T0$ID, seq_len(nrow(d4)))

  X <- as.matrix(ps_obj$X[ord, , drop = FALSE])
  G <- as.character(d4$G[ord])

  phi_hat <- as.numeric(ps_obj$phi_hat)
  pphi <- length(phi_hat)

  P_hat <- multinom_prob_final(phi_hat, X)
  rr_hat <- transport_ratio_final(P_hat, donor_g)

  I1 <- as.numeric(G0 == "1")
  IT <- as.numeric(G0 %in% c("1", "2"))
  Ig <- as.numeric(G0 == as.character(donor_g))

  mu1_hat <- mean(Y0[I1 == 1])

  mbar_hat <- mean(m0[IT == 1])

  resid <- Y0 - m0

  rbar_hat <- sum(
    rr_hat[Ig == 1] * resid[Ig == 1]
  ) / sum(
    rr_hat[Ig == 1]
  )

  mu0_hat <- mbar_hat + rbar_hat

  ate_hat <- mu1_hat - mu0_hat

  # Numerical identity check against the exact old point estimator function.
  old_point <- aipw_var_mu1_mean(
    data_T1 = {
      fold_T1 <- k_fold_predict_yhat_gam(
        df = d4,
        K = K,
        idx_S = which(as.character(d4$G) == "1"),
        idx_U0 = which(as.character(d4$G) %in% c("2", "3", "4")),
        seed = seed
      )
      dt1 <- dplyr::bind_rows(fold_T1)
      dt1[order(dt1$ID), , drop = FALSE]
    },
    data_T0 = data_T0,
    w3 = {
      # Old function indexes the full ratio vector by positions in data_T0.
      # data_T0 contains all rows exactly once, ordered by ID here.
      rr_full <- transport_ratio_final(
        multinom_prob_final(ps_obj$phi_hat, ps_obj$X),
        donor_g
      )
      rr_full[ord]
    },
    cat0 = as.numeric(donor_g),
    id = "ID"
  )

  # If this check fails by more than machine tolerance, stop rather than
  # silently changing the old point estimator.
  if (is.finite(old_point$ATE) &&
      abs(old_point$ATE - ate_hat) > 1e-7) {
    stop(
      "AIPW point-estimate identity check failed. ",
      "Old ATE=", old_point$ATE,
      ", reconstructed ATE=", ate_hat
    )
  }

  N <- nrow(data_T0)

  beta_hat <- c(
    phi_hat,
    mu1_hat,
    mbar_hat,
    rbar_hat
  )

  idx_phi <- seq_len(pphi)
  idx_mu1 <- pphi + 1L
  idx_mbar <- pphi + 2L
  idx_rbar <- pphi + 3L

  psi_fun <- function(beta) {

    phi <- beta[idx_phi]
    mu1 <- beta[idx_mu1]
    mbar <- beta[idx_mbar]
    rbar <- beta[idx_rbar]

    P <- multinom_prob_final(phi, X)
    rr <- transport_ratio_final(P, donor_g)

    cbind(
      multinom_score_final(phi, X, G),
      mu1 = I1 * (Y0 - mu1),
      mbar = IT * (m0 - mbar),
      rbar = Ig * rr * (resid - rbar)
    )
  }

  contrast <- rep(0, length(beta_hat))
  contrast[idx_mu1] <- 1
  contrast[idx_mbar] <- -1
  contrast[idx_rbar] <- -1

  sw <- sandwich_numeric_final(
    beta_hat,
    psi_fun,
    contrast,
    N
  )

  wd <- rr_hat[Ig == 1]

  list(
    ATE = ate_hat,
    old_ATE = old_point$ATE,
    SE = sw$SE,
    old_SE = old_point$se,
    variance = sw$variance,
    CI = transport_ci(ate_hat, sw$SE),
    ESS = transport_ess(wd),
    MaxW = max(wd),
    success = sw$success,
    predictions = data_T0,
    A = sw$A,
    B = sw$B,
    V = sw$V,
    variance_type =
      "joint sandwich: multinomial phi; exact old cross-fitted outcome predictions fixed"
  )
}


###############################################################################
# F. Simulation-style entropy definitions for ET, HD, CE
#
# ET and HD use the simulation inverse-dual maps directly for transport weights.
#
# CE requires a careful transport translation:
#
# simulation raw CE weight:
#   omega = -1/expm1(eta) = 1/pi_D > 1
#
# For target-vs-donor transport:
#   pi_D = p_g / (p_NSW + p_g)
#   1/pi_D = 1 + r_g
#
# The TARGET transport weight is the excess inverse probability:
#   r_g = omega - 1.
#
# Therefore CE calibration is written exactly as the simulation-style
# inverse-probability calibration over the target+donor population:
#
#   sum donor * omega * S = sum (target + donor) * S
#
# which is algebraically equivalent to
#
#   sum donor * (omega - 1) * S = sum target * S.
#
# The transported outcome mean uses (omega - 1), not omega.
###############################################################################

simulation_entropy_transport_maps <- function(
    entropy = c("ET", "HD", "EL")
) {

  entropy <- match.arg(entropy)

  weight_raw <- switch(
    entropy,

    ET = function(eta) {
      exp(eta)
    },

    HD = function(eta) {
      1 / (4 * eta^2)
    },

    EL = function(eta) {
      -1 / eta
    }
  )

  weight_deriv <- switch(
    entropy,

    ET = function(eta) {
      exp(eta)
    },

    HD = function(eta) {
      -1 / (2 * eta^3)
    },

    EL = function(eta) {
      1 / eta^2
    }
  )

  primitive <- switch(
    entropy,

    ET = function(eta) {
      exp(eta)
    },

    HD = function(eta) {
      -1 / (4 * eta)
    },

    EL = function(eta) {
      -log(-eta)
    }
  )

  valid_eta <- switch(
    entropy,

    ET = function(eta) {
      all(is.finite(eta))
    },

    HD = function(eta) {
      all(is.finite(eta) & eta < -1e-8)
    },

    EL = function(eta) {
      all(is.finite(eta) & eta < -1e-8)
    }
  )

  list(
    raw_weight = weight_raw,
    weight_deriv = weight_deriv,
    primitive = primitive,
    valid_eta = valid_eta
  )
}

transport_g_simulation_style <- function(
    r,
    entropy = c("ET", "HD", "EL")
) {
  
  entropy <- match.arg(entropy)
  
  r <- pmax(as.numeric(r), 1e-12)
  
  switch(
    entropy,
    
    # ET: g(r) = log(r)
    ET = log(1+r),
    
    # HD: g(r) = -1 / (2 sqrt(r))
    HD = -1 / (2 * sqrt(1+r)),

    EL = -1 /(1+ r)
  )
}



###############################################################################
# G. Simulation-style G entropy direction adapted to transport
###############################################################################


###############################################################################
# H. Calibration matrix for proposed GEC
#
# Uses EXACT previous GAM predictions.
# theta is NOT included.
###############################################################################

build_GEC_transport_S <- function(
    r,
    yhat0,
    entropy
) {

  cbind(
    intercept = 1,
    yhat = as.numeric(yhat0),
    g_pi =
      transport_g_simulation_style(
        r,
        entropy
      )
  )
}


###############################################################################
# I. Newton + backtracking, same numerical style as simulation
###############################################################################

solve_lambda_transport_simulation_style <- function(
    S,
    donor_ind,
    target_ind,
    entropy = c("ET", "HD", "EL"),
    lambda_start = NULL,
    tol = 1e-8,
    maxit = 2000,
    ridge = 1e-10,
    armijo = 1e-4,
    min_alpha = 1e-12
) {

  entropy <- match.arg(entropy)

  S <- as.matrix(S)
  donor_ind <- as.numeric(donor_ind)
  target_ind <- as.numeric(target_ind)

  N <- nrow(S)
  q <- ncol(S)

  ef <- simulation_entropy_transport_maps(entropy)

  
  #rhs_ind <- target_ind
  rhs_ind <- pmin(target_ind + donor_ind, 1)
  nD <- sum(donor_ind)
  nT <- sum(target_ind)

  if (is.null(lambda_start)) {

    lambda_start <- rep(0, q)

    omega_bar <- 1 + nT / nD
    
    if (entropy == "ET") {
      lambda_start[1] <- log(omega_bar)
    }
    
    if (entropy == "HD") {
      lambda_start[1] <- -1 / (
        2 * sqrt(omega_bar)
      )
    }
    
    if (entropy == "EL") {
      lambda_start[1] <- -1 / omega_bar
    }
  }

  lambda <- as.numeric(lambda_start)

  objective <- function(lam) {

    eta <- as.numeric(S %*% lam)

    if (!ef$valid_eta(eta[donor_ind == 1])) {
      return(Inf)
    }

    (
      sum(
        ef$primitive(
          eta[donor_ind == 1]
        )
      ) -
        sum(
          rhs_ind * eta
        )
    ) / N
  }

  for (iter in seq_len(maxit)) {

    eta <- as.numeric(S %*% lambda)

    if (!ef$valid_eta(eta[donor_ind == 1])) {
      return(list(
        converged = FALSE,
        reason = "dual domain violation",
        lambda = lambda
      ))
    }

    wraw <- numeric(N)
    wraw[donor_ind == 1] <-
      ef$raw_weight(
        eta[donor_ind == 1]
      )

    raw_score <- colSums(
      S *
        as.numeric(
          donor_ind * wraw -
            rhs_ind
        )
    )

    max_raw <- max(abs(raw_score))

    if (is.finite(max_raw) &&
        max_raw < tol) {

      #wtransport <-donor_ind * wraw
      wtransport <- donor_ind * (wraw - 1)
      return(list(
        converged = TRUE,
        reason = "converged",
        lambda = lambda,
        eta = eta,
        raw_weights_all = wraw,
        transport_weights_all = wtransport,
        raw_score = raw_score,
        max_raw_score = max_raw,
        iterations = iter
      ))
    }

    wp <- numeric(N)
    wp[donor_ind == 1] <-
      ef$weight_deriv(
        eta[donor_ind == 1]
      )

    H <- crossprod(
      S,
      S *
        as.numeric(
          donor_ind * wp
        )
    ) / N

    H <- H + diag(ridge, q)

    score <- raw_score / N

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
        raw_score = raw_score,
        max_raw_score = max_raw,
        iterations = iter
      ))
    }

    Q0 <- objective(lambda)

    alpha <- 1
    accepted <- FALSE

    while (alpha >= min_alpha) {

      candidate <- lambda -
        alpha * as.numeric(step)

      eta_new <- as.numeric(
        S %*% candidate
      )

      if (ef$valid_eta(
        eta_new[donor_ind == 1]
      )) {

        Q1 <- objective(candidate)

        if (
          is.finite(Q1) &&
            Q1 <=
              Q0 -
              armijo *
                alpha *
                sum(score * step)
        ) {
          lambda <- candidate
          accepted <- TRUE
          break
        }
      }

      alpha <- alpha / 2
    }

    if (!accepted) {
      return(list(
        converged = FALSE,
        reason = "backtracking failed",
        lambda = lambda,
        raw_score = raw_score,
        max_raw_score = max_raw,
        iterations = iter
      ))
    }
  }

  eta <- as.numeric(S %*% lambda)

  wraw <- numeric(N)

  if (ef$valid_eta(
    eta[donor_ind == 1]
  )) {
    wraw[donor_ind == 1] <-
      ef$raw_weight(
        eta[donor_ind == 1]
      )
  }

  #wtransport <- donor_ind * wraw
  wtransport <- donor_ind * (wraw - 1)
  raw_score <- colSums(
    S *
      as.numeric(
        donor_ind * wraw -
          rhs_ind
      )
  )

  list(
    converged = FALSE,
    reason = "maximum iterations reached",
    lambda = lambda,
    eta = eta,
    raw_weights_all = wraw,
    transport_weights_all = wtransport,
    raw_score = raw_score,
    max_raw_score = max(abs(raw_score)),
    iterations = maxit
  )
}


###############################################################################
# J. Point estimator for ET/HD/EL
#
# EXACT old GAM cross-fitting.
# Same simulation entropy solver.
# No theta inside calibration.
###############################################################################

GEC_transport_point_simulation_style <- function(
    d4,
    donor_g,
    ps_obj,
    entropy = c("ET", "HD", "EL"),
    K = 4,
    seed = 2024202
) {

  entropy <- match.arg(entropy)

  data_T0 <- get_original_crossfit_T0(
    d4 = d4,
    donor_g = donor_g,
    K = K,
    seed = seed,
    model = "gam"
  )

  ord <- match(
    data_T0$ID,
    seq_len(nrow(d4))
  )

  P <- ps_obj$P[
    ord,
    ,
    drop = FALSE
  ]

  rr <- transport_ratio_final(
    P,
    donor_g
  )

  S <- build_GEC_transport_S(
    r = rr,
    yhat0 = data_T0$y.hat,
    entropy = entropy
  )

  G <- as.character(data_T0$G)
  Y <- as.numeric(data_T0$re78)

  ID <- as.numeric(
    G == as.character(donor_g)
  )

  IT <- as.numeric(
    G %in% c("1", "2")
  )

  I1 <- as.numeric(
    G == "1"
  )

  fit <- solve_lambda_transport_simulation_style(
    S = S,
    donor_ind = ID,
    target_ind = IT,
    entropy = entropy
  )

  if (is.null(fit$transport_weights_all) ||
      any(!is.finite(
        fit$transport_weights_all
      ))) {
    stop(
      entropy,
      " calibration failed: ",
      fit$reason
    )
  }

  if (!isTRUE(fit$converged)) {
    warning(
      entropy,
      " calibration warning: ",
      fit$reason,
      "; max raw score=",
      signif(
        fit$max_raw_score,
        5
      )
    )
  }

  wt <- fit$transport_weights_all

  mu1_hat <- mean(
    Y[I1 == 1]
  )

  theta0_hat <- sum(
    ID * wt * Y
  ) / sum(
    ID * wt
  )

  ate_hat <- mu1_hat -
    theta0_hat

  list(
    entropy = entropy,
    ATE = ate_hat,
    mu1 = mu1_hat,
    theta0 = theta0_hat,
    lambda = as.numeric(fit$lambda),
    data_T0 = data_T0,
    ratio = rr,
    S = S,
    raw_weights_all = fit$raw_weights_all,
    transport_weights_all = wt,
    max_calibration_raw =
      fit$max_raw_score,
    converged =
      isTRUE(fit$converged),
    ESS =
      transport_ess(
        wt[ID == 1]
      ),
    MaxW =
      max(
        wt[ID == 1]
      ),
    solver = fit
  )
}


###############################################################################
# K. Joint multinomial sandwich for ET/HD/EL
#
# beta = (phi_multinom, lambda, mu1, theta0)
#
# Cross-fitted GAM prediction is fixed, exactly as in simulation.
###############################################################################

GEC_transport_sandwich_simulation_style <- function(
    d4,
    donor_g,
    ps_obj,
    point_fit
) {

  entropy <- point_fit$entropy

  dat <- point_fit$data_T0
  ord <- match(
    dat$ID,
    seq_len(nrow(d4))
  )

  X <- as.matrix(
    ps_obj$X[
      ord,
      ,
      drop = FALSE
    ]
  )

  G <- as.character(dat$G)
  Y <- as.numeric(dat$re78)

  N <- nrow(dat)

  phi_hat <- as.numeric(ps_obj$phi_hat)
  pphi <- length(phi_hat)

  lambda_hat <- as.numeric(
    point_fit$lambda
  )

  q <- length(lambda_hat)

  mu1_hat <- point_fit$mu1
  theta_hat <- point_fit$theta0

  idx_phi <- seq_len(pphi)
  idx_lambda <- (pphi + 1L):(pphi + q)
  idx_mu1 <- pphi + q + 1L
  idx_theta <- pphi + q + 2L

  ID <- as.numeric(
    G == as.character(donor_g)
  )

  IT <- as.numeric(
    G %in% c("1", "2")
  )

  I1 <- as.numeric(
    G == "1"
  )

  ef <- simulation_entropy_transport_maps(
    entropy
  )

  # Evaluate at the fitted point.
  P <- multinom_prob_final(
    phi_hat,
    X
  )

  rr <- transport_ratio_final(
    P,
    donor_g
  )

  S <- build_GEC_transport_S(
    rr,
    dat$y.hat,
    entropy
  )

  eta <- as.numeric(
    S %*% lambda_hat
  )

  wraw <- numeric(N)
  wraw[ID == 1] <-
    ef$raw_weight(
      eta[ID == 1]
    )

  wdot <- numeric(N)
  wdot[ID == 1] <-
    ef$weight_deriv(
      eta[ID == 1]
    )

  #wtransport <- ID * wraw
  wtransport <- ID * (wraw - 1)
 # rhs_ind <- IT
  rhs_ind <- pmin(IT + ID, 1)
  Psi_phi <- multinom_score_final(
    phi_hat,
    X,
    G
  )

  Psi_lambda <- S *
    as.numeric(
      ID * wraw -
        rhs_ind
    )

  Psi_mu1 <- I1 *
    (Y - mu1_hat)

  Psi_theta <- wtransport *
    (Y - theta_hat)

  Psi <- cbind(
    Psi_phi,
    Psi_lambda,
    mu1 = Psi_mu1,
    theta = Psi_theta
  )

  B <- crossprod(Psi) / N

  # ----------------------
  # Analytic bread
  # ----------------------

  ptotal <- pphi + q + 2L
  A <- matrix(0, ptotal, ptotal)

  A[
    idx_phi,
    idx_phi
  ] <- multinom_information_final(
    phi_hat,
    X
  )

  # Derivative of ONLY the entropy-specific g component wrt multinomial phi.
  # This finite difference does not perturb lambda/eta, so it does not cross
  # HD/EL dual domains.
  Jg <- numDeriv::jacobian(
    func = function(phi) {

      PP <- multinom_prob_final(
        phi,
        X
      )

      rtmp <- transport_ratio_final(
        PP,
        donor_g
      )

      transport_g_simulation_style(
        rtmp,
        entropy
      )
    },
    x = phi_hat
  )

  A_ll <- matrix(0, q, q)
  A_lphi <- matrix(0, q, pphi)
  A_thetal <- matrix(0, 1, q)
  A_thetaphi <- matrix(0, 1, pphi)

  for (i in seq_len(N)) {

    si <- matrix(
      S[i, ],
      ncol = 1
    )

    Ji <- matrix(
      0,
      nrow = q,
      ncol = pphi
    )

    Ji[q, ] <- Jg[i, ]

    A_ll <- A_ll -
      ID[i] *
      wdot[i] *
      (si %*% t(si)) /
      N

    term1 <- (
      ID[i] *
        wraw[i] -
        rhs_ind[i]
    ) * Ji

    term2 <- ID[i] *
      wdot[i] *
      si %*%
      (
        t(lambda_hat) %*%
          Ji
      )

    A_lphi <- A_lphi -
      (term1 + term2) /
      N

    A_thetal <- A_thetal -
      ID[i] *
      (Y[i] - theta_hat) *
      wdot[i] *
      t(si) /
      N

    A_thetaphi <- A_thetaphi -
      ID[i] *
      (Y[i] - theta_hat) *
      wdot[i] *
      (
        t(lambda_hat) %*%
          Ji
      ) /
      N
  }

  A[
    idx_lambda,
    idx_phi
  ] <- A_lphi

  A[
    idx_lambda,
    idx_lambda
  ] <- A_ll

  A[
    idx_mu1,
    idx_mu1
  ] <- sum(I1) / N

  A[
    idx_theta,
    idx_phi
  ] <- A_thetaphi

  A[
    idx_theta,
    idx_lambda
  ] <- A_thetal

  A[
    idx_theta,
    idx_theta
  ] <- sum(wtransport) / N

  Ainv <- tryCatch(
    solve(A),
    error = function(e) NULL
  )

  if (is.null(Ainv) ||
      any(!is.finite(Ainv))) {

    return(list(
      ATE = point_fit$ATE,
      SE = NA_real_,
      variance = NA_real_,
      CI = c(NA_real_, NA_real_),
      success = FALSE,
      A = A,
      B = B,
      V = NULL,
      variance_type =
        "joint sandwich failed: singular bread"
    ))
  }

  V <- (
    Ainv %*%
      B %*%
      t(Ainv)
  ) / N

  contrast <- rep(0, ptotal)
  contrast[idx_mu1] <- 1
  contrast[idx_theta] <- -1

  vv <- as.numeric(
    t(contrast) %*%
      V %*%
      contrast
  )

  se <- if (
    is.finite(vv) &&
      vv >= -1e-8
  ) {
    sqrt(max(vv, 0))
  } else {
    NA_real_
  }

  list(
    ATE = point_fit$ATE,
    SE = se,
    variance = vv,
    CI =
      transport_ci(
        point_fit$ATE,
        se
      ),
    success =
      is.finite(se),
    ESS =
      point_fit$ESS,
    MaxW =
      point_fit$MaxW,
    CalibrationMax =
      point_fit$max_calibration_raw,
    A = A,
    B = B,
    V = V,
    variance_type =
      paste0(
        "simulation-style ",
        entropy,
        " dual + joint sandwich (multinomial phi, lambda, mu1, theta0)"
      )
  )
}


run_GEC_transport_simulation_style <- function(
    d4,
    donor_g,
    ps_obj,
    entropy = c("ET", "HD", "EL"),
    K = 4,
    seed = 2024202
) {

  entropy <- match.arg(entropy)

  pt <- GEC_transport_point_simulation_style(
    d4 = d4,
    donor_g = donor_g,
    ps_obj = ps_obj,
    entropy = entropy,
    K = K,
    seed = seed
  )

  inf <- GEC_transport_sandwich_simulation_style(
    d4 = d4,
    donor_g = donor_g,
    ps_obj = ps_obj,
    point_fit = pt
  )

  inf$point <- pt

  inf
}


###############################################################################
# L. Result row helper
###############################################################################

make_final_real_row <- function(
    donor,
    method,
    est,
    se,
    benchmark,
    ess = NA_real_,
    maxw = NA_real_,
    cal = NA_real_,
    success = TRUE,
    variance_type = ""
) {

  data.frame(
    Donor = donor,
    Method = method,
    Estimate = est,
    EvalDiff = est - benchmark,
    SE = se,
    CI_L = est - qnorm(.975) * se,
    CI_U = est + qnorm(.975) * se,
    ESS0 = ess,
    MaxW0 = maxw,
    CalibrationMax = cal,
    Success = as.integer(isTRUE(success)),
    VarianceType = variance_type,
    stringsAsFactors = FALSE
  )
}



###############################################################################
# PACKAGE COMPARATORS FOR REAL-DATA NSW TRANSPORT
#
# These reproduce the PACKAGE / INFERENCE style used in the simulation:
#
# EBPS  : WeightIt::weightit(method = "ebal") + lm_weightit
# CBPS  : WeightIt::weightit(method = "cbps") + lm_weightit(vcov = "asympt")
# oCBPS : CBPS::CBPS() + CBPS::AsyVar(method = "oCBPS")
# EBCW  : ATE::ATE()
#
# External donors:
#   S = 1 : combined NSW target, G in {1,2}
#   S = 0 : donor, G = 3 or 4
#
# NSW donor:
#   S = 1 : NSW treated, G = 1
#   S = 0 : NSW control, G = 2
###############################################################################


###############################################################################
# 1. Common result helper
###############################################################################

transport_package_result <- function(
    method,
    ate,
    se,
    dat,
    weights = NULL,
    benchmark = NA_real_,
    variance_type = ""
) {
  
  ESS1 <- NA_real_
  ESS0 <- NA_real_
  MaxW1 <- NA_real_
  MaxW0 <- NA_real_
  
  if (!is.null(weights) &&
      length(weights) == nrow(dat) &&
      all(is.finite(weights))) {
    
    w1 <- as.numeric(
      weights[dat$S == 1]
    )
    
    w0 <- as.numeric(
      weights[dat$S == 0]
    )
    
    if (length(w1) > 0 &&
        sum(w1^2) > 0) {
      
      ESS1 <-
        sum(w1)^2 /
        sum(w1^2)
      
      MaxW1 <-
        max(w1)
    }
    
    if (length(w0) > 0 &&
        sum(w0^2) > 0) {
      
      ESS0 <-
        sum(w0)^2 /
        sum(w0^2)
      
      MaxW0 <-
        max(w0)
    }
  }
  
  ci <- if (
    is.finite(ate) &&
    is.finite(se)
  ) {
    
    ate +
      c(-1, 1) *
      qnorm(0.975) *
      se
    
  } else {
    
    c(
      NA_real_,
      NA_real_
    )
  }
  
  list(
    Method = method,
    ATE = as.numeric(ate),
    SE = as.numeric(se),
    variance =
      if (is.finite(se)) {
        se^2
      } else {
        NA_real_
      },
    CI = ci,
    EvalDiff =
      if (is.finite(benchmark)) {
        as.numeric(ate) -
          benchmark
      } else {
        NA_real_
      },
    ESS = ESS0,
    MaxW = MaxW0,
    ESS1 = ESS1,
    ESS0 = ESS0,
    MaxW1 = MaxW1,
    MaxW0 = MaxW0,
    CalibrationMax = NA_real_,
    success =
      is.finite(ate) &&
      is.finite(se),
    variance_type =
      variance_type
  )
}


###############################################################################
# 2. Construct source-transport data
###############################################################################

make_transport_source_data <- function(
    d4,
    donor_g,
    xvars
) {
  
  donor_g <-
    as.character(
      donor_g
    )
  
  Gchar <-
    as.character(
      d4$G
    )
  
  
  ###########################################################################
  # NSW comparison: ordinary randomized ATE
  ###########################################################################
  
  if (donor_g == "2") {
    
    keep <-
      Gchar %in%
      c(
        "1",
        "2"
      )
    
    dat <-
      d4[
        keep,
        ,
        drop = FALSE
      ]
    
    Gk <-
      as.character(
        dat$G
      )
    
    dat$S <-
      as.integer(
        Gk == "1"
      )
    
    dat$ystar <-
      as.numeric(
        dat$re78
      )
    
    return(
      list(
        data = dat,
        xvars = xvars,
        donor_g = donor_g,
        package_estimand = "ATE"
      )
    )
  }
  
  
  ###########################################################################
  # PSID / CPS:
  #
  # target = combined NSW = G 1 + G 2
  # donor  = G 3 or G 4
  ###########################################################################
  
  keep <-
    Gchar %in%
    c(
      "1",
      "2",
      donor_g
    )
  
  dat <-
    d4[
      keep,
      ,
      drop = FALSE
    ]
  
  Gk <-
    as.character(
      dat$G
    )
  
  dat$S <-
    as.integer(
      Gk %in%
        c(
          "1",
          "2"
        )
    )
  
  n_target <-
    sum(
      dat$S == 1
    )
  
  n_treated <-
    sum(
      Gk == "1"
    )
  
  
  ###########################################################################
  # Transformed outcome
  #
  # Mean(Ystar | S=1)
  # =
  # Mean(re78 | G=1)
  ###########################################################################
  
  dat$ystar <- 0
  
  dat$ystar[
    Gk == "1"
  ] <-
    (
      n_target /
        n_treated
    ) *
    dat$re78[
      Gk == "1"
    ]
  
  dat$ystar[
    Gk == donor_g
  ] <-
    dat$re78[
      Gk == donor_g
    ]
  
  
  ###########################################################################
  # Numerical identity check
  ###########################################################################
  
  stopifnot(
    abs(
      mean(
        dat$ystar[
          dat$S == 1
        ]
      ) -
        mean(
          dat$re78[
            Gk == "1"
          ]
        )
    ) <
      1e-8
  )
  
  
  list(
    data = dat,
    xvars = xvars,
    donor_g = donor_g,
    package_estimand = "ATT"
  )
}


###############################################################################
# 3. EBPS
###############################################################################

transport_ebps <- function(
    source_object,
    benchmark = NA_real_
) {
  
  if (!requireNamespace(
    "WeightIt",
    quietly = TRUE
  )) {
    stop(
      "Package 'WeightIt' is required."
    )
  }
  
  dat <-
    source_object$data
  
  xvars <-
    source_object$xvars
  
  form <-
    reformulate(
      xvars,
      response = "S"
    )
  
  estimand_use <-
    source_object$package_estimand
  
  
  fit_w <-
    WeightIt::weightit(
      form,
      data = dat,
      method = "ebal",
      estimand =
        estimand_use
    )
  
  
  fit_y <-
    WeightIt::lm_weightit(
      ystar ~ S,
      data = dat,
      weightit = fit_w
    )
  
  
  ate <-
    as.numeric(
      coef(
        fit_y
      )["S"]
    )
  
  se <-
    sqrt(
      vcov(
        fit_y
      )[
        "S",
        "S"
      ]
    )
  
  
  transport_package_result(
    method = "EBPS",
    ate = ate,
    se = se,
    dat = dat,
    weights =
      as.numeric(
        fit_w$weights
      ),
    benchmark = benchmark,
    variance_type =
      "WeightIt ebal + lm_weightit"
  )
}





###############################################################################
# 4. CBPS
###############################################################################

transport_cbps <- function(
    source_object,
    benchmark = NA_real_
) {
  
  if (!requireNamespace(
    "WeightIt",
    quietly = TRUE
  )) {
    stop(
      "Package 'WeightIt' is required."
    )
  }
  
  dat <-
    source_object$data
  
  xvars <-
    source_object$xvars
  
  form <-
    reformulate(
      xvars,
      response = "S"
    )
  
  estimand_use <-
    source_object$package_estimand
  
  
  fit_w <-
    WeightIt::weightit(
      form,
      data = dat,
      method = "cbps",
      estimand =
        estimand_use,
      over = FALSE
    )
  
  
  fit_y <-
    WeightIt::lm_weightit(
      ystar ~ S,
      data = dat,
      weightit = fit_w,
      vcov = "asympt"
    )
  
  
  ate <-
    as.numeric(
      coef(
        fit_y
      )["S"]
    )
  
  se <-
    sqrt(
      vcov(
        fit_y
      )[
        "S",
        "S"
      ]
    )
  
  
  transport_package_result(
    method = "CBPS",
    ate = ate,
    se = se,
    dat = dat,
    weights =
      as.numeric(
        fit_w$weights
      ),
    benchmark = benchmark,
    variance_type =
      "WeightIt CBPS + lm_weightit(vcov='asympt')"
  )
}







###############################################################################
# 6. EBCW
###############################################################################

transport_ebcw <- function(
    source_object,
    benchmark = NA_real_
) {
  
  if (!requireNamespace(
    "ATE",
    quietly = TRUE
  )) {
    stop(
      "Package 'ATE' is required."
    )
  }
  
  dat <-
    source_object$data
  
  xvars <-
    source_object$xvars
  
  
  X_ebcw <-
    data.frame(
      dat[
        ,
        xvars,
        drop = FALSE
      ]
    )
  
  
  use_ATT <-
    identical(
      source_object$package_estimand,
      "ATT"
    )
  
  
  fit <-
    ATE::ATE(
      as.numeric(
        dat$ystar
      ),
      as.integer(
        dat$S
      ),
      X_ebcw,
      ATT = use_ATT
    )
  
  
  sm <-
    summary(
      fit
    )
  
  
  estmat <-
    as.matrix(
      sm$Estimate
    )
  
  
  ###########################################################################
  # Use ATE row for NSW;
  # ATT row for PSID/CPS
  ###########################################################################
  
  target_row <-
    if (use_ATT) {
      "ATT"
    } else {
      "ATE"
    }
  
  
  if (
    !is.null(
      rownames(
        estmat
      )
    ) &&
    target_row %in%
    rownames(
      estmat
    )
  ) {
    
    row_id <-
      which(
        rownames(
          estmat
        ) ==
          target_row
      )[1]
    
  } else {
    
    stop(
      paste0(
        "Could not find ",
        target_row,
        " row in ATE::ATE summary."
      )
    )
  }
  
  
  ate <-
    as.numeric(
      estmat[
        row_id,
        "Estimate"
      ]
    )
  
  
  if (
    "Std. Error" %in%
    colnames(
      estmat
    )
  ) {
    
    se <-
      as.numeric(
        estmat[
          row_id,
          "Std. Error"
        ]
      )
    
  } else {
    
    # fallback for package versions where SE is column 2
    se <-
      as.numeric(
        estmat[
          row_id,
          2
        ]
      )
  }
  
  
  ###########################################################################
  # Extract package weights where available
  ###########################################################################
  
  w <- NULL
  
  if (
    !is.null(
      sm$weights.treat
    ) &&
    !is.null(
      sm$weights.placebo
    )
  ) {
    
    w <-
      as.numeric(
        sm$weights.treat
      ) +
      as.numeric(
        sm$weights.placebo
      )
  }
  
  
  transport_package_result(
    method = "EBCW",
    ate = ate,
    se = se,
    dat = dat,
    weights = w,
    benchmark = benchmark,
    variance_type =
      if (use_ATT) {
        "ATE::ATE package ATT inference"
      } else {
        "ATE::ATE package ATE inference"
      }
  )
}











transport_ocbps <- function(
    source_object,
    benchmark = NA_real_
) {
  
  if (!requireNamespace(
    "CBPS",
    quietly = TRUE
  )) {
    stop("Package 'CBPS' is required.")
  }
  
  dat <- source_object$data
  xvars <- source_object$xvars
  
  ###########################################################################
  # oCBPS in CBPSOptimal currently supports ATE (ATT = 0),
  # but not ATT = 1.
  #
  # Therefore it can be used for the ordinary NSW randomized comparison,
  # but NOT for PSID/CPS source-transport targeting combined NSW.
  ###########################################################################
  
  if (!identical(
    source_object$package_estimand,
    "ATE"
  )) {
    
    return(
      transport_package_result(
        method = "oCBPS",
        ate = NA_real_,
        se = NA_real_,
        dat = dat,
        weights = NULL,
        benchmark = benchmark,
        variance_type =
          "Not estimated: CBPSOptimal oCBPS does not support ATT=1"
      )
    )
  }
  
  
  ###########################################################################
  # NSW only: ordinary ATE
  ###########################################################################
  
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
  
  
  fit <- CBPS::CBPS(
    dat$S ~ X2,
    ATT = 0,
    method = "exact",
    baseline.formula = ~ X2,
    diff.formula = ~ X2
  )
  
  
  w <- as.numeric(
    fit$weights
  )
  
  
  mu1 <- sum(
    w[dat$S == 1] *
      dat$ystar[dat$S == 1]
  ) /
    sum(
      w[dat$S == 1]
    )
  
  
  mu0 <- sum(
    w[dat$S == 0] *
      dat$ystar[dat$S == 0]
  ) /
    sum(
      w[dat$S == 0]
    )
  
  
  ate <- mu1 - mu0
  
  
  inf <- CBPS::AsyVar(
    Y = dat$ystar,
    CBPS_obj = fit,
    method = "oCBPS",
    CI = 0.95
  )
  
  
  se <- if (
    !is.null(inf$std.err)
  ) {
    as.numeric(
      inf$std.err
    )[1]
  } else {
    NA_real_
  }
  
  
  transport_package_result(
    method = "oCBPS",
    ate = ate,
    se = se,
    dat = dat,
    weights = w,
    benchmark = benchmark,
    variance_type =
      "CBPSOptimal oCBPS ATE + AsyVar"
  )
}





















