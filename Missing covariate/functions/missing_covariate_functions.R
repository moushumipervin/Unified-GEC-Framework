
###############################################################################
# missing_covariate_GEC_clean.R
#
# CLEAN missing-covariate implementation
#
# Keeps only:
#   - data generation
#   - cross-fitting z.hat
#   - IPW point estimate + joint sandwich variance
#   - AIPW point estimate + joint sandwich variance
#   - GEC ET / HD / CE point estimates
#   - GEC joint sandwich variance
#   - simulation-output helper and Monte Carlo summary
#
# Removed on purpose:
#   - variance_theta_diag()
#   - variance_theta_diag1()
#   - old CVXR estimators
#   - nested/no-fold estimators
#   - duplicate ET/HD/CE solvers
#   - IF-based alternative variance
###############################################################################

if (!requireNamespace("numDeriv", quietly = TRUE)) {
  stop("Please install numDeriv: install.packages('numDeriv')")
}

if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("Please install MASS: install.packages('MASS')")
}


###############################################################################
# 1. Generate data
#
# OR1 is aligned with the manuscript:
#   y = 1 + x + 2z + error
#
# OR2 target used in the simulation summary should be
#   beta = (0.500, 0.795, 1.000)
###############################################################################

generate_data <- function(n, OR, PS) {

  x <- rnorm(n, 0, 1)
  z <- rbinom(n, 1, 0.5)

  if (OR == 1) {

    y <- 1 + x + 2 * z + rnorm(n, 0, 1)

  } else {

    y <- 0.5 +
      2 * sin(pi * x) -
      1.5 * cos(2 * pi * x) +
      0.25 * x^3 +
      z +
      rnorm(n, 0, 2)
  }

  if (PS == 1) {

    eta <- -.25 + 0.5 * y + 0.5 * x

  } else {

    eta <- -0.25
  }

  pi_true <- plogis(eta)

  D <- rbinom(
    n,
    1,
    pi_true
  )

  # Same observed-data propensity specification used by variance code.
  ps_fit <- glm(
    D ~ y + x,
    family = binomial()
  )

  D.hat <-fitted(ps_fit)

  dat <- data.frame(
    x = x,
    y = y,
    z = z,
    D = D,
    D.hat = D.hat,
    ID = seq_len(n)
  )

  dat
}

safe_jacobian <- function(func, x) {
  
  numDeriv::jacobian(
    func = func,
    x = x
  )
}
###############################################################################
# 2. Fold construction
###############################################################################

build_SU_folds <- function(
    df,
    K,
    id_col = "ID",
    delta_col = "D",
    seed = 1234) {

  stopifnot(
    all(
      c(
        id_col,
        delta_col
      ) %in%
        names(df)
    )
  )

  set.seed(seed)

  df[[id_col]] <- seq_len(
    nrow(df)
  )

  id <- df[[id_col]]
  delta <- df[[delta_col]]

  idx_S <- which(
    delta == 1
  )

  idx_U0 <- which(
    delta == 0
  )

  if (length(idx_S) < K ||
      length(idx_U0) < K) {

    stop(
      "Need at least K complete and K incomplete observations."
    )
  }

  split_K <- function(idx, K) {

    idx <- sample(idx)

    split(
      idx,
      rep(
        seq_len(K),
        length.out = length(idx)
      )
    )
  }

  S_parts <- split_K(
    idx_S,
    K
  )

  U0_parts <- split_K(
    idx_U0,
    K
  )

  folds <- vector(
    "list",
    K
  )

  for (k in seq_len(K)) {

    S_k_idx <- S_parts[[k]]
    U0_k_idx <- U0_parts[[k]]

    U_k_idx <- c(
      S_k_idx,
      U0_k_idx
    )

    folds[[k]] <- list(
      S_k_idx = S_k_idx,
      U_k_idx = U_k_idx,
      S_k_ids = id[S_k_idx],
      U_k_ids = id[U_k_idx]
    )
  }

  list(
    folds = folds
  )
}


###############################################################################
# 3. Cross-fitting z.hat
#
# z is observed only when D=1.
# Model z | (x,y) is trained only on complete cases.
###############################################################################

k_fold_function <- function(
    df,
    K,
    seed = 1234) {

  folds <- build_SU_folds(
    df = df,
    K = K,
    id_col = "ID",
    delta_col = "D",
    seed = seed
  )

  complete_idx <- which(
    df$D == 1
  )

  out <- vector(
    "list",
    K
  )

  for (k in seq_len(K)) {

    train_idx <- setdiff(
      complete_idx,
      folds$folds[[k]]$S_k_idx
    )

    valid_idx <- folds$folds[[k]]$U_k_idx

    train_dat <- df[
      train_idx,
      ,
      drop = FALSE
    ]

    valid_dat <- df[
      valid_idx,
      ,
      drop = FALSE
    ]

    z_fit <- glm(
      z ~ x + y,
      family = binomial(),
      data = train_dat
    )

    valid_dat$z.hat <- pmin(
      pmax(
        predict(
          z_fit,
          newdata = valid_dat,
          type = "response"
        ),
        1e-8
      ),
      1 - 1e-8
    )

    out[[k]] <- valid_dat
  }

  out
}


###############################################################################
# 4. Prepare one cross-fitted stacked data set
###############################################################################

prepare_missingcov_data <- function(
    data_full,
    K = 3,
    seed = 1234) {

  folds <- k_fold_function(
    df = data_full,
    K = K,
    seed = seed
  )

  data_all <- do.call(
    rbind,
    folds
  )

  rownames(data_all) <- NULL

  # Refit the observed-data propensity model on the same stacked data.
  ps_fit <- glm(
    D ~ x + y,
    family = binomial(),
    data = data_all
  )

  data_all$pi.hat <- pmin(
    pmax(
      fitted(ps_fit),
      1e-8
    ),
    1 - 1e-8
  )

  data_all
}


###############################################################################
# 5. Full-data score U(beta)
#
# TRUE z is used only when D=1.
###############################################################################

missingcov_U <- function(
    beta,
    data_all) {

  beta <- as.numeric(beta)

  Xtrue <- cbind(
    1,
    x = data_all$x,
    z = data_all$z
  )

  U <- matrix(
    0,
    nrow = nrow(data_all),
    ncol = 3
  )

  I1 <- which(
    data_all$D == 1
  )

  resid <- data_all$y[I1] -
    as.numeric(
      Xtrue[
        I1,
        ,
        drop = FALSE
      ] %*%
        beta
    )

  U[I1, ] <-
    Xtrue[
      I1,
      ,
      drop = FALSE
    ] *
    resid

  colnames(U) <- c(
    "beta0",
    "beta1",
    "beta2"
  )

  U
}


###############################################################################
# 6. Plug-in b(beta)
#
# This is the current A4 construction (a):
# U(beta; x, z.hat, y)
###############################################################################

missingcov_b_plugin <- function(
    beta,
    data_all) {

  beta <- as.numeric(beta)

  Xhat <- cbind(
    1,
    x = data_all$x,
    z = data_all$z.hat
  )

  resid_hat <- data_all$y -
    as.numeric(
      Xhat %*%
        beta
    )

  b <- Xhat *
    resid_hat

  colnames(b) <- c(
    "b0",
    "b1",
    "b2"
  )

  b
}


###############################################################################
# 7. Propensity model helper
###############################################################################

fit_ps_missingcov <- function(
    data_all) {

  ps_fit <- glm(
    D ~ x + y,
    family = binomial(),
    data = data_all,
    x = TRUE
  )

  list(
    fit = ps_fit,
    phi = as.numeric(
      coef(ps_fit)
    ),
    O = ps_fit$x,
    pi = pmin(
      pmax(
        fitted(ps_fit),
        1e-8
      ),
      1 - 1e-8
    )
  )
}

estimate_ipw_missingcov <- function(data_all) {
  
  D <- as.numeric(data_all$D)
  y <- as.numeric(data_all$y)
  
 
  
  Xtrue <- cbind(
    1,
    x = data_all$x,
    z = data_all$z
  )
  
  I1 <- which(D == 1)
  
  X1 <- Xtrue[I1, , drop = FALSE]
  Y1 <- y[I1]
  pi_hat <- pmin(
    pmax(data_all$pi.hat, 1e-8),
    1 - 1e-8
  )
  w  <- 1 / pi_hat[I1]
  
  theta_hat <- as.numeric(
    solve(
      crossprod(X1, w * X1),
      crossprod(X1, w * Y1)
    )
  )
  
  list(
    theta = matrix(theta_hat, ncol = 1),
    w = w,
    data_all = data_all,
    converged = TRUE
  )
}

estimate_aipw_missingcov <- function(data_all) {
  
  D <- as.numeric(data_all$D)
  y <- as.numeric(data_all$y)
  
  
  pi_hat <- pmin(
    pmax(data_all$pi.hat, 1e-8),
    1 - 1e-8
  )
  Xhat <- cbind(
    1,
    x = data_all$x,
    z = data_all$z.hat
  )
  
  Xtrue <- cbind(
    1,
    x = data_all$x,
    z = data_all$z
  )
  
  I1 <- which(D == 1)
  
  ## Start with b(beta) contribution from everyone
  A <- crossprod(Xhat, Xhat)
  B <- crossprod(Xhat, y)
  
  ## Add D/pi * {U(beta) - b(beta)} contribution
  for (i in I1) {
    
    wi <- 1 / pi_hat[i]
    
    xt <- Xtrue[i, ]
    xh <- Xhat[i, ]
    
    A <- A +
      wi * tcrossprod(xt) -
      wi * tcrossprod(xh)
    
    B <- B +
      wi * xt * y[i] -
      wi * xh * y[i]
  }
  
  theta_hat <- as.numeric(
    solve(A, B)
  )
  
  list(
    theta = matrix(theta_hat, ncol = 1),
    data_all = data_all,
    converged = TRUE
  )
}
###############################################################################
# 8. Entropy functions used in this missing-covariate analysis
###############################################################################

gec_entropy_missingcov <- function(
    entropy = c(
      "ET",
      "HD",
      "CE"
    )) {

  entropy <- match.arg(
    entropy
  )

  if (entropy == "ET") {

    ginv <- function(eta) {
      exp(eta)
    }

    debias <- function(pi) {
      log(1 / pi)
    }

    valid_eta <- function(eta) {
      rep(
        TRUE,
        length(eta)
      )
    }
  }

  if (entropy == "HD") {

    ginv <- function(eta) {
      1 / (
        4 * eta^2
      )
    }

    debias <- function(pi) {
      -sqrt(pi) / 2
    }

    valid_eta <- function(eta) {
      eta < 0
    }
  }

  if (entropy == "CE") {

    ginv <- function(eta) {
      1 / (
        1 - exp(eta)
      )
    }

    debias <- function(pi) {
      log(1 - pi)
    }

    valid_eta <- function(eta) {
      eta < 0
    }
  }

  list(
    entropy = entropy,
    ginv = ginv,
    debias = debias,
    valid_eta = valid_eta
  )
}


###############################################################################
# 9. ET dual solver
###############################################################################

solve_lambda_ET_dual <- function(
    b_mat,
    pi_hat,
    D,
    lambda_start = NULL,
    maxit = 1000,
    reltol = 1e-10) {

  pi_hat <- pmin(
    pmax(
      pi_hat,
      1e-8
    ),
    1 - 1e-8
  )

  g_pi <- log(
    1 / pi_hat
  )

  S <- cbind(
    b_mat,
    g_pi
  )

  q <- ncol(S)

  if (is.null(lambda_start) ||
      length(lambda_start) != q) {

    lambda_start <- c(
      rep(
        0,
        q - 1
      ),
      1
    )
  }

  objective <- function(lambda) {

    eta <- as.numeric(
      S %*%
        lambda
    )

    if (any(!is.finite(eta)) ||
        max(eta) > 700) {

      return(1e100)
    }

    mean(
      D * exp(eta) -
        eta
    )
  }

  gradient <- function(lambda) {

    eta <- as.numeric(
      S %*%
        lambda
    )

    if (any(!is.finite(eta)) ||
        max(eta) > 700) {

      return(
        rep(
          1e20,
          q
        )
      )
    }

    w_all <- exp(
      eta
    )

    colMeans(
      S *
        as.numeric(
          D * w_all - 1
        )
    )
  }

  fit <- optim(
    par = as.numeric(lambda_start),
    fn = objective,
    gr = gradient,
    method = "BFGS",
    control = list(
      maxit = maxit,
      reltol = reltol
    )
  )

  lambda_hat <- as.numeric(
    fit$par
  )

  eta <- as.numeric(
    S %*%
      lambda_hat
  )

  w_all <- exp(
    eta
  )

  score <- gradient(
    lambda_hat
  )

  list(
    lambda = lambda_hat,
    weights = w_all[D == 1],
    weights_all = w_all,
    eta = eta,
    score = score,
    max_score = max(
      abs(score)
    ),
    converged =
      fit$convergence == 0,
    convergence_code =
      fit$convergence,
    S = S,
    g_pi = g_pi
  )
}


###############################################################################
# 10. HD dual solver
###############################################################################

solve_lambda_HD_dual <- function(
    b_mat,
    pi_hat,
    D,
    lambda_start = NULL,
    margin = 1e-8,
    maxit = 1000,
    reltol = 1e-10) {

  pi_hat <- pmin(
    pmax(
      pi_hat,
      1e-8
    ),
    1 - 1e-8
  )

  g_pi <- -sqrt(
    pi_hat
  ) / 2

  S <- cbind(
    b_mat,
    g_pi
  )

  q <- ncol(S)

  I1 <- which(
    D == 1
  )

  lambda_natural <- c(
    rep(
      0,
      q - 1
    ),
    1
  )

  if (is.null(lambda_start) ||
      length(lambda_start) != q) {

    lambda_start <-
      lambda_natural
  }

  eta_start <- as.numeric(
    S[
      I1,
      ,
      drop = FALSE
    ] %*%
      lambda_start
  )

  if (any(
    eta_start >= -margin
  )) {

    lambda_start <-
      lambda_natural
  }

  objective <- function(lambda) {

    eta <- as.numeric(
      S %*%
        lambda
    )

    eta1 <- eta[I1]

    if (any(!is.finite(eta1)) ||
        any(eta1 >= 0)) {

      return(1e100)
    }

    F_eta1 <- -1 / (
      4 * eta1
    )

    (
      sum(F_eta1) -
        sum(eta)
    ) /
      nrow(S)
  }

  gradient <- function(lambda) {

    eta <- as.numeric(
      S %*%
        lambda
    )

    eta1 <- eta[I1]

    if (any(!is.finite(eta1)) ||
        any(eta1 >= 0)) {

      return(
        rep(
          1e20,
          q
        )
      )
    }

    w1 <- 1 / (
      4 * eta1^2
    )

    multiplier <- rep(
      -1,
      nrow(S)
    )

    multiplier[I1] <-
      w1 - 1

    colMeans(
      S * multiplier
    )
  }

  ui <- -S[
    I1,
    ,
    drop = FALSE
  ]

  ci <- rep(
    margin,
    length(I1)
  )

  fit <- constrOptim(
    theta = as.numeric(
      lambda_start
    ),
    f = objective,
    grad = gradient,
    ui = ui,
    ci = ci,
    method = "BFGS",
    control = list(
      maxit = maxit,
      reltol = reltol
    )
  )

  lambda_hat <- as.numeric(
    fit$par
  )

  eta <- as.numeric(
    S %*%
      lambda_hat
  )

  w_all <- rep(
    NA_real_,
    nrow(S)
  )

  w_all[I1] <- 1 / (
    4 * eta[I1]^2
  )

  score <- gradient(
    lambda_hat
  )

  list(
    lambda = lambda_hat,
    weights = w_all[I1],
    weights_all = w_all,
    eta = eta,
    score = score,
    max_score = max(
      abs(score)
    ),
    converged =
      fit$convergence == 0,
    convergence_code =
      fit$convergence,
    S = S,
    g_pi = g_pi
  )
}

get_b_A4 <- function(
    beta,
    data_all,
    b_method = c("plugin", "exact", "direct"),
    fold_id = NULL,
    K = 3) {
  
  b_method <- match.arg(b_method)
  
  if (b_method == "plugin") {
    
    return(
      b_plugin_A4(
        beta = beta,
        data_all = data_all
      )
    )
  }
  
  if (b_method == "exact") {
    
    return(
      b_exact_A4(
        beta = beta,
        data_all = data_all
      )
    )
  }
  
  if (is.null(fold_id)) {
    stop("fold_id is required for direct-score regression.")
  }
  
  b_direct_A4(
    beta = beta,
    data_all = data_all,
    fold_id = fold_id,
    K = K
  )
}


###############################################################################
# A4-SPECIFIC GEC ESTIMATOR
#
# Separate from original implementation.
#
# Allows comparison of:
#   b_method = "plugin"
#   b_method = "exact"
#   b_method = "direct"
#
# for entropy:
#   "ET", "HD", "CE"
#
# Your original estimate_theta_missingcov_GEC() is NOT changed.
###############################################################################

estimate_theta_missingcov_GEC_A4_compare <- function(
    th,
    data_full,
    K = 3,
    seed = 1234,
    entropy = c("ET", "HD", "CE"),
    b_method = c("plugin", "exact", "direct"),
    max.iter = 100,
    eps = 1e-6,
    lambda_maxit = 1000,
    lambda_tol = 1e-10,
    damping = 1) {
  
  ###########################################################################
  # 1. Match arguments
  ###########################################################################
  
  entropy <- match.arg(entropy)
  b_method <- match.arg(b_method)
  
  
  ###########################################################################
  # 2. Prepare cross-fitted missing-covariate data
  #
  # This gives:
  #   z.hat
  #   pi.hat
  #   D
  #   x
  #   y
  #   z
  ###########################################################################
  
  data_all <- prepare_missingcov_data(
    data_full = data_full,
    K = K,
    seed = seed
  )
  
  
  ###########################################################################
  # 3. Basic objects
  ###########################################################################
  
  D <- as.numeric(
    data_all$D
  )
  
  I1 <- which(
    D == 1
  )
  
  y <- data_all$y
  
  pi_hat <- data_all$pi.hat
  
  
  ###########################################################################
  # 4. True regression design among complete cases
  #
  # True z is used ONLY for complete cases in the beta estimating equation.
  ###########################################################################
  
  Xtrue <- cbind(
    1,
    x = data_all$x,
    z = data_all$z
  )
  
  Q <- Xtrue[
    I1,
    ,
    drop = FALSE
  ]
  
  Y1 <- y[I1]
  
  
  ###########################################################################
  # 5. Fixed folds for direct-score regression
  #
  # IMPORTANT:
  # These folds are created ONCE and kept fixed throughout all theta iterations.
  ###########################################################################
  
  direct_fold <- NULL
  
  if (b_method == "direct") {
    
    direct_fold <- make_A4_fold_id(
      data_all = data_all,
      K = K,
      seed = seed + 50000
    )
  }
  
  
  ###########################################################################
  # 6. Starting theta
  ###########################################################################
  
  theta <- as.numeric(
    th
  )
  
  if (length(theta) != 3) {
    stop("Starting theta must have length 3.")
  }
  
  
  ###########################################################################
  # 7. Current lambda
  ###########################################################################
  
  lambda_current <- NULL
  
  
  ###########################################################################
  # 8. Lambda solver depending on entropy
  ###########################################################################
  
  solve_current_lambda <- function(
    b_mat,
    lambda_start) {
    
    if (entropy == "ET") {
      
      return(
        solve_lambda_ET_dual(
          b_mat = b_mat,
          pi_hat = pi_hat,
          D = D,
          lambda_start = lambda_start,
          maxit = lambda_maxit,
          reltol = lambda_tol
        )
      )
    }
    
    
    if (entropy == "HD") {
      
      return(
        solve_lambda_HD_dual(
          b_mat = b_mat,
          pi_hat = pi_hat,
          D = D,
          lambda_start = lambda_start,
          maxit = lambda_maxit,
          reltol = lambda_tol
        )
      )
    }
    
    
    if (entropy == "CE") {
      
      return(
        solve_lambda_CE_dual(
          b_mat = b_mat,
          pi_hat = pi_hat,
          D = D,
          lambda_start = lambda_start,
          maxit = lambda_maxit,
          reltol = lambda_tol
        )
      )
    }
    
    stop("Unknown entropy.")
  }
  
  
  ###########################################################################
  # 9. Main iteration
  ###########################################################################
  
  diff_theta <- Inf
  iter <- 0
  
  
  for (iter in seq_len(max.iter)) {
    
    #########################################################################
    # 9a. Construct b(beta)
    #
    # plugin:
    #   U(beta; x, z.hat, y)
    #
    # exact:
    #   E{U(beta) | x,y}
    #
    # direct:
    #   cross-fitted direct regression of U(beta) on x,y
    #########################################################################
    
    b_mat <- get_b_A4(
      beta = theta,
      data_all = data_all,
      b_method = b_method,
      fold_id = direct_fold,
      K = K
    )
    
    
    if (any(!is.finite(b_mat))) {
      stop("Non-finite b(beta) encountered.")
    }
    
    
    #########################################################################
    # 9b. Solve calibration equation for lambda
    #########################################################################
    
    lambda_fit <- solve_current_lambda(
      b_mat = b_mat,
      lambda_start = lambda_current
    )
    
    
    if (is.null(lambda_fit) ||
        is.null(lambda_fit$lambda) ||
        any(!is.finite(lambda_fit$lambda))) {
      
      stop("Lambda solver failed.")
    }
    
    
    lambda_current <- as.numeric(
      lambda_fit$lambda
    )
    
    
    #########################################################################
    # 9c. Extract complete-case calibration weights
    #########################################################################
    
    W <- as.numeric(
      lambda_fit$weights
    )
    
    
    if (length(W) != length(I1)) {
      stop("Length of calibration weights does not match complete cases.")
    }
    
    
    if (any(!is.finite(W))) {
      stop("Non-finite calibration weights encountered.")
    }
    
    
    #########################################################################
    # 9d. Weighted regression update for beta
    #
    # Solve:
    #
    # sum_i D_i w_i X_i (Y_i - X_i^T beta) = 0
    #########################################################################
    
    A_beta <- crossprod(
      Q,
      W * Q
    )
    
    B_beta <- crossprod(
      Q,
      W * Y1
    )
    
    
    theta_raw <- tryCatch(
      
      as.numeric(
        solve(
          A_beta,
          B_beta
        )
      ),
      
      error = function(e) {
        
        as.numeric(
          MASS::ginv(
            A_beta
          ) %*%
            B_beta
        )
      }
    )
    
    
    if (any(!is.finite(theta_raw))) {
      stop("Non-finite theta update encountered.")
    }
    
    
    #########################################################################
    # 9e. Convergence difference
    #########################################################################
    
    diff_theta <- max(
      abs(
        theta_raw -
          theta
      )
    )
    
    
    #########################################################################
    # 9f. Damped update
    #########################################################################
    
    theta_new <- as.numeric(
      (1 - damping) *
        theta +
        damping *
        theta_raw
    )
    
    
    theta <- theta_new
    
    
    #########################################################################
    # 9g. Stop if converged
    #########################################################################
    
    if (diff_theta < eps) {
      break
    }
  }
  
  
  ###########################################################################
  # 10. Convergence status
  ###########################################################################
  
  converged <- is.finite(diff_theta) &&
    diff_theta < eps
  
  
  ###########################################################################
  # 11. FINAL SYNCHRONIZATION
  #
  # Recompute b and lambda at final returned theta.
  #
  # This is very important.
  ###########################################################################
  
  b_final <- get_b_A4(
    beta = theta,
    data_all = data_all,
    b_method = b_method,
    fold_id = direct_fold,
    K = K
  )
  
  
  if (any(!is.finite(b_final))) {
    stop("Non-finite final b(beta).")
  }
  
  
  lambda_final_fit <- solve_current_lambda(
    b_mat = b_final,
    lambda_start = lambda_current
  )
  
  
  if (is.null(lambda_final_fit) ||
      is.null(lambda_final_fit$lambda) ||
      any(!is.finite(lambda_final_fit$lambda))) {
    
    stop("Final lambda synchronization failed.")
  }
  
  
  ###########################################################################
  # 12. Final weights
  ###########################################################################
  
  W_final <- as.numeric(
    lambda_final_fit$weights
  )
  
  weights_all_final <- as.numeric(
    lambda_final_fit$weights_all
  )
  
  
  ###########################################################################
  # 13. Check final beta estimating equation
  ###########################################################################
  
  U_final <- missingcov_U(
    beta = theta,
    data_all = data_all
  )
  
  
  beta_score <- colMeans(
    U_final *
      as.numeric(
        D *
          weights_all_final
      ),
    na.rm = TRUE
  )
  
  
  max_beta_score <- max(
    abs(
      beta_score
    )
  )
  
  
  ###########################################################################
  # 14. Calibration residual
  #
  # Use the solver's own residual if available.
  ###########################################################################
  
  calibration_residual <- NA_real_
  
  if (!is.null(
    lambda_final_fit$calibration_residual
  )) {
    
    calibration_residual <- as.numeric(
      lambda_final_fit$calibration_residual
    )
    
  } else if (!is.null(
    lambda_final_fit$residual
  )) {
    
    calibration_residual <- as.numeric(
      lambda_final_fit$residual
    )
  }
  
  
  ###########################################################################
  # 15. Return
  ###########################################################################
  
  list(
    
    theta = matrix(
      theta,
      ncol = 1
    ),
    
    beta = theta,
    
    entropy = entropy,
    
    b_method = b_method,
    
    w = W_final,
    
    weights_all = weights_all_final,
    
    lambda = as.numeric(
      lambda_final_fit$lambda
    ),
    
    b_final = b_final,
    
    data_all = data_all,
    
    direct_fold = direct_fold,
    
    converged = converged,
    
    iterations = iter,
    
    diff_theta = diff_theta,
    
    beta_score = beta_score,
    
    max_beta_score = max_beta_score,
    
    calibration_residual = calibration_residual,
    
    lambda_fit = lambda_final_fit
  )
}

###############################################################################
# A4 ET WRAPPER
###############################################################################

estimate_ET_missingcov_A4_compare <- function(
    th,
    data_full,
    K = 3,
    seed = 1234,
    b_method = c(
      "plugin",
      "exact",
      "direct"
    ),
    max.iter = 100,
    eps = 1e-6,
    lambda_maxit = 1000,
    lambda_tol = 1e-10,
    damping = 1) {
  
  b_method <- match.arg(
    b_method
  )
  
  estimate_theta_missingcov_GEC_A4_compare(
    th = th,
    data_full = data_full,
    K = K,
    seed = seed,
    entropy = "ET",
    b_method = b_method,
    max.iter = max.iter,
    eps = eps,
    lambda_maxit = lambda_maxit,
    lambda_tol = lambda_tol,
    damping = damping
  )
}
###############################################################################
# 11. CE dual solver
###############################################################################

###############################################################################
# CE dual solver -- stable version
###############################################################################

solve_lambda_CE_dual <- function(
    b_mat,
    pi_hat,
    D,
    lambda_start = NULL,
    margin = 1e-8,
    maxit = 1000,
    reltol = 1e-10) {
  
  b_mat <- as.matrix(b_mat)
  storage.mode(b_mat) <- "double"
  
  pi_hat <- as.numeric(pi_hat)
  D      <- as.numeric(D)
  
  pi_hat <- pmin(
    pmax(pi_hat, 1e-8),
    1 - 1e-8
  )
  
  ###########################################################################
  # CE calibration covariate
  ###########################################################################
  
  g_pi <- log(1 - pi_hat)
  
  S <- cbind(
    b_mat,
    g_pi
  )
  
  storage.mode(S) <- "double"
  
  q <- ncol(S)
  
  I1 <- which(
    D == 1
  )
  
  ###########################################################################
  # Natural feasible start
  #
  # eta = log(1-pi) < 0
  ###########################################################################
  
  lambda_natural <- c(
    rep(0, q - 1),
    1
  )
  
  if (is.null(lambda_start) ||
      length(lambda_start) != q ||
      any(!is.finite(lambda_start))) {
    
    lambda_start <- lambda_natural
    
  } else {
    
    lambda_start <- as.numeric(lambda_start)
  }
  
  
  ###########################################################################
  # Check feasibility
  ###########################################################################
  
  eta_start <- as.numeric(
    S[I1, , drop = FALSE] %*%
      lambda_start
  )
  
  if (any(!is.finite(eta_start)) ||
      any(eta_start >= -margin)) {
    
    lambda_start <- lambda_natural
    
    eta_start <- as.numeric(
      S[I1, , drop = FALSE] %*%
        lambda_start
    )
  }
  
  if (any(!is.finite(eta_start)) ||
      any(eta_start >= -margin)) {
    
    stop(
      "Could not construct a feasible CE starting value."
    )
  }
  
  
  ###########################################################################
  # CE objective
  #
  # F(eta) = eta - log(1-exp(eta))
  ###########################################################################
  
  objective <- function(lambda) {
    
    lambda <- as.numeric(lambda)
    
    eta <- as.numeric(
      S %*% lambda
    )
    
    eta1 <- eta[I1]
    
    if (any(!is.finite(eta1)) ||
        any(eta1 >= -margin)) {
      
      return(1e100)
    }
    
    F_eta1 <-
      eta1 -
      log1p(
        -exp(eta1)
      )
    
    value <-
      (
        sum(F_eta1) -
          sum(eta)
      ) /
      nrow(S)
    
    if (!is.finite(value)) {
      return(1e100)
    }
    
    as.numeric(value)
  }
  
  
  ###########################################################################
  # CE calibration gradient
  ###########################################################################
  
  gradient <- function(lambda) {
    
    lambda <- as.numeric(lambda)
    
    eta <- as.numeric(
      S %*% lambda
    )
    
    eta1 <- eta[I1]
    
    if (any(!is.finite(eta1)) ||
        any(eta1 >= -margin)) {
      
      return(
        rep(1e20, q)
      )
    }
    
    w1 <- 1 / (
      1 - exp(eta1)
    )
    
    multiplier <- rep(
      -1,
      nrow(S)
    )
    
    multiplier[I1] <-
      w1 - 1
    
    score <- colMeans(
      S * multiplier
    )
    
    as.numeric(score)
  }
  
  
  ###########################################################################
  # Linear CE-domain constraints
  #
  # S_i lambda <= -margin
  #
  # constrOptim requires
  # ui %*% lambda - ci >= 0
  ###########################################################################
  
  ui <- -S[
    I1,
    ,
    drop = FALSE
  ]
  
  storage.mode(ui) <- "double"
  
  ci_vec <- rep(
    as.numeric(margin),
    length(I1)
  )
  
  
  ###########################################################################
  # Optimization
  ###########################################################################
  
  fit <- constrOptim(
    theta = as.numeric(lambda_start),
    f = objective,
    grad = gradient,
    ui = ui,
    ci = ci_vec,
    method = "BFGS",
    control = list(
      maxit = maxit,
      reltol = reltol
    )
  )
  
  
  ###########################################################################
  # Final quantities
  ###########################################################################
  
  lambda_hat <- as.numeric(
    fit$par
  )
  
  eta <- as.numeric(
    S %*%
      lambda_hat
  )
  
  if (any(
    eta[I1] >= 0
  )) {
    
    stop(
      "Final CE solution violates eta < 0."
    )
  }
  
  w_all <- rep(
    NA_real_,
    nrow(S)
  )
  
  w_all[I1] <-
    1 / (
      1 -
        exp(
          eta[I1]
        )
    )
  
  score <- gradient(
    lambda_hat
  )
  
  
  list(
    lambda =
      lambda_hat,
    
    weights =
      as.numeric(
        w_all[I1]
      ),
    
    weights_all =
      as.numeric(
        w_all
      ),
    
    eta =
      eta,
    
    score =
      score,
    
    max_score =
      max(
        abs(score)
      ),
    
    converged =
      isTRUE(
        fit$convergence == 0
      ),
    
    convergence_code =
      fit$convergence,
    
    S =
      S,
    
    g_pi =
      g_pi
  )
}

###############################################################################
# 12. Single clean GEC point-estimation engine
#
# Public wrappers below preserve ET / HD / CE names.
###############################################################################

estimate_theta_missingcov_GEC <- function(
    th,
    data_full,
    K,
    seed,
    entropy = c(
      "ET",
      "HD",
      "CE"
    ),
    max.iter = 100,
    eps = 1e-6,
    lambda_maxit = 1000,
    lambda_tol = 1e-10,
    damping = 1) {

  entropy <- match.arg(
    entropy
  )

  data_all <- prepare_missingcov_data(
    data_full = data_full,
    K = K,
    seed = seed
  )

  D <- as.numeric(
    data_all$D
  )

  I1 <- which(
    D == 1
  )

  y <- data_all$y

  pi_hat <- data_all$pi.hat

  Xhat <- cbind(
    1,
    x = data_all$x,
    z = data_all$z.hat
  )

  Xtrue <- cbind(
    1,
    x = data_all$x,
    z = data_all$z
  )

  Q <- Xtrue[
    I1,
    ,
    drop = FALSE
  ]

  Y1 <- y[I1]

  theta <- as.numeric(
    th
  )

  lambda_current <- NULL

  solve_current_lambda <- function(
      b_mat,
      lambda_start) {

    if (entropy == "ET") {

      return(
        solve_lambda_ET_dual(
          b_mat = b_mat,
          pi_hat = pi_hat,
          D = D,
          lambda_start = lambda_start,
          maxit = lambda_maxit,
          reltol = lambda_tol
        )
      )
    }

    if (entropy == "HD") {

      return(
        solve_lambda_HD_dual(
          b_mat = b_mat,
          pi_hat = pi_hat,
          D = D,
          lambda_start = lambda_start,
          maxit = lambda_maxit,
          reltol = lambda_tol
        )
      )
    }

    solve_lambda_CE_dual(
      b_mat = b_mat,
      pi_hat = pi_hat,
      D = D,
      lambda_start = lambda_start,
      maxit = lambda_maxit,
      reltol = lambda_tol
    )
  }

  for (iter in seq_len(
    max.iter
  )) {

    b_mat <- missingcov_b_plugin(
      beta = theta,
      data_all = data_all
    )

    lambda_fit <-
      solve_current_lambda(
        b_mat = b_mat,
        lambda_start =
          lambda_current
      )

    lambda_current <-
      lambda_fit$lambda

    W <- lambda_fit$weights

    theta_raw <- as.numeric(
      solve(
        crossprod(
          Q,
          W * Q
        ),
        crossprod(
          Q,
          W * Y1
        )
      )
    )

    diff_theta <- max(
      abs(
        theta_raw -
          theta
      )
    )

    theta_new <-
      as.numeric(
        (
          1 - damping
        ) *
          theta +
          damping *
          theta_raw
      )

    theta <- theta_new

    if (diff_theta < eps) {

      break
    }
  }

  converged <-
    diff_theta < eps

  # Final synchronization:
  # recompute b and lambda at the returned theta.
  b_final <- missingcov_b_plugin(
    beta = theta,
    data_all = data_all
  )

  lambda_final_fit <-
    solve_current_lambda(
      b_mat = b_final,
      lambda_start =
        lambda_current
    )

  W_final <-
    lambda_final_fit$weights

  # Check the weighted beta equation at the synchronized pair.
  beta_score <- colMeans(
    missingcov_U(
      beta = theta,
      data_all = data_all
    ) *
      as.numeric(
        D *
          lambda_final_fit$weights_all
      ),
    na.rm = TRUE
  )

  list(
    theta = matrix(
      theta,
      ncol = 1
    ),
    w = W_final,
    weights_all =
      lambda_final_fit$weights_all,
    lambda =
      lambda_final_fit$lambda,
    data_all =
      data_all,
    b_mat =
      b_final,
    g_pi =
      lambda_final_fit$g_pi,
    eta =
      lambda_final_fit$eta,
    lambda_score =
      lambda_final_fit$score,
    max_lambda_score =
      lambda_final_fit$max_score,
    max_beta_score =
      max(
        abs(beta_score)
      ),
    lambda_converged =
      lambda_final_fit$converged,
    converged =
      converged,
    iterations =
      iter,
    entropy =
      entropy
  )
}


###############################################################################
# 13. Public ET / HD / CE wrappers
###############################################################################

estimate_theta_EM_kfold_dual_ET <- function(
    th,
    data_full,
    K,
    seed = 1234,
    max.iter = 100,
    eps = 1e-6,
    lambda_maxit = 1000,
    lambda_tol = 1e-10,
    damping = 1) {

  estimate_theta_missingcov_GEC(
    th = th,
    data_full = data_full,
    K = K,
    seed = seed,
    entropy = "ET",
    max.iter = max.iter,
    eps = eps,
    lambda_maxit = lambda_maxit,
    lambda_tol = lambda_tol,
    damping = damping
  )
}


estimate_theta_EM_kfold_dual_HD <- function(
    th,
    data_full,
    K,
    seed = 1234,
    max.iter = 100,
    eps = 1e-6,
    lambda_maxit = 1000,
    lambda_tol = 1e-10,
    damping = 1) {

  estimate_theta_missingcov_GEC(
    th = th,
    data_full = data_full,
    K = K,
    seed = seed,
    entropy = "HD",
    max.iter = max.iter,
    eps = eps,
    lambda_maxit = lambda_maxit,
    lambda_tol = lambda_tol,
    damping = damping
  )
}


###############################################################################
# CE entropy helper
###############################################################################

gec_entropy_CE <- function() {
  
  list(
    
    ginv = function(eta) {
      1 / (1 - base::exp(eta))
    },
    
    debias = function(pi) {
      base::log(1 - pi)
    },
    
    valid_eta = function(eta) {
      eta < 0
    }
  )
}


###############################################################################
# CE dual solver
###############################################################################

solve_lambda_CE_dual <- function(
    b_mat,
    pi_hat,
    D,
    lambda_start = NULL,
    margin = 1e-8,
    maxit = 1000,
    reltol = 1e-10) {
  
  b_mat <- as.matrix(b_mat)
  storage.mode(b_mat) <- "double"
  
  pi_hat <- as.numeric(pi_hat)
  D      <- as.numeric(D)
  
  pi_hat <- pmin(
    pmax(pi_hat, 1e-8),
    1 - 1e-8
  )
  
  g_pi <- base::log(
    1 - pi_hat
  )
  
  S <- cbind(
    b_mat,
    g_pi
  )
  
  storage.mode(S) <- "double"
  
  q <- ncol(S)
  
  I1 <- which(
    D == 1
  )
  
  if (length(I1) == 0) {
    stop("No complete cases available for CE.")
  }
  
  ###########################################################################
  # Natural feasible starting value
  ###########################################################################
  
  lambda_natural <- c(
    rep(0, q - 1),
    1
  )
  
  if (
    is.null(lambda_start) ||
    length(lambda_start) != q ||
    any(!is.finite(lambda_start))
  ) {
    
    lambda_start <- lambda_natural
    
  } else {
    
    lambda_start <- as.numeric(lambda_start)
  }
  
  
  ###########################################################################
  # Check feasibility
  ###########################################################################
  
  eta_start <- as.numeric(
    S[I1, , drop = FALSE] %*%
      lambda_start
  )
  
  if (
    any(!is.finite(eta_start)) ||
    any(eta_start >= -margin)
  ) {
    
    lambda_start <- lambda_natural
    
    eta_start <- as.numeric(
      S[I1, , drop = FALSE] %*%
        lambda_start
    )
  }
  
  if (
    any(!is.finite(eta_start)) ||
    any(eta_start >= -margin)
  ) {
    
    stop(
      "Could not construct feasible CE starting value."
    )
  }
  
  
  ###########################################################################
  # CE objective
  #
  # F(eta) = eta - log(1-exp(eta))
  ###########################################################################
  
  objective <- function(lambda) {
    
    lambda <- as.numeric(lambda)
    
    eta <- as.numeric(
      S %*%
        lambda
    )
    
    eta1 <- as.numeric(
      eta[I1]
    )
    
    if (
      any(!is.finite(eta1)) ||
      any(eta1 >= -margin)
    ) {
      
      return(1e100)
    }
    
    exp_eta1 <- base::exp(
      eta1
    )
    
    if (
      any(!is.finite(exp_eta1)) ||
      any(exp_eta1 >= 1)
    ) {
      
      return(1e100)
    }
    
    F_eta1 <-
      eta1 -
      base::log1p(
        -exp_eta1
      )
    
    value <-
      (
        sum(F_eta1) -
          sum(eta)
      ) /
      nrow(S)
    
    if (!is.finite(value)) {
      return(1e100)
    }
    
    as.numeric(value)
  }
  
  
  ###########################################################################
  # CE gradient
  #
  # mean[ s_i {D_i w_i - 1} ]
  ###########################################################################
  
  gradient <- function(lambda) {
    
    lambda <- as.numeric(lambda)
    
    eta <- as.numeric(
      S %*%
        lambda
    )
    
    eta1 <- as.numeric(
      eta[I1]
    )
    
    if (
      any(!is.finite(eta1)) ||
      any(eta1 >= -margin)
    ) {
      
      return(
        rep(
          1e20,
          q
        )
      )
    }
    
    exp_eta1 <- base::exp(
      eta1
    )
    
    if (
      any(!is.finite(exp_eta1)) ||
      any(exp_eta1 >= 1)
    ) {
      
      return(
        rep(
          1e20,
          q
        )
      )
    }
    
    w1 <- 1 / (
      1 - exp_eta1
    )
    
    multiplier <- rep(
      -1,
      nrow(S)
    )
    
    multiplier[I1] <-
      w1 - 1
    
    score <- colMeans(
      S * multiplier
    )
    
    as.numeric(score)
  }
  
  
  ###########################################################################
  # CE domain:
  #
  # S_i lambda <= -margin for complete cases
  #
  # constrOptim uses:
  # ui %*% theta - ci >= 0
  ###########################################################################
  
  ui <- -S[
    I1,
    ,
    drop = FALSE
  ]
  
  storage.mode(ui) <- "double"
  
  ci_vec <- rep(
    margin,
    length(I1)
  )
  
  ci_vec <- as.numeric(
    ci_vec
  )
  
  
  ###########################################################################
  # Optimize
  ###########################################################################
  
  fit <- constrOptim(
    theta = as.numeric(
      lambda_start
    ),
    f = objective,
    grad = gradient,
    ui = ui,
    ci = ci_vec,
    method = "BFGS",
    control = list(
      maxit = maxit,
      reltol = reltol
    )
  )
  
  
  ###########################################################################
  # Final solution
  ###########################################################################
  
  lambda_hat <- as.numeric(
    fit$par
  )
  
  eta <- as.numeric(
    S %*%
      lambda_hat
  )
  
  eta1 <- as.numeric(
    eta[I1]
  )
  
  if (
    any(!is.finite(eta1)) ||
    any(eta1 >= 0)
  ) {
    
    stop(
      "Final CE solution violates eta < 0."
    )
  }
  
  exp_eta1 <- base::exp(
    eta1
  )
  
  w1 <- 1 / (
    1 - exp_eta1
  )
  
  w_all <- rep(
    NA_real_,
    nrow(S)
  )
  
  w_all[I1] <-
    as.numeric(
      w1
    )
  
  score <- gradient(
    lambda_hat
  )
  
  list(
    
    lambda =
      lambda_hat,
    
    weights =
      as.numeric(
        w_all[I1]
      ),
    
    weights_all =
      as.numeric(
        w_all
      ),
    
    eta =
      eta,
    
    score =
      score,
    
    max_score =
      max(
        abs(score)
      ),
    
    converged =
      fit$convergence == 0,
    
    convergence_code =
      fit$convergence,
    
    objective =
      fit$value,
    
    counts =
      fit$counts,
    
    S =
      S,
    
    g_pi =
      g_pi
  )
}


###############################################################################
# CE missing-covariate estimator
###############################################################################

estimate_theta_EM_kfold_dual_CE <- function(
    th,
    data_full,
    K,
    seed = 1234,
    max.iter = 100,
    eps = 1e-6,
    lambda_maxit = 1000,
    lambda_tol = 1e-8,
    damping = 0.1) {
  
  data_all <- prepare_missingcov_data(
    data_full = data_full,
    K = K,
    seed = seed
  )
  
  D <- as.numeric(
    data_all$D
  )
  
  I1 <- which(
    D == 1
  )
  
  y <- as.numeric(
    data_all$y
  )
  
  pi_hat <- as.numeric(
    data_all$pi.hat
  )
  
  Xhat <- cbind(
    1,
    x = data_all$x,
    z = data_all$z.hat
  )
  
  Xtrue <- cbind(
    1,
    x = data_all$x,
    z = data_all$z
  )
  
  storage.mode(Xhat)  <- "double"
  storage.mode(Xtrue) <- "double"
  
  Q <- Xtrue[
    I1,
    ,
    drop = FALSE
  ]
  
  Y1 <- y[I1]
  
  theta <- as.numeric(
    th
  )
  
  lambda_current <- NULL
  
  converged <- FALSE
  
  
  for (iter in seq_len(
    max.iter
  )) {
    
    ###########################################################################
    # b(theta)
    ###########################################################################
    
    b_mat <- missingcov_b_plugin(
      beta = theta,
      data_all = data_all
    )
    
    b_mat <- as.matrix(
      b_mat
    )
    
    storage.mode(b_mat) <- "double"
    
    
    ###########################################################################
    # CE lambda
    ###########################################################################
    
    lambda_fit <- solve_lambda_CE_dual(
      b_mat = b_mat,
      pi_hat = pi_hat,
      D = D,
      lambda_start = lambda_current,
      maxit = lambda_maxit,
      reltol = lambda_tol
    )
    
    lambda_current <-
      as.numeric(
        lambda_fit$lambda
      )
    
    W <- as.numeric(
      lambda_fit$weights
    )
    
    
    ###########################################################################
    # Weighted theta update
    ###########################################################################
    
    XtWX <- crossprod(
      Q,
      W * Q
    )
    
    XtWY <- crossprod(
      Q,
      W * Y1
    )
    
    theta_raw <- as.numeric(
      solve(
        XtWX,
        XtWY
      )
    )
    
    diff_theta <- max(
      abs(
        theta_raw -
          theta
      )
    )
    
    theta_new <-
      (
        1 - damping
      ) *
      theta +
      damping *
      theta_raw
    
    theta_new <- as.numeric(
      theta_new
    )
    
    if (diff_theta < eps) {
      
      theta <- theta_new
      
      converged <- TRUE
      
      break
    }
    
    theta <- theta_new
  }
  
  
  ###########################################################################
  # Final synchronized lambda at returned theta
  ###########################################################################
  
  b_final <- missingcov_b_plugin(
    beta = theta,
    data_all = data_all
  )
  
  b_final <- as.matrix(
    b_final
  )
  
  storage.mode(b_final) <- "double"
  
  lambda_final_fit <- solve_lambda_CE_dual(
    b_mat = b_final,
    pi_hat = pi_hat,
    D = D,
    lambda_start = lambda_current,
    maxit = lambda_maxit,
    reltol = lambda_tol
  )
  
  W_final <- as.numeric(
    lambda_final_fit$weights
  )
  
  
  ###########################################################################
  # Final beta update using final CE weights
  ###########################################################################
  
  theta_final <- as.numeric(
    solve(
      crossprod(
        Q,
        W_final * Q
      ),
      crossprod(
        Q,
        W_final * Y1
      )
    )
  )
  
  
  ###########################################################################
  # Recompute b and lambda one final time at theta_final
  ###########################################################################
  
  b_return <- missingcov_b_plugin(
    beta = theta_final,
    data_all = data_all
  )
  
  b_return <- as.matrix(
    b_return
  )
  
  storage.mode(b_return) <- "double"
  
  lambda_return_fit <- solve_lambda_CE_dual(
    b_mat = b_return,
    pi_hat = pi_hat,
    D = D,
    lambda_start =
      lambda_final_fit$lambda,
    maxit =
      lambda_maxit,
    reltol =
      lambda_tol
  )
  
  
  ###########################################################################
  # Return
  ###########################################################################
  
  list(
    
    theta =
      matrix(
        theta_final,
        ncol = 1
      ),
    
    w =
      as.numeric(
        lambda_return_fit$weights
      ),
    
    weights_all =
      as.numeric(
        lambda_return_fit$weights_all
      ),
    
    lambda =
      as.numeric(
        lambda_return_fit$lambda
      ),
    
    data_all =
      data_all,
    
    b_mat =
      b_return,
    
    g_pi =
      lambda_return_fit$g_pi,
    
    eta =
      lambda_return_fit$eta,
    
    lambda_score =
      lambda_return_fit$score,
    
    max_lambda_score =
      lambda_return_fit$max_score,
    
    lambda_converged =
      lambda_return_fit$converged,
    
    converged =
      converged,
    
    iterations =
      iter,
    
    entropy =
      "CE"
  )
}

###############################################################################
# 17. Joint estimating equations for GEC
###############################################################################

gec_Psi_missingcov <- function(
    par,
    data_all,
    b_fun = missingcov_b_plugin,
    entropy = c(
      "ET",
      "HD",
      "CE"
    ),
    n_phi,
    n_lambda,
    n_beta) {

  entropy <- match.arg(
    entropy
  )

  ENT <- gec_entropy_missingcov(
    entropy
  )

  id_phi <- seq_len(
    n_phi
  )

  id_lambda <-
    n_phi +
    seq_len(
      n_lambda
    )

  id_beta <-
    n_phi +
    n_lambda +
    seq_len(
      n_beta
    )

  phi <- par[
    id_phi
  ]

  lambda <- par[
    id_lambda
  ]

  beta <- par[
    id_beta
  ]

  O <- cbind(
    1,
    x = data_all$x,
    y = data_all$y
  )

  D <- as.numeric(
    data_all$D
  )

  pi <- plogis(
    as.numeric(
      O %*%
        phi
    )
  )

  pi <- pmin(
    pmax(
      pi,
      1e-8
    ),
    1 - 1e-8
  )

  b <- as.matrix(
    b_fun(
      beta = beta,
      data_all = data_all
    )
  )

  g_pi <- ENT$debias(
    pi
  )

  S <- cbind(
    b,
    g_pi
  )

  if (ncol(S) != n_lambda) {

    stop(
      "n_lambda does not match ncol(S)."
    )
  }

  eta_lambda <- as.numeric(
    S %*%
      lambda
  )

  if (!all(
    ENT$valid_eta(
      eta_lambda[
        D == 1
      ]
    )
  )) {

    stop(
      paste0(
        "Invalid dual domain for ",
        entropy
      )
    )
  }

  omega <- ENT$ginv(
    eta_lambda
  )

  # PS score
  Psi_phi <-
    O *
    as.numeric(
      D - pi
    )

  # Calibration equation
  Psi_lambda <-
    S *
    as.numeric(
      D * omega - 1
    )

  # Target score
  U <- missingcov_U(
    beta = beta,
    data_all = data_all
  )

  Psi_beta <-
    U *
    as.numeric(
      D * omega
    )

  cbind(
    Psi_phi,
    Psi_lambda,
    Psi_beta
  )
}


###############################################################################
# 18. Full joint GEC sandwich variance
###############################################################################

gec_sandwich_missingcov <- function(
    data_all,
    lambda_hat,
    beta_hat,
    b_fun = missingcov_b_plugin,
    entropy = c(
      "ET",
      "HD",
      "CE"
    ),
    parameter_names = c(
      "beta0",
      "beta1",
      "beta2"
    )) {

  entropy <- match.arg(
    entropy
  )

  ps <- fit_ps_missingcov(
    data_all
  )

  phi_hat <- ps$phi

  beta_hat <- as.numeric(
    beta_hat
  )

  lambda_hat <- as.numeric(
    lambda_hat
  )

  n_phi <- length(
    phi_hat
  )

  n_lambda <- length(
    lambda_hat
  )

  n_beta <- length(
    beta_hat
  )

  par_hat <- c(
    phi_hat,
    lambda_hat,
    beta_hat
  )

  N <- nrow(
    data_all
  )

  Psi_hat <- gec_Psi_missingcov(
    par = par_hat,
    data_all = data_all,
    b_fun = b_fun,
    entropy = entropy,
    n_phi = n_phi,
    n_lambda = n_lambda,
    n_beta = n_beta
  )

  B_hat <- crossprod(
    Psi_hat
  ) /
    N

  psi_bar <- function(par) {

    colMeans(
      gec_Psi_missingcov(
        par = par,
        data_all = data_all,
        b_fun = b_fun,
        entropy = entropy,
        n_phi = n_phi,
        n_lambda = n_lambda,
        n_beta = n_beta
      )
    )
  }

  if (entropy %in%
      c(
        "HD",
        "CE"
      )) {

    J_hat <- safe_jacobian(
      func = psi_bar,
      x = par_hat
    )

  } else {

    J_hat <- numDeriv::jacobian(
      func = psi_bar,
      x = par_hat
    )
  }

  A_hat <- -J_hat

  A_inv <- tryCatch(
    solve(
      A_hat
    ),
    error =
      function(e)
        MASS::ginv(
          A_hat
        )
  )

  V_full <-
    (
      A_inv %*%
        B_hat %*%
        t(A_inv)
    ) /
    N

  beta_idx <-
    (
      n_phi +
        n_lambda +
        1
    ):(
      n_phi +
        n_lambda +
        n_beta
    )

  V_beta <- V_full[
    beta_idx,
    beta_idx,
    drop = FALSE
  ]

  variance <- diag(
    V_beta
  )

  se <- sqrt(
    pmax(
      variance,
      0
    )
  )

  lower <-
    beta_hat -
    qnorm(0.975) *
    se

  upper <-
    beta_hat +
    qnorm(0.975) *
    se

  list(
    method = entropy,
    table = data.frame(
      Variable =
        parameter_names,
      Estimate =
        beta_hat,
      Variance =
        variance,
      SE =
        se,
      CI_Lower =
        lower,
      CI_Upper =
        upper,
      CI_Width =
        upper - lower,
      row.names =
        NULL
    ),
    V_beta =
      V_beta,
    V_full =
      V_full,
    A_hat =
      A_hat,
    B_hat =
      B_hat,
    Psi_hat =
      Psi_hat,
    phi_hat =
      phi_hat,
    lambda_hat =
      lambda_hat,
    beta_hat =
      beta_hat,
    condition_A =
      tryCatch(
        kappa(A_hat),
        error =
          function(e)
            NA_real_
      ),
    max_mean_score =
      max(
        abs(
          colMeans(
            Psi_hat
          )
        )
      )
  )
}


###############################################################################
# 19. Joint IPW sandwich variance
###############################################################################

ipw_sandwich_missingcov <- function(
    data_all,
    beta_hat,
    parameter_names = c(
      "beta0",
      "beta1",
      "beta2"
    )) {

  ps <- fit_ps_missingcov(
    data_all
  )

  O <- ps$O

  D <- as.numeric(
    data_all$D
  )

  phi_hat <- ps$phi

  beta_hat <- as.numeric(
    beta_hat
  )

  n_phi <- length(
    phi_hat
  )

  n_beta <- length(
    beta_hat
  )

  par_hat <- c(
    phi_hat,
    beta_hat
  )

  N <- nrow(
    data_all
  )

  Psi_fun <- function(par) {

    phi <- par[
      seq_len(
        n_phi
      )
    ]

    beta <- par[
      n_phi +
        seq_len(
          n_beta
        )
    ]

    pi <- plogis(
      as.numeric(
        O %*%
          phi
      )
    )

    pi <- pmin(
      pmax(
        pi,
        1e-8
      ),
      1 - 1e-8
    )

    Psi_phi <-
      O *
      as.numeric(
        D - pi
      )

    U <- missingcov_U(
      beta = beta,
      data_all = data_all
    )

    Psi_beta <-
      U *
      as.numeric(
        D / pi
      )

    cbind(
      Psi_phi,
      Psi_beta
    )
  }

  Psi_hat <- Psi_fun(
    par_hat
  )

  B_hat <- crossprod(
    Psi_hat
  ) /
    N

  psi_bar <- function(par) {

    colMeans(
      Psi_fun(
        par
      )
    )
  }

  J_hat <- numDeriv::jacobian(
    func = psi_bar,
    x = par_hat
  )

  A_hat <- -J_hat

  A_inv <- tryCatch(
    solve(
      A_hat
    ),
    error =
      function(e)
        MASS::ginv(
          A_hat
        )
  )

  V_full <-
    (
      A_inv %*%
        B_hat %*%
        t(A_inv)
    ) /
    N

  beta_idx <-
    n_phi +
    seq_len(
      n_beta
    )

  V_beta <- V_full[
    beta_idx,
    beta_idx,
    drop = FALSE
  ]

  variance <- diag(
    V_beta
  )

  se <- sqrt(
    pmax(
      variance,
      0
    )
  )

  lower <-
    beta_hat -
    qnorm(0.975) *
    se

  upper <-
    beta_hat +
    qnorm(0.975) *
    se

  list(
    method =
      "IPW",
    table = data.frame(
      Variable =
        parameter_names,
      Estimate =
        beta_hat,
      Variance =
        variance,
      SE =
        se,
      CI_Lower =
        lower,
      CI_Upper =
        upper,
      CI_Width =
        upper - lower,
      row.names =
        NULL
    ),
    V_beta =
      V_beta,
    V_full =
      V_full,
    A_hat =
      A_hat,
    B_hat =
      B_hat,
    Psi_hat =
      Psi_hat,
    phi_hat =
      phi_hat,
    beta_hat =
      beta_hat
  )
}


###############################################################################
# 20. Joint AIPW sandwich variance
###############################################################################

aipw_sandwich_missingcov <- function(
    data_all,
    beta_hat,
    b_fun = missingcov_b_plugin,
    parameter_names = c(
      "beta0",
      "beta1",
      "beta2"
    )) {

  ps <- fit_ps_missingcov(
    data_all
  )

  O <- ps$O

  D <- as.numeric(
    data_all$D
  )

  phi_hat <- ps$phi

  beta_hat <- as.numeric(
    beta_hat
  )

  n_phi <- length(
    phi_hat
  )

  n_beta <- length(
    beta_hat
  )

  par_hat <- c(
    phi_hat,
    beta_hat
  )

  N <- nrow(
    data_all
  )

  Psi_fun <- function(par) {

    phi <- par[
      seq_len(
        n_phi
      )
    ]

    beta <- par[
      n_phi +
        seq_len(
          n_beta
        )
    ]

    pi <- plogis(
      as.numeric(
        O %*%
          phi
      )
    )

    pi <- pmin(
      pmax(
        pi,
        1e-8
      ),
      1 - 1e-8
    )

    Psi_phi <-
      O *
      as.numeric(
        D - pi
      )

    U <- missingcov_U(
      beta = beta,
      data_all = data_all
    )

    b <- as.matrix(
      b_fun(
        beta = beta,
        data_all = data_all
      )
    )

    Psi_beta <- b

    I1 <- which(
      D == 1
    )

    Psi_beta[I1, ] <-
      b[
        I1,
        ,
        drop = FALSE
      ] +
      (
        U[
          I1,
          ,
          drop = FALSE
        ] -
          b[
            I1,
            ,
            drop = FALSE
          ]
      ) /
      pi[I1]

    cbind(
      Psi_phi,
      Psi_beta
    )
  }

  Psi_hat <- Psi_fun(
    par_hat
  )

  B_hat <- crossprod(
    Psi_hat
  ) /
    N

  psi_bar <- function(par) {

    colMeans(
      Psi_fun(
        par
      )
    )
  }

  J_hat <- numDeriv::jacobian(
    func = psi_bar,
    x = par_hat
  )

  A_hat <- -J_hat

  A_inv <- tryCatch(
    solve(
      A_hat
    ),
    error =
      function(e)
        MASS::ginv(
          A_hat
        )
  )

  V_full <-
    (
      A_inv %*%
        B_hat %*%
        t(A_inv)
    ) /
    N

  beta_idx <-
    n_phi +
    seq_len(
      n_beta
    )

  V_beta <- V_full[
    beta_idx,
    beta_idx,
    drop = FALSE
  ]

  variance <- diag(
    V_beta
  )

  se <- sqrt(
    pmax(
      variance,
      0
    )
  )

  lower <-
    beta_hat -
    qnorm(0.975) *
    se

  upper <-
    beta_hat +
    qnorm(0.975) *
    se

  list(
    method =
      "AIPW",
    table = data.frame(
      Variable =
        parameter_names,
      Estimate =
        beta_hat,
      Variance =
        variance,
      SE =
        se,
      CI_Lower =
        lower,
      CI_Upper =
        upper,
      CI_Width =
        upper - lower,
      row.names =
        NULL
    ),
    V_beta =
      V_beta,
    V_full =
      V_full,
    A_hat =
      A_hat,
    B_hat =
      B_hat,
    Psi_hat =
      Psi_hat,
    phi_hat =
      phi_hat,
    beta_hat =
      beta_hat
  )
}


###############################################################################
# 21. Weight diagnostics
###############################################################################

weight_diagnostics <- function(
    w) {

  w <- as.numeric(
    w
  )

  w <- w[
    is.finite(w) &
      w > 0
  ]

  if (length(w) == 0) {

    return(
      list(
        ESS = NA_real_,
        max_weight = NA_real_
      )
    )
  }

  list(
    ESS =
      sum(w)^2 /
      sum(w^2),
    max_weight =
      max(w)
  )
}


###############################################################################
# 22. Convert one fitted method to simulation rows
###############################################################################

process_missingcov_fit <- function(
    fit,
    method,
    truth,
    entropy = NULL,
    parameter_names = c(
      "beta0",
      "beta1",
      "beta2"
    ),
    replication = NA_integer_) {

  beta_hat <- as.numeric(
    fit$theta
  )

  data_all <- fit$data_all

  if (method == "IPW") {

    var_fit <-
      ipw_sandwich_missingcov(
        data_all = data_all,
        beta_hat = beta_hat,
        parameter_names =
          parameter_names
      )

    wd <- weight_diagnostics(
      fit$w
    )

  } else if (method == "AIPW") {

    var_fit <-
      aipw_sandwich_missingcov(
        data_all = data_all,
        beta_hat = beta_hat,
        parameter_names =
          parameter_names
      )

    # AIPW is not fundamentally a weighting estimator,
    # but report IPW weights for the requested diagnostic.
    ps <- fit_ps_missingcov(
      data_all
    )

    wd <- weight_diagnostics(
      1 /
        ps$pi[
          data_all$D == 1
        ]
    )

  } else {

    var_fit <-
      gec_sandwich_missingcov(
        data_all = data_all,
        lambda_hat = fit$lambda,
        beta_hat = beta_hat,
        entropy = entropy,
        parameter_names =
          parameter_names
      )

    wd <- weight_diagnostics(
      fit$w
    )
  }

  se <- as.numeric(
    var_fit$table$SE
  )

  lower <- beta_hat -
    qnorm(0.975) *
    se

  upper <- beta_hat +
    qnorm(0.975) *
    se

  data.frame(
    replication =
      replication,
    method =
      method,
    parameter =
      parameter_names,
    estimate =
      beta_hat,
    truth =
      as.numeric(
        truth
      ),
    error =
      beta_hat -
      as.numeric(
        truth
      ),
    squared_error =
      (
        beta_hat -
          as.numeric(
            truth
          )
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
      wd$ESS,
    max_weight =
      wd$max_weight,
    estimator_success =
      as.integer(
        isTRUE(
          fit$converged
        )
      ),
    lambda_success =
      if (!is.null(
        fit$lambda_converged
      )) {
        as.integer(
          isTRUE(
            fit$lambda_converged
          )
        )
      } else {
        NA_integer_
      },
    calibration_residual =
      if (!is.null(
        fit$max_lambda_score
      )) {
        fit$max_lambda_score
      } else {
        NA_real_
      },
    stringsAsFactors =
      FALSE
  )
}


###############################################################################
# 23. Monte Carlo summary
#
# Main table:
# Bias
# MC_SD
# Mean analytic SE
# analytic SE / MC SD
# RMSE
# 95% coverage
# CI width
# ESS
# max weight
# failure rates
###############################################################################

summarize_missingcov_results <- function(
    results) {

  if (!requireNamespace(
    "dplyr",
    quietly = TRUE
  )) {

    stop(
      "Please install dplyr."
    )
  }

  results |>
    dplyr::group_by(
      method,
      parameter
    ) |>
    dplyr::summarise(
      N =
        sum(
          is.finite(
            estimate
          )
        ),
      Bias =
        mean(
          estimate -
            truth,
          na.rm = TRUE
        ),
      MC_SD =
        sd(
          estimate,
          na.rm = TRUE
        ),
      Mean_Analytic_SE =
        mean(
          analytic_se,
          na.rm = TRUE
        ),
      SE_to_MCSD =
        Mean_Analytic_SE /
        MC_SD,
      RMSE =
        sqrt(
          mean(
            squared_error,
            na.rm = TRUE
          )
        ),
      Coverage =
        mean(
          covered,
          na.rm = TRUE
        ),
      Mean_CI_Width =
        mean(
          ci_width,
          na.rm = TRUE
        ),
      Mean_ESS =
        mean(
          ESS,
          na.rm = TRUE
        ),
      Mean_MaxWeight =
        mean(
          max_weight,
          na.rm = TRUE
        ),
      Estimator_Failure_Rate =
        mean(
          estimator_success == 0,
          na.rm = TRUE
        ),
      Lambda_Failure_Rate =
        if (all(
          is.na(
            lambda_success
          )
        )) {
          NA_real_
        } else {
          mean(
            lambda_success == 0,
            na.rm = TRUE
          )
        },
      Mean_Calibration_Residual =
        if (all(
          is.na(
            calibration_residual
          )
        )) {
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
}




###############################################################################
# 24. Example: one replication
###############################################################################

# dat <- generate_data(
#   n = 1000,
#   OR = 1,
#   PS = 1
# )
#
# data_all <- prepare_missingcov_data(
#   data_full = dat,
#   K = 3,
#   seed = 1234
# )
#
# theta_start <- c(
#   0,
#   0,
#   0
# )
#
# fit_IPW <- estimate_ipw_missingcov(
#   data_all
# )
#
# fit_AIPW <- estimate_aipw_missingcov(
#   data_all
# )
#
# fit_ET <- estimate_theta_EM_kfold_dual_ET(
#   th = theta_start,
#   data_full = dat,
#   K = 3,
#   seed = 1234
# )
#
# fit_HD <- estimate_theta_EM_kfold_dual_HD(
#   th = theta_start,
#   data_full = dat,
#   K = 3,
#   seed = 1234
# )
#
# fit_CE <- estimate_theta_EM_kfold_dual_CE(
#   th = theta_start,
#   data_full = dat,
#   K = 3,
#   seed = 1234
# )
#
# truth <- c(
#   1,
#   1,
#   2
# )
#
# out_IPW <- process_missingcov_fit(
#   fit = fit_IPW,
#   method = "IPW",
#   truth = truth,
#   replication = 1
# )
#
# out_AIPW <- process_missingcov_fit(
#   fit = fit_AIPW,
#   method = "AIPW",
#   truth = truth,
#   replication = 1
# )
#
# out_ET <- process_missingcov_fit(
#   fit = fit_ET,
#   method = "ET",
#   entropy = "ET",
#   truth = truth,
#   replication = 1
# )
#
# out_HD <- process_missingcov_fit(
#   fit = fit_HD,
#   method = "HD",
#   entropy = "HD",
#   truth = truth,
#   replication = 1
# )
#
# out_CE <- process_missingcov_fit(
#   fit = fit_CE,
#   method = "CE",
#   entropy = "CE",
#   truth = truth,
#   replication = 1
# )
#
# results_one <-
#   dplyr::bind_rows(
#     out_IPW,
#     out_AIPW,
#     out_ET,
#     out_HD,
#     out_CE
#   )
#
# print(
#   results_one
# )
###############################################################################
