###############################################################################
# NHANES FUNCTIONS
# Semi-supervised NHANES application
#
# Repository location:
#   Semi-supervised/functions/nhanes_functions.R
#
# Used by:
#   Semi-supervised/application/nhanes/nhanes_semi_supervised_analysis.R
#
# Provides only the functions required to reproduce:
#   Table 6  - compact repeated-split summary
#   Table S8 - single-split NHANES results
#   Table S9 - repeated-split coefficient-level results
#
# Required packages are loaded by the application script.
# The application script also defines these objects before run_one_split() is
# called: encode_data(), theta_full, and parameter_names.
###############################################################################


###############################################################################
# PSSE AND DRESS COMPARATORS
###############################################################################

polynomial<- function(data,order)
{
  polynomial_Z<- rep(1,nrow(data))
  for (i in 1:order)
  {
    polynomial_i <- data^i
    polynomial_Z <- cbind(polynomial_Z,polynomial_i)
  }
  return(polynomial_Z)
}
PSSE1 <- function(labelled_data,
                  unlabelled_data,
                  c1 = NULL,
                  type = "linear",
                  tau = 0.5,
                  alpha = 1,
                  gamma = 10,
                  sd = TRUE,
                  Kfolds = 5) {
  if (type != "linear") {
    stop("This NHANES implementation of PSSE1 supports type = 'linear' only.")
  }
  if (is.null(alpha)) {
    stop("For the NHANES application, alpha must be supplied; alpha = 1 is used.")
  }

  n <- nrow(labelled_data)
  N <- nrow(unlabelled_data)
  p <- ncol(labelled_data) - 1

  if (p != ncol(unlabelled_data)) {
    stop("The labelled and unlabelled covariate dimensions do not match.")
  }
  if (is.null(c1)) {
    c1 <- n / (n + N)
  }

  labelled_Z <- polynomial(labelled_data[, -1, drop = FALSE], alpha)
  unlabelled_Z <- polynomial(unlabelled_data, alpha)

  unlabelled_Z_mean <- colMeans(unlabelled_Z)
  labelled_Z_secondmoment <- crossprod(labelled_Z) / n
  weights_loss <- as.vector(
    c1 +
      (1 - c1) *
        t(unlabelled_Z_mean) %*%
        solve(labelled_Z_secondmoment) %*%
        t(labelled_Z)
  )

  X_lab <- cbind(1, labelled_data[, -1, drop = FALSE])
  theta_hat <- as.vector(
    solve(
      crossprod(X_lab, weights_loss * X_lab),
      crossprod(X_lab, weights_loss * labelled_data[, 1])
    )
  )

  out <- list(Hattheta = theta_hat)
  if (!sd) {
    out$alpha <- alpha
    return(out)
  }

  X_unlab <- cbind(1, unlabelled_data)
  X_secondmoment_inverse <- solve(crossprod(X_unlab) / N)

  set.seed(20218080)
  index <- caret::createFolds(seq_len(n), k = Kfolds)
  W1_test <- NULL

  for (k in seq_len(Kfolds)) {
    index_k <- as.vector(index[[k]])

    Yt_train <- labelled_data[-index_k, 1]
    Xt_train <- labelled_data[-index_k, -1, drop = FALSE]
    Zt_train <- labelled_Z[-index_k, , drop = FALSE]

    Yt_test <- labelled_data[index_k, 1]
    Xt_test <- labelled_data[index_k, -1, drop = FALSE]
    Zt_test <- labelled_Z[index_k, , drop = FALSE]

    n_train <- length(Yt_train)
    n_test <- length(Yt_test)

    X_train <- cbind(1, Xt_train)
    X_test <- cbind(1, Xt_test)

    L_firstder_train <-
      as.vector(Yt_train - X_train %*% theta_hat) * X_train

    projection_coef <-
      solve(crossprod(Zt_train) / n_train) %*%
      crossprod(Zt_train, L_firstder_train) / n_train

    L_firstder_test <-
      as.vector(Yt_test - X_test %*% theta_hat) * X_test

    W1_test_k <-
      L_firstder_test +
      (c1 - 1) * Zt_test %*% projection_coef

    W1_test <- rbind(W1_test, W1_test_k)
  }

  L_firstder_total <-
    as.vector(labelled_data[, 1] - X_lab %*% theta_hat) * X_lab

  projection_coef_total <-
    solve(crossprod(labelled_Z) / n) %*%
    crossprod(labelled_Z, L_firstder_total) / n

  W2_total <-
    (1 - c1) * unlabelled_Z %*% projection_coef_total

  W1_covariance <- crossprod(W1_test) / n
  W2_covariance <- crossprod(W2_total) / N

  Vc_hat_semi <-
    W1_covariance + (n / N) * W2_covariance

  Var_matrix_hat <-
    X_secondmoment_inverse %*%
    Vc_hat_semi %*%
    X_secondmoment_inverse

  out$sd.of.hattheta <- sqrt(diag(Var_matrix_hat / n))
  out$alpha <- alpha
  out
}

DRESS<- function(labelled_data,unlabelled_data,L,Kfolds)
{

  n=nrow(labelled_data)
  p=ncol(labelled_data)-1
  N=nrow(unlabelled_data)

  base_labelled<- polynomial(labelled_data[,-1],L)
  base_unlabelled <- polynomial(unlabelled_data,L)
  alpha_first_derivative_unlabelled <- colMeans(base_unlabelled)

  tol=10^(-5)
  max_i = 50

  alpha_initial <- rep(0,L*p+1)
  alpha_distance<-10
  i=1
  error_svd=FALSE
  while(alpha_distance >tol& i<=max_i){
    exponential_phi_labelled <- as.vector(exp(base_labelled%*%alpha_initial))
    exponential_phi_labelled[which(exponential_phi_labelled==Inf)]=10^(-7)

    alpha_first_derivative <- colMeans(exponential_phi_labelled*base_labelled)-
      alpha_first_derivative_unlabelled

    alpha_second_derivative <- t(exponential_phi_labelled*base_labelled)%*%base_labelled/n

    if (max(abs(svd(alpha_second_derivative)$d))>99999|min(abs(svd(alpha_second_derivative)$d))<tol)
    {
      error_svd = TRUE
      alpha_initial<- rep(0,L*p+1)
      break
    }

    alpha_new <- as.vector(alpha_initial-solve(alpha_second_derivative)%*%alpha_first_derivative)

    alpha_distance <- sqrt(sum((alpha_new-alpha_initial)^2))
    alpha_initial<- alpha_new
    i=i+1
  }

  exponential_phi_labelled <- as.vector(exp(base_labelled%*%alpha_initial))

  hattheta_DRESS<- lm(labelled_data[,1]~labelled_data[,-1],weights=exponential_phi_labelled)$coefficients

  hattheta_DRESS = append(list("Hattheta"=hattheta_DRESS),list("error"=error_svd))

  X_secondmoment_inverse=solve(t(cbind(rep(1,N),unlabelled_data))%*%cbind(rep(1,N),unlabelled_data)/N)

  set.seed(20218080)
  index=caret::createFolds(1:n, k = Kfolds) # data splitting
  W1_test = vector()

  for(k in 1:Kfolds){

    index_k=as.vector(index[[k]])
    Yt_train = labelled_data[-index_k,1]
    Xt_train = labelled_data[-index_k,-1]
    Zt_train = base_labelled[-index_k,]
    Yt_test = labelled_data[index_k,1]
    Xt_test = labelled_data[index_k,-1]
    Zt_test = base_labelled[index_k,]

    nrow_Xt_train = n-length(index_k)
    nrow_Xt_test = length(index_k)

    hattheta_supervised_train_cv <-  hattheta_DRESS[[1]]#lm(Yt_train~Xt_train)$coefficients
    L_firstder_train_cv <- as.vector(Yt_train-cbind(rep(1,nrow_Xt_train),Xt_train)%*%hattheta_supervised_train_cv)*
      cbind(rep(1,nrow_Xt_train),Xt_train)
    L_firstder_projection_cof <- solve(t(Zt_train)%*%Zt_train/nrow_Xt_train)%*%
      t(Zt_train)%*%L_firstder_train_cv/nrow_Xt_train

    L_firstder_test_cv <- as.vector(Yt_test-cbind(rep(1,nrow_Xt_test),Xt_test)%*%hattheta_supervised_train_cv)*
      cbind(rep(1,nrow_Xt_test),Xt_test)

    W1_test_k <- L_firstder_test_cv-Zt_test%*%L_firstder_projection_cof
    W1_test<- rbind(W1_test,W1_test_k)

  }

  L_firstder_total <- as.vector(labelled_data[,1]-(cbind(rep(1,n),labelled_data[,-1]))%*%hattheta_DRESS[[1]])*
    (cbind(rep(1,n),labelled_data[,-1]))
  L_firstder_projection_cof_total <- solve(t(base_labelled)%*%base_labelled/n)%*%t(base_labelled)%*%L_firstder_total/n
  W2_total <- base_unlabelled%*%L_firstder_projection_cof_total

  W1_covariance <- t(W1_test)%*%W1_test/n
  W2_covariance <- t(W2_total)%*%W2_total/N
  Vc_hat_semi = W1_covariance+(n/N)*W2_covariance

  Var_matrix_hat = X_secondmoment_inverse%*%Vc_hat_semi%*%X_secondmoment_inverse
  sd_hattheta_DRESS = sqrt(diag(Var_matrix_hat/n))

  hattheta_DRESS <- append(hattheta_DRESS,list("sd.of.hattheta"=sd_hattheta_DRESS))

  return(hattheta_DRESS)
}


###############################################################################
# CROSS-FITTING
###############################################################################

build_SU_folds <- function(df, K , id_col = "ID", delta_col = "D", seed = seed) {
  set.seed(seed)

  id <- df[[id_col]]
  delta <- df[[delta_col]]

  idx_S  <- which(delta == 1)  # labeled indices
  idx_U0 <- which(delta == 0)  # unlabeled-only indices

  nS  <- length(idx_S)
  nU0 <- length(idx_U0)
  if (nS < K || nU0 < K) stop("Need at least K labeled and K unlabeled observations.")

  split_K <- function(idx, K) {
    if (length(idx) == 0) return(rep(list(integer(0)), K))
    idx <- sample(idx)  # shuffle
    split(idx, rep(1:K, length.out = length(idx)))
  }

  S_parts  <- split_K(idx_S,  K)  # disjoint labeled parts S_k
  U0_parts <- split_K(idx_U0, K)  # disjoint unlabeled-only parts U0_k

  folds <- vector("list", K)
  fold_id_S  <- integer(nrow(df))
  fold_id_U  <- integer(nrow(df))

  for (k in 1:K) {
    S_k_idx  <- S_parts[[k]]
    U0_k_idx <- U0_parts[[k]]
    U_k_idx  <- c(S_k_idx, U0_k_idx)  # ensure S_k ⊆ U_k

    folds[[k]] <- list(
      S_k_ids  = id[S_k_idx],
      U_k_ids  = id[U_k_idx],
      S_k_idx  = S_k_idx,
      U_k_idx  = U_k_idx
    )

    fold_id_S[S_k_idx] <- k
    fold_id_U[U_k_idx] <- k
  }

  list(
    folds = folds,
    fold_id_S = fold_id_S,  # for labeled rows (0 if unlabeled)
    fold_id_U = fold_id_U   # for all rows in U_k (0 if outside that fold)
  )
}

k_fold_function<-function(df,K,seed){
  res <- build_SU_folds(df, K , id_col = "ID", delta_col = "D", seed = seed)

  lab_idx_all <- which(df$D == 1)

  data_unlabeled<-list()
  for(k in 1:K){
    train_idx_lab <- setdiff(lab_idx_all, res$folds[[k]]$S_k_idx)
    validation_data_labeled <- df[res$folds[[k]]$S_k_idx,]
    validation_data_unlabeled<-df[res$folds[[k]]$U_k_idx,]

    train_data_labeled <-df[train_idx_lab,] ###get the data frame for trained data

    gam_model <- mgcv::gam(Y ~ s(age) + s(BMI) + s(SBP) + s(DBP) + sex + race,data = as.data.frame(
        subset(
          train_data_labeled,
          select = -c(D, ID, pi.hat)
        )
      ),
      method = "REML"
    )

    y_hat <- predict(gam_model,newdata = as.data.frame(validation_data_unlabeled))
    validation_data_unlabeled$y.hat<-y_hat

    data_unlabeled[[k]]<-validation_data_unlabeled
  }
  return(data_unlabeled)

}


###############################################################################
# GEC ENTROPY AND STACKED SANDWICH
###############################################################################

gec_entropy <- function(entropy = c("SL", "EL", "ET", "HD", "CE")) {

  entropy <- match.arg(entropy)

  if (entropy == "SL") {

    g <- function(w) w

    ginv <- function(eta) eta

    gp <- function(w) rep(1, length(w))

    debias <- function(pi) 1 / pi

    valid_eta <- function(eta) rep(TRUE, length(eta))
  }

  if (entropy == "EL") {

    g <- function(w) -1 / w

    ginv <- function(eta) -1 / eta

    gp <- function(w) 1 / w^2

    debias <- function(pi) -pi

    valid_eta <- function(eta) eta < 0
  }

  if (entropy == "ET") {

    g <- function(w) log(w)

    ginv <- function(eta) exp(eta)

    gp <- function(w) 1 / w

    debias <- function(pi) log(1 / pi)

    valid_eta <- function(eta) rep(TRUE, length(eta))
  }

  if (entropy == "HD") {

    g <- function(w) -1 / (2 * sqrt(w))

    ginv <- function(eta) 1 / (4 * eta^2)

    gp <- function(w) 1 / (4 * w^(3/2))

    debias <- function(pi) -sqrt(pi) / 2

    valid_eta <- function(eta) eta < 0
  }

  if (entropy == "CE") {

    g <- function(w) log((w - 1) / w)

    ginv <- function(eta) 1 / (1 - exp(eta))

    gp <- function(w) 1 / (w * (w - 1))

    debias <- function(pi) log(1 - pi)

    valid_eta <- function(eta) eta < 0
  }

  list(
    name = entropy,
    g = g,
    ginv = ginv,
    gp = gp,
    debias = debias,
    valid_eta = valid_eta
  )
}

pi_logistic <- function(phi, O) {

  eta <- as.numeric(O %*% phi)

  plogis(eta)
}

h_logistic <- function(phi, O) {

  pi <- pi_logistic(phi, O)

  O * pi
}

gec_Psi <- function(beta,
                    D,
                    O,
                    U_fun,
                    b_fun,
                    n_phi,
                    n_lambda,
                    n_theta,
                    entropy) {

  ENT <- gec_entropy(entropy)

  id_phi <- seq_len(n_phi)

  id_lambda <- n_phi + seq_len(n_lambda)

  id_theta <- n_phi + n_lambda + seq_len(n_theta)

  phi <- beta[id_phi]
  lambda <- beta[id_lambda]
  theta <- beta[id_theta]

  pi <- pi_logistic(phi, O)

  h <- h_logistic(phi, O)

  b <- as.matrix(
    b_fun(theta)
  )

  gpi <- ENT$debias(pi)

  S <- cbind(
    b,
    gpi
  )

  if (ncol(S) != n_lambda) {
    stop("n_lambda does not match number of calibration covariates.")
  }

  eta <- as.numeric(
    S %*% lambda
  )

  if (!all(ENT$valid_eta(eta[D == 1]))) {
    stop(
      paste0(
        "Invalid dual domain for entropy ",
        entropy
      )
    )
  }

  omega <- ENT$ginv(eta)

  U <- as.matrix(
    U_fun(theta)
  )

  Psi_phi <- h *
    as.numeric(D / pi - 1)

  Psi_lambda <- S *
    as.numeric(D * omega - 1)

  U[D == 0, ] <- 0

  Psi_theta <- U *
    as.numeric(D * omega)

  cbind(
    Psi_phi,
    Psi_lambda,
    Psi_theta
  )
}

safe_jacobian <- function(
    func,
    x,
    eps = 1e-6,
    min_eps = 1e-10,
    max_shrink = 30) {

  x <- as.numeric(x)
  f0 <- func(x)

  m <- length(f0)
  p <- length(x)

  J <- matrix(
    NA_real_,
    nrow = m,
    ncol = p
  )

  for (j in seq_len(p)) {

    h <- eps * max(1, abs(x[j]))

    success <- FALSE

    for (k in seq_len(max_shrink)) {

      x_plus  <- x
      x_minus <- x

      x_plus[j]  <- x_plus[j]  + h
      x_minus[j] <- x_minus[j] - h

      f_plus <- tryCatch(
        func(x_plus),
        error = function(e) NULL
      )

      f_minus <- tryCatch(
        func(x_minus),
        error = function(e) NULL
      )

      if (!is.null(f_plus) &&
          !is.null(f_minus) &&
          all(is.finite(f_plus)) &&
          all(is.finite(f_minus))) {

        J[, j] <-
          (f_plus - f_minus) / (2 * h)

        success <- TRUE
        break
      }

      if (!is.null(f_plus) &&
          all(is.finite(f_plus))) {

        J[, j] <-
          (f_plus - f0) / h

        success <- TRUE
        break
      }

      if (!is.null(f_minus) &&
          all(is.finite(f_minus))) {

        J[, j] <-
          (f0 - f_minus) / h

        success <- TRUE
        break
      }

      h <- h / 2

      if (h < min_eps) {
        break
      }
    }

    if (!success) {

      stop(
        paste(
          "Unable to calculate safe Jacobian",
          "for parameter", j
        )
      )
    }
  }

  J
}

gec_sandwich <- function(phi_hat,
                         lambda_hat,
                         theta_hat,
                         D,
                         O,
                         U_fun,
                         b_fun,
                         entropy = c("SL", "EL", "ET", "HD", "CE"),
                         parameter_names = NULL) {

  entropy <- match.arg(entropy)

  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("Please install package 'numDeriv'.")
  }

  phi_hat <- as.numeric(phi_hat)
  lambda_hat <- as.numeric(lambda_hat)
  theta_hat <- as.numeric(theta_hat)

  n_phi <- length(phi_hat)
  n_lambda <- length(lambda_hat)
  n_theta <- length(theta_hat)

  beta_hat <- c(
    phi_hat,
    lambda_hat,
    theta_hat
  )

  N <- length(D)

  Psi_hat <- gec_Psi(
    beta = beta_hat,
    D = D,
    O = O,
    U_fun = U_fun,
    b_fun = b_fun,
    n_phi = n_phi,
    n_lambda = n_lambda,
    n_theta = n_theta,
    entropy = entropy
  )

  B_hat <- crossprod(Psi_hat) / N

  psi_bar <- function(beta) {

    Psi <- gec_Psi(
      beta = beta,
      D = D,
      O = O,
      U_fun = U_fun,
      b_fun = b_fun,
      n_phi = n_phi,
      n_lambda = n_lambda,
      n_theta = n_theta,
      entropy = entropy
    )

    colMeans(Psi)
  }

  if (entropy %in% c("HD", "CE")) {

    J_hat <- safe_jacobian(
      func = psi_bar,
      x = beta_hat,
      eps = 1e-6,
      min_eps = 1e-10,
      max_shrink = 30
    )

  } else {

    J_hat <- numDeriv::jacobian(
      func = psi_bar,
      x = beta_hat
    )
  }

  A_hat <- -J_hat
  theta_idx <- (
    n_phi + n_lambda + 1
  ):(
    n_phi + n_lambda + n_theta
  )

  A_inv <- solve(A_hat)

  V_beta <- (
    A_inv %*%
      B_hat %*%
      t(A_inv)
  ) / N

  A_tt_numeric <-
    A_hat[
      theta_idx,
      theta_idx,
      drop = FALSE
    ]

  V_theta <- V_beta[
    theta_idx,
    theta_idx,
    drop = FALSE
  ]

  se <- sqrt(
    pmax(diag(V_theta), 0)
  )

  lower <- theta_hat -
    qnorm(0.975) * se

  upper <- theta_hat +
    qnorm(0.975) * se

  if (is.null(parameter_names)) {
    parameter_names <- paste0(
      "theta",
      seq_len(n_theta)
    )
  }

  table <- data.frame(
    Variable = parameter_names,
    Estimate = theta_hat,
    Variance = diag(V_theta),
    SE = se,
    CI_Lower = lower,
    CI_Upper = upper,
    CI_Width = upper - lower,
    row.names = NULL
  )

  list(
    entropy = entropy,
    table = table,
    V_theta = V_theta,
    V_beta = V_beta,
    A_hat = A_hat,
    B_hat = B_hat,
    Psi_hat = Psi_hat,
    beta_hat = beta_hat
  )
}


###############################################################################
# ET ESTIMATOR
###############################################################################

solve_lambda_ET_dual <- function(
    b_mat,
    pi_hat,
    D,
    lambda_start = NULL,
    maxit = 1000,
    reltol = 1e-10) {

  g_pi <- log(1 / pi_hat)

  S <- cbind(
    b_mat,
    g_pi
  )

  q <- ncol(S)

  if (is.null(lambda_start)) {

    lambda_start <- c(
      rep(0, q - 1),
      1
    )
  }

  objective <- function(lambda) {

    eta <- as.numeric(
      S %*% lambda
    )

    mean(
      D * exp(eta) - eta
    )
  }

  gradient <- function(lambda) {

    eta <- as.numeric(
      S %*% lambda
    )

    w <- exp(eta)

    colMeans(
      S *
        as.numeric(
          D * w - 1
        )
    )
  }

  fit <- optim(
    par = lambda_start,
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
    S %*% lambda_hat
  )

  w_all <- exp(eta)

  score <- gradient(
    lambda_hat
  )

  list(
    lambda = lambda_hat,

    weights = w_all[D == 1],

    weights_all = w_all,

    eta = eta,

    score = score,

    max_score = max(abs(score)),

    converged = fit$convergence == 0,

    convergence_code = fit$convergence,

    iterations = fit$counts,

    S = S,

    g_pi = g_pi
  )
}

estimate_theta_EM_kfold_dual_ET <- function(
    th,
    data_full,
    K,
    seed=2025,
    max.iter = 500,
    eps = 1e-4,
    lambda_tol = 1e-10,
    lambda_maxit = 1000,
    damping = 1) {

  fold_hat <- k_fold_function(
    df = data_full,
    K = K,
    seed = seed
  )

  data_all <- do.call(
    rbind,
    fold_hat
  )

  x_cols <- subset(
    data_all,
    select = -c(
      Y,
      D,
      pi.hat,
      ID,
      y.hat
    )
  )

  model.formula <- as.formula(
    paste(
      "~",
      paste(
        colnames(x_cols),
        collapse = "+"
      )
    )
  )

  X <- model.matrix(
    model.formula,
    data = x_cols
  )

  y.hat <- data_all$y.hat
  D     <- data_all$D
  I1    <- which(D == 1)

  theta <- as.numeric(th)

  iter <- 0

  lambda_current <- NULL

  repeat {

    iter <- iter + 1

    error_hat <- y.hat -
      as.numeric(
        X %*% theta
      )

    b_mat <- X * error_hat

    lambda_fit <- solve_lambda_ET_dual(
      b_mat = b_mat,
      pi_hat = data_all$pi.hat,
      D = D,
      lambda_start = lambda_current,
      maxit = lambda_maxit,
      reltol = lambda_tol
    )

    lambda_hat <- lambda_fit$lambda

    lambda_current <- lambda_hat

    W <- lambda_fit$weights

    Q <- X[
      I1,
      ,
      drop = FALSE
    ]

    Y1 <- data_all$Y[I1]

    th.raw <- as.numeric(
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
        th.raw - theta
      )
    )

    if (diff_theta < eps) {

      return(
        list(
          theta = matrix(
            th.raw,
            ncol = 1
          ),

          w = W,

          lambda = lambda_hat,

          data_all = data_all,

          b_mat = b_mat,

          g_pi = lambda_fit$g_pi,

          eta = lambda_fit$eta,

          lambda_score = lambda_fit$score,

          max_lambda_score =
            lambda_fit$max_score,

          lambda_converged =
            lambda_fit$converged,

          converged = TRUE,

          iterations = iter
        )
      )
    }

    if (iter >= max.iter) {

      return(
        list(
          theta = matrix(
            th.raw,
            ncol = 1
          ),

          w = W,

          lambda = lambda_hat,

          data_all = data_all,

          b_mat = b_mat,

          g_pi = lambda_fit$g_pi,

          eta = lambda_fit$eta,

          lambda_score = lambda_fit$score,

          max_lambda_score =
            lambda_fit$max_score,

          lambda_converged =
            lambda_fit$converged,

          converged = TRUE,

          iterations = iter
        )
      )
    }

    theta <-
      (1 - damping) * theta +
      damping * th.raw
  }
}


###############################################################################
# HD ESTIMATOR
###############################################################################

solve_lambda_HD_dual <- function(
    b_mat,
    pi_hat,
    D,
    lambda_start = NULL,
    margin = 1e-8,
    maxit = 1000,
    reltol = 1e-10) {

  g_pi <- -sqrt(pi_hat) / 2

  S <- cbind(
    b_mat,
    g_pi
  )

  q <- ncol(S)

  I1 <- which(D == 1)

  if (is.null(lambda_start)) {

    lambda_start <- c(
      rep(0, q - 1),
      1
    )

  } else {

    lambda_start <- as.numeric(
      lambda_start
    )
  }

  eta_start <- as.numeric(
    S[I1, , drop = FALSE] %*%
      lambda_start
  )

  if (any(eta_start >= -margin)) {

    lambda_start <- c(
      rep(0, q - 1),
      1
    )

    eta_start <- as.numeric(
      S[I1, , drop = FALSE] %*%
        lambda_start
    )
  }

  if (any(eta_start >= -margin)) {
    stop(
      "Could not construct a feasible HD lambda starting value."
    )
  }

  objective <- function(lambda) {

    eta <- as.numeric(
      S %*% lambda
    )

    eta1 <- eta[I1]

    if (any(eta1 >= 0)) {
      return(1e100)
    }

    F_eta1 <- -1 / (
      4 * eta1
    )

    (
      sum(F_eta1) -
        sum(eta)
    ) / nrow(S)
  }

  gradient <- function(lambda) {

    eta <- as.numeric(
      S %*% lambda
    )

    eta1 <- eta[I1]

    if (any(eta1 >= 0)) {
      return(
        rep(1e20, q)
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
      S *
        multiplier
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
    theta = lambda_start,
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
    S %*% lambda_hat
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

  return(
    list(

      lambda = lambda_hat,

      weights = w_all[I1],

      weights_all = w_all,

      eta = eta,

      score = score,

      max_score =
        max(abs(score)),

      converged =
        fit$convergence == 0,

      convergence_code =
        fit$convergence,

      objective =
        fit$value,

      counts =
        fit$counts,

      S = S,

      g_pi = g_pi
    )
  )
}

estimate_theta_EM_kfold_dual_HD <- function(
    th,
    data_full,
    K,
    seed=2025,
    max.iter = 500,
    eps = 1e-4,
    lambda_tol = 1e-10,
    lambda_maxit = 1000,
    damping ) {

  fold_hat <- k_fold_function(
    df = data_full,
    K = K,
    seed = seed
  )

  data_all <- do.call(
    rbind,
    fold_hat
  )

  x_cols <- subset(
    data_all,
    select = -c(
      Y,
      D,
      pi.hat,
      ID,
      y.hat
    )
  )

  model.formula <- as.formula(
    paste(
      "~",
      paste(
        colnames(x_cols),
        collapse = "+"
      )
    )
  )

  X <- model.matrix(
    model.formula,
    data = as.data.frame(
      data_all[
        ,
        colnames(x_cols),
        drop = FALSE
      ]
    )
  )

  y.hat <- data_all$y.hat

  D <- data_all$D

  I1 <- which(
    D == 1
  )

  theta <- as.numeric(
    th
  )

  iter <- 0

  lambda_current <- NULL

  repeat {

    iter <- iter + 1

    error_hat <- y.hat -
      as.numeric(
        X %*% theta
      )

    b_mat <- X *
      error_hat

    lambda_fit <- solve_lambda_HD_dual(
      b_mat = b_mat,
      pi_hat = data_all$pi.hat,
      D = D,
      lambda_start = lambda_current,
      maxit = lambda_maxit,
      reltol = lambda_tol
    )

    if (!lambda_fit$converged) {

      warning(
        paste(
          "HD lambda optimization did not fully converge",
          "at outer iteration",
          iter,
          ". convergence code =",
          lambda_fit$convergence_code
        )
      )
    }

    lambda_hat <-
      lambda_fit$lambda

    lambda_current <-
      lambda_hat

    W <-
      lambda_fit$weights

    Q <- X[
      I1,
      ,
      drop = FALSE
    ]

    Y1 <- data_all$Y[I1]

    th.raw <- as.numeric(
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
        th.raw -
          theta
      )
    )

    ESS <- sum(W)^2 /
      sum(W^2)

    ESS_ratio <-
      ESS / length(W)

    if (diff_theta < eps) {

      theta_final <-
        th.raw

      error_final <- y.hat -
        as.numeric(
          X %*% theta_final
        )

      b_final <- X *
        error_final

      lambda_final_fit <-
        solve_lambda_HD_dual(
          b_mat = b_final,
          pi_hat = data_all$pi.hat,
          D = D,
          lambda_start = lambda_hat,
          maxit = lambda_maxit,
          reltol = lambda_tol
        )

      lambda_final <-
        lambda_final_fit$lambda

      W_final <-
        lambda_final_fit$weights

      theta_final2 <- as.numeric(
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

      y <- data_all$Y[I1]

      return(
        list(

          theta = matrix(
            theta_final2,
            ncol = 1
          ),

          w =
            W_final,

          lambda =
            lambda_final,

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

          lambda_converged =
            lambda_final_fit$converged,

          converged =
            TRUE,

          iterations =
            iter
        )
      )
    }

    if (iter >= max.iter) {

      y <- data_all$Y[I1]

      return(
        list(

          theta =
            matrix(
              th.raw,
              ncol = 1
            ),

          w =
            W,

          lambda =
            lambda_hat,

          data_all =
            data_all,

          b_mat =
            b_mat,

          g_pi =
            lambda_fit$g_pi,

          eta =
            lambda_fit$eta,

          lambda_score =
            lambda_fit$score,

          max_lambda_score =
            lambda_fit$max_score,

          lambda_converged =
            lambda_fit$converged,

          converged =
            FALSE,

          iterations =
            iter
        )
      )
    }

    theta <- as.numeric(
      (1 - damping) * theta +
        damping * th.raw
    )
  }
}


###############################################################################
# CE ESTIMATOR
###############################################################################

solve_lambda_CE_dual <- function(
    b_mat,
    pi_hat,
    D,
    lambda_start = NULL,
    margin = 1e-8,
    maxit = 1000,
    reltol = 1e-10) {

  g_pi <- log(1 - pi_hat)
  S <- cbind(
    b_mat,
    g_pi
  )

  q <- ncol(S)

  I1 <- which(D == 1)

  lambda_natural <- c(
    rep(0, q - 1),
    1
  )

  if (is.null(lambda_start)) {

    lambda_start <- lambda_natural

  } else {

    lambda_start <- as.numeric(lambda_start)
  }

  eta_start <- as.numeric(
    S[I1, , drop = FALSE] %*%
      lambda_start
  )

  if (any(eta_start >= -margin)) {

    lambda_start <- lambda_natural

    eta_start <- as.numeric(
      S[I1, , drop = FALSE] %*%
        lambda_start
    )
  }

  if (any(eta_start >= -margin)) {

    stop(
      "Could not construct a feasible CE lambda starting value."
    )
  }

  objective <- function(lambda) {

    eta <- as.numeric(
      S %*% lambda
    )

    eta1 <- eta[I1]

    if (any(eta1 >= 0)) {
      return(1e100)
    }

    F_eta1 <-
      eta1 -
      log1p(
        -exp(eta1)
      )

    mean(
      D * 0 # placeholder only to preserve N scaling
    ) +
      (
        sum(F_eta1) -
          sum(eta)
      ) / nrow(S)
  }

  gradient <- function(lambda) {

    eta <- as.numeric(
      S %*% lambda
    )

    eta1 <- eta[I1]

    if (any(eta1 >= 0)) {
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

  eta_start <- as.numeric(
    S[I1, , drop = FALSE] %*%
      lambda_start
  )

  if (any(eta_start >= -margin)) {

    stop(
      "Initial lambda is not strictly inside the CE dual domain."
    )
  }

  fit <- constrOptim(
    theta = lambda_start,
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
    S %*% lambda_hat
  )

  w_all <- rep(
    NA_real_,
    nrow(S)
  )

  w_all[I1] <- 1 / (
    1 - exp(
      eta[I1]
    )
  )

  score <- gradient(
    lambda_hat
  )

  list(

    lambda = lambda_hat,

    weights = w_all[I1],

    eta = eta,

    score = score,

    max_score =
      max(abs(score)),

    converged =
      fit$convergence == 0,

    convergence_code =
      fit$convergence,

    objective =
      fit$value,

    counts =
      fit$counts,

    S = S,

    g_pi = g_pi
  )
}

estimate_theta_EM_kfold_CE <- function(
    th,
    data_full,
    K,
    seed,
    max.iter = 500,
    eps = 1e-4,
    lambda_tol = 1e-8,
    lambda_maxit = 1000,
    damping =0.1) {

  fold_hat <- k_fold_function(
    df = data_full,
    K = K,
    seed = seed
  )

  data_all <- do.call(
    rbind,
    fold_hat
  )

  x_cols <- subset(
    data_all,
    select = -c(
      Y,
      D,
      pi.hat,
      ID,
      y.hat
    )
  )

  model.formula <- as.formula(
    paste(
      "~",
      paste(
        colnames(x_cols),
        collapse = "+"
      )
    )
  )

  X <- model.matrix(
    model.formula,
    data = as.data.frame(
      data_all[
        ,
        colnames(x_cols),
        drop = FALSE
      ]
    )
  )

  y.hat <- data_all$y.hat
  D     <- data_all$D
  I1    <- which(D == 1)

  theta <- as.numeric(th)

  iter <- 0

  lambda_current <- NULL

  repeat {

    iter <- iter + 1

    error_hat <- y.hat -
      as.numeric(
        X %*% theta
      )

    b_mat <- X * error_hat

    lambda_fit <- solve_lambda_CE_dual(
      b_mat = b_mat,
      pi_hat = data_all$pi.hat,
      D = D,
      lambda_start = lambda_current,
      maxit = lambda_maxit,
      reltol = lambda_tol
    )

    if (!lambda_fit$converged) {

      warning(
        paste(
          "Lambda optimization did not fully converge",
          "at outer iteration",
          iter
        )
      )
    }

    lambda_hat <- lambda_fit$lambda

    lambda_current <- lambda_hat

    W <- lambda_fit$weights

    Q <- X[
      I1,
      ,
      drop = FALSE
    ]

    Y1 <- data_all$Y[I1]

    th.raw <- as.numeric(
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
        th.raw - theta
      )
    )

    if (diff_theta < eps) {

      theta_final <- th.raw

      error_final <- y.hat -
        as.numeric(
          X %*% theta_final
        )

      b_final <- X * error_final

      lambda_final_fit <- solve_lambda_CE_dual(
        b_mat = b_final,
        pi_hat = data_all$pi.hat,
        D = D,
        lambda_start = lambda_hat,
        maxit = lambda_maxit,
        reltol = lambda_tol
      )

      lambda_final <- lambda_final_fit$lambda
      W_final      <- lambda_final_fit$weights
      g_pi_final   <- lambda_final_fit$g_pi

      theta_final2 <- as.numeric(
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

      y <- data_all$Y[I1]

      return(
        list(

          theta = matrix(
            theta_final2,
            ncol = 1
          ),

          w = W_final,

          lambda = lambda_final,

          data_all = data_all,

          b_mat = b_final,

          g_pi = g_pi_final,

          eta = lambda_final_fit$eta,

          lambda_score =
            lambda_final_fit$score,

          max_lambda_score =
            lambda_final_fit$max_score,

          lambda_converged =
            lambda_final_fit$converged,

          converged = TRUE,

          iterations = iter
        )
      )
    }

    if (iter >= max.iter) {

      y <- data_all$Y[I1]

      return(
        list(

          theta = matrix(
            th.raw,
            ncol = 1
          ),

          w = W,

          lambda = lambda_hat,

          data_all = data_all,

          b_mat = b_mat,

          g_pi = lambda_fit$g_pi,

          eta = lambda_fit$eta,

          lambda_score =
            lambda_fit$score,

          max_lambda_score =
            lambda_fit$max_score,

          lambda_converged =
            lambda_fit$converged,

          converged = TRUE,

          iterations = iter
        )
      )
    }

    theta <- as.numeric(
      (1 - damping) * theta +
        damping * th.raw
    )
  }
}


###############################################################################
# NHANES SPLIT RUNNER
###############################################################################

run_one_split <- function(split_seed,datX) {

  tryCatch({

    set.seed(split_seed)

    original_labeled <- datX %>%
      dplyr::filter(!is.na(Y))

    original_unlabeled <- datX %>%
      dplyr::filter(is.na(Y))

    labeled <- original_labeled %>%
      dplyr::slice_sample(prop = 0.5)

    hidden_labeled <- original_labeled %>%
      dplyr::filter(!SEQN %in% labeled$SEQN) %>%
      dplyr::mutate(Y = NA_real_)

    unlabeled <- dplyr::bind_rows(
      original_unlabeled,
      hidden_labeled
    )

    labeled_final <- cbind(
      Y = labeled$Y,
      encode_data(labeled)
    )

    unlabeled_final <- encode_data(
      unlabeled
    )

    model.formula <- as.formula(

      paste(
        "Y ~",
        paste(
          sprintf(
            "`%s`",
            colnames(labeled_final)[-1]
          ),
          collapse = " + "
        )
      )
    )

    model.fit <- lm(
      model.formula,
      data = as.data.frame(labeled_final)
    )

    theta_start <- as.numeric(
      coef(model.fit)
    )

    S1 <- summary(model.fit)$coefficients

    SUP_est <- as.numeric(S1[, "Estimate"])
    SUP_se  <-sqrt(diag(sandwich::vcovHC(model.fit, type = "HC0")))
    SUP_lower <- SUP_est - 1.96 * SUP_se
    SUP_upper <- SUP_est + 1.96 * SUP_se
    SUP_width <- SUP_upper - SUP_lower

    out_SUP <- data.frame(
      split = split_seed,
      method = "Supervised",
      parameter = parameter_names,
      estimate = SUP_est,
      analytic_se = SUP_se,
      ci_lower = SUP_lower,
      ci_upper = SUP_upper,
      ci_width = SUP_width,
      benchmark = as.numeric(theta_full),
      covered = as.numeric(
        SUP_lower <= theta_full &
          SUP_upper >= theta_full
      ),
      se_type = "OLS",
      stringsAsFactors = FALSE
    )

    PSSE_fit <- PSSE1(
      labelled_data = as.matrix(labeled_final),
      unlabelled_data = as.matrix(unlabeled_final),
      c1 = NULL,
      type = "linear",
      tau = 0,
      alpha = 1,
      gamma = 10,
      sd = TRUE,
      Kfolds = 5
    )

    PSSE_est <- as.numeric(
      PSSE_fit$Hattheta
    )

    PSSE_se <- as.numeric(
      PSSE_fit$sd.of.hattheta
    )

    PSSE_lower <- PSSE_est - 1.96 * PSSE_se
    PSSE_upper <- PSSE_est + 1.96 * PSSE_se
    PSSE_width <- PSSE_upper - PSSE_lower

    out_PSSE <- data.frame(
      split = split_seed,
      method = "PSSE",
      parameter = parameter_names,
      estimate = PSSE_est,
      analytic_se = PSSE_se,
      ci_lower = PSSE_lower,
      ci_upper = PSSE_upper,
      ci_width = PSSE_width,
      benchmark = as.numeric(theta_full),
      covered = as.numeric(
        PSSE_lower <= theta_full &
          PSSE_upper >= theta_full
      ),
      se_type = "PSSE",
      stringsAsFactors = FALSE
    )

    DRESS_fit <- DRESS(
      labelled_data = as.matrix(labeled_final),
      unlabelled_data = as.matrix(unlabeled_final),
      L = 1,
      Kfolds = 5
    )

    DRESS_est <- as.numeric(
      DRESS_fit$Hattheta
    )

    DRESS_se <- as.numeric(
      DRESS_fit$sd.of.hattheta
    )

    DRESS_lower <- DRESS_est - 1.96 * DRESS_se
    DRESS_upper <- DRESS_est + 1.96 * DRESS_se
    DRESS_width <- DRESS_upper - DRESS_lower

    out_DRESS <- data.frame(
      split = split_seed,
      method = "DRESS",
      parameter = parameter_names,
      estimate = DRESS_est,
      analytic_se = DRESS_se,
      ci_lower = DRESS_lower,
      ci_upper = DRESS_upper,
      ci_width = DRESS_width,
      benchmark = as.numeric(theta_full),
      covered = as.numeric(
        DRESS_lower <= theta_full &
          DRESS_upper >= theta_full
      ),
      se_type = "DRESS",
      stringsAsFactors = FALSE
    )

    labeled_noseqn <- labeled %>%
      dplyr::select(-SEQN)

    unlabeled_noseqn <- unlabeled %>%
      dplyr::select(-SEQN)

    data_full_real <- as.data.frame(

      rbind(

        cbind(
          labeled_noseqn,
          D = 1
        ),

        cbind(
          unlabeled_noseqn,
          D = 0
        )
      )
    )

    data_full_real <- data_full_real %>%
      dplyr::select(
        Y,
        everything()
      )

    glm.model.formula <- as.formula(

      paste(
        "D ~",
        paste(
          colnames(data_full_real)[
            !colnames(data_full_real) %in%
              c("Y", "D")
          ],
          collapse = " + "
        )
      )
    )

    ps_initial <- glm(
      glm.model.formula,
      family = binomial(link = "logit"),
      data = as.data.frame(
        data_full_real[, -1]
      )
    )

    data_full_real$pi.hat <-
      fitted(ps_initial)

    data_full_real$ID <-
      seq_len(
        nrow(data_full_real)
      )

    ET <- estimate_theta_EM_kfold_dual_ET(

      th =
        theta_start,

      data_full =
        data_full_real,

      K =
        4,

      seed =2025,

      max.iter =
        500,

      eps =
        1e-4,

      damping = 1
    )

    data_all_ET <- ET$data_all

    x_cols_ET <- subset(
      data_all_ET,
      select =
        -c(
          Y,
          D,
          pi.hat,
          ID,
          y.hat
        )
    )

    X_ET <- model.matrix(
      ~ .,
      data = x_cols_ET
    )

    D_ET <- data_all_ET$D

    U_fun_ET <- function(theta) {

      U <- matrix(
        0,
        nrow = nrow(X_ET),
        ncol = ncol(X_ET)
      )

      id <- which(
        D_ET == 1
      )

      resid <- data_all_ET$Y[id] -
        as.numeric(
          X_ET[id, , drop = FALSE] %*%
            theta
        )

      U[id, ] <-
        X_ET[id, , drop = FALSE] *
        resid

      U
    }

    b_fun_ET <- function(theta) {

      resid_hat <-
        data_all_ET$y.hat -
        as.numeric(
          X_ET %*%
            theta
        )

      X_ET *
        resid_hat
    }

    ps_fit_ET <- glm(

      D_ET ~ .,

      data = data.frame(
        D_ET = D_ET,
        x_cols_ET
      ),

      family =
        binomial(),

      x =
        TRUE
    )

    phi_hat_ET <-
      coef(ps_fit_ET)

    O_ET <-
      ps_fit_ET$x

    ET_var <- gec_sandwich(

      phi_hat =
        phi_hat_ET,

      lambda_hat =
        ET$lambda,

      theta_hat =
        as.numeric(
          ET$theta
        ),

      D =
        D_ET,

      O =
        O_ET,

      U_fun =
        U_fun_ET,

      b_fun =
        b_fun_ET,

      entropy =
        "ET",

      parameter_names =
        colnames(X_ET)
    )

    ET_est <-
      as.numeric(
        ET_var$table$Estimate
      )

    ET_se <-
      as.numeric(
        ET_var$table$SE
      )

    ET_lower <-
      ET_est -
      1.96 * ET_se

    ET_upper <-
      ET_est +
      1.96 * ET_se

    ET_width <-
      ET_upper -
      ET_lower

    HD <- estimate_theta_EM_kfold_dual_HD(

      th =
        theta_start,

      data_full =
        data_full_real,

      K =
        4,

      seed =2025,

      max.iter =
        500,

      eps =
        1e-4,

      damping = 0.1
    )
    omega_hat<-HD$w
    eta1 <- HD$eta[HD$data_all$D == 1]

    data_all_HD <- HD$data_all

    x_cols_HD <- subset(
      data_all_HD,
      select =
        -c(
          Y,
          D,
          pi.hat,
          ID,
          y.hat
        )
    )

    X_HD <- model.matrix(
      ~ .,
      data = x_cols_HD
    )

    D_HD <- data_all_HD$D

    U_fun_HD <- function(theta) {

      U <- matrix(
        0,
        nrow = nrow(X_HD),
        ncol = ncol(X_HD)
      )

      id <- which(
        D_HD == 1
      )

      resid <- data_all_HD$Y[id] -
        as.numeric(
          X_HD[id, , drop = FALSE] %*%
            theta
        )

      U[id, ] <-
        X_HD[id, , drop = FALSE] *
        resid

      U
    }

    b_fun_HD <- function(theta) {

      resid_hat <-
        data_all_HD$y.hat -
        as.numeric(
          X_HD %*%
            theta
        )

      X_HD *
        resid_hat
    }

    ps_fit_HD <- glm(

      D_HD ~ .,

      data = data.frame(
        D_HD = D_HD,
        x_cols_HD
      ),

      family =
        binomial(),

      x =
        TRUE
    )

    phi_hat_HD <-
      coef(ps_fit_HD)

    O_HD <-
      ps_fit_HD$x

    HD_var <- gec_sandwich(

      phi_hat =
        phi_hat_HD,

      lambda_hat =
        HD$lambda,

      theta_hat =
        as.numeric(
          HD$theta
        ),

      D =
        D_HD,

      O =
        O_HD,

      U_fun =
        U_fun_HD,

      b_fun =
        b_fun_HD,

      entropy =
        "HD",

      parameter_names =
        colnames(X_HD)
    )

    HD_est <-
      as.numeric(
        HD_var$table$Estimate
      )

    HD_se <-
      as.numeric(
        HD_var$table$SE
      )

    HD_lower <-
      HD_est -
      1.96 * HD_se

    HD_upper <-
      HD_est +
      1.96 * HD_se

    HD_width <-
      HD_upper -
      HD_lower

    CE <- estimate_theta_EM_kfold_CE(
      th = as.numeric(coef(model.fit)),
      data_full = data_full_real,
      K = 4,

      seed = 2025,

      max.iter = 500,
      eps = 1e-3,
      lambda_tol = 1e-4,
      lambda_maxit = 1000,
      damping = 0.1
    )

    data_all_CE <- CE$data_all

    x_cols_CE <- subset(
      data_all_CE,
      select = -c(
        Y,
        D,
        pi.hat,
        ID,
        y.hat
      )
    )

    model.formula_CE <- as.formula(
      paste(
        "~",
        paste(
          colnames(x_cols_CE),
          collapse = "+"
        )
      )
    )

    X_CE <- model.matrix(
      model.formula_CE,
      data = as.data.frame(
        data_all_CE[
          ,
          colnames(x_cols_CE),
          drop = FALSE
        ]
      )
    )

    D_CE <- data_all_CE$D

    ps_formula_CE <- as.formula(
      paste(
        "D ~",
        paste(
          colnames(x_cols_CE),
          collapse = "+"
        )
      )
    )

    ps_fit_CE <- glm(
      ps_formula_CE,
      family = binomial(),
      data = data_all_CE,
      x = TRUE
    )

    O_CE <-
      ps_fit_CE$x

    U_fun_CE <- function(theta) {

      resid <- data_all_CE$Y -
        as.numeric(X_CE %*% theta)

      X_CE * resid
    }

    b_fun_CE <- function(theta) {

      resid_hat <- data_all_CE$y.hat -
        as.numeric(X_CE %*% theta)

      X_CE * resid_hat
    }

    CE_var <- gec_sandwich(
      phi_hat = coef(ps_fit_CE),
      lambda_hat = CE$lambda,
      theta_hat = as.numeric(CE$theta),

      D = D_CE,
      O = O_CE,
      U_fun = U_fun_CE,
      b_fun = b_fun_CE,
      entropy = "CE"
    )

    CE_est <- as.numeric(
      CE$theta
    )

    CE_se <- sqrt(
      diag(
        CE_var$V_theta
      )
    )

    CE_lower <-
      CE_est - 1.96 * CE_se

    CE_upper <-
      CE_est + 1.96 * CE_se

    CE_width <-
      2 * 1.96 * CE_se

    out_CE <- data.frame(

      split =
        split_seed,

      method =
        "CE",

      parameter =
        parameter_names,

      estimate =
        CE_est,

      analytic_se =
        CE_se,

      ci_lower =
        CE_lower,

      ci_upper =
        CE_upper,

      ci_width =
        CE_width,

      benchmark =
        as.numeric(
          theta_full
        ),

      covered =
        as.numeric(
          CE_lower <= theta_full &
            CE_upper >= theta_full
        ),
      se_type = "GEC sandwich"
    )

    out_ET <- data.frame(

      split =
        split_seed,

      method =
        "ET",

      parameter =
        parameter_names,

      estimate =
        ET_est,

      analytic_se = ET_se,

      ci_lower =
        ET_lower,

      ci_upper =
        ET_upper,

      ci_width =
        ET_width,

      benchmark =
        as.numeric(
          theta_full
        ),

      covered =
        as.numeric(
          ET_lower <= theta_full &
            ET_upper >= theta_full
        ),
      se_type = "GEC sandwich"
    )

    out_HD <- data.frame(

      split =
        split_seed,

      method =
        "HD",

      parameter =
        parameter_names,

      estimate =
        HD_est,

      analytic_se =
        HD_se,

      ci_lower =
        HD_lower,

      ci_upper =
        HD_upper,

      ci_width =
        HD_width,

      benchmark =
        as.numeric(
          theta_full
        ),

      covered =
        as.numeric(
          HD_lower <= theta_full &
            HD_upper >= theta_full
        ),
      se_type = "GEC sandwich"
    )

    bind_rows(
      out_SUP,
      out_DRESS,
      out_PSSE,
      out_ET,
      out_HD,
      out_CE
    )

  }, error = function(e) {

    message(
      "Split ",
      split_seed,
      " failed: ",
      conditionMessage(e)
    )

    NULL
  })
}
