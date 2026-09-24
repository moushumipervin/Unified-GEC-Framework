###############################################################################
# REQUIRED PACKAGES
###############################################################################

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



GenerateData<- function(n,N,p,OR, MAR){
  sigma= 2  # the sd of the error term eta   
  
  alpha0=1
  alpha1=rep(1,p)
  alpha2=rep(1,p)
  
  X=mvrnorm(n+N,rep(0,p),diag(rep(1,p)))
  
  if(OR==1){
    
    Y=alpha0+X%*%alpha1+rnorm(n+N,0,1) ####OR1
    
  }else{
    Y=alpha0+X%*%alpha1+(X^3-X^2+exp(X))%*%alpha2+rnorm(n+N,0,2) ###OR2 
  }
  
  
  if(MAR==1){

    px <- 1 / (1 + exp(1 + X[, 1] + 0.5 * X[, 2] - 0.5 * X[, 3] - 0.1 * X[, 4])) ###MAR mechanism
    D=rbinom(n+N,1,px)

  }else{
    D=rbinom(n+N,1,n/(n+N)) ###MCAR mechanism
  }
  
  data_full<-cbind(D,Y,X)
  
  
  colnames(data_full)<-c("D","Y",paste0("X",1:p))
  data_full<-as.data.frame(data_full)
  data_full$Y[data_full$D == 0] <- NA
  data_full$ID <- seq_len(nrow(data_full))
  pi.hat1<-fitted(glm(D~.,data=subset(data_full, select = -c(Y, ID)),family=binomial()))
  
  data_full$pi.hat<-  pi.hat1
  
  
  
  data_labelled<-as.matrix(data_full[D==1,c("Y",paste0("X",1:p))])
  data_unlabelled<-as.matrix(data_full[D==0,paste0("X",1:p)])
  
  return(list(Data.labelled =data_labelled,Data.unlabelled=data_unlabelled,data_full=data_full))
}


#---------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------#
#
#
#                                   Supervised coefficients estimation functions
#          (mainly for linear working model;logistic working model;quantile working model)
#
#---------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------#


SupervisedEst<- function(Data){
  #-----------------------------------------Arguments--------------------------------------------------------#
  # Purpose: This function is to calculate the supervised coefficients estimate based on the given labelled data. 
  #
  # Input: 
  #      Data: A matrix, whose each row is an observation of the response and the predictor vector and first
  #              column is the observation vector of the response. 
  #      option: The type of data settings. We offer nine choices: "i","ii","iii","W1","W2","W3",S1","S2","S3".
  #               These nine choices correspond to settings (i) - (iii) in Section 4.1 of the paper, 
  #               settings (W1) - (W2), settings (S1) - (S3) in Section 3.2 of the supplementary material 
  #               in sequence. 
  #      tau: For "option" equal to "iii", "W3" or "S3", the quantile level; otherwise, this quantity is useless. 
  #
  # Output: 
  #       Est.coef: the estimated value of the coefficients based on the given data 
  #----------------------------------------------------------------------------------------------------------#
  
  n=nrow(Data)
  if (n>2500){
    print(paste("NOTE: The coefficients estimate is computed only using labelled data of size ",n,".",sep=""))
  }
  
  
  
    Target<- lm(as.vector(Data[,1])~Data[,-1],data=data.frame(Data))$coefficients
    
  return(list("Est.coef"=as.vector(Target)))
  
  
}

###############################################################################
# Entropy machinery
###############################################################################

U_fun1 <- function(theta) {
  
  U <- matrix(
    0,
    nrow = nrow(X),
    ncol = ncol(X)
  )
  
  id <- which(D == 1)
  
  resid <- data_all$Y[id] -
    as.numeric(
      X[id, , drop = FALSE] %*% theta
    )
  
  U[id, ] <- X[id, , drop = FALSE] * resid
  
  U
}

b_fun1 <- function(theta) {
  
  resid_hat <- data_all$y.hat -
    as.numeric(
      X %*% theta
    )
  
  X * resid_hat
}
gec_entropy <- function(entropy = c("SL", "EL", "ET", "HD", "CE")) {
  
  entropy <- match.arg(entropy)
  
  if (entropy == "SL") {
    
    # G(w) = w^2 / 2
    # g(w) = w
    # g^{-1}(eta) = eta
    
    g <- function(w) w
    
    ginv <- function(eta) eta
    
    gp <- function(w) rep(1, length(w))
    
    debias <- function(pi) 1 / pi
    
    valid_eta <- function(eta) rep(TRUE, length(eta))
  }
  
  
  if (entropy == "EL") {
    
    # G(w) = -log(w)
    # g(w) = -1/w
    # g^{-1}(eta) = -1/eta
    # eta < 0
    
    g <- function(w) -1 / w
    
    ginv <- function(eta) -1 / eta
    
    gp <- function(w) 1 / w^2
    
    debias <- function(pi) -pi
    
    valid_eta <- function(eta) eta < 0
  }
  
  
  if (entropy == "ET") {
    
    # G(w) = w log(w) - w
    # g(w) = log(w)
    # g^{-1}(eta) = exp(eta)
    
    g <- function(w) log(w)
    
    ginv <- function(eta) exp(eta)
    
    gp <- function(w) 1 / w
    
    debias <- function(pi) log(1 / pi)
    
    valid_eta <- function(eta) rep(TRUE, length(eta))
  }
  
  
  if (entropy == "HD") {
    
    # G(w) = -sqrt(w)
    # g(w) = -1/(2 sqrt(w))
    # g^{-1}(eta) = 1/(4 eta^2)
    # eta < 0
    
    g <- function(w) -1 / (2 * sqrt(w))
    
    ginv <- function(eta) 1 / (4 * eta^2)
    
    gp <- function(w) 1 / (4 * w^(3/2))
    
    debias <- function(pi) -sqrt(pi) / 2
    
    valid_eta <- function(eta) eta < 0
  }
  
  
  if (entropy == "CE") {
    
    # G(w) = (w-1)log(w-1) - w log(w)
    #
    # g(w) = log((w-1)/w)
    #      = log(1 - 1/w)
    #
    # g^{-1}(eta) = 1/(1-exp(eta))
    # eta < 0
    
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





###############################################################################
# Logistic PS
###############################################################################

pi_logistic <- function(phi, O) {
  
  eta <- as.numeric(O %*% phi)
  
  plogis(eta)
}


h_logistic <- function(phi, O) {
  
  pi <- pi_logistic(phi, O)
  
  O * pi
}









###############################################################################
# Build individual stacked estimating equations Psi_i(beta)
###############################################################################

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
  
  # ---------------------------------------------------------------
  # Split beta = (phi, lambda, theta)
  # ---------------------------------------------------------------
  
  id_phi <- seq_len(n_phi)
  
  id_lambda <- n_phi + seq_len(n_lambda)
  
  id_theta <- n_phi + n_lambda + seq_len(n_theta)
  
  phi <- beta[id_phi]
  lambda <- beta[id_lambda]
  theta <- beta[id_theta]
  
  
  # ---------------------------------------------------------------
  # Propensity
  # ---------------------------------------------------------------
  
  pi <- pi_logistic(phi, O)
  
 # #pi <- pmin(
   # pmax(pi, 1e-8),
    #1 - 1e-8
  #)
  
  h <- h_logistic(phi, O)
  
  
  # ---------------------------------------------------------------
  # b_i(theta)
  #
  # User-supplied function returning N x q_b matrix
  # ---------------------------------------------------------------
  
  b <- as.matrix(
    b_fun(theta)
  )
  
  
  # ---------------------------------------------------------------
  # Debiasing covariate g(pi^{-1})
  # ---------------------------------------------------------------
  
  gpi <- ENT$debias(pi)
  
  
  # ---------------------------------------------------------------
  # s_i = ( b_i^T, g(pi^{-1}) )^T
  # ---------------------------------------------------------------
  
  S <- cbind(
    b,
    gpi
  )
  
  if (ncol(S) != n_lambda) {
    stop("n_lambda does not match number of calibration covariates.")
  }
  
  
  # ---------------------------------------------------------------
  # omega_i = g^{-1}(lambda^T s_i)
  # ---------------------------------------------------------------
  
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
  
  
  # ---------------------------------------------------------------
  # U_i(theta)
  #
  # N x q_theta
  # ---------------------------------------------------------------
  
  U <- as.matrix(
    U_fun(theta)
  )
  
  
  # ---------------------------------------------------------------
  # Block 1:
  # (D/pi - 1) h_i(phi)
  # ---------------------------------------------------------------
  
  Psi_phi <- h *
    as.numeric(D / pi - 1)
  
  
  # ---------------------------------------------------------------
  # Block 2:
  # D omega_i s_i - s_i
  # ---------------------------------------------------------------
  
  Psi_lambda <- S *
    as.numeric(D * omega - 1)
  
  
  # ---------------------------------------------------------------
  # Block 3:
  # D omega_i U_i(theta)
  #
  # Important:
  # U may be unavailable when D=0.
  # Set those rows to zero before multiplying.
  # ---------------------------------------------------------------
  
  U[D == 0, ] <- 0
  
  Psi_theta <- U *
    as.numeric(D * omega)
  
  
  cbind(
    Psi_phi,
    Psi_lambda,
    Psi_theta
  )
}


###############################################################################
# Full joint sandwich
###############################################################################

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
  
  
  #################################################################
  # Psi_i at beta_hat
  #################################################################
  
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
  
  
  #################################################################
  # B_hat
  #################################################################
  
  B_hat <- crossprod(Psi_hat) / N
  
  
  #################################################################
  # Mean estimating equation
  #################################################################
  
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
  
  
  #################################################################
  # Numerical Jacobian
  #################################################################
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
  
  
  #################################################################
  # Sandwich
  #################################################################
  
  A_inv <- solve(A_hat)
  
  V_beta <- (
    A_inv %*%
      B_hat %*%
      t(A_inv)
  ) / N
  
  
  #################################################################
  # Extract theta block
  #################################################################
  
 
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





solve_lambda_ET_dual <- function(
    b_mat,
    pi_hat,
    D,
    lambda_start = NULL,
    maxit = 1000,
    reltol = 1e-10) {
  
  # ------------------------------------------------------------
  # s_i = [ b_i , log(pi_i^{-1}) ]
  # ------------------------------------------------------------
  
  g_pi <- log(1 / pi_hat)
  
  S <- cbind(
    b_mat,
    g_pi
  )
  
  q <- ncol(S)
  
  if (is.null(lambda_start)) {
    
    # Correct-PS motivated initial value
    lambda_start <- c(
      rep(0, q - 1),
      1
    )
  }
  
  
  # ------------------------------------------------------------
  # ET dual objective
  #
  # F(eta) = exp(eta)
  #
  # Q(lambda)
  # = 1/N sum [ D_i exp(eta_i) - eta_i ]
  # ------------------------------------------------------------
  
  objective <- function(lambda) {
    
    eta <- as.numeric(
      S %*% lambda
    )
    
    mean(
      D * exp(eta) - eta
    )
  }
  
  
  # ------------------------------------------------------------
  # Gradient = calibration equation
  # ------------------------------------------------------------
  
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
  
  # ============================================================
  # Cross-fitting
  # ============================================================
  
  fold_hat <- k_fold_function(
    df = data_full,
    K = K,
    seed = seed
  )
  
  data_all <- do.call(
    rbind,
    fold_hat
  )
  
  
  # ============================================================
  # Design matrix
  # ============================================================
  
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
  
  
  # ============================================================
  # Alternating updates
  # ============================================================
  
  repeat {
    
    iter <- iter + 1
    
    
    # ----------------------------------------------------------
    # b_i(theta)
    # ----------------------------------------------------------
    
    error_hat <- y.hat -
      as.numeric(
        X %*% theta
      )
    
    b_mat <- X * error_hat
    
    
    # ----------------------------------------------------------
    # Direct ET lambda estimation
    # ----------------------------------------------------------
    
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
    
    
    # ----------------------------------------------------------
    # Weighted regression
    # ----------------------------------------------------------
    
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
    
    
   
    # ----------------------------------------------------------
    # Convergence
    # ----------------------------------------------------------
    
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
      
      message(
        "Maximum iterations reached."
      )
      
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
          
          converged = FALSE,
          
          iterations = iter
        )
      )
    }
    
    
    # damping = 1 means identical full theta update
    theta <-
      (1 - damping) * theta +
      damping * th.raw
  }
}

solve_lambda_CE_dual <- function(
    b_mat,
    pi_hat,
    D,
    lambda_start = NULL,
    margin = 1e-8,
    maxit = 1000,
    reltol = 1e-10) {
  
  # ------------------------------------------------------------
  # s_i = [ b_i , g(pi_i^{-1}) ]
  #
  # CE:
  # g(pi^{-1}) = log(1-pi)
  # ------------------------------------------------------------
  
  g_pi <- log(1 - pi_hat)
  S <- cbind(
    b_mat,
    g_pi
  )
  
  q <- ncol(S)
  
  I1 <- which(D == 1)
  
  
  # ============================================================
  # ADD THE FEASIBILITY-RESET BLOCK HERE
  # ============================================================
  
  lambda_natural <- c(
    rep(0, q - 1),
    1
  )
  
  
  # ------------------------------------------------------------
  # Starting value
  #
  # lambda = (0,...,0,1)
  #
  # gives eta = log(1-pi) < 0
  # so it is naturally feasible.
  # ------------------------------------------------------------
  if (is.null(lambda_start)) {
    
    lambda_start <- lambda_natural
    
  } else {
    
    lambda_start <- as.numeric(lambda_start)
  }
  
  
  eta_start <- as.numeric(
    S[I1, , drop = FALSE] %*%
      lambda_start
  )
  
  
  # Previous warm start may no longer be valid
  # after theta changes and therefore S changes
  if (any(eta_start >= -margin)) {
    
    lambda_start <- lambda_natural
    
    eta_start <- as.numeric(
      S[I1, , drop = FALSE] %*%
        lambda_start
    )
  }
  
  # Make sure even the natural start is feasible
  if (any(eta_start >= -margin)) {
    
    stop(
      "Could not construct a feasible CE lambda starting value."
    )
  }
  
  # ------------------------------------------------------------
  # Dual objective
  #
  # F(eta)
  # =
  # eta - log(1-exp(eta))
  #
  # domain eta < 0 for respondents
  # ------------------------------------------------------------
  
  objective <- function(lambda) {
    
    eta <- as.numeric(
      S %*% lambda
    )
    
    eta1 <- eta[I1]
    
    # safety
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
  
  
  # ------------------------------------------------------------
  # Gradient
  #
  # exactly:
  #
  # mean[
  #   s_i {D_i w_i - 1}
  # ]
  # ------------------------------------------------------------
  
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
  
  
  # ------------------------------------------------------------
  # Linear constraints:
  #
  # eta_i = S_i lambda <= -margin
  #
  # constrOptim uses
  #
  # ui %*% lambda - ci >= 0
  #
  # Therefore:
  #
  # -S_i lambda >= margin
  # ------------------------------------------------------------
  
  ui <- -S[
    I1,
    ,
    drop = FALSE
  ]
  
  ci <- rep(
    margin,
    length(I1)
  )
  
  
  # ------------------------------------------------------------
  # Check initial feasibility
  # ------------------------------------------------------------
  
  eta_start <- as.numeric(
    S[I1, , drop = FALSE] %*%
      lambda_start
  )
  
  if (any(eta_start >= -margin)) {
    
    stop(
      "Initial lambda is not strictly inside the CE dual domain."
    )
  }
  
  
  # ------------------------------------------------------------
  # Optimize dual
  # ------------------------------------------------------------
  
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
  
  
  # ------------------------------------------------------------
  # Final weights
  # ------------------------------------------------------------
  
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
  
  
  # ------------------------------------------------------------
  # Calibration score
  # ------------------------------------------------------------
  
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
  
  # ============================================================
  # 1. Cross-fitted predictions
  # ============================================================
  
  fold_hat <- k_fold_function(
    df = data_full,
    K = K,
    seed = seed
  )
  
  data_all <- do.call(
    rbind,
    fold_hat
  )
  
  
  # ============================================================
  # 2. Design matrix
  # ============================================================
  
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
  
  # ------------------------------------------------------------
  # Warm start for lambda
  # ------------------------------------------------------------
  
  lambda_current <- NULL
  
  
  # ============================================================
  # 3. Alternating theta/lambda iterations
  # ============================================================
  
  repeat {
    
    iter <- iter + 1
    
    
    # ----------------------------------------------------------
    # b_i(theta)
    #
    # b_i(theta)
    # = X_i { yhat_i - X_i' theta }
    # ----------------------------------------------------------
    
    error_hat <- y.hat -
      as.numeric(
        X %*% theta
      )
    
    b_mat <- X * error_hat
    
    
    # ----------------------------------------------------------
    # Estimate lambda directly from CE dual
    # ----------------------------------------------------------
    
    lambda_fit <- solve_lambda_CE_dual(
      b_mat = b_mat,
      pi_hat = data_all$pi.hat,
      D = D,
      lambda_start = lambda_current,
      maxit = lambda_maxit,
      reltol = lambda_tol
    )
    
    
    # ----------------------------------------------------------
    # Check lambda solution
    # ----------------------------------------------------------
    
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
    
    # warm start next iteration
    lambda_current <- lambda_hat
    
    W <- lambda_fit$weights
    
    
    # ----------------------------------------------------------
    # Weighted regression update
    # ----------------------------------------------------------
    
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
    
    
    # ----------------------------------------------------------
    # Difference before damping
    # ----------------------------------------------------------
    
    diff_theta <- max(
      abs(
        th.raw - theta
      )
    )
    
    
    #cat(
    #  "iter =", iter,
    #  " theta diff =", diff_theta
    #)
    
    
    # ==========================================================
    # 4. Check convergence
    # ==========================================================
    
    if (diff_theta < eps) {
      
      theta_final <- th.raw
      
      
      # --------------------------------------------------------
      # IMPORTANT:
      #
      # Recompute b(theta_final) and lambda one final time
      # so returned lambda corresponds exactly to returned theta.
      # --------------------------------------------------------
      
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
      
      
      # --------------------------------------------------------
      # One final theta update using final weights
      # --------------------------------------------------------
      
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
    
    
    # ==========================================================
    # 5. Maximum iterations
    # ==========================================================
    
    if (iter >= max.iter) {
      
      message(
        "Maximum outer iterations reached."
      )
      
      
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
          
         
          converged = FALSE,
          
          iterations = iter
        )
      )
    }
    
    
    # ==========================================================
    # 6. Damped theta update
    # ==========================================================
    
    theta <- as.numeric(
      (1 - damping) * theta +
        damping * th.raw
    )
  }
}


solve_lambda_HD_dual <- function(
    b_mat,
    pi_hat,
    D,
    lambda_start = NULL,
    margin = 1e-8,
    maxit = 1000,
    reltol = 1e-10) {
  
  # ------------------------------------------------------------
  # s_i = [ b_i , g(pi_i^{-1}) ]
  #
  # HD:
  # g(pi^{-1}) = -sqrt(pi)/2
  # ------------------------------------------------------------
  
  g_pi <- -sqrt(pi_hat) / 2
  
  S <- cbind(
    b_mat,
    g_pi
  )
  
  q <- ncol(S)
  
  I1 <- which(D == 1)
  
  
  # ------------------------------------------------------------
  # Starting value
  #
  # lambda = (0,...,0,1)
  #
  # then eta_i = g(pi_i^{-1}) < 0
  # ------------------------------------------------------------
  
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
  
  
  # ------------------------------------------------------------
  # Check feasibility of starting lambda
  # ------------------------------------------------------------
  
  eta_start <- as.numeric(
    S[I1, , drop = FALSE] %*%
      lambda_start
  )
  
  # If warm start is no longer feasible,
  # reset to natural starting value
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
  
  
  # ============================================================
  # HD dual objective
  #
  # w(eta) = 1/(4 eta^2)
  #
  # F'(eta) = w(eta)
  #
  # therefore
  #
  # F(eta) = -1/(4 eta)
  #
  # Dual objective:
  #
  # mean[ D_i F(eta_i) - eta_i ]
  # ============================================================
  
  objective <- function(lambda) {
    
    eta <- as.numeric(
      S %*% lambda
    )
    
    eta1 <- eta[I1]
    
    
    # HD domain
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
  
  
  # ============================================================
  # Gradient
  #
  # mean[
  #   s_i {D_i w_i - 1}
  # ]
  #
  # This is exactly the calibration estimating equation.
  # ============================================================
  
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
  
  
  # ============================================================
  # Domain constraints:
  #
  # eta_i = S_i lambda <= -margin
  #
  # constrOptim requires
  #
  # ui %*% lambda - ci >= 0
  #
  # therefore:
  #
  # -S_i lambda >= margin
  # ============================================================
  
  ui <- -S[
    I1,
    ,
    drop = FALSE
  ]
  
  ci <- rep(
    margin,
    length(I1)
  )
  
  
  # ============================================================
  # Optimize
  # ============================================================
  
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
  
  
  # ============================================================
  # Final eta and weights
  # ============================================================
  
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
  
  
  # ------------------------------------------------------------
  # Calibration score
  # ------------------------------------------------------------
  
  score <- gradient(
    lambda_hat
  )
  
  
  return(
    list(
      
      # PAPER convention lambda
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
  
  # ============================================================
  # 1. Cross-fitted predictions
  # ============================================================
  
  fold_hat <- k_fold_function(
    df = data_full,
    K = K,
    seed = seed
  )
  
  data_all <- do.call(
    rbind,
    fold_hat
  )
  
  
  # ============================================================
  # 2. Design matrix
  # ============================================================
  
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
  
  
  # ------------------------------------------------------------
  # Warm start for lambda
  # ------------------------------------------------------------
  
  lambda_current <- NULL
  
  
  # ============================================================
  # 3. Alternating lambda / theta updates
  # ============================================================
  
  repeat {
    
    iter <- iter + 1
    
    
    # ----------------------------------------------------------
    # b_i(theta)
    #
    # b_i(theta)
    # =
    # X_i { yhat_i - X_i' theta }
    # ----------------------------------------------------------
    
    error_hat <- y.hat -
      as.numeric(
        X %*% theta
      )
    
    b_mat <- X *
      error_hat
    
    
    # ----------------------------------------------------------
    # Estimate HD lambda directly
    # ----------------------------------------------------------
    
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
    
    
    # Use current solution as starting value
    # for next outer iteration
    lambda_current <-
      lambda_hat
    
    
    W <-
      lambda_fit$weights
    
    
    # ==========================================================
    # 4. Weighted regression theta update
    # ==========================================================
    
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
    
    
    # ----------------------------------------------------------
    # Difference from current theta
    # ----------------------------------------------------------
    
    diff_theta <- max(
      abs(
        th.raw -
          theta
      )
    )
    
    
    # ----------------------------------------------------------
    # Diagnostics
    # ----------------------------------------------------------
    
    ESS <- sum(W)^2 /
      sum(W^2)
    
    ESS_ratio <-
      ESS / length(W)
    
    
    
    # ==========================================================
    # 5. Convergence
    # ==========================================================
    
    if (diff_theta < eps) {
      
      theta_final <-
        th.raw
      
      
      # --------------------------------------------------------
      # Recalculate b at final theta
      # --------------------------------------------------------
      
      error_final <- y.hat -
        as.numeric(
          X %*% theta_final
        )
      
      b_final <- X *
        error_final
      
      
      # --------------------------------------------------------
      # Re-estimate lambda at final theta
      # --------------------------------------------------------
      
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
      
      
      # --------------------------------------------------------
      # One final theta update
      # --------------------------------------------------------
      
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
          
          # IMPORTANT:
          # already PAPER convention
          # NO minus sign needed
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
    
    
    # ==========================================================
    # 6. Maximum iterations
    # ==========================================================
    
    if (iter >= max.iter) {
      
      message(
        "Maximum outer iterations reached."
      )
      
      
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
    
    
    # ==========================================================
    # 7. Theta update
    #
    # damping = 1:
    # same full update as your original CVXR algorithm.
    # ==========================================================
    
    theta <- as.numeric(
      (1 - damping) * theta +
        damping * th.raw
    )
  }
}



##############Using efficient influence function approach to estimate the variance

## ============================================================
##  IF-based variance for GEC
##  Compare directly with your joint sandwich V_theta
## ============================================================

if (!requireNamespace("numDeriv", quietly = TRUE)) {
  stop("Please install numDeriv: install.packages('numDeriv')")
}


## ------------------------------------------------------------
## 1. Entropy-specific functions
## ------------------------------------------------------------

gec_entropy_IF <- function(entropy = c("SL", "EL", "ET", "HD", "CE")) {
  
  entropy <- match.arg(entropy)
  
  if (entropy == "SL") {
    
    ginv <- function(eta) eta
    
    ## derivative of g^{-1}(eta)
    fprime <- function(eta) rep(1, length(eta))
    
    debias <- function(pi) 1 / pi
  }
  
  if (entropy == "EL") {
    
    ginv <- function(eta) -1 / eta
    
    fprime <- function(eta) 1 / eta^2
    
    debias <- function(pi) -pi
  }
  
  if (entropy == "ET") {
    
    ginv <- function(eta) exp(eta)
    
    fprime <- function(eta) exp(eta)
    
    debias <- function(pi) log(1 / pi)
  }
  
  if (entropy == "HD") {
    
    ginv <- function(eta) 1 / (4 * eta^2)
    
    fprime <- function(eta) -1 / (2 * eta^3)
    
    ## IMPORTANT: corrected HD convention
    debias <- function(pi) -sqrt(pi) / 2
  }
  
  if (entropy == "CE") {
    
    ginv <- function(eta) 1 / (1 - exp(eta))
    
    fprime <- function(eta) {
      exp(eta) / (1 - exp(eta))^2
    }
    
    debias <- function(pi) log(1 - pi)
  }
  
  list(
    ginv = ginv,
    fprime = fprime,
    debias = debias
  )
}


## ------------------------------------------------------------
## 2. Logistic propensity score pieces
## ------------------------------------------------------------

pi_logistic_IF <- function(phi, O) {
  
  plogis(as.numeric(O %*% phi))
}


## In your notation:
##
## h(phi) = {1-pi(phi)}^{-1} d pi(phi)/d phi
##
## For logistic PS:
##
## d pi/d phi = pi(1-pi)O
##
## therefore h = pi O
##

h_logistic_IF <- function(phi, O) {
  
  pi <- pi_logistic_IF(phi, O)
  
  O * pi
}


## ------------------------------------------------------------
## 3. Main IF variance function
## ------------------------------------------------------------

gec_IF_variance <- function(
    phi_hat,
    lambda_hat,
    theta_hat,
    D,
    O,
    U_fun,
    b_fun,
    entropy = c("SL", "EL", "ET", "HD", "CE"),
    parameter_names = NULL
) {
  
  entropy <- match.arg(entropy)
  
  ENT <- gec_entropy_IF(entropy)
  
  phi_hat    <- as.numeric(phi_hat)
  lambda_hat <- as.numeric(lambda_hat)
  theta_hat  <- as.numeric(theta_hat)
  
  D <- as.numeric(D)
  
  N <- length(D)
  
  ## ----------------------------------------------------------
  ## A. propensity score
  ## ----------------------------------------------------------
  
  pi_hat <- pi_logistic_IF(phi_hat, O)
  
  pi_hat <- pmin(
    pmax(pi_hat, 1e-8),
    1 - 1e-8
  )
  
  h_hat <- h_logistic_IF(phi_hat, O)
  
  
  ## ----------------------------------------------------------
  ## B. Calibration functions S
  ## ----------------------------------------------------------
  
  b_hat <- as.matrix(
    b_fun(theta_hat)
  )
  
  gpi_hat <- ENT$debias(pi_hat)
  
  S_hat <- cbind(
    b_hat,
    gpi_hat
  )
  
  
  if (ncol(S_hat) != length(lambda_hat)) {
    stop(
      "length(lambda_hat) must equal ncol(S_hat)"
    )
  }
  
  
  ## ----------------------------------------------------------
  ## C. GEC weights
  ## ----------------------------------------------------------
  
  eta_hat <- as.numeric(
    S_hat %*% lambda_hat
  )
  
  omega_hat <- ENT$ginv(eta_hat)
  
  
  ## ----------------------------------------------------------
  ## D. Outcome score U(theta)
  ## ----------------------------------------------------------
  
  U_hat <- as.matrix(
    U_fun(theta_hat)
  )
  
  q <- ncol(U_hat)
  r <- ncol(S_hat)
  
  
  ## ----------------------------------------------------------
  ## E. gamma_hat
  ##
  ## gamma =
  ## E[D f'(lambda'S) U S']
  ## [E[D f'(lambda'S) S S']]^{-1}
  ## ----------------------------------------------------------
  
  fp_hat <- ENT$fprime(eta_hat)
  
  a_gamma <- D * fp_hat
  
  
  ## numerator: q x r
  Gamma_num <-
    crossprod(
      U_hat,
      S_hat * a_gamma
    ) / N
  
  
  ## denominator: r x r
  Gamma_den <-
    crossprod(
      S_hat,
      S_hat * a_gamma
    ) / N
  
  
  gamma_hat <-
    Gamma_num %*%
    solve(Gamma_den)
  
  
  ## gamma S_i for each observation
  ##
  ## gamma is q x r
  ## S is N x r
  ##
  ## result N x q
  gammaS_hat <-
    S_hat %*% t(gamma_hat)
  
  
  ## ----------------------------------------------------------
  ## F. kappa_hat
  ##
  ## Rather than manually differentiating the complicated
  ## expression, calculate the derivative numerically.
  ##
  ## This is also useful for checking A9(i).
  ## ----------------------------------------------------------
  
  linearized_mean_phi <- function(phi) {
    
    pi <- pi_logistic_IF(phi, O)
    
    pi <- pmin(
      pmax(pi, 1e-8),
      1 - 1e-8
    )
    
    gpi <- ENT$debias(pi)
    
    S <- cbind(
      b_hat,
      gpi
    )
    
    eta <- as.numeric(
      S %*% lambda_hat
    )
    
    omega <- ENT$ginv(eta)
    
    
    ## gamma S
    gammaS <-
      S %*% t(gamma_hat)
    
    
    ## quantity inside derivative:
    ##
    ## gamma S +
    ## D omega {U - gamma S}
    ##
    
    M <-
      gammaS +
      (U_hat - gammaS) *
      as.numeric(D * omega)
    
    
    colMeans(M)
  }
  
  
  ## q x p_phi derivative
  Dphi_hat <-
    numDeriv::jacobian(
      func = linearized_mean_phi,
      x = phi_hat
    )
  
  
  ## ----------------------------------------------------------
  ## weighted h inner-product matrix
  ##
  ## E[(1-pi)/pi h h']
  ## ----------------------------------------------------------
  
  h_weight <-
    (1 - pi_hat) / pi_hat
  
  Hh_hat <-
    crossprod(
      h_hat,
      h_hat * h_weight
    ) / N
  
  
  ## ----------------------------------------------------------
  ## IMPORTANT:
  ##
  ## corrected A9 sign convention
  ##
  ## kappa = - derivative * Hh^{-1}
  ## ----------------------------------------------------------
  
  kappa_hat <-
    -Dphi_hat %*%
    solve(Hh_hat)
  
  
  ## ----------------------------------------------------------
  ## G. Construct d_i
  ## ----------------------------------------------------------
  
  ## component 1
  part1 <- gammaS_hat
  
  
  ## component 2
  part2 <-
    (U_hat - gammaS_hat) *
    as.numeric(D * omega_hat)
  
  
  ## kappa h_i
  ##
  ## h_hat = N x p_phi
  ## kappa = q x p_phi
  ##
  ## result N x q
  kappa_h_hat <-
    h_hat %*% t(kappa_hat)
  
  
  ## component 3
  part3 <-
    kappa_h_hat *
    as.numeric(
      1 - D / pi_hat
    )
  
  
  ## influence contribution d_i
  d_hat <-
    part1 +
    part2 +
    part3
  
  
  ## ----------------------------------------------------------
  ## H. tau_hat
  ##
  ## tau =
  ## E[d/dtheta {D omega U(theta)}]
  ##
  ## Compute numerically so that theta-dependence of
  ## b(theta), S(theta), and omega(theta) is included.
  ## ----------------------------------------------------------
  
  weighted_score_mean <- function(theta) {
    
    b <- as.matrix(
      b_fun(theta)
    )
    
    S <- cbind(
      b,
      gpi_hat
    )
    
    eta <- as.numeric(
      S %*% lambda_hat
    )
    
    omega <- ENT$ginv(eta)
    
    U <- as.matrix(
      U_fun(theta)
    )
    
    colMeans(
      U * as.numeric(D * omega)
    )
  }
  
  
  tau_hat <-
    numDeriv::jacobian(
      func = weighted_score_mean,
      x = theta_hat
    )
  
  
  ## ----------------------------------------------------------
  ## I. empirical Var(d_i)
  ##
  ## use population-style 1/N scaling to correspond directly
  ## with the asymptotic formula
  ## ----------------------------------------------------------
  
  d_centered <-
    sweep(
      d_hat,
      2,
      colMeans(d_hat),
      "-"
    )
  
  
  Vd_hat <-
    crossprod(d_centered) / N
  
  
  ## ----------------------------------------------------------
  ## J. IF variance of theta_hat
  ##
  ## Var(theta_hat)
  ## =
  ## (1/N) tau^{-1} V(d) tau^{-T}
  ## ----------------------------------------------------------
  
  tau_inv <- solve(tau_hat)
  
  
  V_IF <-
    tau_inv %*%
    Vd_hat %*%
    t(tau_inv) / N
  
  
  SE_IF <-
    sqrt(
      pmax(diag(V_IF), 0)
    )
  
  
  if (is.null(parameter_names)) {
    parameter_names <-
      paste0(
        "theta",
        seq_along(theta_hat)
      )
  }
  
  
  table_IF <- data.frame(
    
    Variable = parameter_names,
    
    Estimate = theta_hat,
    
    Variance_IF = diag(V_IF),
    
    SE_IF = SE_IF,
    
    CI_Lower_IF =
      theta_hat -
      qnorm(0.975) * SE_IF,
    
    CI_Upper_IF =
      theta_hat +
      qnorm(0.975) * SE_IF
  )
  
  
  return(
    list(
      
      table = table_IF,
      
      V_IF = V_IF,
      
      SE_IF = SE_IF,
      
      d_hat = d_hat,
      
      gamma_hat = gamma_hat,
      
      kappa_hat = kappa_hat,
      
      tau_hat = tau_hat,
      
      Vd_hat = Vd_hat,
      
      omega_hat = omega_hat,
      
      pi_hat = pi_hat,
      
      S_hat = S_hat,
      
      h_hat = h_hat,
      
      parts = list(
        gammaS = part1,
        weighted_residual = part2,
        propensity_correction = part3
      )
    )
  )
}









run_one_split <- function(split_seed,damping) {
  
  tryCatch({
    
    ###########################################################################
    # Random 50% label deletion
    ###########################################################################
    
    set.seed(split_seed)
    
    
    labeled <- datX %>%
      filter(!is.na(Y)) %>%
      slice_sample(prop = 0.5)
    
    
    unlabeled <- datX %>%
      filter(is.na(Y))
    ###############################################################################
    # A7 CHECK 1: NHANES sample counts
    #
    # Excluded originally labeled observations are dropped entirely.
    # Expected analysis sample:
    #   1266 retained labeled
    # + 3697 originally unlabeled
    # = 4963 total
    ###############################################################################
    
    if (nrow(labeled) != 1266) {
      stop(
        paste0(
          "A7 count check failed: retained labeled n = ",
          nrow(labeled),
          ", expected 1266."
        )
      )
    }
    
    if (nrow(unlabeled) != 3697) {
      stop(
        paste0(
          "A7 count check failed: original unlabeled n = ",
          nrow(unlabeled),
          ", expected 3697."
        )
      )
    }
    
    if (nrow(labeled) + nrow(unlabeled) != 4963) {
      stop(
        paste0(
          "A7 total sample check failed: N = ",
          nrow(labeled) + nrow(unlabeled),
          ", expected 4963."
        )
      )
    }
    
    ###########################################################################
    # Encode data
    ###########################################################################
    
    labeled_final <- cbind(
      Y = labeled$Y,
      encode_data(labeled)
    )
    
    
    unlabeled_final <- encode_data(
      unlabeled
    )
    
    
    ###########################################################################
    # Make sure labeled/unlabeled X columns are identical
    ###########################################################################
    
    if (!identical(
      colnames(labeled_final)[-1],
      colnames(unlabeled_final)
    )) {
      
      stop(
        "Labeled and unlabeled encoded columns do not match."
      )
    }
    
    
    ###########################################################################
    # OLS on labeled subset
    #
    # Used as starting theta for ET/HD
    ###########################################################################
    
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
    
    ###############################################################################
    # SUPERVISED OLS
    ###############################################################################
    
    S1 <- summary(model.fit)$coefficients
    
    SUP_est <- as.numeric(S1[, "Estimate"])
    SUP_se  <- as.numeric(S1[, "Std. Error"])
    
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
    
    
    ###############################################################################
    # PSSE
    ###############################################################################
    
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
    
    
    ###############################################################################
    # DRESS
    ###############################################################################
    
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
    
    ###########################################################################
    # Build semi-supervised dataset
    ###########################################################################
    
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
    
    
    ###########################################################################
    # Initial propensity-score fit
    #
    # Same construction as your current script
    ###########################################################################
    
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
    
    
    ###########################################################################
    # ============================================================
    # ET
    # ============================================================
    ###########################################################################
    
    ET <- estimate_theta_EM_kfold_dual_ET(
      
      th =
        theta_start,
      
      data_full =
        data_full_real,
      
      K =
        3,
      
      # IMPORTANT:
      #
      # Make cross-fitting seed depend on repetition.
      seed =2025,
        
      
      max.iter =
        500,
      
      eps =
        1e-4,
      
      damping = damping
    )
    
    
    ###########################################################################
    # ET sandwich setup
    ###########################################################################
    
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
    
    
    ###########################################################################
    # U function for ET
    ###########################################################################
    
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
    
    
    ###########################################################################
    # b function for ET
    ###########################################################################
    
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
    
    
    ###########################################################################
    # Exact propensity design matrix used in glm
    ###########################################################################
    
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
    
    
    ###########################################################################
    # Joint sandwich for ET
    ###########################################################################
    
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
    
    
    ###########################################################################
    # Save ET quantities
    ###########################################################################
    
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
    
    
    ###########################################################################
    # ============================================================
    # HD
    # ============================================================
    ###########################################################################
    
    HD <- estimate_theta_EM_kfold_dual_HD(
      
      th =
        theta_start,
      
      data_full =
        data_full_real,
      
      K =
        3,
      
      seed =2025,
      
      max.iter =
        500,
      
      eps =
        1e-4,
      
      damping = damping
    )
    omega_hat<-HD$w
    # Check how close HD is to the dual-domain boundary eta = 0
    eta1 <- HD$eta[HD$data_all$D == 1]
    
    
    ###########################################################################
    # HD sandwich setup
    ###########################################################################
    
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
    
    # ================================================================
    # CE
    # ================================================================
    
    CE <- estimate_theta_EM_kfold_CE(
      th = as.numeric(coef(model.fit)),
      data_full = data_full_real,
      K = 3,
      
      # If you want ONLY label-split variability in A5,
      # keep this fixed at 2025 for every split.
      seed = 2025,
      
      max.iter = 500,
      eps = 1e-3,
      lambda_tol = 1e-8,
      lambda_maxit = 1000,
      damping = damping
    )
    
    
    # ---------------------------------------------------------------
    # CE data returned from estimator
    # ---------------------------------------------------------------
    
    data_all_CE <- CE$data_all
    
    
    # ---------------------------------------------------------------
    # Design matrix
    # ---------------------------------------------------------------
    
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
    
    
    # ---------------------------------------------------------------
    # Propensity model
    # ---------------------------------------------------------------
    
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
    
    # ---------------------------------------------------------------
    # U_i(theta)
    #
    # U_i(theta) = X_i(Y_i - X_i' theta)
    # ---------------------------------------------------------------
    
    U_fun_CE <- function(theta) {
      
      resid <- data_all_CE$Y -
        as.numeric(X_CE %*% theta)
      
      X_CE * resid
    }
    
    
    # ---------------------------------------------------------------
    # b_i(theta)
    #
    # b_i(theta) = X_i(yhat_i - X_i' theta)
    # ---------------------------------------------------------------
    
    b_fun_CE <- function(theta) {
      
      resid_hat <- data_all_CE$y.hat -
        as.numeric(X_CE %*% theta)
      
      X_CE * resid_hat
    }
    
    
    # ================================================================
    # CE sandwich variance
    # ================================================================
    
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
    ###########################################################################
    # Combine ET + HD results for this repetition
    ###########################################################################
    
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
    
    ###########################################################################
    # If one repetition fails, don't stop all 200.
    #
    # Record the failure and continue.
    ###########################################################################
    
    message(
      "Split ",
      split_seed,
      " failed: ",
      conditionMessage(e)
    )
    
    
    NULL
  })
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
      
      # Central difference if both perturbations are valid
      if (!is.null(f_plus) &&
          !is.null(f_minus) &&
          all(is.finite(f_plus)) &&
          all(is.finite(f_minus))) {
        
        J[, j] <-
          (f_plus - f_minus) / (2 * h)
        
        success <- TRUE
        break
      }
      
      # Forward difference
      if (!is.null(f_plus) &&
          all(is.finite(f_plus))) {
        
        J[, j] <-
          (f_plus - f0) / h
        
        success <- TRUE
        break
      }
      
      # Backward difference
      if (!is.null(f_minus) &&
          all(is.finite(f_minus))) {
        
        J[, j] <-
          (f0 - f_minus) / h
        
        success <- TRUE
        break
      }
      
      # Shrink perturbation if neither side is valid
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



###############################################################################
# Helper 1: effective sample size and max weight
###############################################################################

weight_diagnostics <- function(w) {
  
  w <- as.numeric(w)
  
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
  
  ESS <-
    sum(w)^2 /
    sum(w^2)
  
  list(
    ESS = ESS,
    max_weight = max(w)
  )
}


###############################################################################
# Helper 2: construct one method's simulation output
###############################################################################

make_simulation_rows <- function(
    rep_id,
    method,
    estimate,
    truth,
    se = NULL,
    ESS = NA_real_,
    max_weight = NA_real_,
    success = TRUE,
    parameter_names = NULL
) {
  
  estimate <- as.numeric(estimate)
  truth    <- as.numeric(truth)
  
  q <- length(truth)
  
  if (is.null(parameter_names)) {
    parameter_names <- paste0("theta", seq_len(q))
  }
  
  if (is.null(se)) {
    se <- rep(NA_real_, q)
  }
  
  se <- as.numeric(se)
  
  lower <-
    estimate -
    qnorm(0.975) * se
  
  upper <-
    estimate +
    qnorm(0.975) * se
  
  covered <-
    ifelse(
      is.finite(se),
      as.numeric(
        lower <= truth &
          truth <= upper
      ),
      NA_real_
    )
  
  data.frame(
    
    replication = rep_id,
    
    method = method,
    
    parameter = parameter_names,
    
    estimate = estimate,
    
    truth = truth,
    
    error =
      estimate - truth,
    
    squared_error =
      (estimate - truth)^2,
    
    analytic_se = se,
    
    ci_lower = lower,
    
    ci_upper = upper,
    
    ci_width =
      upper - lower,
    
    covered = covered,
    
    ESS = ESS,
    
    max_weight = max_weight,
    
    success =
      as.integer(success),
    
    stringsAsFactors = FALSE
  )
}


###############################################################################
# Wrapper: GEC estimator -> sandwich SE -> simulation output
#
# Works for ET / HD / CE / EL / SL
###############################################################################

process_gec_simulation <- function(
    fit_gec,
    entropy,
    rep_id,
    target_parameter,
    parameter_names
) {
  
  # -------------------------------------------------------------------------
  # Basic failure check
  # -------------------------------------------------------------------------
  
  if (is.null(fit_gec)) {
    return(NULL)
  }
  
  if (is.null(fit_gec$theta) ||
      is.null(fit_gec$lambda) ||
      is.null(fit_gec$data_all)) {
    
    message(
      entropy,
      " incomplete fit in replication ",
      rep_id
    )
    
    return(NULL)
  }
  
  
  # -------------------------------------------------------------------------
  # 1. Cross-fitted data
  # -------------------------------------------------------------------------
  
  data_all <- fit_gec$data_all
  
  
  # -------------------------------------------------------------------------
  # 2. Covariates
  # -------------------------------------------------------------------------
  
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
  
  
  # -------------------------------------------------------------------------
  # 3. Design matrix
  # -------------------------------------------------------------------------
  
  X <- model.matrix(
    ~ .,
    data = as.data.frame(x_cols)
  )
  
  D <- as.numeric(
    data_all$D
  )
  
  
  # -------------------------------------------------------------------------
  # 4. Observed estimating function U_i(theta)
  #
  # U_i(theta) = X_i {Y_i - X_i' theta}
  #
  # Only labeled observations contribute.
  # -------------------------------------------------------------------------
  
  U_fun <- function(theta) {
    
    U <- matrix(
      0,
      nrow = nrow(X),
      ncol = ncol(X)
    )
    
    id <- which(
      D == 1
    )
    
    resid <-
      data_all$Y[id] -
      as.numeric(
        X[
          id,
          ,
          drop = FALSE
        ] %*%
          theta
      )
    
    U[
      id,
    ] <-
      X[
        id,
        ,
        drop = FALSE
      ] *
      resid
    
    U
  }
  
  
  # -------------------------------------------------------------------------
  # 5. Cross-fitted calibration function
  #
  # b_i(theta) = X_i {yhat_i - X_i' theta}
  # -------------------------------------------------------------------------
  
  b_fun <- function(theta) {
    
    resid_hat <-
      data_all$y.hat -
      as.numeric(
        X %*%
          theta
      )
    
    X *
      resid_hat
  }
  
  
  # -------------------------------------------------------------------------
  # 6. Propensity model
  # -------------------------------------------------------------------------
  
  ps_fit <- glm(
    D ~ .,
    data = data.frame(
      D = D,
      x_cols
    ),
    family = binomial(),
    x = TRUE
  )
  
  
  phi_hat <-
    coef(
      ps_fit
    )
  
  O <-
    ps_fit$x
  
  
  # -------------------------------------------------------------------------
  # 7. Joint GEC sandwich
  # -------------------------------------------------------------------------
  
  var_fit <- tryCatch(
    
    gec_sandwich(
      
      phi_hat =
        phi_hat,
      
      lambda_hat =
        as.numeric(
          fit_gec$lambda
        ),
      
      theta_hat =
        as.numeric(
          fit_gec$theta
        ),
      
      D =
        D,
      
      O =
        O,
      
      U_fun =
        U_fun,
      
      b_fun =
        b_fun,
      
      entropy =
        entropy,
      
      parameter_names =
        colnames(X)
    ),
    
    error = function(e) {
      
      message(
        entropy,
        " sandwich failed in replication ",
        rep_id,
        ": ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  
  if (is.null(var_fit)) {
    return(NULL)
  }
  
  
  # -------------------------------------------------------------------------
  # 8. Point estimate
  # -------------------------------------------------------------------------
  
  est <-
    as.numeric(
      fit_gec$theta
    )
  
  
  # -------------------------------------------------------------------------
  # 9. Sandwich SE
  # -------------------------------------------------------------------------
  
  se <-
    sqrt(
      diag(
        var_fit$V_theta
      )
    )
  
  
  # -------------------------------------------------------------------------
  # 10. Weight diagnostics
  # -------------------------------------------------------------------------
  
  w <-
    as.numeric(
      fit_gec$w
    )
  
  w <-
    w[
      is.finite(w) &
        w > 0
    ]
  
  
  if (length(w) > 0) {
    
    ESS <-
      sum(w)^2 /
      sum(w^2)
    
    max_weight <-
      max(w)
    
  } else {
    
    ESS <-
      NA_real_
    
    max_weight <-
      NA_real_
  }
  
  
  # -------------------------------------------------------------------------
  # 11. Return simulation-format rows
  # -------------------------------------------------------------------------
  
  make_simulation_rows(
    
    rep_id =
      rep_id,
    
    method =
      entropy,
    
    estimate =
      est,
    
    truth =
      target_parameter,
    
    se =
      se,
    
    ESS =
      ESS,
    
    max_weight =
      max_weight,
    
    success =
      TRUE,
    
    parameter_names =
      parameter_names
  )
}




#-----------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------#
#
#
#                  Various semi-supervised estimation methods for M-estimation
#          (mainly for linear working model;logistic working model;quantile working model)
#
#-----------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------#


#################### (1) proposed method (projection-based semi-supervised estimator): PSSE################# 
PSSE <- function(labelled_data,unlabelled_data,c1=NULL,
                 type="linear",alpha,gamma=10,Kfolds=5)
{
  #-----------------------------------------Arguments-----------------------------------------------------#
  # Purpose: This function is to implement our proposed semi-supervised method. 
  #           This function is the main function, and it also involves several sub-functions defined in
  #           (1.1) - (1.6) below.
  #
  # Input: 
  #      labelled_data: A matrix, whose each row is an observation of predictor vector and the response
  #                      variable and its first column is the observation vector of the response variable.  
  #      unlabelled_data: A matrix, whose each row is an observation of predictor vector. 
  #      c1: Weight. It is the weight to balance the contribution of labelled data and unlabelled data.
  #           Default value is nrow(labelled_data)/(nrow(labelled_data)+nrow(unlabelled_data)).
  #      type: A specification for the type of loss function $L$. We offer three choices: 
  #             "linear", "logistic", "quantile". Default value is "linear". 
  #      tau: The quantile level, which is only useful for "type=quantile". Default value is 0.5.
  #      alpha: The polynomial order. If it is NULL, then we use the data-driven selector GBIC_ppo to 
  #              determine the polynomial order. Default is NULL. 
  #      gamma: The maximum of possible polynomial order. Default is 10. Only if alpha is NULL, gamma will 
  #              be usful for the selection of polynomial order in the construction of Z. 
  #      sd: Return the pointwise standard deviation estimate based on the available samples? Default 
  #           value is FALSE.  
  #      Kfolds: K-folds cross validation method to avoiding the over-fitting during estimating variance, 
  #          which is only useful When "sd=TRUE".
  #
  # Output: 
  #       Hattheta: The estimate for the target parameter. 
  #       sd.of.hattheta: If "sd=TRUE", the elementwise standard deviation estimate for the estimator; 
  #                       otherwise, NULL.
  #       alpha: The selected polynomial order. 
  #------------------------------------------------------------------------------------------------------#
  
  # some basic parameters 
  n=nrow(labelled_data)
  N=nrow(unlabelled_data)
  p=ncol(labelled_data)-1
  

  if (p!=ncol(unlabelled_data)){
    print("the dimension of the labelled and unlabelled data is unmatched")
    
    return(NULL)
  }
  
  if(is.null(c1)==TRUE){
    c1=n/(n+N) ###For semi-paramteric one
    #c1=1 ####For supervised one 
  }
  
  # whether the polynomial order needs to be selected or not 
  if (is.null(alpha)==TRUE){ 
    
    # data splitting 
    labelled_data_part1<- labelled_data[1:round(n/2),]
    labelled_data_part2<- labelled_data[(round(n/2)+1):n,]
    
    
    # the supervised estimate
    hattheta_supervised_part1<- supervised_one(labelled_data_part1)
    hattheta_supervised_part2<- supervised_one(labelled_data_part2)
    
      
    # the first derivative of the loss function L using the supervised estimate 
    L_first_derivative_part1<- L_first_derivative(labelled_data_part1,hattheta_supervised_part2)
    L_first_derivative_part2<- L_first_derivative(labelled_data_part2,hattheta_supervised_part1)
    
    
    # new data integrating the first derivatives of loss function
    Ynew<- rbind(L_first_derivative_part1,L_first_derivative_part2)
    covariates_labelled<- labelled_data[,-1]
    
    
    # GBIC_ppo to select the polynomial order 
    GBICscrores<-apply(as.matrix(1:gamma,1,gamma), 1, function(t) GBIC_ppo(Ynew,covariates_labelled,t))
    alpha<- which.min(GBICscrores)
    
  }
  
  # determine Z
  labelled_Z <- polynomial(labelled_data[,-1],alpha)
  unlabelled_Z <- polynomial(unlabelled_data,alpha)
  
  
  
  # weights (w_i) using in the our optimization problem to get the proposed semi-supervised estimate 
  unlabelled_Z_mean <- colMeans(unlabelled_Z)
  labelled_Z_seondmoment <- t(labelled_Z)%*%labelled_Z/n
  weights_loss <- as.vector(c1+(1-c1)*t(unlabelled_Z_mean)%*%solve(labelled_Z_seondmoment)%*%t(labelled_Z))
####This function needs to be changed
  
  
  
  
 
    hattheta_PSSE <- as.vector(solve(t(cbind(rep(1,n),labelled_data[,-1]))%*%
                                       (weights_loss*cbind(rep(1,n),labelled_data[,-1])))%*%t(cbind(rep(1,n),labelled_data[,-1]))%*%
                                 (weights_loss*labelled_data[,1]))

  hattheta_PSSE = list("Hattheta"=hattheta_PSSE)
  
  
  return(hattheta_PSSE)
}


###############(1.1) polynomial generation function 

polynomial<- function(data,order) 
{
  #---------------------------Arguments----------------------------------------#
  # Purpose: This function is to produce the polynomial basis of X including 
  #           the intercept vector. 
  #
  # Input: 
  #       data: A matrix, whose each row is an observation of predictor vector.
  #       order: The polynomial order.
  #----------------------------------------------------------------------------#
  polynomial_Z<- rep(1,nrow(data))
  for (i in 1:order)
  {
    polynomial_i <- data^i
    polynomial_Z <- cbind(polynomial_Z,polynomial_i)
  }
  return(polynomial_Z)
}




###############(1.2) the data driven selector GBIC_ppo  
GBIC_ppo <- function (Y_new,covariates_matrix,order_poly){
  #-----------------------Arguments--------------------------------------------------# 
  # Purpose: GBIC_ppo criterion function 
  #
  # Input:
  #       Y_new: The first derivatives of the loss function L.
  #       covariates_matrix: A matrix, each row is an observation of predictor vector.
  #       order_poly: The polynomial order.
  #
  #---------------------------------------------------------------------------------#
  d <- ncol(Y_new)
  n <- nrow(Y_new)
  p <- ncol(covariates_matrix)
  
  polynomial_matrix <- polynomial(covariates_matrix, order_poly)
  
  gammahat <- lm(Y_new~polynomial_matrix-1)$coefficients
  residual_square <- (Y_new - polynomial_matrix %*% gammahat)^2
  sigmahat_square <- colSums(residual_square)/(n-p*order_poly-1)
  design_matrix <- t(polynomial_matrix)%*%polynomial_matrix/n
  
  
  if(class(try(solve(design_matrix),silent=T))[1]=="try-error"){
    
    GBIC<- 9999999
    
  }else {
    trace_det_AB = numeric()
    for (j in 1:d){
      Bhat_j<- t(polynomial_matrix*residual_square[,j])%*%polynomial_matrix/n
      cova_constrast <- 1/sigmahat_square[j]*solve(design_matrix)%*%Bhat_j
      trace_AB <- sum(diag(cova_constrast))
      det_AB <- det(cova_constrast)
      trace_det_AB=c(trace_det_AB,trace_AB-log(det_AB))
    }
    GBIC = 1/n*(d*(n-p*order_poly-1)+ n*sum(log(sigmahat_square)) + d*log(n)*p*order_poly+sum(trace_det_AB))
  }
  
  return(GBIC)
  
}


#################(1.3) the supervised estimate based on the given samples
supervised_one<- function(labelled_data){
  #---------------------Arguments-------------------------------------------#
  # Purpose: According the type of loss function, this function is used to  
  #         get the estimate for the target parameter only based on the given 
  #         labelled samples.
  #------------------------------------------------------------------------#
 
    hattheta_supervised <- lm(labelled_data[,1]~labelled_data[,-1])$coefficients
  
  return(hattheta_supervised)
}




#################(1.4)the first derivative of loss function L
L_first_derivative<- function(labelled_data,hattheta_supervised)
{
  #---------------------Arguments-------------------------------------------#
  # Purpose: According the type of loss function, this function is used to  
  #         determine the first-order derivatives of the loss function based on 
  #         an estimate of the target parameter and a labelled data set.
  #------------------------------------------------------------------------#
  n=nrow(labelled_data)
  
    
    residuals<- as.vector(labelled_data[,1]-cbind(rep(1,n),labelled_data[,-1])%*%hattheta_supervised)
    L_first_derivative <- residuals*cbind(rep(1,n),labelled_data[,-1])
    
  
  return(L_first_derivative)
}



###########################(2) PI proposed by Azriel et al. (2021) #################################
PI <- function(labelled_data,unlabelled_data)
{
  #---------------------------------Arguments------------------------------------------#
  # Purpose: This function is to implement the method proposed by Azriel et al.(2021)
  #          for linear working model. 
  #
  # Input:
  #       labelled_data: Same as in the function "PSSE". 
  #       unlabelled_data: Same as in the function "PSSE".     
  #
  # Output: 
  #        Hattheta: The estimate for the target parameter. 
  #------------------------------------------------------------------------------------#
  n=nrow(labelled_data)
  p=ncol(labelled_data)-1
  N=nrow(unlabelled_data)
  
  # combine all the covairates 
  X_combined <- rbind(labelled_data[,-1],unlabelled_data)
  X_labelled <- labelled_data[,-1]
  
  hat_beta_initial <- numeric()
  X_dot <- matrix(rep(0,n*p),n,p) 
  delta_tilde <- matrix(rep(0,n*p),n,p)
  
  for (j in 1:p)
  {
    ## first step 
    coefficients_negtive_j = lm(X_combined[,j]~X_combined[,-j])$coefficients #unlabeled data
    X_j_dot = X_labelled[,j] - cbind(rep(1,n),X_labelled[,-j]) %*% coefficients_negtive_j #labeled data projection error
    X_j_dot_total = X_combined[,j] - cbind(rep(1,n+N),X_combined[,-j]) %*% coefficients_negtive_j
    X_j_dot_square_sampleaverage=mean(X_j_dot_total^2)
    
    
    ## step step 
    W_j = as.vector(labelled_data[,1])*X_j_dot/X_j_dot_square_sampleaverage
    U1 = X_j_dot/X_j_dot_square_sampleaverage
    X_dot[,j] = as.vector(U1) 
    U = matrix(rep(0,n*p),n,p) 
    for (jj in 1:p)
    {
      if (jj==j)
      {
        U[,jj] = X_labelled[,jj]*X_j_dot/X_j_dot_square_sampleaverage-1
      }
      
      else
      {
        U[,jj] = X_labelled[,jj]*X_j_dot/X_j_dot_square_sampleaverage
      }
      
    }
    
    delta_tilde[,j] = W_j-as.matrix(cbind(rep(1,n),U1,U)) %*% (lm(W_j~U1+U)$coefficients)
    hat_beta_j = lm(W_j~U1+U)$coefficients[1]
    hat_beta_initial = c(hat_beta_initial,hat_beta_j)
  }
  hat_alpha = mean(labelled_data[,1])-t(hat_beta_initial)%*%colMeans(X_labelled)
  hat_theta_Azriel = as.vector(c(hat_alpha,hat_beta_initial))
  
  hat_theta_Azriel = list("Hattheta"=hat_theta_Azriel)
  

  
  return(hat_theta_Azriel)
}
#PI(data_labelled,data_unlabelled)



##################(3) EASE proposed by Chakrabortty and Cai (2018) ##########################
EASE<-function(labelled_data,unlabelled_data,K,H,r)
{
  #---------------------------------Arguments-------------------------------------------------#
  # Purpose: This function is to implement the method proposed by Chakrabortty and Cai (2018) 
  #          for linear working model. This is the main function, and it also involves a 
  #          subfunction defined in (3.1). 
  #   
  # Input: 
  #       labelled_data: Same as in the function "PSSE". 
  #       unlabelled_data: Same as in the function "PSSE".
  #       K: K-folds for sliced inverse regression.
  #       H: H slices in sliced inverse regression.
  #       r: The reduced dimension. 
  #
  # Output: 
  #        Hattheta: The estimate for the target parameter. 
  #------------------------------------------------------------------------------------------#
  
  n=nrow(labelled_data)
  p=ncol(labelled_data)-1
  N=nrow(unlabelled_data)
  
  X_combined<- rbind(labelled_data[,-1],unlabelled_data)
  hat_Gamma_semisupervised= t(cbind(rep(1,n+N),X_combined))%*%cbind(rep(1,n+N),X_combined)/(n+N)
  coef_supervised = lm(labelled_data[,1]~labelled_data[,-1])$coefficients
  set.seed(20210907)
  index=createFolds(1:n, k = K) # data splitting
  
  
  #####estimate phi_k(k=1,...,K) and parameters
  hat_nonpara_test_vector <- numeric()
  hat_nonpara_unlabelled_matrix <- vector()
  covariates_test = vector()
  response_test = numeric()
  
  for (k in 1:K)
  {
    
    index_k=as.vector(index[[k]])
    Yt_train = labelled_data[-index_k,1]
    Xt_train = labelled_data[-index_k,-1]
    Yt_test = labelled_data[index_k,1]
    Xt_test = labelled_data[index_k,-1]
    X_combined_train = rbind(Xt_train, unlabelled_data)
    nrow_Xt_train = n-length(index_k)
    nrow_Xt_test = length(index_k)
    
    covariates_test = rbind(covariates_test,Xt_test) #save the covariates for test 
    response_test <- c(response_test,Yt_test)
    
    
    ###SS-STR
    
    ##step(0):estimate mean of covariance of predictors X and standardized 
    X_combined_train_mean = colMeans(X_combined_train)
    X_train_covariance = cov(X_combined_train)
    eigen_X_train_covariance = eigen(X_train_covariance)
    lamda_sqrt <- diag(sqrt(eigen_X_train_covariance$values))
    Sigma_sqrt <- (eigen_X_train_covariance$vectors)%*%(lamda_sqrt*solve(eigen_X_train_covariance$vectors))
    Sigma_negative_sqrt <- solve(Sigma_sqrt)  # calculate the Sigma^{-1/2}
    X_train_standard = t(Sigma_negative_sqrt %*% t(X_combined_train - X_combined_train_mean))
    Xt_train_standard = X_train_standard[1:nrow_Xt_train,]
    Xv_standard = X_train_standard[(nrow_Xt_train+1):(nrow_Xt_train+N),]
    
    
    ##step(i)
    range_Yt_train = range(Yt_train)
    intrval_length = (range_Yt_train[2]-range_Yt_train[1])/H
    indicator_Yt_train = floor(Yt_train/intrval_length)+1 # the interval indicator
    group_with_train = sort(unique(indicator_Yt_train))
    prob_h_k_train = prop.table(table(group_with_train))
    
    
    ##step(ii)
    distance_Xv_Xt_train = matrix(rep(0,N*nrow_Xt_train),N,nrow_Xt_train)
    for(i in 1:N)
    {
      
      for(j in 1:nrow_Xt_train)
      {
        distance_Xv_Xt_train[i,j] = sum((Xv_standard[i,]-Xt_train_standard[j,])^2)
      }
      
    }
    index_Xv_imputed = apply(distance_Xv_Xt_train,1,which.min)
    #index_Xv_tmp_imputed
    indicator_Xv_imputed = indicator_Yt_train[index_Xv_imputed]
    #indicator_Xv_tmp_imputed
    indicator_X_train = c(indicator_Yt_train,indicator_Xv_imputed)
    
    
    ##step(iii)
    M_matrix_train_semi = matrix(rep(0,p*p),p,p)
    for(kk in 1: length(group_with_train))
    {         
      index_train_standard_kk<- which(indicator_X_train==group_with_train[kk])
      if(length(index_train_standard_kk)==1)
      {
        M_bar_sampleaverage = X_train_standard[index_train_standard_kk,]
        M_matrix_train_semi = M_matrix_train_semi + prob_h_k_train[kk]*M_bar_sampleaverage%*%t(M_bar_sampleaverage)
      }
      if (length(index_train_standard_kk)>1)
      {
        M_bar_sampleaverage = colMeans(X_train_standard[index_train_standard_kk,])
        M_matrix_train_semi = M_matrix_train_semi + prob_h_k_train[kk]*M_bar_sampleaverage%*%t(M_bar_sampleaverage)
      }
      
    }
    
    hat_pstar_h_k = eigen(M_matrix_train_semi)$vectors[,1:r]
    hat_Zero_pstar_h_k = Sigma_negative_sqrt %*% hat_pstar_h_k
    
    Xt_train_SIR = Xt_train %*% hat_Zero_pstar_h_k
    Xt_test_SIR = Xt_test %*% hat_Zero_pstar_h_k
    Xt_unlabelled_SIR = unlabelled_data%*%hat_Zero_pstar_h_k
    
    
    ###calculate the nonparametric estimator in test and unlabelled data
    
    h_opt_chark = optimize(bandwidth.cv,c(0.1,2),Y=Yt_train,X=Xt_train_SIR,d=(K-1),order=r)$minimum
    
    
    hat_nonpara_test = rep(0,nrow_Xt_test)
    for (kkk in 1:nrow_Xt_test)
    {
      kernel = apply(Xt_train_SIR,1, function(t){kernel_guassian(sqrt(sum((Xt_test_SIR[kkk,]-t)^2))/h_opt_chark^r)}) 
      kernel_mean = sum(kernel)
      kernel_Y_mean = t(kernel)%*%Yt_train
      hat_nonpara_test[kkk] = kernel_Y_mean/(kernel_mean+10^(-6))
      
    }
    #hat_nonpara_test  #the nonparametric estimator in test data
    
    hat_nonpara_test_vector <- c(hat_nonpara_test_vector,hat_nonpara_test)
    
    
    
    hat_nonpara_unlabelled = rep(0,N)
    for (kkkk in 1:N)
    {
      kernel = apply(Xt_train_SIR,1, function(t){kernel_guassian(sqrt(sum((Xt_unlabelled_SIR[kkkk,]-t)^2))/h_opt_chark^r)}) 
      kernel_mean = sum(kernel)
      kernel_Y_mean = t(kernel)%*%Yt_train
      hat_nonpara_unlabelled[kkkk] = kernel_Y_mean/(kernel_mean+10^(-6))
      
    }
    #hat_nonpara_unlabelled #the nonparametric estimator in unlabelled data
    
    hat_nonpara_unlabelled_matrix<- cbind(hat_nonpara_unlabelled_matrix,hat_nonpara_unlabelled)
    
  }
  
  
  
  ###estimate eta
  
  hat_eta_k_semisupervised = lm((response_test-hat_nonpara_test_vector)~covariates_test)$coefficients
  #hat_eta_k_semisupervised
  
  ###estimate mu for the unlabelled data
  
  hat_mu_unlabelled = rowMeans(hat_nonpara_unlabelled_matrix)+
    as.vector(cbind(rep(1,N),unlabelled_data) %*% hat_eta_k_semisupervised)
  
  
  ###(SNP estimator) estimate theta by unlabelled data
  
  hattheta_K <- lm(hat_mu_unlabelled~unlabelled_data)$coefficients
  
  
  
  ### EASE estimator
  
  
  ##step(i) estimate hat_eta_k, hat_mu_k,hat_phi_k,hat_phi_0
  
  length_intial <- 0
  hat_Gamma_semisupervised_inverse <- solve(hat_Gamma_semisupervised)
  
  hat_phi_k_cv_matrix<- vector()
  hat_phi_0_cv_matrix <- vector()
  for (k in 1:K)
  {
    length_k = length(index[[k]])
    length_intial = length_intial+length_k
    index_test_k = (1+length_intial-length_k):length_intial
    index_exclude_k = (1:n)[-index_test_k]
    hat_eta_k_cv <- lm((response_test[index_exclude_k]-hat_nonpara_test_vector[index_exclude_k])~
                         covariates_test[index_exclude_k,])$coefficients
    
    hat_mu_k_cv <- hat_nonpara_test_vector[index_test_k]+
      as.vector(cbind(rep(1,length_k),covariates_test[index_test_k,])%*%hat_eta_k_cv)
    
    hat_phi_k_cv <- as.vector(response_test[index_test_k]-hat_mu_k_cv)*cbind(rep(1,length_k),
                                                                             covariates_test[index_test_k,])%*% hat_Gamma_semisupervised_inverse
    
    hat_phi_0_cv <- as.vector(response_test[index_test_k]-cbind(rep(1,length_k),
                                                                covariates_test[index_test_k,])%*%coef_supervised)*cbind(rep(1,length_k),
                                                                                                                         covariates_test[index_test_k,])%*% hat_Gamma_semisupervised_inverse
    
    hat_phi_k_cv_matrix <- rbind(hat_phi_k_cv_matrix,hat_phi_k_cv)
    hat_phi_0_cv_matrix <- rbind(hat_phi_0_cv_matrix,hat_phi_0_cv)
    
  }
  
  ##step(ii) calculate delta_lk
  
  hat_delta_12 = colMeans(hat_phi_0_cv_matrix*(hat_phi_k_cv_matrix-hat_phi_0_cv_matrix))
  hat_delta_22 =colMeans((hat_phi_k_cv_matrix-hat_phi_0_cv_matrix)^2)
  epsilon_n = 1/n^(3/8)
  delta_lk = hat_delta_12/(hat_delta_22+epsilon_n)
  
  
  
  ##step(iii) calculate EASE estimator
  
  hat_theta_semi_Chark<- coef_supervised+delta_lk*(hattheta_K-coef_supervised)
  
  
  
  hat_theta_semi_Chark = list("Hattheta"=hat_theta_semi_Chark)
  
  

  
  return(hat_theta_semi_Chark)
}
#EASE(data_labelled,data_unlabelled,K=5,H=40,r=2)[[1]]



#############(3.1) bandwidth selection
kernel_guassian <- function(a) #guassian kernel
{
  1/sqrt(2*pi)*exp(-a^2/2)
}



bandwidth.cv <- function(Y,X,d,h,order) #bandwidth selection
{
  #------------------------Arguments---------------------------------#
  # Purpose: This function is used for the bandwidth selection
  #          in the function "EASE". After SIR, we need to estimate
  #          the $E[Y|P_rX]$ using kernel smoothing. 
  #
  # Input:
  #       Y: Response observation vector.
  #       X: Covariate matrix.
  #       d: d-folds.
  #       h: bandwidth
  #       order: The power.
  #------------------------------------------------------------------#
  
  n = length(Y)
  p = ncol(X)
  
  set.seed(20220907)
  index=createFolds(1:n, k = d)
  
  err_test=0
  for (i in 1:d)
  {
    index_i=as.vector(index[[i]])
    Yt_train = Y[-index_i]
    Xt_train = X[-index_i,]
    Yt_test = Y[index_i]
    Xt_test = X[index_i,]
    n_positive=length(index_i)
    
    
    for(j in 1:n_positive)
    {
      kernel = apply(Xt_train,1, function(t){kernel_guassian(sqrt(sum((Xt_test[j,]-t)^2))/h^order)}) 
      kernel_mean = sum(kernel)
      kernel_Y_mean = t(kernel)%*%Yt_train
      bias_nonpara_estimator = Yt_test[j]-kernel_Y_mean/kernel_mean
      err_test = err_test + bias_nonpara_estimator^2
    }
    
  } 
  return(err_test)
  
}




#################(4) DRESS proposed by Kawakita and Kanamori (2013) ########################

DRESS<- function(labelled_data,unlabelled_data,type="linear",tau=0.5,L,sd,Kfolds=5)
{
  #-------------------------Arguments-------------------------------------------#
  # Purpose: This function is to implement the method proposed by
  #         Kawakita and Kanamori (2013) for M-estimation,
  #         which uses the desity-ratio to improve the estimation
  #         efficiency.
  #
  # Input:
  #       labelled_data: Same as in the function "PSSE". 
  #       unlabelled_data: Same as in the function "PSSE".     
  #       type: Same as in the function "PSSE".
  #       tau: Same as in the function "PSSE".
  #        L: The polynomial order when we construct the polynomial function as
  #            base function for density ratio estimation.
  #       sd: Return the pointwise standard deviation estimate based on the available
  #           samples? Default value is FALSE.
  #       Kfolds: Same as in the function "PSSE". 
  #
  # Output: 
  #       Hattheta: The estimate for the target parameter. 
  #       sd.of.hattheta: If "sd=TRUE", the elementwise standard deviation 
  #                        estimate for the estimator; otherwise, NULL. 
  #       error: indicator of whether the parameters in density-ratio is properly
  #              estimated. "TRUE" represents "NOT"; "FALSE" represents "YES".
  #-----------------------------------------------------------------------------#
  
  n=nrow(labelled_data)
  p=ncol(labelled_data)-1
  N=nrow(unlabelled_data)
  
  
  base_labelled<- polynomial(labelled_data[,-1],L)
  base_unlabelled <- polynomial(unlabelled_data,L)
  alpha_first_derivative_unlabelled <- colMeans(base_unlabelled)
  
  ## estimating alpha (the parameter in density-ratio)
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
    #print(i)
  }


  exponential_phi_labelled <- as.vector(exp(base_labelled%*%alpha_initial))
  
  
  ## estimating theta*(the target parameter)
  if (type=="linear"){
    hattheta_DRESS<- lm(labelled_data[,1]~labelled_data[,-1],weights=exponential_phi_labelled)$coefficients
    
  }else if (type=="logistic"){
    hattheta_DRESS <- glm(labelled_data[,1]~labelled_data[,-1],family="quasibinomial",
                          weights=exponential_phi_labelled)$coefficients
  }
  else if (type=="quantile"){
    hattheta_DRESS <- rq(labelled_data[,1]~labelled_data[,-1],tau=tau,weights=exponential_phi_labelled)$coefficients
  }
  
  
  hattheta_DRESS = append(list("Hattheta"=hattheta_DRESS),list("error"=error_svd))
  
  
  ### estimating the sd
  if(sd==TRUE){
    
    if(type=="linear"){
      
      
      X_secondmoment_inverse=solve(t(cbind(rep(1,N),unlabelled_data))%*%cbind(rep(1,N),unlabelled_data)/N) 
      
      set.seed(20218080)
      index=createFolds(1:n, k = Kfolds) # data splitting
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
        
        
        
        # projection matrix A(theta) estimation 
        
        hattheta_supervised_train_cv <-  hattheta_DRESS[[1]]#lm(Yt_train~Xt_train)$coefficients
        L_firstder_train_cv <- as.vector(Yt_train-cbind(rep(1,nrow_Xt_train),Xt_train)%*%hattheta_supervised_train_cv)*
          cbind(rep(1,nrow_Xt_train),Xt_train)
        L_firstder_projection_cof <- solve(t(Zt_train)%*%Zt_train/nrow_Xt_train)%*%
          t(Zt_train)%*%L_firstder_train_cv/nrow_Xt_train
        
        
        # estimate variance by test data
        
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
    }
    
    if(type=="logistic"){
      
      exp_linear_combined <- as.vector(cbind(rep(1,N),unlabelled_data)%*%hattheta_DRESS[[1]])
      weights_second_derivative <- 1/(1+exp(-exp_linear_combined))^2*exp(-exp_linear_combined)
      X_secondmoment_inverse=solve(t(weights_second_derivative*cbind(rep(1,N),unlabelled_data))%*%cbind(rep(1,N),unlabelled_data)/N) 
      
      set.seed(20218080)
      index=createFolds(1:n, k = Kfolds) # data splitting
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
        
        
        
        # projection matrix A(theta) estimation 
        
        hattheta_supervised_train_cv <- hattheta_DRESS[[1]]#lm(Yt_train~Xt_train)$coefficients
        exp_linear_combined_train <- as.vector(cbind(rep(1,nrow_Xt_train),Xt_train)%*%hattheta_supervised_train_cv)
        L_firstder_train_cv <- as.vector(1/(1+exp(-exp_linear_combined_train))-Yt_train)*
          cbind(rep(1,nrow_Xt_train),Xt_train)
        L_firstder_projection_cof <- solve(t(Zt_train)%*%Zt_train/nrow_Xt_train)%*%
          t(Zt_train)%*%L_firstder_train_cv/nrow_Xt_train
        
        
        # estimate variance by test data
        
        exp_linear_combined_test <- as.vector(cbind(rep(1,nrow_Xt_test),Xt_test)%*%hattheta_supervised_train_cv)
        L_firstder_test_cv <- as.vector(1/(1+exp(-exp_linear_combined_test))-Yt_test)*
          cbind(rep(1,nrow_Xt_test),Xt_test)
        
        W1_test_k <- L_firstder_test_cv-Zt_test%*%L_firstder_projection_cof
        W1_test<- rbind(W1_test,W1_test_k)
        
      }
      
      
      exp_linear_combined_total <- as.vector(cbind(rep(1,n),labelled_data[,-1])%*%hattheta_DRESS[[1]])
      L_firstder_total <- (1/(1+exp(-exp_linear_combined_total))-labelled_data[,1])*
        (cbind(rep(1,n),labelled_data[,-1]))
      L_firstder_projection_cof_total <- solve(t(base_labelled)%*%base_labelled/n)%*%t(base_labelled)%*%L_firstder_total/n
      W2_total <-  base_unlabelled%*%L_firstder_projection_cof_total
      
      W1_covariance <- t(W1_test)%*%W1_test/n
      W2_covariance <- t(W2_total)%*%W2_total/N
      Vc_hat_semi = W1_covariance+(n/N)*W2_covariance
      
      Var_matrix_hat = X_secondmoment_inverse%*%Vc_hat_semi%*%X_secondmoment_inverse 
      sd_hattheta_DRESS = sqrt(diag(Var_matrix_hat/n))
      
    }
    
    if(type=="quantile"){
      
      
      # estimating $V_c$ by K-folds CV
      
      set.seed(20218080)
      index=createFolds(1:n, k = Kfolds) # data splitting
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
        
        
        residual_train_indicator <- as.vector(Yt_train-cbind(rep(1,nrow_Xt_train),Xt_train)%*%hattheta_DRESS[[1]])<=0
        L_firstder_train_cv <- (residual_train_indicator-tau)*cbind(rep(1,nrow_Xt_train),Xt_train)
        L_firstder_projection_cof <- solve(t(Zt_train)%*%Zt_train)%*%
          t(Zt_train)%*%L_firstder_train_cv
        
        
        # estimate variance by test data
        
        residual_test_indicator <- as.vector(Yt_test-cbind(rep(1,nrow_Xt_test),Xt_test)%*%hattheta_DRESS[[1]])<=0
        L_firstder_test_cv <- (residual_test_indicator-tau)*cbind(rep(1,nrow_Xt_test),Xt_test)
        
        W1_test_k <- L_firstder_test_cv-Zt_test%*%L_firstder_projection_cof
        W1_test<- rbind(W1_test,W1_test_k)
        
      }
      
      
      residual_total_indicator <- as.vector(labelled_data[,1]-cbind(rep(1,n),labelled_data[,-1])%*%hattheta_DRESS[[1]])<=0
      L_firstder_total <- (residual_total_indicator-tau)*(cbind(rep(1,n),labelled_data[,-1]))
      L_firstder_projection_cof_total <- solve(t(base_labelled)%*%base_labelled)%*%t(base_labelled)%*%L_firstder_total
      W2_total <-  base_unlabelled%*%L_firstder_projection_cof_total
      
      W1_covariance <- t(W1_test)%*%W1_test/n
      W2_covariance <- t(W2_total)%*%W2_total/N
      Vc_hat_semi = W1_covariance+(n/N)*W2_covariance
      
      
      
      ##estimating the second derivatives (M)
      B=2000
      G=mvrnorm(B,rep(0,(p+1)),diag(rep(1,(p+1))))
      theta_check_semi <- apply(1/sqrt(n)*G,1,function(t) hattheta_DRESS[[1]]+t)
      residual_semi_indicator= apply(cbind(rep(1,n),labelled_data[,-1])%*%theta_check_semi,2, function(t) as.vector(labelled_data[,1]-t)<=0)
      U_check_semi <- t(residual_semi_indicator-tau)%*%cbind(rep(1,n),labelled_data[,-1])/sqrt(n)
      
      hat_M_semi <- t(apply(U_check_semi,2, function(t) lm(t~G-1)$coefficients))
      #X_secondmoment_inverse <- solve((hat_M_semi+t(hat_M_semi))/2)
      X_secondmoment_inverse<- solve(hat_M_semi)
      
      
      
      Var_matrix_hat = X_secondmoment_inverse%*%Vc_hat_semi%*%t(X_secondmoment_inverse)
      sd_hattheta_DRESS = sqrt(diag(Var_matrix_hat)/n)
      
      
      
      
    }
    
    hattheta_DRESS <- append(hattheta_DRESS,list("sd.of.hattheta"=sd_hattheta_DRESS))
  }
  

  
  return(hattheta_DRESS)
}






PSSE1 <- function(labelled_data,unlabelled_data,c1=NULL,
                 type="linear",tau=0.5,alpha=NULL,gamma=10,sd=FALSE,Kfolds=5)
{
  #-----------------------------------------Arguments-----------------------------------------------------#
  # Purpose: This function is to implement our proposed semi-supervised method. 
  #           This function is the main function, and it also involves several sub-functions defined in
  #           (1.1) - (1.6) below.
  #
  # Input: 
  #      labelled_data: A matrix, whose each row is an observation of predictor vector and the response
  #                      variable and its first column is the observation vector of the response variable.  
  #      unlabelled_data: A matrix, whose each row is an observation of predictor vector. 
  #      c1: Weight. It is the weight to balance the contribution of labelled data and unlabelled data.
  #           Default value is nrow(labelled_data)/(nrow(labelled_data)+nrow(unlabelled_data)).
  #      type: A specification for the type of loss function $L$. We offer three choices: 
  #             "linear", "logistic", "quantile". Default value is "linear". 
  #      tau: The quantile level, which is only useful for "type=quantile". Default value is 0.5.
  #      alpha: The polynomial order. If it is NULL, then we use the data-driven selector GBIC_ppo to 
  #              determine the polynomial order. Default is NULL. 
  #      gamma: The maximum of possible polynomial order. Default is 10. Only if alpha is NULL, gamma will 
  #              be usful for the selection of polynomial order in the construction of Z. 
  #      sd: Return the pointwise standard deviation estimate based on the available samples? Default 
  #           value is FALSE.  
  #      Kfolds: K-folds cross validation method to avoiding the over-fitting during estimating variance, 
  #          which is only useful When "sd=TRUE".
  #
  # Output: 
  #       Hattheta: The estimate for the target parameter. 
  #       sd.of.hattheta: If "sd=TRUE", the elementwise standard deviation estimate for the estimator; 
  #                       otherwise, NULL.
  #       alpha: The selected polynomial order. 
  #------------------------------------------------------------------------------------------------------#
  
  # some basic parameters 
  n=nrow(labelled_data)
  N=nrow(unlabelled_data)
  p=ncol(labelled_data)-1
  
  
  if (p!=ncol(unlabelled_data)){
    print("the dimension of the labelled and unlabelled data is unmatched")
    
    return(NULL)
  }
  
  if(is.null(c1)==TRUE){
    c1=n/(n+N)
  }
  
  # whether the polynomial order needs to be selected or not 
  if (is.null(alpha)==TRUE){ 
    
    # data splitting 
    labelled_data_part1<- labelled_data[1:round(n/2),]
    labelled_data_part2<- labelled_data[(round(n/2)+1):n,]
    
    
    # the supervised estimate
    hattheta_supervised_part1<- supervised_one(labelled_data_part1)
    hattheta_supervised_part2<- supervised_one(labelled_data_part2)
    
    
    # the first derivative of the loss function L using the supervised estimate 
    L_first_derivative_part1<- L_first_derivative(labelled_data_part1,hattheta_supervised_part2)
    L_first_derivative_part2<- L_first_derivative(labelled_data_part2,hattheta_supervised_part1)
    
    
    # new data integrating the first derivatives of loss function
    Ynew<- rbind(L_first_derivative_part1,L_first_derivative_part2)
    covariates_labelled<- labelled_data[,-1]
    
    
    # GBIC_ppo to select the polynomial order 
    GBICscrores<-apply(as.matrix(1:gamma,1,gamma), 1, function(t) GBIC_ppo(Ynew,covariates_labelled,t))
    alpha<- which.min(GBICscrores)
    
  }
  
  # determine Z
  labelled_Z <- polynomial(labelled_data[,-1],alpha)
  unlabelled_Z <- polynomial(unlabelled_data,alpha)
  
  
  
  # weights (w_i) using in the our optimization problem to get the proposed semi-supervised estimate 
  unlabelled_Z_mean <- colMeans(unlabelled_Z)
  labelled_Z_seondmoment <- t(labelled_Z)%*%labelled_Z/n
  weights_loss <- as.vector(c1+(1-c1)*t(unlabelled_Z_mean)%*%solve(labelled_Z_seondmoment)%*%t(labelled_Z))
  
  
  
  # estimate of the target parameter
  if (type=="linear"){
    hattheta_PSSE <- as.vector(solve(t(cbind(rep(1,n),labelled_data[,-1]))%*%
                                       (weights_loss*cbind(rep(1,n),labelled_data[,-1])))%*%t(cbind(rep(1,n),labelled_data[,-1]))%*%
                                 (weights_loss*labelled_data[,1]))
  }
  
  if (type=="logistic"){
    hattheta_PSSE <- newton_iteration(labelled_data,tol=10^(-5),max_i=50,weights=weights_loss)
  }
  
  else if(type=="quantile"){
    
    beta_initial<- rq(as.vector(labelled_data[,1]) ~ labelled_data[,-1], tau=tau, weights=weights_loss*(weights_loss>=0), 
                      data = as.data.frame(labelled_data))$coefficients
    if(sum(weights_loss<0)>0)
    {
      tol=10^(-5)
      max_i=100
      rate=10^(-3)
      beta_initial<- gradient_descent(labelled_data,beta_initial=beta_initial,rate=rate,
                                      tol,max_i,weights=weights_loss,type="quantile",tau=tau)
      
    }
    hattheta_PSSE<- beta_initial
  }
  hattheta_PSSE = list("Hattheta"=hattheta_PSSE)
  
  ### estimating the sd
  if(sd==TRUE){
    
    if(type=="linear"){
      
      
      X_secondmoment_inverse=solve(t(cbind(rep(1,N),unlabelled_data))%*%cbind(rep(1,N),unlabelled_data)/N) 
      
      set.seed(20218080)
      index=createFolds(1:n, k = Kfolds) # data splitting
      W1_test = vector()
      
      
      for(k in 1:Kfolds){
        
        index_k=as.vector(index[[k]])
        Yt_train = labelled_data[-index_k,1]
        Xt_train = labelled_data[-index_k,-1]
        Zt_train = labelled_Z[-index_k,]
        Yt_test = labelled_data[index_k,1]
        Xt_test = labelled_data[index_k,-1]
        Zt_test = labelled_Z[index_k,]
        
        nrow_Xt_train = n-length(index_k)
        nrow_Xt_test = length(index_k)
        
        
        
        # projection matrix A(theta) estimation 
        
        hattheta_supervised_train_cv <- hattheta_PSSE[[1]]#lm(Yt_train~Xt_train)$coefficients
        L_firstder_train_cv <- as.vector(Yt_train-cbind(rep(1,nrow_Xt_train),Xt_train)%*%hattheta_supervised_train_cv)*
          cbind(rep(1,nrow_Xt_train),Xt_train)
        L_firstder_projection_cof <- solve(t(Zt_train)%*%Zt_train/nrow_Xt_train)%*%
          t(Zt_train)%*%L_firstder_train_cv/nrow_Xt_train
        
        
        # estimate variance by test data
        
        L_firstder_test_cv <- as.vector(Yt_test-cbind(rep(1,nrow_Xt_test),Xt_test)%*%hattheta_supervised_train_cv)*
          cbind(rep(1,nrow_Xt_test),Xt_test)
        
        W1_test_k <- L_firstder_test_cv+(c1-1)*Zt_test%*%L_firstder_projection_cof
        W1_test<- rbind(W1_test,W1_test_k)
        
        
        
      }
      
      
      
      L_firstder_total <- as.vector(labelled_data[,1]-(cbind(rep(1,n),labelled_data[,-1]))%*%hattheta_PSSE[[1]])*
        (cbind(rep(1,n),labelled_data[,-1]))
      L_firstder_projection_cof_total <- solve(t(labelled_Z)%*%labelled_Z/n)%*%t(labelled_Z)%*%L_firstder_total/n
      W2_total <- (1-c1)*unlabelled_Z%*%L_firstder_projection_cof_total
      
      W1_covariance <- t(W1_test)%*%W1_test/n
      W2_covariance <- t(W2_total)%*%W2_total/N
      Vc_hat_semi = W1_covariance+(n/N)*W2_covariance
      
      Var_matrix_hat = X_secondmoment_inverse%*%Vc_hat_semi%*%X_secondmoment_inverse 
      sd_hattheta_PSSE = sqrt(diag(Var_matrix_hat/n))
    }
    
    if(type=="logistic"){
      
      exp_linear_combined <- as.vector(cbind(rep(1,N),unlabelled_data)%*%hattheta_PSSE[[1]])
      weights_second_derivative <- 1/(1+exp(-exp_linear_combined))^2*exp(-exp_linear_combined)
      X_secondmoment_inverse=solve(t(weights_second_derivative*cbind(rep(1,N),unlabelled_data))%*%cbind(rep(1,N),unlabelled_data)/N) 
      
      set.seed(20218080)
      index=createFolds(1:n, k = Kfolds) # data splitting
      W1_test = vector()
      
      
      for(k in 1:Kfolds){
        
        index_k=as.vector(index[[k]])
        Yt_train = labelled_data[-index_k,1]
        Xt_train = labelled_data[-index_k,-1]
        Zt_train = labelled_Z[-index_k,]
        Yt_test = labelled_data[index_k,1]
        Xt_test = labelled_data[index_k,-1]
        Zt_test = labelled_Z[index_k,]
        
        nrow_Xt_train = n-length(index_k)
        nrow_Xt_test = length(index_k)
        
        
        
        # projection matrix A(theta) estimation 
        
        hattheta_supervised_train_cv <- hattheta_PSSE[[1]]#lm(Yt_train~Xt_train)$coefficients
        exp_linear_combined_train <- as.vector(cbind(rep(1,nrow_Xt_train),Xt_train)%*%hattheta_supervised_train_cv)
        L_firstder_train_cv <- as.vector(1/(1+exp(-exp_linear_combined_train))-Yt_train)*
          cbind(rep(1,nrow_Xt_train),Xt_train)
        L_firstder_projection_cof <- solve(t(Zt_train)%*%Zt_train/nrow_Xt_train)%*%
          t(Zt_train)%*%L_firstder_train_cv/nrow_Xt_train
        
        
        # estimate variance by test data
        
        exp_linear_combined_test <- as.vector(cbind(rep(1,nrow_Xt_test),Xt_test)%*%hattheta_supervised_train_cv)
        L_firstder_test_cv <- as.vector(1/(1+exp(-exp_linear_combined_test))-Yt_test)*
          cbind(rep(1,nrow_Xt_test),Xt_test)
        
        W1_test_k <- L_firstder_test_cv+(c1-1)*Zt_test%*%L_firstder_projection_cof
        W1_test<- rbind(W1_test,W1_test_k)
        
      }
      
      
      exp_linear_combined_total <- as.vector(cbind(rep(1,n),labelled_data[,-1])%*%hattheta_PSSE[[1]])
      L_firstder_total <- (1/(1+exp(-exp_linear_combined_total))-labelled_data[,1])*
        (cbind(rep(1,n),labelled_data[,-1]))
      L_firstder_projection_cof_total <- solve(t(labelled_Z)%*%labelled_Z/n)%*%t(labelled_Z)%*%L_firstder_total/n
      W2_total <- (1-c1)*unlabelled_Z%*%L_firstder_projection_cof_total
      
      W1_covariance <- t(W1_test)%*%W1_test/n
      W2_covariance <- t(W2_total)%*%W2_total/N
      Vc_hat_semi = W1_covariance+(n/N)*W2_covariance
      
      Var_matrix_hat = X_secondmoment_inverse%*%Vc_hat_semi%*%X_secondmoment_inverse 
      sd_hattheta_PSSE = sqrt(diag(Var_matrix_hat/n))
      
    }
    
    if(type=="quantile"){
      
      set.seed(20218080)
      index=createFolds(1:n, k = Kfolds) # data splitting
      W1_test = vector()
      
      
      for(k in 1:Kfolds){
        
        index_k=as.vector(index[[k]])
        Yt_train = labelled_data[-index_k,1]
        Xt_train = labelled_data[-index_k,-1]
        Zt_train = labelled_Z[-index_k,]
        Yt_test = labelled_data[index_k,1]
        Xt_test = labelled_data[index_k,-1]
        Zt_test = labelled_Z[index_k,]
        
        nrow_Xt_train = n-length(index_k)
        nrow_Xt_test = length(index_k)
        
        
        residual_train_indicator <- as.vector(Yt_train-cbind(rep(1,nrow_Xt_train),Xt_train)%*%hattheta_PSSE[[1]])<=0
        L_firstder_train_cv <- (residual_train_indicator-tau)*cbind(rep(1,nrow_Xt_train),Xt_train)
        L_firstder_projection_cof <- solve(t(Zt_train)%*%Zt_train)%*%
          t(Zt_train)%*%L_firstder_train_cv
        
        
        # estimate variance by test data
        
        residual_test_indicator <- as.vector(Yt_test-cbind(rep(1,nrow_Xt_test),Xt_test)%*%hattheta_PSSE[[1]])<=0
        L_firstder_test_cv <- (residual_test_indicator-tau)*cbind(rep(1,nrow_Xt_test),Xt_test)
        
        W1_test_k <- L_firstder_test_cv+(c1-1)*Zt_test%*%L_firstder_projection_cof
        W1_test<- rbind(W1_test,W1_test_k)
        
      }
      
      
      residual_total_indicator <- as.vector(labelled_data[,1]-cbind(rep(1,n),labelled_data[,-1])%*%hattheta_PSSE[[1]])<=0
      L_firstder_total <- (residual_total_indicator-tau)*(cbind(rep(1,n),labelled_data[,-1]))
      L_firstder_projection_cof_total <- solve(t(labelled_Z)%*%labelled_Z)%*%t(labelled_Z)%*%L_firstder_total
      W2_total <- (1-c1)*unlabelled_Z%*%L_firstder_projection_cof_total
      
      W1_covariance <- t(W1_test)%*%W1_test/n
      W2_covariance <- t(W2_total)%*%W2_total/N
      Vc_hat_semi = W1_covariance+(n/N)*W2_covariance
      
      
      
      ##estimating the second derivatives (M)
      B=2000
      G=mvrnorm(B,rep(0,(p+1)),diag(rep(1,(p+1))))
      theta_check_semi <- apply(1/sqrt(n)*G,1,function(t) hattheta_PSSE[[1]]+t)
      residual_semi_indicator= apply(cbind(rep(1,n),labelled_data[,-1])%*%theta_check_semi,2, function(t) as.vector(labelled_data[,1]-t)<=0)
      U_check_semi <- t(residual_semi_indicator-tau)%*%cbind(rep(1,n),labelled_data[,-1])/sqrt(n)
      
      hat_M_semi <- t(apply(U_check_semi,2, function(t) lm(t~G-1)$coefficients))
      
      
      if(class(try(solve(hat_M_semi),silent=T))[1]=="try-error"){
        X_secondmoment_inverse <- solve((hat_M_semi+t(hat_M_semi))/2)
      } else{
        X_secondmoment_inverse<- solve(hat_M_semi)
      }
      
      
      
      
      Var_matrix_hat = X_secondmoment_inverse%*%Vc_hat_semi%*%t(X_secondmoment_inverse)
      sd_hattheta_PSSE = sqrt(diag(Var_matrix_hat)/n)
    }
    
    
    hattheta_PSSE = append(hattheta_PSSE,list("sd.of.hattheta"=sd_hattheta_PSSE))
  }
  
  hattheta_PSSE = append(hattheta_PSSE,list("alpha"=alpha))
  
  return(hattheta_PSSE)
}


##################(3) EASE proposed by Chakrabortty and Cai (2018) ##########################
EASE<-function(labelled_data,unlabelled_data,K,H,r)
{
  #---------------------------------Arguments-------------------------------------------------#
  # Purpose: This function is to implement the method proposed by Chakrabortty and Cai (2018) 
  #          for linear working model. This is the main function, and it also involves a 
  #          subfunction defined in (3.1). 
  #   
  # Input: 
  #       labelled_data: Same as in the function "PSSE". 
  #       unlabelled_data: Same as in the function "PSSE".
  #       K: K-folds for sliced inverse regression.
  #       H: H slices in sliced inverse regression.
  #       r: The reduced dimension. 
  #
  # Output: 
  #        Hattheta: The estimate for the target parameter. 
  #------------------------------------------------------------------------------------------#
  
  n=nrow(labelled_data)
  p=ncol(labelled_data)-1
  N=nrow(unlabelled_data)
  
  X_combined<- rbind(labelled_data[,-1],unlabelled_data)
  hat_Gamma_semisupervised= t(cbind(rep(1,n+N),X_combined))%*%cbind(rep(1,n+N),X_combined)/(n+N)
  coef_supervised = lm(labelled_data[,1]~labelled_data[,-1])$coefficients
  set.seed(20210907)
  index=createFolds(1:n, k = K) # data splitting
  
  
  #####estimate phi_k(k=1,...,K) and parameters
  hat_nonpara_test_vector <- numeric()
  hat_nonpara_unlabelled_matrix <- vector()
  covariates_test = vector()
  response_test = numeric()
  
  for (k in 1:K)
  {
    
    index_k=as.vector(index[[k]])
    Yt_train = labelled_data[-index_k,1]
    Xt_train = labelled_data[-index_k,-1]
    Yt_test = labelled_data[index_k,1]
    Xt_test = labelled_data[index_k,-1]
    X_combined_train = rbind(Xt_train, unlabelled_data)
    nrow_Xt_train = n-length(index_k)
    nrow_Xt_test = length(index_k)
    
    covariates_test = rbind(covariates_test,Xt_test) #save the covariates for test 
    response_test <- c(response_test,Yt_test)
    
    
    ###SS-STR
    
    ##step(0):estimate mean of covariance of predictors X and standardized 
    X_combined_train_mean = colMeans(X_combined_train)
    X_train_covariance = cov(X_combined_train)
    eigen_X_train_covariance = eigen(X_train_covariance)
    lamda_sqrt <- diag(sqrt(eigen_X_train_covariance$values))
    Sigma_sqrt <- (eigen_X_train_covariance$vectors)%*%(lamda_sqrt*solve(eigen_X_train_covariance$vectors))
    Sigma_negative_sqrt <- solve(Sigma_sqrt)  # calculate the Sigma^{-1/2}
    X_train_standard = t(Sigma_negative_sqrt %*% t(X_combined_train - X_combined_train_mean))
    Xt_train_standard = X_train_standard[1:nrow_Xt_train,]
    Xv_standard = X_train_standard[(nrow_Xt_train+1):(nrow_Xt_train+N),]
    
    
    ##step(i)
    range_Yt_train = range(Yt_train)
    intrval_length = (range_Yt_train[2]-range_Yt_train[1])/H
    indicator_Yt_train = floor(Yt_train/intrval_length)+1 # the interval indicator
    group_with_train = sort(unique(indicator_Yt_train))
    prob_h_k_train = prop.table(table(group_with_train))
    
    
    ##step(ii)
    distance_Xv_Xt_train = matrix(rep(0,N*nrow_Xt_train),N,nrow_Xt_train)
    for(i in 1:N)
    {
      
      for(j in 1:nrow_Xt_train)
      {
        distance_Xv_Xt_train[i,j] = sum((Xv_standard[i,]-Xt_train_standard[j,])^2)
      }
      
    }
    index_Xv_imputed = apply(distance_Xv_Xt_train,1,which.min)
    #index_Xv_tmp_imputed
    indicator_Xv_imputed = indicator_Yt_train[index_Xv_imputed]
    #indicator_Xv_tmp_imputed
    indicator_X_train = c(indicator_Yt_train,indicator_Xv_imputed)
    
    
    ##step(iii)
    M_matrix_train_semi = matrix(rep(0,p*p),p,p)
    for(kk in 1: length(group_with_train))
    {         
      index_train_standard_kk<- which(indicator_X_train==group_with_train[kk])
      if(length(index_train_standard_kk)==1)
      {
        M_bar_sampleaverage = X_train_standard[index_train_standard_kk,]
        M_matrix_train_semi = M_matrix_train_semi + prob_h_k_train[kk]*M_bar_sampleaverage%*%t(M_bar_sampleaverage)
      }
      if (length(index_train_standard_kk)>1)
      {
        M_bar_sampleaverage = colMeans(X_train_standard[index_train_standard_kk,])
        M_matrix_train_semi = M_matrix_train_semi + prob_h_k_train[kk]*M_bar_sampleaverage%*%t(M_bar_sampleaverage)
      }
      
    }
    
    hat_pstar_h_k = eigen(M_matrix_train_semi)$vectors[,1:r]
    hat_Zero_pstar_h_k = Sigma_negative_sqrt %*% hat_pstar_h_k
    
    Xt_train_SIR = Xt_train %*% hat_Zero_pstar_h_k
    Xt_test_SIR = Xt_test %*% hat_Zero_pstar_h_k
    Xt_unlabelled_SIR = unlabelled_data%*%hat_Zero_pstar_h_k
    
    
    ###calculate the nonparametric estimator in test and unlabelled data
    
    h_opt_chark = optimize(bandwidth.cv,c(0.1,2),Y=Yt_train,X=Xt_train_SIR,d=(K-1),order=r)$minimum
    
    
    hat_nonpara_test = rep(0,nrow_Xt_test)
    for (kkk in 1:nrow_Xt_test)
    {
      kernel = apply(Xt_train_SIR,1, function(t){kernel_guassian(sqrt(sum((Xt_test_SIR[kkk,]-t)^2))/h_opt_chark^r)}) 
      kernel_mean = sum(kernel)
      kernel_Y_mean = t(kernel)%*%Yt_train
      hat_nonpara_test[kkk] = kernel_Y_mean/(kernel_mean+10^(-6))
      
    }
    #hat_nonpara_test  #the nonparametric estimator in test data
    
    hat_nonpara_test_vector <- c(hat_nonpara_test_vector,hat_nonpara_test)
    
    
    
    hat_nonpara_unlabelled = rep(0,N)
    for (kkkk in 1:N)
    {
      kernel = apply(Xt_train_SIR,1, function(t){kernel_guassian(sqrt(sum((Xt_unlabelled_SIR[kkkk,]-t)^2))/h_opt_chark^r)}) 
      kernel_mean = sum(kernel)
      kernel_Y_mean = t(kernel)%*%Yt_train
      hat_nonpara_unlabelled[kkkk] = kernel_Y_mean/(kernel_mean+10^(-6))
      
    }
    #hat_nonpara_unlabelled #the nonparametric estimator in unlabelled data
    
    hat_nonpara_unlabelled_matrix<- cbind(hat_nonpara_unlabelled_matrix,hat_nonpara_unlabelled)
    
  }
  
  
  
  ###estimate eta
  
  hat_eta_k_semisupervised = lm((response_test-hat_nonpara_test_vector)~covariates_test)$coefficients
  #hat_eta_k_semisupervised
  
  ###estimate mu for the unlabelled data
  
  hat_mu_unlabelled = rowMeans(hat_nonpara_unlabelled_matrix)+
    as.vector(cbind(rep(1,N),unlabelled_data) %*% hat_eta_k_semisupervised)
  
  
  ###(SNP estimator) estimate theta by unlabelled data
  
  hattheta_K <- lm(hat_mu_unlabelled~unlabelled_data)$coefficients
  
  
  
  ### EASE estimator
  
  
  ##step(i) estimate hat_eta_k, hat_mu_k,hat_phi_k,hat_phi_0
  
  length_intial <- 0
  hat_Gamma_semisupervised_inverse <- solve(hat_Gamma_semisupervised)
  
  hat_phi_k_cv_matrix<- vector()
  hat_phi_0_cv_matrix <- vector()
  for (k in 1:K)
  {
    length_k = length(index[[k]])
    length_intial = length_intial+length_k
    index_test_k = (1+length_intial-length_k):length_intial
    index_exclude_k = (1:n)[-index_test_k]
    hat_eta_k_cv <- lm((response_test[index_exclude_k]-hat_nonpara_test_vector[index_exclude_k])~
                         covariates_test[index_exclude_k,])$coefficients
    
    hat_mu_k_cv <- hat_nonpara_test_vector[index_test_k]+
      as.vector(cbind(rep(1,length_k),covariates_test[index_test_k,])%*%hat_eta_k_cv)
    
    hat_phi_k_cv <- as.vector(response_test[index_test_k]-hat_mu_k_cv)*cbind(rep(1,length_k),
                                                                             covariates_test[index_test_k,])%*% hat_Gamma_semisupervised_inverse
    
    hat_phi_0_cv <- as.vector(response_test[index_test_k]-cbind(rep(1,length_k),
                                                                covariates_test[index_test_k,])%*%coef_supervised)*cbind(rep(1,length_k),
                                                                                                                         covariates_test[index_test_k,])%*% hat_Gamma_semisupervised_inverse
    
    hat_phi_k_cv_matrix <- rbind(hat_phi_k_cv_matrix,hat_phi_k_cv)
    hat_phi_0_cv_matrix <- rbind(hat_phi_0_cv_matrix,hat_phi_0_cv)
    
  }
  
  ##step(ii) calculate delta_lk
  
  hat_delta_12 = colMeans(hat_phi_0_cv_matrix*(hat_phi_k_cv_matrix-hat_phi_0_cv_matrix))
  hat_delta_22 =colMeans((hat_phi_k_cv_matrix-hat_phi_0_cv_matrix)^2)
  epsilon_n = 1/n^(3/8)
  delta_lk = hat_delta_12/(hat_delta_22+epsilon_n)
  
  
  
  ##step(iii) calculate EASE estimator
  
  hat_theta_semi_Chark<- coef_supervised+delta_lk*(hattheta_K-coef_supervised)
  
  
  
  hat_theta_semi_Chark = list("Hattheta"=hat_theta_semi_Chark)
  
  

  
  return(hat_theta_semi_Chark)
}
#EASE(data_labelled,data_unlabelled,K=5,H=40,r=2)[[1]]



#############(3.1) bandwidth selection
kernel_guassian <- function(a) #guassian kernel
{
  1/sqrt(2*pi)*exp(-a^2/2)
}



bandwidth.cv <- function(Y,X,d,h,order) #bandwidth selection
{
  #------------------------Arguments---------------------------------#
  # Purpose: This function is used for the bandwidth selection
  #          in the function "EASE". After SIR, we need to estimate
  #          the $E[Y|P_rX]$ using kernel smoothing. 
  #
  # Input:
  #       Y: Response observation vector.
  #       X: Covariate matrix.
  #       d: d-folds.
  #       h: bandwidth
  #       order: The power.
  #------------------------------------------------------------------#
  
  n = length(Y)
  p = ncol(X)
  
  set.seed(20220907)
  index=createFolds(1:n, k = d)
  
  err_test=0
  for (i in 1:d)
  {
    index_i=as.vector(index[[i]])
    Yt_train = Y[-index_i]
    Xt_train = X[-index_i,]
    Yt_test = Y[index_i]
    Xt_test = X[index_i,]
    n_positive=length(index_i)
    
    
    for(j in 1:n_positive)
    {
      kernel = apply(Xt_train,1, function(t){kernel_guassian(sqrt(sum((Xt_test[j,]-t)^2))/h^order)}) 
      kernel_mean = sum(kernel)
      kernel_Y_mean = t(kernel)%*%Yt_train
      bias_nonpara_estimator = Yt_test[j]-kernel_Y_mean/kernel_mean
      err_test = err_test + bias_nonpara_estimator^2
    }
    
  } 
  return(err_test)
  
}


polynomial<- function(data,order) 
{
  #---------------------------Arguments----------------------------------------#
  # Purpose: This function is to produce the polynomial basis of X including 
  #           the intercept vector. 
  #
  # Input: 
  #       data: A matrix, whose each row is an observation of predictor vector.
  #       order: The polynomial order.
  #----------------------------------------------------------------------------#
  polynomial_Z<- rep(1,nrow(data))
  for (i in 1:order)
  {
    polynomial_i <- data^i
    polynomial_Z <- cbind(polynomial_Z,polynomial_i)
  }
  return(polynomial_Z)
}


###############(1.2) the data driven selector GBIC_ppo  
GBIC_ppo <- function (Y_new,covariates_matrix,order_poly){
  #-----------------------Arguments--------------------------------------------------# 
  # Purpose: GBIC_ppo criterion function 
  #
  # Input:
  #       Y_new: The first derivatives of the loss function L.
  #       covariates_matrix: A matrix, each row is an observation of predictor vector.
  #       order_poly: The polynomial order.
  #
  #---------------------------------------------------------------------------------#
  d <- ncol(Y_new)
  n <- nrow(Y_new)
  p <- ncol(covariates_matrix)
  
  polynomial_matrix <- polynomial(covariates_matrix, order_poly)
  
  gammahat <- lm(Y_new~polynomial_matrix-1)$coefficients
  residual_square <- (Y_new - polynomial_matrix %*% gammahat)^2
  sigmahat_square <- colSums(residual_square)/(n-p*order_poly-1)
  design_matrix <- t(polynomial_matrix)%*%polynomial_matrix/n
  
  
  if(class(try(solve(design_matrix),silent=T))[1]=="try-error"){
    
    GBIC<- 9999999
    
  }else {
    trace_det_AB = numeric()
    for (j in 1:d){
      Bhat_j<- t(polynomial_matrix*residual_square[,j])%*%polynomial_matrix/n
      cova_constrast <- 1/sigmahat_square[j]*solve(design_matrix)%*%Bhat_j
      trace_AB <- sum(diag(cova_constrast))
      det_AB <- det(cova_constrast)
      trace_det_AB=c(trace_det_AB,trace_AB-log(det_AB))
    }
    GBIC = 1/n*(d*(n-p*order_poly-1)+ n*sum(log(sigmahat_square)) + d*log(n)*p*order_poly+sum(trace_det_AB))
  }
  
  return(GBIC)
  
}


#################(1.3) the supervised estimate based on the given samples
supervised_one<- function(labelled_data){
  #---------------------Arguments-------------------------------------------#
  # Purpose: According the type of loss function, this function is used to  
  #         get the estimate for the target parameter only based on the given 
  #         labelled samples.
  #------------------------------------------------------------------------#
  
  hattheta_supervised <- lm(labelled_data[,1]~labelled_data[,-1])$coefficients
  
  return(hattheta_supervised)
}




#################(1.4)the first derivative of loss function L
L_first_derivative<- function(labelled_data,hattheta_supervised)
{
  #---------------------Arguments-------------------------------------------#
  # Purpose: According the type of loss function, this function is used to  
  #         determine the first-order derivatives of the loss function based on 
  #         an estimate of the target parameter and a labelled data set.
  #------------------------------------------------------------------------#
  n=nrow(labelled_data)
  
  
  residuals<- as.vector(labelled_data[,1]-cbind(rep(1,n),labelled_data[,-1])%*%hattheta_supervised)
  L_first_derivative <- residuals*cbind(rep(1,n),labelled_data[,-1])
  
  
  return(L_first_derivative)
}



PSSE1 <- function(labelled_data,unlabelled_data,c1=NULL,
                  type="linear",tau=0.5,alpha=NULL,gamma=10,sd=FALSE,Kfolds=5)
{
  #-----------------------------------------Arguments-----------------------------------------------------#
  # Purpose: This function is to implement our proposed semi-supervised method. 
  #           This function is the main function, and it also involves several sub-functions defined in
  #           (1.1) - (1.6) below.
  #
  # Input: 
  #      labelled_data: A matrix, whose each row is an observation of predictor vector and the response
  #                      variable and its first column is the observation vector of the response variable.  
  #      unlabelled_data: A matrix, whose each row is an observation of predictor vector. 
  #      c1: Weight. It is the weight to balance the contribution of labelled data and unlabelled data.
  #           Default value is nrow(labelled_data)/(nrow(labelled_data)+nrow(unlabelled_data)).
  #      type: A specification for the type of loss function $L$. We offer three choices: 
  #             "linear", "logistic", "quantile". Default value is "linear". 
  #      tau: The quantile level, which is only useful for "type=quantile". Default value is 0.5.
  #      alpha: The polynomial order. If it is NULL, then we use the data-driven selector GBIC_ppo to 
  #              determine the polynomial order. Default is NULL. 
  #      gamma: The maximum of possible polynomial order. Default is 10. Only if alpha is NULL, gamma will 
  #              be usful for the selection of polynomial order in the construction of Z. 
  #      sd: Return the pointwise standard deviation estimate based on the available samples? Default 
  #           value is FALSE.  
  #      Kfolds: K-folds cross validation method to avoiding the over-fitting during estimating variance, 
  #          which is only useful When "sd=TRUE".
  #
  # Output: 
  #       Hattheta: The estimate for the target parameter. 
  #       sd.of.hattheta: If "sd=TRUE", the elementwise standard deviation estimate for the estimator; 
  #                       otherwise, NULL.
  #       alpha: The selected polynomial order. 
  #------------------------------------------------------------------------------------------------------#
  
  # some basic parameters 
  n=nrow(labelled_data)
  N=nrow(unlabelled_data)
  p=ncol(labelled_data)-1
  
  
  if (p!=ncol(unlabelled_data)){
    print("the dimension of the labelled and unlabelled data is unmatched")
    
    return(NULL)
  }
  
  if(is.null(c1)==TRUE){
    c1=n/(n+N)
  }
  
  # whether the polynomial order needs to be selected or not 
  if (is.null(alpha)==TRUE){ 
    
    # data splitting 
    labelled_data_part1<- labelled_data[1:round(n/2),]
    labelled_data_part2<- labelled_data[(round(n/2)+1):n,]
    
    
    # the supervised estimate
    hattheta_supervised_part1<- supervised_one(labelled_data_part1)
    hattheta_supervised_part2<- supervised_one(labelled_data_part2)
    
    
    # the first derivative of the loss function L using the supervised estimate 
    L_first_derivative_part1<- L_first_derivative(labelled_data_part1,hattheta_supervised_part2)
    L_first_derivative_part2<- L_first_derivative(labelled_data_part2,hattheta_supervised_part1)
    
    
    # new data integrating the first derivatives of loss function
    Ynew<- rbind(L_first_derivative_part1,L_first_derivative_part2)
    covariates_labelled<- labelled_data[,-1]
    
    
    # GBIC_ppo to select the polynomial order 
    GBICscrores<-apply(as.matrix(1:gamma,1,gamma), 1, function(t) GBIC_ppo(Ynew,covariates_labelled,t))
    alpha<- which.min(GBICscrores)
    
  }
  
  # determine Z
  labelled_Z <- polynomial(labelled_data[,-1],alpha)
  unlabelled_Z <- polynomial(unlabelled_data,alpha)
  
  
  
  # weights (w_i) using in the our optimization problem to get the proposed semi-supervised estimate 
  unlabelled_Z_mean <- colMeans(unlabelled_Z)
  labelled_Z_seondmoment <- t(labelled_Z)%*%labelled_Z/n
  weights_loss <- as.vector(c1+(1-c1)*t(unlabelled_Z_mean)%*%solve(labelled_Z_seondmoment)%*%t(labelled_Z))
  
  
  
  # estimate of the target parameter
  if (type=="linear"){
    hattheta_PSSE <- as.vector(solve(t(cbind(rep(1,n),labelled_data[,-1]))%*%
                                       (weights_loss*cbind(rep(1,n),labelled_data[,-1])))%*%t(cbind(rep(1,n),labelled_data[,-1]))%*%
                                 (weights_loss*labelled_data[,1]))
  }
  
  if (type=="logistic"){
    hattheta_PSSE <- newton_iteration(labelled_data,tol=10^(-5),max_i=50,weights=weights_loss)
  }
  
  else if(type=="quantile"){
    
    beta_initial<- rq(as.vector(labelled_data[,1]) ~ labelled_data[,-1], tau=tau, weights=weights_loss*(weights_loss>=0), 
                      data = as.data.frame(labelled_data))$coefficients
    if(sum(weights_loss<0)>0)
    {
      tol=10^(-5)
      max_i=100
      rate=10^(-3)
      beta_initial<- gradient_descent(labelled_data,beta_initial=beta_initial,rate=rate,
                                      tol,max_i,weights=weights_loss,type="quantile",tau=tau)
      
    }
    hattheta_PSSE<- beta_initial
  }
  hattheta_PSSE = list("Hattheta"=hattheta_PSSE)
  
  ### estimating the sd
  if(sd==TRUE){
    
    if(type=="linear"){
      
      
      X_secondmoment_inverse=solve(t(cbind(rep(1,N),unlabelled_data))%*%cbind(rep(1,N),unlabelled_data)/N) 
      
      set.seed(20218080)
      index=createFolds(1:n, k = Kfolds) # data splitting
      W1_test = vector()
      
      
      for(k in 1:Kfolds){
        
        index_k=as.vector(index[[k]])
        Yt_train = labelled_data[-index_k,1]
        Xt_train = labelled_data[-index_k,-1]
        Zt_train = labelled_Z[-index_k,]
        Yt_test = labelled_data[index_k,1]
        Xt_test = labelled_data[index_k,-1]
        Zt_test = labelled_Z[index_k,]
        
        nrow_Xt_train = n-length(index_k)
        nrow_Xt_test = length(index_k)
        
        
        
        # projection matrix A(theta) estimation 
        
        hattheta_supervised_train_cv <- hattheta_PSSE[[1]]#lm(Yt_train~Xt_train)$coefficients
        L_firstder_train_cv <- as.vector(Yt_train-cbind(rep(1,nrow_Xt_train),Xt_train)%*%hattheta_supervised_train_cv)*
          cbind(rep(1,nrow_Xt_train),Xt_train)
        L_firstder_projection_cof <- solve(t(Zt_train)%*%Zt_train/nrow_Xt_train)%*%
          t(Zt_train)%*%L_firstder_train_cv/nrow_Xt_train
        
        
        # estimate variance by test data
        
        L_firstder_test_cv <- as.vector(Yt_test-cbind(rep(1,nrow_Xt_test),Xt_test)%*%hattheta_supervised_train_cv)*
          cbind(rep(1,nrow_Xt_test),Xt_test)
        
        W1_test_k <- L_firstder_test_cv+(c1-1)*Zt_test%*%L_firstder_projection_cof
        W1_test<- rbind(W1_test,W1_test_k)
        
        
        
      }
      
      
      
      L_firstder_total <- as.vector(labelled_data[,1]-(cbind(rep(1,n),labelled_data[,-1]))%*%hattheta_PSSE[[1]])*
        (cbind(rep(1,n),labelled_data[,-1]))
      L_firstder_projection_cof_total <- solve(t(labelled_Z)%*%labelled_Z/n)%*%t(labelled_Z)%*%L_firstder_total/n
      W2_total <- (1-c1)*unlabelled_Z%*%L_firstder_projection_cof_total
      
      W1_covariance <- t(W1_test)%*%W1_test/n
      W2_covariance <- t(W2_total)%*%W2_total/N
      Vc_hat_semi = W1_covariance+(n/N)*W2_covariance
      
      Var_matrix_hat = X_secondmoment_inverse%*%Vc_hat_semi%*%X_secondmoment_inverse 
      sd_hattheta_PSSE = sqrt(diag(Var_matrix_hat/n))
    }
    
    if(type=="logistic"){
      
      exp_linear_combined <- as.vector(cbind(rep(1,N),unlabelled_data)%*%hattheta_PSSE[[1]])
      weights_second_derivative <- 1/(1+exp(-exp_linear_combined))^2*exp(-exp_linear_combined)
      X_secondmoment_inverse=solve(t(weights_second_derivative*cbind(rep(1,N),unlabelled_data))%*%cbind(rep(1,N),unlabelled_data)/N) 
      
      set.seed(20218080)
      index=createFolds(1:n, k = Kfolds) # data splitting
      W1_test = vector()
      
      
      for(k in 1:Kfolds){
        
        index_k=as.vector(index[[k]])
        Yt_train = labelled_data[-index_k,1]
        Xt_train = labelled_data[-index_k,-1]
        Zt_train = labelled_Z[-index_k,]
        Yt_test = labelled_data[index_k,1]
        Xt_test = labelled_data[index_k,-1]
        Zt_test = labelled_Z[index_k,]
        
        nrow_Xt_train = n-length(index_k)
        nrow_Xt_test = length(index_k)
        
        
        
        # projection matrix A(theta) estimation 
        
        hattheta_supervised_train_cv <- hattheta_PSSE[[1]]#lm(Yt_train~Xt_train)$coefficients
        exp_linear_combined_train <- as.vector(cbind(rep(1,nrow_Xt_train),Xt_train)%*%hattheta_supervised_train_cv)
        L_firstder_train_cv <- as.vector(1/(1+exp(-exp_linear_combined_train))-Yt_train)*
          cbind(rep(1,nrow_Xt_train),Xt_train)
        L_firstder_projection_cof <- solve(t(Zt_train)%*%Zt_train/nrow_Xt_train)%*%
          t(Zt_train)%*%L_firstder_train_cv/nrow_Xt_train
        
        
        # estimate variance by test data
        
        exp_linear_combined_test <- as.vector(cbind(rep(1,nrow_Xt_test),Xt_test)%*%hattheta_supervised_train_cv)
        L_firstder_test_cv <- as.vector(1/(1+exp(-exp_linear_combined_test))-Yt_test)*
          cbind(rep(1,nrow_Xt_test),Xt_test)
        
        W1_test_k <- L_firstder_test_cv+(c1-1)*Zt_test%*%L_firstder_projection_cof
        W1_test<- rbind(W1_test,W1_test_k)
        
      }
      
      
      exp_linear_combined_total <- as.vector(cbind(rep(1,n),labelled_data[,-1])%*%hattheta_PSSE[[1]])
      L_firstder_total <- (1/(1+exp(-exp_linear_combined_total))-labelled_data[,1])*
        (cbind(rep(1,n),labelled_data[,-1]))
      L_firstder_projection_cof_total <- solve(t(labelled_Z)%*%labelled_Z/n)%*%t(labelled_Z)%*%L_firstder_total/n
      W2_total <- (1-c1)*unlabelled_Z%*%L_firstder_projection_cof_total
      
      W1_covariance <- t(W1_test)%*%W1_test/n
      W2_covariance <- t(W2_total)%*%W2_total/N
      Vc_hat_semi = W1_covariance+(n/N)*W2_covariance
      
      Var_matrix_hat = X_secondmoment_inverse%*%Vc_hat_semi%*%X_secondmoment_inverse 
      sd_hattheta_PSSE = sqrt(diag(Var_matrix_hat/n))
      
    }
    
    if(type=="quantile"){
      
      set.seed(20218080)
      index=createFolds(1:n, k = Kfolds) # data splitting
      W1_test = vector()
      
      
      for(k in 1:Kfolds){
        
        index_k=as.vector(index[[k]])
        Yt_train = labelled_data[-index_k,1]
        Xt_train = labelled_data[-index_k,-1]
        Zt_train = labelled_Z[-index_k,]
        Yt_test = labelled_data[index_k,1]
        Xt_test = labelled_data[index_k,-1]
        Zt_test = labelled_Z[index_k,]
        
        nrow_Xt_train = n-length(index_k)
        nrow_Xt_test = length(index_k)
        
        
        residual_train_indicator <- as.vector(Yt_train-cbind(rep(1,nrow_Xt_train),Xt_train)%*%hattheta_PSSE[[1]])<=0
        L_firstder_train_cv <- (residual_train_indicator-tau)*cbind(rep(1,nrow_Xt_train),Xt_train)
        L_firstder_projection_cof <- solve(t(Zt_train)%*%Zt_train)%*%
          t(Zt_train)%*%L_firstder_train_cv
        
        
        # estimate variance by test data
        
        residual_test_indicator <- as.vector(Yt_test-cbind(rep(1,nrow_Xt_test),Xt_test)%*%hattheta_PSSE[[1]])<=0
        L_firstder_test_cv <- (residual_test_indicator-tau)*cbind(rep(1,nrow_Xt_test),Xt_test)
        
        W1_test_k <- L_firstder_test_cv+(c1-1)*Zt_test%*%L_firstder_projection_cof
        W1_test<- rbind(W1_test,W1_test_k)
        
      }
      
      
      residual_total_indicator <- as.vector(labelled_data[,1]-cbind(rep(1,n),labelled_data[,-1])%*%hattheta_PSSE[[1]])<=0
      L_firstder_total <- (residual_total_indicator-tau)*(cbind(rep(1,n),labelled_data[,-1]))
      L_firstder_projection_cof_total <- solve(t(labelled_Z)%*%labelled_Z)%*%t(labelled_Z)%*%L_firstder_total
      W2_total <- (1-c1)*unlabelled_Z%*%L_firstder_projection_cof_total
      
      W1_covariance <- t(W1_test)%*%W1_test/n
      W2_covariance <- t(W2_total)%*%W2_total/N
      Vc_hat_semi = W1_covariance+(n/N)*W2_covariance
      
      
      
      ##estimating the second derivatives (M)
      B=2000
      G=mvrnorm(B,rep(0,(p+1)),diag(rep(1,(p+1))))
      theta_check_semi <- apply(1/sqrt(n)*G,1,function(t) hattheta_PSSE[[1]]+t)
      residual_semi_indicator= apply(cbind(rep(1,n),labelled_data[,-1])%*%theta_check_semi,2, function(t) as.vector(labelled_data[,1]-t)<=0)
      U_check_semi <- t(residual_semi_indicator-tau)%*%cbind(rep(1,n),labelled_data[,-1])/sqrt(n)
      
      hat_M_semi <- t(apply(U_check_semi,2, function(t) lm(t~G-1)$coefficients))
      
      
      if(class(try(solve(hat_M_semi),silent=T))[1]=="try-error"){
        X_secondmoment_inverse <- solve((hat_M_semi+t(hat_M_semi))/2)
      } else{
        X_secondmoment_inverse<- solve(hat_M_semi)
      }
      
      
      
      
      Var_matrix_hat = X_secondmoment_inverse%*%Vc_hat_semi%*%t(X_secondmoment_inverse)
      sd_hattheta_PSSE = sqrt(diag(Var_matrix_hat)/n)
    }
    
    
    hattheta_PSSE = append(hattheta_PSSE,list("sd.of.hattheta"=sd_hattheta_PSSE))
  }
  
  hattheta_PSSE = append(hattheta_PSSE,list("alpha"=alpha))
  
  return(hattheta_PSSE)
}





#####################DRESS###########################
DRESS1<- function(labelled_data,unlabelled_data,L,Kfolds)
{
  #-------------------------Arguments-------------------------------------------#
  # Purpose: This function is to implement the method proposed by
  #         Kawakita and Kanamori (2013) for M-estimation,
  #         which uses the desity-ratio to improve the estimation
  #         efficiency.
  #
  # Input:
  #       labelled_data: Same as in the function "PSSE". 
  #       unlabelled_data: Same as in the function "PSSE".     
  #       type: Same as in the function "PSSE".
  #       tau: Same as in the function "PSSE".
  #        L: The polynomial order when we construct the polynomial function as
  #            base function for density ratio estimation.
  #       sd: Return the pointwise standard deviation estimate based on the available
  #           samples? Default value is FALSE.
  #       Kfolds: Same as in the function "PSSE". 
  #
  # Output: 
  #       Hattheta: The estimate for the target parameter. 
  #       sd.of.hattheta: If "sd=TRUE", the elementwise standard deviation 
  #                        estimate for the estimator; otherwise, NULL. 
  #       error: indicator of whether the parameters in density-ratio is properly
  #              estimated. "TRUE" represents "NOT"; "FALSE" represents "YES".
  #-----------------------------------------------------------------------------#
  
  n=nrow(labelled_data)
  p=ncol(labelled_data)-1
  N=nrow(unlabelled_data)
  
  
  base_labelled<- polynomial(labelled_data[,-1],L)
  base_unlabelled <- polynomial(unlabelled_data,L)
  alpha_first_derivative_unlabelled <- colMeans(base_unlabelled)
  
  ## estimating alpha (the parameter in density-ratio)
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
    #print(i)
  }
  
  
  exponential_phi_labelled <- as.vector(exp(base_labelled%*%alpha_initial))
  
  
  ## estimating theta*(the target parameter)

    hattheta_DRESS<- lm(labelled_data[,1]~labelled_data[,-1],weights=exponential_phi_labelled)$coefficients
    
  
  
  hattheta_DRESS = append(list("Hattheta"=hattheta_DRESS),list("error"=error_svd))
  
  
  ### estimating the sd

    
   
      
      
      X_secondmoment_inverse=solve(t(cbind(rep(1,N),unlabelled_data))%*%cbind(rep(1,N),unlabelled_data)/N) 
      
      set.seed(20218080)
      index=createFolds(1:n, k = Kfolds) # data splitting
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
        
        
        
        # projection matrix A(theta) estimation 
        
        hattheta_supervised_train_cv <-  hattheta_DRESS[[1]]#lm(Yt_train~Xt_train)$coefficients
        L_firstder_train_cv <- as.vector(Yt_train-cbind(rep(1,nrow_Xt_train),Xt_train)%*%hattheta_supervised_train_cv)*
          cbind(rep(1,nrow_Xt_train),Xt_train)
        L_firstder_projection_cof <- solve(t(Zt_train)%*%Zt_train/nrow_Xt_train)%*%
          t(Zt_train)%*%L_firstder_train_cv/nrow_Xt_train
        
        
        # estimate variance by test data
        
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



###########################(2) PI proposed by Azriel et al. (2021) #################################
PI <- function(labelled_data,unlabelled_data)
{
  #---------------------------------Arguments------------------------------------------#
  # Purpose: This function is to implement the method proposed by Azriel et al.(2021)
  #          for linear working model. 
  #
  # Input:
  #       labelled_data: Same as in the function "PSSE". 
  #       unlabelled_data: Same as in the function "PSSE".     
  #
  # Output: 
  #        Hattheta: The estimate for the target parameter. 
  #------------------------------------------------------------------------------------#
  n=nrow(labelled_data)
  p=ncol(labelled_data)-1
  N=nrow(unlabelled_data)
  
  # combine all the covairates 
  X_combined <- rbind(labelled_data[,-1],unlabelled_data)
  X_labelled <- labelled_data[,-1]
  
  hat_beta_initial <- numeric()
  X_dot <- matrix(rep(0,n*p),n,p) 
  delta_tilde <- matrix(rep(0,n*p),n,p)
  
  for (j in 1:p)
  {
    ## first step 
    coefficients_negtive_j = lm(X_combined[,j]~X_combined[,-j])$coefficients #unlabeled data
    X_j_dot = X_labelled[,j] - cbind(rep(1,n),X_labelled[,-j]) %*% coefficients_negtive_j #labeled data projection error
    X_j_dot_total = X_combined[,j] - cbind(rep(1,n+N),X_combined[,-j]) %*% coefficients_negtive_j
    X_j_dot_square_sampleaverage=mean(X_j_dot_total^2)
    
    
    ## step step 
    W_j = as.vector(labelled_data[,1])*X_j_dot/X_j_dot_square_sampleaverage
    U1 = X_j_dot/X_j_dot_square_sampleaverage
    X_dot[,j] = as.vector(U1) 
    U = matrix(rep(0,n*p),n,p) 
    for (jj in 1:p)
    {
      if (jj==j)
      {
        U[,jj] = X_labelled[,jj]*X_j_dot/X_j_dot_square_sampleaverage-1
      }
      
      else
      {
        U[,jj] = X_labelled[,jj]*X_j_dot/X_j_dot_square_sampleaverage
      }
      
    }
    
    delta_tilde[,j] = W_j-as.matrix(cbind(rep(1,n),U1,U)) %*% (lm(W_j~U1+U)$coefficients)
    hat_beta_j = lm(W_j~U1+U)$coefficients[1]
    hat_beta_initial = c(hat_beta_initial,hat_beta_j)
  }
  hat_alpha = mean(labelled_data[,1])-t(hat_beta_initial)%*%colMeans(X_labelled)
  hat_theta_Azriel = as.vector(c(hat_alpha,hat_beta_initial))
  
  hat_theta_Azriel = list("Hattheta"=hat_theta_Azriel)
  
  
  
  return(hat_theta_Azriel)
}
#PI(data_labelled,data_unlabelled)








#############################################################################################################################
##############################################ET functions###############################################################
variance_theta_diag <- function(y, X, w_hat, theta_hat, ci_level = 0.95) {
  N <- length(y)
  p <- length(theta_hat)
  
  # Ensure column names exist
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", 1:p)
  }
  
  # Residuals
  resid <- as.numeric(y - X %*% theta_hat)
  
  # Compute d_i
  D <- matrix(NA, nrow = N, ncol = p)
  for (i in 1:N) {
    x_i <- as.numeric(X[i, ])
    D[i, ] <- w_hat[i] * resid[i] * x_i
  }
  d_bar <- colMeans(D)
  
  # τ_hat
  tau_hat <- matrix(0, p, p)
  for (i in 1:N) {
    x_i <- as.numeric(X[i, ])
    tau_hat <- tau_hat + w_hat[i] * (x_i %o% x_i)
  }
  tau_hat <- tau_hat / N
  
  # Middle term
  middle <- matrix(0, p, p)
  for (i in 1:N) {
    diff_i <- D[i, ] - d_bar
    middle <- middle + diff_i %o% diff_i
  }
  middle <- middle / (N * (N - 1))
  
  # Sandwich variance
  tau_inv <- solve(tau_hat)
  var_theta <- tau_inv %*% middle %*% t(tau_inv)
  
  # Extract diagonals (variances)
  variances <- diag(var_theta)
  names(variances) <- colnames(X)
  
  # Standard errors and confidence intervals
  se <- sqrt(variances)
  alpha <- 1 - ci_level
  z_val <- qnorm(1 - alpha / 2)
  lower <- theta_hat - z_val * se
  upper <- theta_hat + z_val * se
  width=upper-lower
  # Output table
  result <- data.frame(
    Variable = colnames(X),
    Estimate = theta_hat,
    Variance = variances,
    SE = se,
    CI_Lower = lower,
    CI_Upper = upper,
    CI_Width = width,
    row.names = NULL
  )
  
  return(result)
}


variance_theta_diag_scaled <- function(y, X, w_hat, theta_hat, ci_level = 0.95,
                                       scale_continuous = TRUE) {
  N <- length(y)
  p <- length(theta_hat)
  X <- as.matrix(X)
  
  # Ensure column names
  if (is.null(colnames(X))) colnames(X) <- paste0("X", 1:p)
  
  # Detect binary (0/1) columns (don't scale these)
  is_binary <- apply(X, 2, function(col) {
    vals <- unique(col)
    all(vals %in% c(0,1)) && length(vals) <= 2
  })
  
  # Choose columns to scale = non-binary (continuous) columns
  scale_idx <- if (scale_continuous) which(!is_binary) else integer(0)
  
  # SDs for columns we will scale (avoid divide-by-zero)
  sd_vec <- rep(1, p)
  if (length(scale_idx) > 0) {
    sd_tmp <- apply(X[, scale_idx, drop = FALSE], 2, sd, na.rm = TRUE)
    sd_tmp[!is.finite(sd_tmp) | sd_tmp == 0] <- 1
    sd_vec[scale_idx] <- sd_tmp
  }
  
  # Scale design WITHOUT centering: x* = x / sd
  Xs <- X
  if (length(scale_idx) > 0) {
    Xs[, scale_idx] <- sweep(X[, scale_idx, drop = FALSE], 2, sd_vec[scale_idx], "/")
  }
  
  # Transform coefficients to the scaled parameterization:
  # beta*_j = beta_j * sd_j  (intercept and dummies unchanged since sd=1)
  thetas <- theta_hat * sd_vec
  
  # Residuals computed on the scaled system give the same y - X beta
  resid <- as.numeric(y - Xs %*% thetas)
  
  # d_i rows (vectorized)
  D <- resid * w_hat * Xs
  d_bar <- colMeans(D)
  
  # τ_hat = (1/N) Σ w_i x_i x_i^T
  Xw <- sqrt(w_hat) * Xs
  tau_hat <- crossprod(Xw) / N
  
  # Middle term: (1/[N(N-1)]) Σ (d_i - d̄)(d_i - d̄)^T
  Dc <- sweep(D, 2, d_bar, "-")
  middle <- crossprod(Dc) / (N * (N - 1))
  
  # Sandwich variance on the scaled system
  tau_inv <- solve(tau_hat)
  var_theta_scaled <- tau_inv %*% middle %*% t(tau_inv)
  
  # Backscale variances to original units:
  # Var(beta_j) = Var(beta*_j) / sd_j^2
  variances <- diag(var_theta_scaled) / (sd_vec^2)
  names(variances) <- colnames(X)
  
  # SE and CI in original units
  se <- sqrt(variances)
  alpha <- 1 - ci_level
  z_val <- qnorm(1 - alpha / 2)
  lower <- theta_hat - z_val * se
  upper <- theta_hat + z_val * se
  
  # Result table in original units
  result <- data.frame(
    Variable = colnames(X),
    Estimate = theta_hat,
    Variance = variances,
    SE = se,
    CI_Lower = lower,
    CI_Upper = upper,
    row.names = NULL
  )
  
  return(result)
}


build_SU_folds <- function(df, K , id_col = "ID", delta_col = "D", seed = seed) {
  #stopifnot(all(c(id_col, delta_col) %in% names(df)))
  set.seed(seed)
  #df$ID <- seq_len(nrow(df))
  
  id <- df[[id_col]]
  delta <- df[[delta_col]]
  
  idx_S  <- which(delta == 1)  # labeled indices
  idx_U0 <- which(delta == 0)  # unlabeled-only indices
  
  nS  <- length(idx_S)
  nU0 <- length(idx_U0)
  if (nS < K || nU0 < K) stop("Need at least K labeled and K unlabeled observations.")
  
  # helper: random, roughly-equal K-way split
  split_K <- function(idx, K) {
    if (length(idx) == 0) return(rep(list(integer(0)), K))
    idx <- sample(idx)  # shuffle
    split(idx, rep(1:K, length.out = length(idx)))
  }
  
  S_parts  <- split_K(idx_S,  K)  # disjoint labeled parts S_k
  U0_parts <- split_K(idx_U0, K)  # disjoint unlabeled-only parts U0_k
  
  folds <- vector("list", K)
  # fold assignment vectors (same length as df rows), 0 = not in that split
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
    
    
    model.formula <- as.formula(paste("Y ~", paste(sprintf("`%s`", colnames( subset(train_data_labeled, select = -c(D, Y, ID,pi.hat)))), collapse = " + ")))
    
    #lm_model<-lm(model.formula,data=as.data.frame( subset(train_data_labeled, select = -c(D, ID,pi.hat))))
    # gam_model <- gam(Y ~ s(age) + s(BMI) + s(SBP) + s(DBP) + sex + race,
    #             data=as.data.frame( subset(train_data_labeled, select = -c(D, ID,pi.hat))), method = "REML")
    
    #y_hat<-predict( lm_model, newdata=as.data.frame(validation_data_unlabeled))
    
    gam_model <-
      mgcv::gam(
        Y ~
          s(X1) +
          s(X2) +
          s(X3) +
          s(X4),
        data =
          as.data.frame(
            subset(
              train_data_labeled,
              select = -c(
                D,
                ID,
                pi.hat
              )
            )
          ),
        method = "REML"
      )
    
    y_hat <-
      predict(
        gam_model,
        newdata =
          as.data.frame(
            validation_data_unlabeled
          )
      )
    
    #y_hat<-predict( gam_model, newdata=as.data.frame(validation_data_unlabeled))
    #rf_model <- randomForest(model.formula, data = train_data_labeled)
    #y_hat <- predict(rf_model, newdata = validation_data_unlabeled)
    validation_data_unlabeled$y.hat<-y_hat
    
    
    ##############################partition for pi.hat
    #train_data_labeled_pi<-df[-res$folds[[k]]$U_k_idx,]
    #train_data_pi<-as.data.frame(subset(train_data_labeled_pi, select = -c(Y,ID,pi.hat)))
    #validation_data_pi<-as.data.frame(subset(df[res$folds[[k]]$U_k_idx,], select = -c(Y,ID,pi.hat)))
    #formula.glm<-as.formula(paste0("D~",paste0(colnames( train_data_pi[,-1]),collapse = "+")))
    #glm_model<-glm(formula.glm,data=train_data_pi,family=binomial())
    #pi_hat_cv<-predict( glm_model, newdata=as.data.frame(validation_data_pi[,-1]))
    #validation_data_unlabeled$pi.hat.cv<-pi_hat_cv
    data_unlabeled[[k]]<-validation_data_unlabeled
  }
  return(data_unlabeled)
  
}

estimate_theta_EM_kfold_CVXR_ET <- function(th, data_full, K, seed, max.iter,eps) {
  
  # --- prepare cross-fitted data, stacked in fold order ---
  fold_hat <- k_fold_function(df = data_full, K = K, seed = seed)
  data_all <- do.call(rbind, fold_hat)
  n_k      <- sapply(fold_hat, nrow)
  ofs      <- c(0, cumsum(n_k))
  
  # features
  # define p properly:
  x_cols <- subset(data_all,select=-c(Y,D,pi.hat,ID,y.hat))
  p      <- length(x_cols)
  model.formula<-as.formula(paste("~",paste(colnames(x_cols),collapse="+")))
  new.data  <- as.data.frame(data_all[, colnames(x_cols), drop=FALSE])
  new.data1 <- (model.matrix(  model.formula,
                               data = as.data.frame(new.data)))
  
  
  y.hat <- data_all$y.hat
  D     <- data_all$D
  I1    <- which(D == 1)
  
  theta <- as.matrix(th)
  iter  <- 0
  
  repeat {
    iter <- iter + 1
    
    
    
    # residual-based moments
    error <- y.hat - as.numeric(new.data1 %*% theta)
    # H: N x (p+1)
    H <- new.data1 * error  # vectorized (faster & clearer than apply)
    
    # Moment matrix for calibration (choose what you want to balance)
    A <- H # N x r, r = p+1
    
    y<-data_all$Y[I1]
    resid <- y -  as.numeric(new.data1[I1,] %*% theta)
    U_mat1 <-  new.data1[I1,] * as.vector(resid)  # N x p matrix
    
    
    # --- CVXR problem on treated only ---
    r <- ncol(A)
    w <- CVXR::Variable(length(I1), pos = TRUE)
    a <- rep(1/length(I1), length(I1))  # base (uniform)
    g<-log(data_all$pi.hat)
    
    A_treated <- A[I1, , drop = FALSE]
    
    constr <- list(
      sum(w) == 1,
      sum(w * g[I1]) == mean(g),
      t(A_treated) %*% w ==
        matrix(colMeans(A), ncol = 1)
    )
    
    ones <- rep(1,length(I1))
    #objective <- CVXR::Minimize( sum(CVXR::kl_div(w, a)) )
    objective <- CVXR::Minimize( sum(CVXR::kl_div(w, ones)) )
    
    
    prob <- CVXR::Problem(objective, constr)
    sol  <- CVXR::psolve(prob,solver = "ECOS")
    
    
    if (CVXR::status(prob) != "optimal") {
      return(matrix(NA_real_, nrow=1, ncol=ncol(new.data1)))
    }
    
    W  <- as.numeric(CVXR::value(w))
    Q  <- new.data1[I1, , drop=FALSE]
    Y1 <- data_all$Y[I1]
    A<-t(Q) %*% diag(W) %*% Q
    B<-t(Q) %*% diag(W) %*% Y1
    th.new <- solve(A, B)
    #print(th.new)
    #print(max(abs(th.new - as.numeric(theta))))
    if (max(abs(th.new - as.numeric(theta))) < eps) {
      
      
      return(list(theta=matrix(th.new, ncol=1),w=W,estimate.table=variance_theta_diag(
        y = y,
        X = new.data1[I1, ],
        w_hat = W,
        theta_hat = th.new
      )))
    }
    if (iter >= max.iter) {
      message("Maximum iterations reached.")
      return(list(theta=matrix(th.new, ncol=1),w=W,estimate.table=variance_theta_diag(
        y = y,
        X = new.data1[I1, ],
        w_hat = W,
        theta_hat = th.new
      )))
    }
    theta <- as.matrix(th.new)
  }
}



estimate_theta_EM_kfold_CVXR_HD <- function(th, data_full, K, seed, max.iter,eps) {
  
  # --- prepare cross-fitted data, stacked in fold order ---
  fold_hat <- k_fold_function(df = data_full, K = K, seed = seed)
  data_all <- do.call(rbind, fold_hat)
  n_k      <- sapply(fold_hat, nrow)
  ofs      <- c(0, cumsum(n_k))
  
  # features
  # define p properly:
  x_cols <- subset(data_all,select=-c(Y,D,pi.hat,ID,y.hat))
  p      <- length(x_cols)
  model.formula<-as.formula(paste("~",paste(colnames(x_cols),collapse="+")))
  new.data  <- as.data.frame(data_all[, colnames(x_cols), drop=FALSE])
  new.data1 <- (model.matrix(  model.formula,
                               data = as.data.frame(new.data)))
  
  
  y.hat <- data_all$y.hat
  D     <- data_all$D
  I1    <- which(D == 1)
  
  theta <- as.matrix(th)
  iter  <- 0
  
  repeat {
    iter <- iter + 1
    
    
    
    # residual-based moments
    error <- y.hat - as.numeric(new.data1 %*% theta)
    # H: N x (p+1)
    H <- new.data1 * error  # vectorized (faster & clearer than apply)
    
    # Moment matrix for calibration (choose what you want to balance)
    A <- H # N x r, r = p+1
    
    y<-data_all$Y[I1]
    resid <- y -  as.numeric(new.data1[I1,] %*% theta)
    U_mat1 <-  new.data1[I1,] * as.vector(resid)  # N x p matrix
    
    
    # --- CVXR problem on treated only ---
    r <- ncol(A)
    w <- CVXR::Variable(length(I1), pos = TRUE)
    a <- rep(1/length(I1), length(I1))  # base (uniform)
    g<--sqrt(data_all$pi.hat)/2
    
    constr <- list(sum(w) == 1,
                   sum(w*g[I1]) == mean(g))
    
    
    # per-fold constraints on treated rows (means match)
    for (k in seq_len(K)) {
      idx_k_all <- (ofs[k]+1):ofs[k+1]
      idx_k     <- intersect(idx_k_all, I1)
      if (length(idx_k) == 0) next
      
      A_k <- A[idx_k, , drop=FALSE]
      map <- match(idx_k, I1)
      w_k <- w[map]
      A_k_all<-A[ idx_k_all, , drop=FALSE]
      constr <- c(constr, list(
        t(A_k) %*% w_k ==  matrix(colSums(A_k_all)/nrow(A), ncol=1)
      ))
    }
    
    objective <- CVXR::Minimize(-sum(sqrt(w)))
    
    prob <- CVXR::Problem(objective, constr)
    sol  <- CVXR::psolve(prob,solver = "ECOS")
    
    
    if (CVXR::status(prob) != "optimal") {
      return(matrix(NA_real_, nrow=1, ncol=ncol(new.data1)))
    }
    
    W  <- as.numeric(CVXR::value(w))
    Q  <- new.data1[I1, , drop=FALSE]
    Y1 <- data_all$Y[I1]
    A<-t(Q) %*% diag(W) %*% Q
    B<-t(Q) %*% diag(W) %*% Y1
    th.new <- solve(A, B)
    #print(th.new)
    print(max(abs(th.new - as.numeric(theta))))
    if (max(abs(th.new - as.numeric(theta))) < eps) {
      
      
      return(list(theta=matrix(th.new, ncol=1),w=W,estimate.table=variance_theta_diag(
        y = y,
        X = new.data1[I1, ],
        w_hat = W,
        theta_hat = th.new
      )))
    }
    if (iter >= max.iter) {
      message("Maximum iterations reached.")
      return(list(theta=matrix(th.new, ncol=1),w=W,estimate.table=variance_theta_diag(
        y = y,
        X = new.data1[I1, ],
        w_hat = W,
        theta_hat = th.new
      )))
    }
    theta <- as.matrix(th.new)
  }
}




estimate_theta_nested_kfold_CVXR_optim <- function(theta_init, data_full, K , seed ) {
  # ------------------------------
  # Define the objective function L(theta)
  # ------------------------------
  # --- prepare cross-fitted data, stacked in fold order ---
  fold_hat <- k_fold_function(df = data_full, K = K, seed = seed)
  data_all <- do.call(rbind, fold_hat)
  n_k      <- sapply(fold_hat, nrow)
  ofs      <- c(0, cumsum(n_k))
  
  # features
  # define p properly:
  x_cols <- subset(data_all,select=-c(Y,D,pi.hat,ID,y.hat))
  p      <- length(x_cols)
  model.formula<-as.formula(paste("~",paste(colnames(x_cols),collapse="+")))
  new.data  <- as.data.frame(data_all[, colnames(x_cols), drop=FALSE])
  new.data1 <- (model.matrix(  model.formula,
                               data = as.data.frame(new.data)))
  
  
  y.hat <- data_all$y.hat
  D     <- data_all$D
  I1    <- which(D == 1)
  
  g <- log(data_all$pi.hat)
  y <- data_all[I1, "Y"]
  
  L_fn <- function(theta) {
    
    
    # residual-based moments
    error <- y.hat - as.numeric(new.data1 %*% theta)
    H <- new.data1 * error
    A <- H
    r <- ncol(A)
    
    # --- CVXR setup ---
    w <- CVXR::Variable(length(I1), pos = TRUE)
    
    resid <- y - as.matrix(new.data1[I1, ]) %*% theta
    U_mat1 <- as.matrix(new.data1[I1, ]) * as.vector(resid)
    
    # Constraints: global + per-fold
    constr <- list(
      sum(w) == 1,
      sum(w * g[I1]) == mean(g),
      t(U_mat1) %*% w == rep(0, ncol(U_mat1))
    )
    
    mu <- c()
    for (k in seq_len(K)) {
      idx_k_all <- (ofs[k] + 1):ofs[k + 1]
      idx_k     <- intersect(idx_k_all, I1)
      if (length(idx_k) == 0) next
      
      A_k <- A[idx_k, , drop = FALSE]
      map <- match(idx_k, I1)
      w_k <- w[map]
      A_k_all <- A[idx_k_all, , drop = FALSE]
      constr <- c(constr, list(
        t(A_k) %*% w_k == matrix(colSums(A_k_all) / nrow(A), ncol = 1)
      ))
      mu <- c(mu, colSums(A_k_all) / nrow(A))
    }
    
    ones <- rep(1, length(I1))
    U1 <- c(1, mean(g), mu)
    
    objective <- CVXR::Minimize(sum(CVXR::kl_div(w, ones)))
    prob <- CVXR::Problem(objective, constr)
    sol  <- CVXR::psolve(prob,solver = "ECOS")
    
    if (CVXR::status(prob) != "optimal") {
      return(NA_real_)
    }
    
    # Extract duals (lambda)
    lamda <- lapply(constr, sol$getDualValue)
    lamda <- lamda[-3]                    # drop constraint #3 if needed
    lambda_hat <- -unlist(lamda)
    
    W <- as.numeric(CVXR::value(w))
    val <- -(sum(W)) + as.vector(U1 %*% lambda_hat)
    return(val)
  }
  
  # ------------------------------
  # Run optimization over theta
  # ------------------------------
  optim_result <- optim(
    par = theta_init,
    fn = L_fn,
    method = "BFGS",
    control = list(reltol = 1e-3, maxit = 300)
  )
  
  theta_hat = optim_result$par
  
  # residual-based moments
  error <- y.hat - as.numeric(new.data1 %*%   theta_hat)
  H <- new.data1 * error
  A <- H
  r <- ncol(A)
  
  # --- CVXR setup ---
  w <- CVXR::Variable(length(I1), pos = TRUE)
  
  resid <- y - as.matrix(new.data1[I1, ]) %*%  theta_hat
  U_mat1 <- as.matrix(new.data1[I1, ]) * as.vector(resid)
  
  # Constraints: global + per-fold
  constr <- list(
    sum(w) == 1,
    sum(w * g[I1]) == mean(g),
    t(U_mat1) %*% w == rep(0, ncol(U_mat1))
  )
  
  mu <- c()
  for (k in seq_len(K)) {
    idx_k_all <- (ofs[k] + 1):ofs[k + 1]
    idx_k     <- intersect(idx_k_all, I1)
    if (length(idx_k) == 0) next
    
    A_k <- A[idx_k, , drop = FALSE]
    map <- match(idx_k, I1)
    w_k <- w[map]
    A_k_all <- A[idx_k_all, , drop = FALSE]
    constr <- c(constr, list(
      t(A_k) %*% w_k == matrix(colSums(A_k_all) / nrow(A), ncol = 1)
    ))
    mu <- c(mu, colSums(A_k_all) / nrow(A))
  }
  
  ones <- rep(1, length(I1))
  U1 <- c(1, mean(g), mu)
  
  objective <- CVXR::Minimize(sum(CVXR::kl_div(w, ones)))
  prob <- CVXR::Problem(objective, constr)
  sol  <- CVXR::psolve(prob,solver = "ECOS")
  
  if (CVXR::status(prob) != "optimal") {
    return(NA_real_)
  }
  
  
  
  W <- as.numeric(CVXR::value(w))
  
  # Return estimated theta and value
  
  list(
    theta_hat = optim_result$par,
    value = optim_result$value,
    convergence = optim_result$convergence,
    message = optim_result$message,
    estimate.table=variance_theta_diag(
      y = y,
      X = new.data1[I1, ],
      w_hat = W,
      theta_hat = optim_result$par
    )
  )
}











estimate_theta_nested_kfold_CVXR_nlm <- function(theta_init, data_full, K, seed) {
  
  # --- 1. Prepare data (unchanged) ---
  fold_hat <- k_fold_function(df = data_full, K = K, seed = seed)
  data_all <- do.call(rbind, fold_hat)
  n_k      <- sapply(fold_hat, nrow)
  ofs      <- c(0, cumsum(n_k))
  
  x_cols <- subset(data_all, select = -c(Y, D, pi.hat, ID, y.hat))
  p      <- length(x_cols)
  model.formula <- as.formula(paste("~", paste(colnames(x_cols), collapse = "+")))
  new.data  <- as.data.frame(data_all[, colnames(x_cols), drop = FALSE])
  new.data1 <- model.matrix(model.formula, data = new.data)
  
  y.hat <- data_all$y.hat
  D     <- data_all$D
  I1    <- which(D == 1)
  
  g <- log(data_all$pi.hat)
  y <- data_all[I1, "Y"]
  
  # --- 2. Define objective L(theta) (unchanged) ---
  L_fn <- function(theta) {
    
    error <- y.hat - as.numeric(new.data1 %*% theta)
    H <- new.data1 * error
    A <- H
    r <- ncol(A)
    
    w <- CVXR::Variable(length(I1), pos = TRUE)
    
    resid <- y - as.matrix(new.data1[I1, ]) %*% theta
    U_mat1 <- as.matrix(new.data1[I1, ]) * as.vector(resid)
    
    constr <- list(
      sum(w) == 1,
      sum(w * g[I1]) == mean(g),
      t(U_mat1) %*% w == rep(0, ncol(U_mat1))
    )
    
    mu <- c()
    for (k in seq_len(K)) {
      idx_k_all <- (ofs[k] + 1):ofs[k + 1]
      idx_k     <- intersect(idx_k_all, I1)
      if (length(idx_k) == 0) next
      
      A_k <- A[idx_k, , drop = FALSE]
      map <- match(idx_k, I1)
      w_k <- w[map]
      A_k_all <- A[idx_k_all, , drop = FALSE]
      constr <- c(constr, list(
        t(A_k) %*% w_k == matrix(colSums(A_k_all) / nrow(A), ncol = 1)
      ))
      mu <- c(mu, colSums(A_k_all) / nrow(A))
    }
    
    ones <- rep(1, length(I1))
    U1 <- c(1, mean(g), mu)
    
    objective <- CVXR::Minimize(sum(CVXR::kl_div(w, ones)))
    prob <- CVXR::Problem(objective, constr)
    sol  <- CVXR::psolve(prob)
    
    if (CVXR::status(prob) != "optimal") {
      return(NA_real_)
    }
    
    lamda <- lapply(constr, sol$getDualValue)
    lamda <- lamda[-3]
    lambda_hat <- -unlist(lamda)
    
    W <- as.numeric(CVXR::value(w))
    val <- -(sum(W)) + as.vector(U1 %*% lambda_hat)
    return(val)
  }
  
  # --- 3. Outer optimization replaced with nlminb() ---
  nlm_result <- tryCatch({
    nlminb(
      start = theta_init,
      objective = L_fn,
      control = list(rel.tol = 1e-2, x.tol = 1e-2, eval.max = 300, iter.max = 300)
    )
  }, error = function(e) {
    message("[nlminb] error: ", e$message)
    list(par = theta_init, objective = NA, convergence = 1, message = e$message)
  })
  
  theta_hat <- nlm_result$par
  
  # --- 4. Final CVXR solve for w_hat using estimated theta_hat ---
  error <- y.hat - as.numeric(new.data1 %*% theta_hat)
  H <- new.data1 * error
  A <- H
  r <- ncol(A)
  
  w <- CVXR::Variable(length(I1), pos = TRUE)
  resid <- y - as.matrix(new.data1[I1, ]) %*% theta_hat
  U_mat1 <- as.matrix(new.data1[I1, ]) * as.vector(resid)
  
  constr <- list(
    sum(w) == 1,
    sum(w * g[I1]) == mean(g),
    t(U_mat1) %*% w == rep(0, ncol(U_mat1))
  )
  
  mu <- c()
  for (k in seq_len(K)) {
    idx_k_all <- (ofs[k] + 1):ofs[k + 1]
    idx_k     <- intersect(idx_k_all, I1)
    if (length(idx_k) == 0) next
    
    A_k <- A[idx_k, , drop = FALSE]
    map <- match(idx_k, I1)
    w_k <- w[map]
    A_k_all <- A[idx_k_all, , drop = FALSE]
    constr <- c(constr, list(
      t(A_k) %*% w_k == matrix(colSums(A_k_all) / nrow(A), ncol = 1)
    ))
    mu <- c(mu, colSums(A_k_all) / nrow(A))
  }
  
  ones <- rep(1, length(I1))
  U1 <- c(1, mean(g), mu)
  
  objective <- CVXR::Minimize(sum(CVXR::kl_div(w, ones)))
  prob <- CVXR::Problem(objective, constr)
  sol  <- CVXR::psolve(prob)
  
  if (CVXR::status(prob) != "optimal") {
    message("[Final CVXR] solver did not converge, returning NA for w_hat.")
    W <- rep(NA_real_, length(I1))
  } else {
    W <- as.numeric(CVXR::value(w))
  }
  
  # --- 5. Return ---
  list(
    theta_hat = theta_hat,
    value = nlm_result$objective,
    convergence = nlm_result$convergence,
    message = nlm_result$message,
    w_hat = W,
    estimate.table = variance_theta_diag(
      y = y,
      X = new.data1[I1, ],
      w_hat = W,
      theta_hat = theta_hat
    )
  )
}





























