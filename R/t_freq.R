#' @title Frequentist estimation for BRMA model parameters
#' @description Computes point estimates and confidence intervals for the Bivariate Random-Effects Meta-Analysis (BRMA) model using frequentist methods.
#' @importFrom nleqslv nleqslv
#' @importFrom MASS mvrnorm
#' @param data A list of numeric vectors containing transformed treatment effect estimates
#'   and corresponding standard errors for two endpoints (e.g., potential surrogate and true endpoint).
#'   Each vector must have length \eqn{n}, with elements ordered by study:
#'   \describe{
#'     \item{y1}{Log-transformed effect estimates for endpoint 1.}
#'     \item{se1}{Standard errors corresponding to \code{y1}.}
#'     \item{y2}{Log-transformed effect estimates for endpoint 2.}
#'     \item{se2}{Standard errors corresponding to \code{y2}.}
#'   }
#'   The four vectors must be aligned such that the \eqn{i}th elements of
#'   \code{y1}, \code{se1}, \code{y2}, and \code{se2} all refer to the same study.
#' @param method Character string specifying the frequentist estimation method. Supported options: \code{c("REML" , "ML")}. Default is \code{"REML"}.
#' @param interval.method Character string specifying the method to construct confidence intervals. Supported options: \code{c("normal", "bootstrap", "auto")}.
#'                        Default is \code{"auto"}, which chooses appropriate method based on sample size.
#' @param nleqslv.param A list of control settings passed directly to
#'   \code{\link[nleqslv]{nleqslv}} for solving the system of nonlinear equations
#'   used to obtain point estimates.
#'
#'   The default is:
#'   \preformatted{
#'   list(
#'     method = "Broyden",
#'     control = list(maxit = 200, ftol = 1e-6, xtol = 1e-8)
#'   )
#'   }
#'
#'   The list typically includes:
#'   \describe{
#'     \item{\code{method}}{Character string specifying the numerical method used.
#'       Common options include \code{"Broyden"}, \code{"Newton"}, or \code{"dfp"}.}
#'     \item{\code{control}}{A list of algorithmic control parameters, such as
#'       \code{maxit} (maximum number of iterations),
#'       \code{ftol} (function tolerance), and
#'       \code{xtol} (solution tolerance).}
#'   }
#'
#'   For a complete description of available methods and control options,
#'   see \code{\link[nleqslv]{nleqslv}}.
#' @param bootstrap.times Integer. Number of bootstrap replications if \code{interval.method = "bootstrap"}. Default is 1000.
#' @param pars.start Numeric vector of initial values for iterative parameter estimation. If \code{NULL}, data-driven initial values are used.
#' @param verbose Logical. If \code{TRUE}, shows progress messages. Default is \code{TRUE}.
#' @param alpha Numeric. Significance level for confidence intervals. Default is 0.05.
#' @param seed Optional positive integer used to fix the random number generator when performing bootstrap sampling. If NULL, no seed is set.
#' @return A list with two elements:
#' \describe{
#'   \item{Results}{A data frame containing point estimates and confidence intervals of model parameters.
#'     Output details depend on \code{method} and \code{interval.method}.}
#'   \item{Summary}{A character string summarizing the estimation method and how confidence intervals were computed.}
#' }
#' @examples
#' data <- list(y1 = c(0.198, -0.091, 0.345, 1.193, 0.828, 1.671, 0.152, 1.416),
#'              se1 = c(0.076, 0.081, 0.062, 0.387, 0.497, 0.195, 0.432, 0.254),
#'              y2 = c(-0.234, -0.614, -0.054, -0.670, 0.105, 1.045, -0.910, 0.485),
#'              se2 = c(0.465, 0.408, 0.240, 0.384, 0.239, 0.434, 0.141, 0.159))
#' t_freq(data) # parameters estimated by REML;
#'              # confidence intervals constructed by bootstrap (sample size <10)
#' t_freq(data, method = "ML", interval.method = "normal", verbose = TRUE)
#' @export

t_freq <- function(data,
                   method = c("REML", "ML"),
                   interval.method = c("auto", "normal", "bootstrap"),
                   bootstrap.times = 1000,
                   pars.start = NULL,
                   nleqslv.param = list(
                     method = "Broyden",
                     control = list(maxit = 200,
                                    ftol = 1e-6,
                                    xtol = 1e-8)),
                   verbose = TRUE,
                   alpha = 0.05,
                   seed = NULL) {
  # Input checks
  stopifnot(is.list(data))
  required <- c("y1", "y2", "se1", "se2")
  if (!all(required %in% names(data))) stop("data must contain y1, y2, se1, se2")
  if (length(data$y1) != length(data$y2)) stop("y1 and y2 must have the same length")

  n <- length(data$y1)
  method <- match.arg(method)
  interval.method <- match.arg(interval.method)

  # Determine interval method
  if (interval.method == "auto") {
    interval.method <- if (n < 10) "bootstrap" else "normal"
    if (verbose) message(sprintf("Sample size = %d -> using %s CI by default.", n, interval.method))
  }

  # CI name
  ci_label <- sprintf("%d", round((1 - alpha) * 100))
  ci_colname <- if (interval.method == "bootstrap") {
    paste0("bootstrapped ", ci_label, "% CI")
  } else {
    paste0("normal ", ci_label, "% CI")
  }
  param.names <- c("beta1", "beta2", "psi1_2", "psi2_2", "rho")
  if(!is.null(seed)){set.seed(seed)}
  # Run selected method
  if (method == "ML") {
    if (verbose) message("Running ML estimation...")
    old_warn <- getOption("warn"); options(warn = -1)
    mlsolve.res <- ml_solve(data,
                            pars.start = pars.start,
                            verbose = verbose,
                            nleqslv.param = nleqslv.param)
    if (is.null(mlsolve.res)) return(NULL)
    ml.root <- mlsolve.res$x
    # Compute Hessian & SE
    ml.Htrans_Htheta <- ml_hessian(data, ml.root)
    ml.setrans_setheta <- freq_se(ml.Htrans_Htheta)
    if(sum(sapply(ml.setrans_setheta, function(x) is.nan(x)))>0) return(NULL)
    # Compute normal CI
    ml.est_normal.ci <- freq_normalCI(ml.root, ml.setrans_setheta$se.trans, alpha = alpha)
    df.mlnormal <- data.frame(
      Parameter = param.names,
      Estimates = ml.est_normal.ci$Estimates,
      SE = ml.setrans_setheta$se.orig
    )
    df.mlnormal[[ci_colname]] <- ml.est_normal.ci$CI
    options(warn = old_warn)
    if (interval.method == "bootstrap") {
      if (verbose) message(sprintf("Running ML bootstrap (%d resamples)...", bootstrap.times))
      mlboot.res <- ml_paramBoot(data,
                                 ml.results = df.mlnormal,
                                 bootstrap.times = bootstrap.times,
                                 pars.start = pars.start,
                                 nleqslv.param = nleqslv.param,
                                 alpha = alpha)
      df.mlboot <- data.frame(
        Parameter = param.names,
        Estimates = ml.est_normal.ci$Estimates
      )
      df.mlboot[[ci_colname]] <- sprintf("(%.3f, %.3f)", mlboot.res$bootCI_lower, mlboot.res$bootCI_upper)
      ML <- df.mlboot
    } else if (interval.method == "normal") {
      if (verbose) message("Computing ML normal CIs...")
      ML <- df.mlnormal
    }
    ML = ML[c(5,1:4),];row.names(ML)<-1:5
    text_summary <- paste0(
      "Point estimates obtained using Maximum Likelihood (ML). ",
      sprintf("%s%% confidence intervals computed using the %s method%s.",
              round((1 - alpha) * 100),
              interval.method,
              if (interval.method == "bootstrap") sprintf(" with %d bootstrap replicates", bootstrap.times) else "")
    )
    return(list(Results = ML, Summary = text_summary))
  }
  if (method == "REML") {
    if (verbose) message("Running REML estimation...")
    old_warn <- getOption("warn"); options(warn = -1)
    remlsolve.res <- reml_solve(data,
                                pars.start = pars.start,
                                verbose = verbose,
                                nleqslv.param = nleqslv.param)
    if (is.null(remlsolve.res)) return(NULL)
    reml.root <- remlsolve.res$x
    reml.beta <- solve_remlBeta(reml.root, data)
    # Compute Hessian & SE
    reml.Htrans_Htheta <- reml_hessian(data, reml.root)
    reml.setrans_setheta <- freq_se(reml.Htrans_Htheta)
    # Compute normal CI
    reml.est_normal.ci <- freq_normalCI(c(reml.beta$beta_hat, reml.root),
                                    c(reml.beta$se_beta_hat, reml.setrans_setheta$se.trans), alpha = alpha)

    df.remlnormal <- data.frame(
      Parameter = param.names,
      Estimates = c(reml.est_normal.ci$Estimates),
      SE = c(reml.beta$se_beta_hat, reml.setrans_setheta$se.orig)
    )
    df.remlnormal[[ci_colname]] <- reml.est_normal.ci$CI
    options(warn = old_warn)
    if (interval.method == "bootstrap") {
      if (verbose) message(sprintf("Running REML bootstrap (%d resamples)...", bootstrap.times))
      remlboot.res <- reml_paramBoot(data,
                                     reml.results = df.remlnormal,
                                     bootstrap.times = bootstrap.times,
                                     pars.start = pars.start,
                                     nleqslv.param = nleqslv.param,
                                     alpha = alpha)
      df.remlboot <- data.frame(
        Parameter = param.names,
        Estimates = reml.est_normal.ci$Estimates
      )
      df.remlboot[[ci_colname]] <- sprintf("(%.3f, %.3f)", remlboot.res$bootCI_lower, remlboot.res$bootCI_upper)
      REML <- df.remlboot
    } else if (interval.method == "normal") {
      if (verbose) message("Computing REML normal CIs...")
      REML <- df.remlnormal
    }
    REML = REML[c(5,1:4),];row.names(REML)<-1:5
    text_summary <- paste0(
      "Point estimates obtained using Restricted Maximum Likelihood (REML). ",
      sprintf("%s%% confidence intervals computed using the %s method%s.",
              round((1 - alpha) * 100),
              interval.method,
              if (interval.method == "bootstrap") sprintf(" with %d bootstrap replicates", bootstrap.times) else "")
    )
    return(list(Results = REML, Summary = text_summary))
  }
}


### Not export
generate_data <- function(n, beta1, beta2, rho, psi1_2, psi2_2, se1, se2) {
  stopifnot(length(se1) == n, length(se2) == n)
  y1 <- numeric(n)
  y2 <- numeric(n)
  X <- matrix(1, nrow = n, ncol = 1)
  s1_2 <- se1^2
  s2_2 <- se2^2
  for (i in 1:n) {
    mu_i1 <- sum(X[i, ] * beta1)
    mu_i2 <- sum(X[i, ] * beta2)

    Sigma_i <- matrix(c(
      psi1_2 + s1_2[i],
      rho * sqrt((psi1_2 + s1_2[i]) * (psi2_2 + s2_2[i])),
      rho * sqrt((psi1_2 + s1_2[i]) * (psi2_2 + s2_2[i])),
      psi2_2 + s2_2[i]
    ), nrow = 2)

    obs <- MASS::mvrnorm(n = 1, mu = c(mu_i1, mu_i2), Sigma = Sigma_i)

    y1[i] <- obs[1]
    y2[i] <- obs[2]
  }

  return(list(y1 = y1, y2 = y2, se1 = se1, se2 = se2))
}

## ML
ml_scores <- function(data, x) {
  p <- 1
  # Parameters
  beta1 <- x[1:p]
  beta2 <- x[(p+1):(2*p)]
  psi1_2 <- exp(x[2*p + 1]) + 1e-3
  psi2_2 <- exp(x[2*p + 2]) + 1e-3
  rho <- (2 * plogis(x[2*p + 3])) - 1  # inverse logit scaled to (-1, 1)

  n <- length(data$y1)
  data$X <- matrix(1, nrow = n, ncol = 1)

  y1 <- data$y1
  y2 <- data$y2
  s1_2 <- (data$se1)^2
  s2_2 <- (data$se2)^2
  X1 <- X2 <- X <- data$X

  mu1 <- X1 %*% beta1
  mu2 <- X2 %*% beta2
  e1 <- y1 - mu1
  e2 <- y2 - mu2
  v1 <- psi1_2 + s1_2
  v2 <- psi2_2 + s2_2
  sqrt_v1 <- sqrt(v1)
  sqrt_v2 <- sqrt(v2)

  # Scores
  score1 <- colSums(1/(1 - rho^2) * as.vector(
    (e1 / v1) - rho * (e2 / (sqrt_v1 * sqrt_v2))
  ) * X1)

  score2 <- colSums(1/(1 - rho^2) * as.vector(
    (e2 / v2) - rho * (e1 / (sqrt_v1 * sqrt_v2))
  ) * X2)

  score3 <- sum(
    -1 / (2 * v1) + 1 / (2 * (1 - rho^2)) * (
      (e1^2) / (v1^2) - rho * e1 * e2 / (v1^(3/2) * v2^(1/2))
    )
  )

  score4 <- sum(
    -1 / (2 * v2) + 1 / (2 * (1 - rho^2)) * (
      (e2^2) / (v2^2) - rho * e1 * e2 / (v1^(1/2) * v2^(3/2))
    )
  )

  score5 <- n * rho / (1 - rho^2) -
    rho / (1 - rho^2)^2 * sum(e1^2 / v1) -
    rho / (1 - rho^2)^2 * sum(e2^2 / v2) +
    (rho^2 + 1) / (1 - rho^2)^2 * sum(e1 * e2 / (sqrt_v1 * sqrt_v2))
  d_logpsi1 <- psi1_2 * score3      # *d/d log(psi1^2) by chain rule
  d_logpsi2 <- psi2_2 * score4      # *d/d log(psi2^2) by chain rule
  d_zrho    <- ((1 - rho^2) / 2) * score5   # *d/d z_rho by chain rule
  res <- c(score1, score2, d_logpsi1, d_logpsi2, d_zrho)
  return(res)
}

##### Hessian
ml_hessian <- function(data, ml.root){
  n <- length(data$y1)
  p <- 1

  # Parameters
  beta1.hat <- ml.root[1:p]
  beta2.hat <- ml.root[(p+1):(2*p)]
  psi1_2.hat <- exp(ml.root[2*p + 1])+1e-3
  psi2_2.hat <- exp(ml.root[2*p + 2])+1e-3
  rho.hat <- (2 * plogis(ml.root[2*p + 3])) - 1  # inverse logit scaled to (-1, 1)

  data$X <- matrix(1, nrow = n, ncol = 1)

  y1 <- data$y1
  y2 <- data$y2
  s1_2 <- (data$se1)^2
  s2_2 <- (data$se2)^2
  X1 <- X2 <- X <- data$X

  xtx <- list()
  for (i in 1:n) { xtx[[i]] <- (X[i,]) %*% t(X[i,]) }

  mu1 <- as.vector(X1 %*% beta1.hat)
  mu2 <- as.vector(X2 %*% beta2.hat)
  e1 <- y1 - mu1
  e2 <- y2 - mu2
  v1 <- psi1_2.hat + s1_2
  v2 <- psi2_2.hat + s2_2
  sqrt_v1 <- sqrt(v1)
  sqrt_v2 <- sqrt(v2)

  w11 <- -1 / (1 - rho.hat^2) / v1
  h11 <- Reduce("+", Map("*", w11, xtx))
  w12 <-  rho.hat * 1 / (1 - rho.hat^2) / (sqrt(v1) * sqrt(v2))
  h12 <- Reduce("+", Map("*", w12, xtx))
  h13 <- colSums((1/(1-rho.hat^2))*(-(e1)/(v1)^2 + rho.hat*(e2)/2/(v1)^(3/2)/(v2)^(1/2)) * X)
  h14 <- colSums((1/(1-rho.hat^2))*(rho.hat*(e2)/2/(v1)^(1/2)/(v2)^(3/2)) * X)
  h15 <- colSums((2*rho.hat*(e1)/(1-rho.hat^2)^2/(v1) - (1+rho.hat^2)*(e2)/(1-rho.hat^2)^2/sqrt_v1/sqrt_v2)*X)
  w22 <- -1 / (1 - rho.hat^2) / v2
  h22 <- Reduce("+", Map("*", w22, xtx))
  h23 <- (1/(1-rho.hat^2))*colSums(rho.hat*(e1)/2/(v1)^(3/2)/(v2)^(1/2) * X)
  h24 <- colSums((1/(1-rho.hat^2))*(-(e2)/(v2)^2 + rho.hat*(e1)/2/(v1)^(1/2)/(v2)^(3/2) * X))

  h25 <- colSums((2*rho.hat*(e2)/(1-rho.hat^2)^2/(v2) - (1+rho.hat^2)*(e1)/(1-rho.hat^2)^2/sqrt_v1/sqrt_v2)*X)
  h33 <- sum(1/2/(v1)^2)-sum((e1)^2/(1-rho.hat^2)/(v1)^3) + sum(3*rho.hat*(e1)*(e2)/4/(1-rho.hat^2)/(v1)^(5/2)/(v2)^(1/2))
  h34 <- sum(rho.hat*(e1)*(e2)/4/(1-rho.hat^2)/(v1)^(3/2)/(v2)^(3/2))
  h35 <- sum(rho.hat*(e1)^2/(1-rho.hat^2)^2/(v1)^2) - sum((1+rho.hat^2)*(e1)*(e2)/2/(1-rho.hat^2)^2/(v1)^(3/2)/(v2)^(1/2))
  h44 <- sum(1/2/(v2)^2)-sum((e2)^2/(1-rho.hat^2)/(v2)^3) + sum(3*rho.hat*(e1)*(e2)/4/(1-rho.hat^2)/(v1)^(1/2)/(v2)^(5/2))
  h45 <- sum(rho.hat*(e2)^2/(1-rho.hat^2)^2/(v2)^2) - sum((1+rho.hat^2)*(e1)*(e2)/2/(1-rho.hat^2)^2/(v1)^(1/2)/(v2)^(3/2))
  h55 <- n*(1+rho.hat^2)/(1-rho.hat^2)^2-(1+3*rho.hat^2)/(1-rho.hat^2)^3*sum((e1)^2/(v1)) -
    (1+3*rho.hat^2)/(1-rho.hat^2)^3*sum((e2)^2/(v2)) +
    2*rho.hat*(rho.hat^2+3)/(1-rho.hat^2)^3*sum((e1)*(e2)/sqrt((v1)*(v2)))
  H <- rbind(cbind(h11, h12, h13, h14, h15),
             cbind(h12, h22, h23, h24, h25),
             c(h13, h23, h33, h34, h35),
             c(h14, h24, h34, h44, h45),
             c(h15, h25, h35, h45, h55))
  J <- diag(c(rep(1, p*2), psi1_2.hat, psi2_2.hat, (1-rho.hat^2)/2 ))
  H.trans <- t(J)%*%H%*%J
  return(list(H.orig = H,
              H.trans = H.trans)) # covariance matrix in the transformed scale
}
freq_se <- function(hessians, ridge = 1e-6) {
  H.orig.stable <- hessians$H.orig - diag(ridge, nrow(hessians$H.orig))
  H.trans.stable <- hessians$H.trans - diag(ridge, nrow(hessians$H.trans))

  se.original <- sqrt(diag(solve(-H.orig.stable)))
  se.trans <- sqrt(diag(solve(-H.trans.stable)))

  list(se.orig = se.original, se.trans = se.trans)
}

freq_normalCI <- function(root, se.trans, alpha) {
  p <- 1
  z <- qnorm(1 - alpha/2)

  # Parameter estimates
  beta1 <- root[1]
  beta2 <- root[2]
  logpsi1 <- root[3]
  logpsi2 <- root[4]
  z_rho <- root[5]

  # CI in transformed scale
  lower <- c(beta1, beta2, logpsi1, logpsi2, z_rho) - z*se.trans
  upper <- c(beta1, beta2, logpsi1, logpsi2, z_rho) + z*se.trans

  # Back-transform
  lower[3:4] <- exp(lower[3:4])
  upper[3:4] <- exp(upper[3:4])
  lower[5] <- 2*plogis(lower[5]) - 1
  upper[5] <- 2*plogis(upper[5]) - 1

  CI <- paste0("(", round(lower,3), ", ", round(upper,3), ")")

  Estimates <- c(beta1, beta2, exp(logpsi1), exp(logpsi2), 2*plogis(z_rho)-1)

  return(list(Estimates = Estimates, CI = CI))
}

ml_solve <- function(data, pars.start,
                     verbose, nleqslv.param) {
  # Default starting values
  if (is.null(pars.start)) {
    # Estimate initial psi variances by subtracting mean squared measurement error
    v1_init <- var(data$y1) - mean(data$se1^2)
    v2_init <- var(data$y2) - mean(data$se2^2)

    # Make sure variances are positive
    v1_init <- max(v1_init, 1e-3)
    v2_init <- max(v2_init, 1e-3)

    # Initial values for root solver
    pars.start <- c(
      mean(data$y1),   # beta1
      mean(data$y2),   # beta2
      log(v1_init),    # log(psi1^2)
      log(v2_init),    # log(psi2^2)
      0                # rho
    )

  }

  if (!is.numeric(pars.start) || length(pars.start) != 5) {
    stop("pars.start must be a numeric vector of length 5: (beta1, beta2, (log(psi1), log(psi2), tanh(rho/2))")
  }

  # Solve equations
  res <- tryCatch({
    nleqslv::nleqslv(
      x = pars.start,
      fn = ml_scores,
      data = data,
      method = nleqslv.param$method,
      global = nleqslv.param$global,
      control = nleqslv.param$control
    )
  }, error = function(e) {
    if (verbose) message("nleqslv failed: ", e$message)
    return(NULL)
  })

  # --- check convergence safely ---
  if (!is.null(res)) {
    converged <- res$termcd %in% c(1, 2) &&
      all(is.finite(res$fvec)) &&
      max(abs(res$fvec)) < 1e-5
  } else {
    converged <- FALSE
  }

  # --- verbose reporting ---
  if (verbose) {
    if (is.null(res)) {
      message("ml_solve(): error when applying nleqslv.")
    } else if (!converged) {
      msg <- sprintf(
        "ml_solve(): non-convergence [termcd=%d]. Message: %s",
        res$termcd, res$message
      )
      msg <- paste0(msg, ". May consider adjusting the initial values.")
      message(msg)
    } else {
      msg <- sprintf(
        "ml_solve(): converged successfully [termcd=%d].",
        res$termcd)
      message(msg)
    }
  }
  if (!converged) {
    return(NULL)
  }
  return(res)
}
ml_paramBoot <- function(data, ml.results,
                         bootstrap.times,
                         pars.start,
                         nleqslv.param,
                         alpha) {
  n <- length(data$y1)
  p <- 1
  if(is.null(pars.start)){
    v1_init <- var(data$y1) - mean(data$se1^2)
    v2_init <- var(data$y2) - mean(data$se2^2)
    v1_init <- max(v1_init, 1e-3)
    v2_init <- max(v2_init, 1e-3)
    pars.start <- c(
      mean(data$y1),   # beta1
      mean(data$y2),   # beta2
      log(v1_init),    # log(psi1^2)
      log(v2_init),    # log(psi2^2)
      0                # rho
    )
  }

  beta1.est <- ml.results$Estimates[grepl("^beta1", ml.results$Parameter)]
  beta2.est <- ml.results$Estimates[grepl("^beta2", ml.results$Parameter)]
  psi1_2.est <- ml.results$Estimates[grepl("^psi1_2$", ml.results$Parameter)]
  psi2_2.est <- ml.results$Estimates[grepl("^psi2_2$", ml.results$Parameter)]
  rho.est <- ml.results$Estimates[grepl("^rho$", ml.results$Parameter)]

  se1 <- (data$se1)
  se2 <- (data$se2)
  bootstrap.results <- list()
  n.success <- 0
  n.attempts <- 0

  while (n.success < bootstrap.times) {
    n.attempts <- n.attempts + 1

    boot.data <- generate_data(n = n,
                               beta1 = beta1.est, beta2 = beta2.est, rho = rho.est,
                               psi1_2 = psi1_2.est, psi2_2 = psi2_2.est,
                               se1 = se1, se2 = se2)

    result <- tryCatch({
      suppressWarnings({
        invisible(capture.output({
          out <- nleqslv::nleqslv(
            x = pars.start,
            fn = ml_scores,
            data = boot.data,
            method = nleqslv.param$method,
            global = nleqslv.param$global,
            control = nleqslv.param$control
          )
        }, type = "output"))

        root <- out$x
        if (!is.null(root) &&
            all(out$fvec < 1e-6) &&
            all(is.finite(root)) &&
            all(abs(root) < 10) &&
            out$termcd %in% c(1, 2)) {
          root
        } else {
          NULL
        }
      })
    }, error = function(e) NULL)

    if (!is.null(result)) {
      n.success <- n.success + 1
      bootstrap.results[[n.success]] <- result
    }
  }

  bootstrap.matrix <- do.call(rbind, bootstrap.results)

  beta1.mat <- bootstrap.matrix[, 1:p, drop = FALSE]
  beta2.mat <- bootstrap.matrix[, (p + 1):(2 * p), drop = FALSE]
  psi1_2.vec <- exp(bootstrap.matrix[, 2 * p + 1])
  psi2_2.vec <- exp(bootstrap.matrix[, 2 * p + 2])
  rho.vec    <- 2 * plogis(bootstrap.matrix[, 2 * p + 3]) - 1

  bootstrap.transformed <- cbind(beta1.mat, beta2.mat, psi1_2 = psi1_2.vec, psi2_2 = psi2_2.vec, rho = rho.vec)

  bootstrap.mean <- colMeans(bootstrap.transformed)
  bootstrap.se <- apply(bootstrap.transformed, 2, sd)
  bootstrap.ci <- t(apply(bootstrap.transformed, 2, quantile, probs = c(alpha/2, 1 - alpha/2)))
  colnames(bootstrap.ci) <- c("CI_lower", "CI_upper")

  param.names <- c(paste0("beta1", 1:p),
                   paste0("beta2", 1:p),
                   "psi1_2", "psi2_2", "rho")

  data.frame(
    Parameter = param.names,
    bootMean = bootstrap.mean,
    bootSE = bootstrap.se,
    bootCI_lower = bootstrap.ci[, 1],
    bootCI_upper = bootstrap.ci[, 2]
  )
}

## REML
reml_scores <- function(data, x) {
  n <- length(data$y1)

  # Transform parameters
  psi1_2 <- exp(x[1]) + 1e-3
  psi2_2 <- exp(x[2]) + 1e-3
  rho <- 2 * plogis(x[3]) - 1  # scale to (-1,1)

  # Data
  y1 <- data$y1
  y2 <- data$y2
  s1_2 <- data$se1^2
  s2_2 <- data$se2^2
  X <- matrix(1, nrow = n, ncol = 1)  # scalar intercept per i

  # Initialize sums for beta
  sum_Xt_PhiInv_X <- matrix(0, nrow = 2, ncol = 2)
  sum_Xt_PhiInv_y <- matrix(0, nrow = 2, ncol = 1)

  # Store Phi_i, Phi_inv, derivatives for later
  Phi_list <- vector("list", n)
  Phi_inv_list <- vector("list", n)
  dPhi1_list <- vector("list", n)
  dPhi2_list <- vector("list", n)
  dPhiR_list <- vector("list", n)

  for (i in 1:n) {
    psi1s <- psi1_2 + s1_2[i]
    psi2s <- psi2_2 + s2_2[i]
    sqrt_psi1s <- sqrt(psi1s)
    sqrt_psi2s <- sqrt(psi2s)

    # Phi_i
    Phi_i <- matrix(c(psi1s, rho * sqrt_psi1s * sqrt_psi2s,
                      rho * sqrt_psi1s * sqrt_psi2s, psi2s), nrow = 2)
    Phi_inv <- solve(Phi_i)

    # Store
    Phi_list[[i]] <- Phi_i
    Phi_inv_list[[i]] <- Phi_inv

    # Derivatives
    dPhi1_list[[i]] <- matrix(c(
      1, rho * psi2s / (2 * sqrt_psi1s * sqrt_psi2s),
      rho * psi2s / (2 * sqrt_psi1s * sqrt_psi2s), 0
    ), 2, 2)

    dPhi2_list[[i]] <- matrix(c(
      0, rho * psi1s / (2 * sqrt_psi1s * sqrt_psi2s),
      rho * psi1s / (2 * sqrt_psi1s * sqrt_psi2s), 1
    ), 2, 2)

    dPhiR_list[[i]] <- matrix(c(
      0, sqrt_psi1s * sqrt_psi2s,
      sqrt_psi1s * sqrt_psi2s, 0
    ), 2, 2)

    # Contributions to Xt_PhiInv_X and Xt_PhiInv_y
    X_i <- t(kronecker(diag(2), X[i, ]))  # intercept
    y_i <- matrix(c(y1[i], y2[i]), nrow = 2)
    sum_Xt_PhiInv_X <- sum_Xt_PhiInv_X + t(X_i) %*% Phi_inv %*% X_i
    sum_Xt_PhiInv_y <- sum_Xt_PhiInv_y + t(X_i) %*% Phi_inv %*% y_i
  }

  # Estimate beta_hat
  beta_hat <- solve(sum_Xt_PhiInv_X) %*% sum_Xt_PhiInv_y

  # Initialize scores
  score_psi1 <- 0
  score_psi2 <- 0
  score_rho  <- 0

  # Compute score contributions per i
  for (i in 1:n) {
    X_i <- t(kronecker(diag(2), X[i, ]))
    y_i <- matrix(c(y1[i], y2[i]), nrow = 2)
    Phi_i <- Phi_list[[i]]
    Phi_inv <- Phi_inv_list[[i]]
    dPhi1 <- dPhi1_list[[i]]
    dPhi2 <- dPhi2_list[[i]]
    dPhiR <- dPhiR_list[[i]]

    resid_i <- y_i - X_i %*% beta_hat

    # term1: trace(Phi_i^{-1} dPhi_i/dtheta)
    term1_psi1 <- sum(diag(Phi_inv %*% dPhi1))
    term1_psi2 <- sum(diag(Phi_inv %*% dPhi2))
    term1_rho  <- sum(diag(Phi_inv %*% dPhiR))

    # term2: trace((X'Phi^{-1}X)^{-1} X' Phi^{-1} dPhi Phi^{-1} X)
    term2_psi1 <- sum(diag(solve(sum_Xt_PhiInv_X) %*% t(X_i) %*% Phi_inv %*% dPhi1 %*% Phi_inv %*% X_i))
    term2_psi2 <- sum(diag(solve(sum_Xt_PhiInv_X) %*% t(X_i) %*% Phi_inv %*% dPhi2 %*% Phi_inv %*% X_i))
    term2_rho  <- sum(diag(solve(sum_Xt_PhiInv_X) %*% t(X_i) %*% Phi_inv %*% dPhiR %*% Phi_inv %*% X_i))

    # term3: resid' Phi^{-1} dPhi Phi^{-1} resid
    term3_psi1 <- t(resid_i) %*% Phi_inv %*% dPhi1 %*% Phi_inv %*% resid_i
    term3_psi2 <- t(resid_i) %*% Phi_inv %*% dPhi2 %*% Phi_inv %*% resid_i
    term3_rho  <- t(resid_i) %*% Phi_inv %*% dPhiR %*% Phi_inv %*% resid_i

    # add contribution to scores
    score_psi1 <- score_psi1 - 0.5 * (term1_psi1 - term2_psi1 - term3_psi1)
    score_psi2 <- score_psi2 - 0.5 * (term1_psi2 - term2_psi2 - term3_psi2)
    score_rho  <- score_rho  - 0.5 * (term1_rho  - term2_rho  - term3_rho)
  }

  # chain rule for transformed parameters
  d_logpsi1 <- psi1_2 * score_psi1
  d_logpsi2 <- psi2_2 * score_psi2
  d_zrho    <- 0.5 * (1 - rho^2) * score_rho

  return(c(d_logpsi1, d_logpsi2, d_zrho))
}

reml_hessian <- function(data, reml.root) {
  # Extract data
  y1 <- data$y1
  y2 <- data$y2
  s1_2 <- data$se1^2
  s2_2 <- data$se2^2
  n <- length(y1)

  # Design matrix (intercept-only)
  X <- matrix(1, n, 1)
  p <- ncol(X)

  # Parameters
  psi1_2 <- exp(reml.root[1])+1e-3
  psi2_2 <- exp(reml.root[2])+1e-3
  rho <- (2 * plogis(reml.root[3])) - 1
  beta_hat <- solve_remlBeta(reml.root = reml.root,
                             data = data)$beta_hat

  compute_matrices.hessian <- function(i) {
    # Subject-specific variances
    s1 <- s1_2[i]
    s2 <- s2_2[i]
    psi1s <- psi1_2 + s1
    psi2s <- psi2_2 + s2
    sqrt_psi1s <- sqrt(psi1s)
    sqrt_psi2s <- sqrt(psi2s)

    # Covariance matrix
    Phi <- matrix(c(psi1s, rho*sqrt_psi1s*sqrt_psi2s,
                    rho*sqrt_psi1s*sqrt_psi2s, psi2s), nrow = 2)
    Phi.inv <- solve(Phi)

    # First derivatives
    Phi.psi1 <- matrix(c(1, 0.5*rho*sqrt_psi2s/sqrt_psi1s,
                         0.5*rho*sqrt_psi2s/sqrt_psi1s, 0), 2)
    Phi.psi2 <- matrix(c(0, 0.5*rho*sqrt_psi1s/sqrt_psi2s,
                         0.5*rho*sqrt_psi1s/sqrt_psi2s, 1), 2)
    Phi.rho  <- matrix(c(0, sqrt_psi1s*sqrt_psi2s,
                         sqrt_psi1s*sqrt_psi2s, 0), 2)

    # Second derivatives
    Phi.psi1psi1 <- matrix(c(0, -0.25*rho*sqrt_psi2s/(sqrt_psi1s^3),
                             -0.25*rho*sqrt_psi2s/(sqrt_psi1s^3), 0), 2)
    Phi.psi2psi2 <- matrix(c(0, -0.25*rho*sqrt_psi1s/(sqrt_psi2s^3),
                             -0.25*rho*sqrt_psi1s/(sqrt_psi2s^3), 0), 2)
    Phi.psi1psi2 <- matrix(c(0, 0.25*rho/(sqrt_psi1s*sqrt_psi2s),
                             0.25*rho/(sqrt_psi1s*sqrt_psi2s), 0), 2)
    Phi.psi1rho  <- matrix(c(0, 0.5*sqrt_psi2s/sqrt_psi1s,
                             0.5*sqrt_psi2s/sqrt_psi1s, 0), 2)
    Phi.psi2rho  <- matrix(c(0, 0.5*sqrt_psi1s/sqrt_psi2s,
                             0.5*sqrt_psi1s/sqrt_psi2s, 0), 2)
    Phi.rhorho   <- matrix(0, 2, 2)

    # Stack y vector
    y <- matrix(c(y1[i], y2[i]), nrow = 2)
    x <- X[i,,drop=FALSE]  # p x 1
    x <- kronecker(diag(2), x)

    # Precompute repeated terms for Hessian
    xPhiInv <- t(x) %*% Phi.inv
    yres <- y - x %*% matrix(beta_hat, ncol=1)

    list(
      x = x,
      y = y,
      Phi = Phi,
      Phi.inv = Phi.inv,
      Phi.derivatives = list(psi1 = Phi.psi1, psi2 = Phi.psi2, rho = Phi.rho),
      Phi.second = list(
        psi1psi1 = Phi.psi1psi1, psi2psi2 = Phi.psi2psi2, psi1psi2 = Phi.psi1psi2,
        psi1rho = Phi.psi1rho, psi2rho = Phi.psi2rho, rhorho = Phi.rhorho
      ),
      xPhiInv = xPhiInv,
      yres = yres
    )
  }

  compute_matrices <- lapply(1:n, compute_matrices.hessian)
  sum_xPhiInvX <- matrix(0, p*2, p*2)
  sum_xPhiInvY <- matrix(0, p*2, 1)
  sum_xPhipsi1Phix <- matrix(0, p*2, p*2)
  sum_xPhipsi2Phix <- matrix(0, p*2, p*2)
  sum_xPhirhoPhix <- matrix(0, p*2, p*2)
  sum_dbetapsi1term2 <- 0
  sum_dbetapsi2term2 <- 0
  sum_dbetarhoterm2 <- 0
  for(i in 1:n) {
    res <- compute_matrices[[i]]
    sum_xPhiInvX <- sum_xPhiInvX + t(res$x) %*% res$Phi.inv %*% res$x
    sum_xPhiInvY <- sum_xPhiInvY + t(res$x) %*% res$Phi.inv %*% res$y
    sum_xPhipsi1Phix <- sum_xPhipsi1Phix + t(res$x) %*% res$Phi.inv %*% res$Phi.derivatives$psi1 %*% res$Phi.inv %*% res$x
    sum_xPhipsi2Phix <- sum_xPhipsi2Phix + t(res$x) %*% res$Phi.inv %*% res$Phi.derivatives$psi2 %*% res$Phi.inv %*% res$x
    sum_xPhirhoPhix <- sum_xPhirhoPhix + t(res$x) %*% res$Phi.inv %*% res$Phi.derivatives$rho %*% res$Phi.inv %*% res$x
    sum_dbetapsi1term2 <- sum_dbetapsi1term2 + t(res$x) %*% res$Phi.inv %*% res$Phi.derivatives$psi1 %*% res$Phi.inv %*% res$yres
    sum_dbetapsi2term2 <- sum_dbetapsi2term2 + t(res$x) %*% res$Phi.inv %*% res$Phi.derivatives$psi2 %*% res$Phi.inv %*% res$yres
    sum_dbetarhoterm2 <- sum_dbetarhoterm2 + t(res$x) %*% res$Phi.inv %*% res$Phi.derivatives$rho %*% res$Phi.inv %*% res$yres
  }
  dbeta_dpsi1 <- -solve(sum_xPhiInvX) %*% sum_dbetapsi1term2
  dbeta_dpsi2 <- -solve(sum_xPhiInvX) %*% sum_dbetapsi2term2
  dbeta_drho <- -solve(sum_xPhiInvX) %*% sum_dbetarhoterm2

  Hpsi1psi1 <- 0
  Hpsi1psi2 <- 0
  Hpsi1rho <- 0
  Hpsi2psi2 <- 0
  Hpsi2rho <- 0
  Hrhorho <- 0
  for (i in 1:n) {
    res <- compute_matrices[[i]]
    Hpsi1psi1 <- Hpsi1psi1 - 0.5*(sum(diag(res$Phi.inv%*%res$Phi.second$psi1psi1)) - sum(diag(res$Phi.inv%*%res$Phi.derivatives$psi1%*%res$Phi.inv%*%res$Phi.derivatives$psi1))-
                                    sum(diag(solve(sum_xPhiInvX)%*%sum_xPhipsi1Phix%*%solve(sum_xPhiInvX)%*%res$x%*%res$Phi.inv%*%res$Phi.derivatives$psi1%*%res$Phi.inv%*%res$x))-
                                    sum(diag(solve(sum_xPhiInvX)%*%res$x%*%res$Phi.inv%*%res$Phi.second$psi1psi1%*%res$Phi.inv%*%res$x))+
                                    sum(diag(solve(sum_xPhiInvX)%*%res$x%*%res$Phi.inv%*%res$Phi.derivatives$psi1%*%res$Phi.inv%*%res$Phi.derivatives$psi1%*%res$Phi.inv%*%res$x))+
                                    sum(diag(solve(sum_xPhiInvX)%*%res$x%*%res$Phi.inv%*%res$Phi.derivatives$psi1%*%res$Phi.inv%*%res$Phi.derivatives$psi1%*%res$Phi.inv%*%res$x))+
                                    2*t(res$yres)%*%res$Phi.inv%*%res$Phi.derivatives$psi1%*%res$Phi.inv%*%res$x%*%dbeta_dpsi1-
                                    t(res$yres)%*%(-res$Phi.inv%*%res$Phi.derivatives$psi1%*%res$Phi.inv%*%res$Phi.derivatives$psi1%*%res$Phi.inv-
                                                     res$Phi.inv%*%res$Phi.derivatives$psi1%*%res$Phi.inv%*%res$Phi.derivatives$psi1%*%res$Phi.inv+
                                                     res$Phi.inv%*%res$Phi.second$psi1psi1%*%res$Phi.inv)%*%res$yres)
    Hpsi1psi2 <- Hpsi1psi2 - 0.5*(sum(diag(res$Phi.inv%*%res$Phi.second$psi1psi2)) - sum(diag(res$Phi.inv%*%res$Phi.derivatives$psi1%*%res$Phi.inv%*%res$Phi.derivatives$psi2))-
                                    sum(diag(solve(sum_xPhiInvX)%*%sum_xPhipsi1Phix%*%solve(sum_xPhiInvX)%*%res$x%*%res$Phi.inv%*%res$Phi.derivatives$psi2%*%res$Phi.inv%*%res$x))-
                                    sum(diag(solve(sum_xPhiInvX)%*%res$x%*%res$Phi.inv%*%res$Phi.second$psi1psi2%*%res$Phi.inv%*%res$x))+
                                    sum(diag(solve(sum_xPhiInvX)%*%res$x%*%res$Phi.inv%*%res$Phi.derivatives$psi1%*%res$Phi.inv%*%res$Phi.derivatives$psi2%*%res$Phi.inv%*%res$x))+
                                    sum(diag(solve(sum_xPhiInvX)%*%res$x%*%res$Phi.inv%*%res$Phi.derivatives$psi2%*%res$Phi.inv%*%res$Phi.derivatives$psi1%*%res$Phi.inv%*%res$x))+
                                    2*t(res$yres)%*%res$Phi.inv%*%res$Phi.derivatives$psi2%*%res$Phi.inv%*%res$x%*%dbeta_dpsi1-
                                    t(res$yres)%*%(-res$Phi.inv%*%res$Phi.derivatives$psi1%*%res$Phi.inv%*%res$Phi.derivatives$psi2%*%res$Phi.inv-
                                                     res$Phi.inv%*%res$Phi.derivatives$psi2%*%res$Phi.inv%*%res$Phi.derivatives$psi1%*%res$Phi.inv+
                                                     res$Phi.inv%*%res$Phi.second$psi1psi2%*%res$Phi.inv)%*%res$yres)
    Hpsi1rho <- Hpsi1rho - 0.5*(sum(diag(res$Phi.inv%*%res$Phi.second$psi1rho)) - sum(diag(res$Phi.inv%*%res$Phi.derivatives$psi1%*%res$Phi.inv%*%res$Phi.derivatives$rho))-
                                  sum(diag(solve(sum_xPhiInvX)%*%sum_xPhipsi1Phix%*%solve(sum_xPhiInvX)%*%res$x%*%res$Phi.inv%*%res$Phi.derivatives$rho%*%res$Phi.inv%*%res$x))-
                                  sum(diag(solve(sum_xPhiInvX)%*%res$x%*%res$Phi.inv%*%res$Phi.second$psi1rho%*%res$Phi.inv%*%res$x))+
                                  sum(diag(solve(sum_xPhiInvX)%*%res$x%*%res$Phi.inv%*%res$Phi.derivatives$psi1%*%res$Phi.inv%*%res$Phi.derivatives$rho%*%res$Phi.inv%*%res$x))+
                                  sum(diag(solve(sum_xPhiInvX)%*%res$x%*%res$Phi.inv%*%res$Phi.derivatives$rho%*%res$Phi.inv%*%res$Phi.derivatives$psi1%*%res$Phi.inv%*%res$x))+
                                  2*t(res$yres)%*%res$Phi.inv%*%res$Phi.derivatives$rho%*%res$Phi.inv%*%res$x%*%dbeta_dpsi1-
                                  t(res$yres)%*%(-res$Phi.inv%*%res$Phi.derivatives$psi1%*%res$Phi.inv%*%res$Phi.derivatives$rho%*%res$Phi.inv-
                                                   res$Phi.inv%*%res$Phi.derivatives$rho%*%res$Phi.inv%*%res$Phi.derivatives$psi1%*%res$Phi.inv+
                                                   res$Phi.inv%*%res$Phi.second$psi1rho%*%res$Phi.inv)%*%res$yres)
    Hpsi2psi2 <- Hpsi2psi2 - 0.5*(sum(diag(res$Phi.inv%*%res$Phi.second$psi2psi2)) - sum(diag(res$Phi.inv%*%res$Phi.derivatives$psi2%*%res$Phi.inv%*%res$Phi.derivatives$psi2))-
                                    sum(diag(solve(sum_xPhiInvX)%*%sum_xPhipsi2Phix%*%solve(sum_xPhiInvX)%*%res$x%*%res$Phi.inv%*%res$Phi.derivatives$psi2%*%res$Phi.inv%*%res$x))-
                                    sum(diag(solve(sum_xPhiInvX)%*%res$x%*%res$Phi.inv%*%res$Phi.second$psi2psi2%*%res$Phi.inv%*%res$x))+
                                    sum(diag(solve(sum_xPhiInvX)%*%res$x%*%res$Phi.inv%*%res$Phi.derivatives$psi2%*%res$Phi.inv%*%res$Phi.derivatives$psi2%*%res$Phi.inv%*%res$x))+
                                    sum(diag(solve(sum_xPhiInvX)%*%res$x%*%res$Phi.inv%*%res$Phi.derivatives$psi2%*%res$Phi.inv%*%res$Phi.derivatives$psi2%*%res$Phi.inv%*%res$x))+
                                    2*t(res$yres)%*%res$Phi.inv%*%res$Phi.derivatives$psi2%*%res$Phi.inv%*%res$x%*%dbeta_dpsi2-
                                    t(res$yres)%*%(-res$Phi.inv%*%res$Phi.derivatives$psi2%*%res$Phi.inv%*%res$Phi.derivatives$psi2%*%res$Phi.inv-
                                                     res$Phi.inv%*%res$Phi.derivatives$psi2%*%res$Phi.inv%*%res$Phi.derivatives$psi2%*%res$Phi.inv+
                                                     res$Phi.inv%*%res$Phi.second$psi2psi2%*%res$Phi.inv)%*%res$yres)
    Hpsi2rho <- Hpsi2rho - 0.5*(sum(diag(res$Phi.inv%*%res$Phi.second$psi2rho)) - sum(diag(res$Phi.inv%*%res$Phi.derivatives$psi2%*%res$Phi.inv%*%res$Phi.derivatives$rho))-
                                  sum(diag(solve(sum_xPhiInvX)%*%sum_xPhipsi2Phix%*%solve(sum_xPhiInvX)%*%res$x%*%res$Phi.inv%*%res$Phi.derivatives$rho%*%res$Phi.inv%*%res$x))-
                                  sum(diag(solve(sum_xPhiInvX)%*%res$x%*%res$Phi.inv%*%res$Phi.second$psi2rho%*%res$Phi.inv%*%res$x))+
                                  sum(diag(solve(sum_xPhiInvX)%*%res$x%*%res$Phi.inv%*%res$Phi.derivatives$psi2%*%res$Phi.inv%*%res$Phi.derivatives$rho%*%res$Phi.inv%*%res$x))+
                                  sum(diag(solve(sum_xPhiInvX)%*%res$x%*%res$Phi.inv%*%res$Phi.derivatives$rho%*%res$Phi.inv%*%res$Phi.derivatives$psi2%*%res$Phi.inv%*%res$x))+
                                  2*t(res$yres)%*%res$Phi.inv%*%res$Phi.derivatives$rho%*%res$Phi.inv%*%res$x%*%dbeta_dpsi2-
                                  t(res$yres)%*%(-res$Phi.inv%*%res$Phi.derivatives$psi2%*%res$Phi.inv%*%res$Phi.derivatives$rho%*%res$Phi.inv-
                                                   res$Phi.inv%*%res$Phi.derivatives$rho%*%res$Phi.inv%*%res$Phi.derivatives$psi2%*%res$Phi.inv+
                                                   res$Phi.inv%*%res$Phi.second$psi2rho%*%res$Phi.inv)%*%res$yres)
    Hrhorho <- Hrhorho - 0.5*(sum(diag(res$Phi.inv%*%res$Phi.second$rhorho)) - sum(diag(res$Phi.inv%*%res$Phi.derivatives$rho%*%res$Phi.inv%*%res$Phi.derivatives$rho))-
                                sum(diag(solve(sum_xPhiInvX)%*%sum_xPhirhoPhix%*%solve(sum_xPhiInvX)%*%res$x%*%res$Phi.inv%*%res$Phi.derivatives$rho%*%res$Phi.inv%*%res$x))-
                                sum(diag(solve(sum_xPhiInvX)%*%res$x%*%res$Phi.inv%*%res$Phi.second$rhorho%*%res$Phi.inv%*%res$x))+
                                sum(diag(solve(sum_xPhiInvX)%*%res$x%*%res$Phi.inv%*%res$Phi.derivatives$rho%*%res$Phi.inv%*%res$Phi.derivatives$rho%*%res$Phi.inv%*%res$x))+
                                sum(diag(solve(sum_xPhiInvX)%*%res$x%*%res$Phi.inv%*%res$Phi.derivatives$rho%*%res$Phi.inv%*%res$Phi.derivatives$rho%*%res$Phi.inv%*%res$x))+
                                2*t(res$yres)%*%res$Phi.inv%*%res$Phi.derivatives$rho%*%res$Phi.inv%*%res$x%*%dbeta_drho-
                                t(res$yres)%*%(-res$Phi.inv%*%res$Phi.derivatives$rho%*%res$Phi.inv%*%res$Phi.derivatives$rho%*%res$Phi.inv-
                                                 res$Phi.inv%*%res$Phi.derivatives$rho%*%res$Phi.inv%*%res$Phi.derivatives$rho%*%res$Phi.inv+
                                                 res$Phi.inv%*%res$Phi.second$rhorho%*%res$Phi.inv)%*%res$yres)

  }
  H <- matrix(0, 3, 3)
  H[1, 1] <- Hpsi1psi1
  H[2, 2] <- Hpsi2psi2
  H[3, 3] <- Hrhorho
  H[1, 2] <- H[2, 1] <- Hpsi1psi2
  H[1, 3] <- H[3, 1] <- Hpsi1rho
  H[2, 3] <- H[3, 2] <- Hpsi2rho
  J <- diag(c(psi1_2, psi2_2, (1 - rho^2)/2))
  H_trans <- t(J) %*% H %*% J
  return(list(H.orig = H,
              H.trans = H_trans))
}

solve_remlBeta <- function(reml.root, data){
  n <- length(data$y1)
  p <- 1
  X <- matrix(1, nrow = n, ncol = 1)
  # Parameters
  psi1_2 <- exp(reml.root[1])
  psi2_2 <- exp(reml.root[2])
  rho <- (2 * plogis(reml.root[3])) - 1  # inverse logit scaled to (-1, 1)

  # Data
  y1 <- data$y1
  y2 <- data$y2
  s1_2 <- (data$se1)^2
  s2_2 <- (data$se2)^2
  compute_reml.betas <- function(i) {
    s1 <- s1_2[i]
    s2 <- s2_2[i]
    psi1s <- psi1_2 + s1
    psi2s <- psi2_2 + s2
    sqrt_psi1s <- sqrt(psi1s)
    sqrt_psi2s <- sqrt(psi2s)
    rho_sqrt <- rho * sqrt_psi1s * sqrt_psi2s

    # Covariance matrix and its inverse
    Phi <- matrix(c(psi1s, rho_sqrt, rho_sqrt, psi2s), nrow = 2)
    Phi.inv <- solve(Phi)

    # Kronecker product and related matrix operations
    x <- t(kronecker(diag(2), X[i, ]))
    y <- matrix(c(y1[i], y2[i]), nrow = 2)
    xPhiInv <- t(x) %*% Phi.inv

    list(
      xphix       = xPhiInv %*% x,
      xphiy       = xPhiInv %*% y
    )
  }
  # Compute all results
  results <- lapply(1:n, compute_reml.betas)

  # Reduce sums
  sum_xphix      <- Reduce("+", lapply(results, `[[`, "xphix"))
  sum_xphiy      <- Reduce("+", lapply(results, `[[`, "xphiy"))
  inv_sum_xphix <- solve(sum_xphix)
  beta_hat <- inv_sum_xphix %*% sum_xphiy

  return(list(beta_hat = as.vector(beta_hat),
              se_beta_hat = sqrt(diag(inv_sum_xphix))))
}
reml_solve <- function(data,
                       pars.start,
                       verbose,
                       nleqslv.param){
  # Default starting values
  if (is.null(pars.start)) {
    v1_init <- var(data$y1) - mean(data$se1^2)
    v2_init <- var(data$y2) - mean(data$se2^2)
    v1_init <- max(v1_init, 1e-3)
    v2_init <- max(v2_init, 1e-3)
    pars.start <- c(
      log(v1_init),    # log(psi1^2)
      log(v2_init),    # log(psi2^2)
      0                # rho
    )

  }

  if (!is.numeric(pars.start) || length(pars.start) != 3) {
    stop("pars.start must be a numeric vector of length 3: (log(psi1), log(psi2), tanh(rho/2))")
  }

  # Solve equations
  res <- tryCatch({
    nleqslv::nleqslv(
      x = pars.start,
      fn = reml_scores,
      data = data,
      method = nleqslv.param$method,
      global = nleqslv.param$global,
      control = nleqslv.param$control
    )
  }, error = function(e) {
    if (verbose) message("nleqslv failed: ", e$message)
    return(NULL)
  })

  # --- check convergence safely ---
  if (!is.null(res)) {
    converged <- res$termcd %in% c(1, 2) &&
      all(is.finite(res$fvec)) &&
      max(abs(res$fvec)) < 1e-5
  } else {
    converged <- FALSE
  }

  # --- verbose reporting ---
  if (verbose) {
    if (is.null(res)) {
      message("reml_solve(): error when applying nleqslv.")
    } else if (!converged) {
      msg <- sprintf(
        "reml_solve(): non-convergence [termcd=%d]. Message: %s",
        res$termcd, res$message
      )
      msg <- paste0(msg, ". May consider adjusting the initial values.")
      message(msg)
    } else {
      msg <- sprintf(
        "reml_solve(): converged successfully [termcd=%d].",
        res$termcd
      )
      message(msg)
    }
  }
  if (!converged) {
    return(NULL)
  }
  return(res)
}

reml_paramBoot <- function(data, reml.results,
                           bootstrap.times,
                           pars.start,
                           nleqslv.param,
                           alpha) {
  n <- length(data$y1)
  p <- 1
  if(is.null(pars.start)){
    v1_init <- var(data$y1) - mean(data$se1^2)
    v2_init <- var(data$y2) - mean(data$se2^2)
    v1_init <- max(v1_init, 1e-3)
    v2_init <- max(v2_init, 1e-3)
    pars.start <- c(
      log(v1_init),    # log(psi1^2)
      log(v2_init),    # log(psi2^2)
      0                # rho
    )
  }
  beta1.est <- reml.results$Estimates[grepl("^beta1", reml.results$Parameter)]
  beta2.est <- reml.results$Estimates[grepl("^beta2", reml.results$Parameter)]
  psi1_2.est <- reml.results$Estimates[grepl("^psi1_2$", reml.results$Parameter)]
  psi2_2.est <- reml.results$Estimates[grepl("^psi2_2$", reml.results$Parameter)]
  rho.est <- reml.results$Estimates[grepl("^rho$", reml.results$Parameter)]

  se1 <- (data$se1)
  se2 <- (data$se2)

  bootstrap.results <- list()
  n.success <- 0
  n.attempts <- 0

  while (n.success < bootstrap.times) {
    n.attempts <- n.attempts + 1
    boot.data <- generate_data(
      n = n,
      beta1 = beta1.est, beta2 = beta2.est, rho = rho.est,
      psi1_2 = psi1_2.est, psi2_2 = psi2_2.est,
      se1 = se1, se2 = se2
    )

    result <- tryCatch({
      suppressWarnings({
        invisible(capture.output({
          out <- nleqslv::nleqslv(
            x = pars.start,
            fn = reml_scores,
            data = boot.data,
            method = nleqslv.param$method,
            global = nleqslv.param$global,
            control = nleqslv.param$control
          )
        }, type = "output"))

        root <- out$x
        if (!is.null(root) &&
            all(out$fvec < 1e-6) &&
            all(is.finite(root)) &&
            all(abs(root) < 10) &&
            out$termcd %in% c(1, 2)) {
          reml.beta <- solve_remlBeta(root, boot.data)
          c(reml.beta$beta_hat, root)
        } else {
          NULL
        }
      })
    }, error = function(e) NULL)

    if (!is.null(result)) {
      n.success <- n.success + 1
      bootstrap.results[[n.success]] <- result
    }
  }

  bootstrap.matrix <- do.call(rbind, bootstrap.results)

  beta1.mat <- bootstrap.matrix[, 1:p, drop = FALSE]
  beta2.mat <- bootstrap.matrix[, (p + 1):(2 * p), drop = FALSE]
  psi1_2.vec <- exp(bootstrap.matrix[, 2 * p + 1])
  psi2_2.vec <- exp(bootstrap.matrix[, 2 * p + 2])
  rho.vec <- 2 * plogis(bootstrap.matrix[, 2 * p + 3]) - 1

  bootstrap.transformed <- cbind(beta1.mat, beta2.mat,
                                 psi1_2 = psi1_2.vec,
                                 psi2_2 = psi2_2.vec,
                                 rho = rho.vec)

  bootstrap.mean <- colMeans(bootstrap.transformed)
  bootstrap.se <- apply(bootstrap.transformed, 2, sd)
  bootstrap.ci <- t(apply(bootstrap.transformed, 2, quantile, probs = c(alpha/2, 1 - alpha/2)))
  colnames(bootstrap.ci) <- c("CI_lower", "CI_upper")

  param.names <- c(paste0("beta1", 1:p),
                   paste0("beta2", 1:p),
                   "psi1_2", "psi2_2", "rho")

  data.frame(
    Parameter = param.names,
    bootMean = bootstrap.mean,
    bootSE = bootstrap.se,
    bootCI_lower = bootstrap.ci[, 1],
    bootCI_upper = bootstrap.ci[, 2]
  )
}
