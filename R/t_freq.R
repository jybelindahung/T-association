#' @title Frequentist estimation for BRMA model parameters
#' @description Computes point estimates and confidence intervals for the Bivariate Random-Effects Meta-Analysis (BRMA) model using frequentist methods.
#' @importFrom nleqslv nleqslv
#' @importFrom MASS mvrnorm

#' @param data A list containing transformed treatment effect estimates and corresponding standard errors for both endpoints.
#'             Required columns: `y1`, `se1`, `y2`, `se2`.
#' @param method Character string specifying the frequentist estimation method. Supported options: \code{c("REML" , "ML")}. Default is \code{"REML"}.
#' @param interval.method Character string specifying the method to construct confidence intervals. Supported options: \code{c("wald", "bootstrap", "auto")}.
#'                        Default is \code{"auto"}, which chooses appropriate method based on sample size.
#' @param nleqslv.param A list of control settings passed directly to
#'   \code{\link[nleqslv]{nleqslv}} for solving the system of nonlinear equations
#'   used to obtain point estimates.
#'
#'   The default is:
#'   \preformatted{
#'   list(
#'     method = "Broyden",
#'     control = list(maxit = 200, ftol = 1e-6, xtol = 1e-6)
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
#' t_freq(data, method = "ML", interval.method = "Wald", verbose = TRUE)
#' @export

t_freq <- function(data,
                   method = c("REML", "ML"),
                   interval.method = c("auto", "Wald", "bootstrap"),
                   bootstrap.times = 1000,
                   pars.start = NULL,
                   nleqslv.param = list(
                     method = "Broyden",
                     control = list(maxit = 200,
                                    ftol = 1e-6,
                                    xtol = 1e-6)),
                   verbose = TRUE,
                   alpha = 0.05) {
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
    interval.method <- if (n < 10) "bootstrap" else "Wald"
    if (verbose) message(sprintf("Sample size = %d -> using %s CI by default.", n, interval.method))
  }

  # CI name
  ci_label <- sprintf("%d", round((1 - alpha) * 100))
  ci_colname <- if (interval.method == "bootstrap") {
    paste0("bootstrapped ", ci_label, "% CI")
  } else {
    paste0("Wald ", ci_label, "% CI")
  }
  param.names <- c("beta1", "beta2", "psi1_2", "psi2_2", "rho")
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
    # Compute Wald CI
    ml.est_wald.ci <- freq_WaldCI(ml.root, ml.setrans_setheta$se.trans, alpha = alpha)
    df.mlwald <- data.frame(
      Parameter = param.names,
      Estimates = ml.est_wald.ci$Estimates,
      SE = ml.setrans_setheta$se.orig
    )
    df.mlwald[[ci_colname]] <- ml.est_wald.ci$CI
    options(warn = old_warn)
    if (interval.method == "bootstrap") {
      if (verbose) message(sprintf("Running ML bootstrap (%d resamples)...", bootstrap.times))
      mlboot.res <- ml_paramBoot(data,
                                 ml.results = df.mlwald,
                                 bootstrap.times = bootstrap.times,
                                 pars.start = pars.start,
                                 nleqslv.param = nleqslv.param,
                                 alpha = alpha)
      df.mlboot <- data.frame(
        Parameter = param.names,
        Estimates = ml.est_wald.ci$Estimates
      )
      df.mlboot[[ci_colname]] <- sprintf("(%.3f, %.3f)", mlboot.res$bootCI_lower, mlboot.res$bootCI_upper)
      ML <- df.mlboot
    } else if (interval.method == "Wald") {
      if (verbose) message("Computing ML Wald CIs...")
      ML <- df.mlwald
    }

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
    # Compute Wald CI
    reml.est_wald.ci <- freq_WaldCI(c(reml.beta$beta_hat, reml.root),
                                    c(reml.beta$se_beta_hat, reml.setrans_setheta$se.trans), alpha = alpha)

    df.remlwald <- data.frame(
      Parameter = param.names,
      Estimates = c(reml.est_wald.ci$Estimates),
      SE = c(reml.beta$se_beta_hat, reml.setrans_setheta$se.orig)
    )
    df.remlwald[[ci_colname]] <- reml.est_wald.ci$CI
    options(warn = old_warn)
    if (interval.method == "bootstrap") {
      if (verbose) message(sprintf("Running REML bootstrap (%d resamples)...", bootstrap.times))
      remlboot.res <- reml_paramBoot(data,
                                     reml.results = df.remlwald,
                                     bootstrap.times = bootstrap.times,
                                     pars.start = pars.start,
                                     nleqslv.param = nleqslv.param,
                                     alpha = alpha)
      df.remlboot <- data.frame(
        Parameter = param.names,
        Estimates = reml.est_wald.ci$Estimates
      )
      df.remlboot[[ci_colname]] <- sprintf("(%.3f, %.3f)", remlboot.res$bootCI_lower, remlboot.res$bootCI_upper)
      REML <- df.remlboot
    } else if (interval.method == "Wald") {
      if (verbose) message("Computing REML Wald CIs...")
      REML <- df.remlwald
    }

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
  psi1_2 <- exp(x[2*p + 1])
  psi2_2 <- exp(x[2*p + 2])
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
  psi1_2.hat <- exp(ml.root[2*p + 1])
  psi2_2.hat <- exp(ml.root[2*p + 2])
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

freq_WaldCI <- function(root, se.trans, alpha) {
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

  # Check convergence
  if (is.null(res) || res$termcd > 2 || any(abs(res$fvec) > 1e-5)) {  # 1 or 2 mean success
    if (verbose) message("ml_solve(): did not converge. termcd = ", res$termcd)
    return(NULL)
  }

  if (verbose) message("ml_solve(): converged successfully.")
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
  p <- 1

  # Parameters
  psi1_2 <- exp(x[1])
  psi2_2 <- exp(x[2])
  rho <- (2 * plogis(x[3])) - 1  # inverse logit scaled to (-1, 1)

  # Data
  y1 <- data$y1
  y2 <- data$y2
  s1_2 <- (data$se1)^2
  s2_2 <- (data$se2)^2
  X1 <- X2 <- X <- data$X <- matrix(1, ncol = 1, nrow = n)

  # function for a single observation
  compute_matrices <- function(i) {
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

    # Derivatives of Phi
    Phi.psi1 <- matrix(c(
      1, 0.5 * rho * sqrt_psi2s / sqrt_psi1s,
      0.5 * rho * sqrt_psi2s / sqrt_psi1s, 0
    ), nrow = 2)

    Phi.psi2 <- matrix(c(
      0, 0.5 * rho * sqrt_psi1s / sqrt_psi2s,
      0.5 * rho * sqrt_psi1s / sqrt_psi2s, 1
    ), nrow = 2)

    Phi.rho <- matrix(c(
      0, sqrt_psi1s * sqrt_psi2s,
      sqrt_psi1s * sqrt_psi2s, 0
    ), nrow = 2)

    # Kronecker product and related matrix operations
    x <- t(kronecker(diag(2), X[i, ]))
    y <- matrix(c(y1[i], y2[i]), nrow = 2)
    xPhiInv <- t(x) %*% Phi.inv

    list(
      xphix       = xPhiInv %*% x,
      xphiy       = xPhiInv %*% y,
      xphip1phix  = xPhiInv %*% Phi.psi1 %*% Phi.inv %*% x,
      xphip2phix  = xPhiInv %*% Phi.psi2 %*% Phi.inv %*% x,
      xphiprphix  = xPhiInv %*% Phi.rho %*% Phi.inv %*% x
    )
  }

  # Compute all results
  results <- lapply(1:n, compute_matrices)

  # Reduce sums
  sum_xphix      <- Reduce("+", lapply(results, `[[`, "xphix"))
  sum_xphiy      <- Reduce("+", lapply(results, `[[`, "xphiy"))
  sum_xphip1phix <- Reduce("+", lapply(results, `[[`, "xphip1phix"))
  sum_xphip2phix <- Reduce("+", lapply(results, `[[`, "xphip2phix"))
  sum_xphiprphix <- Reduce("+", lapply(results, `[[`, "xphiprphix"))

  inv_sum_xphix <- solve(sum_xphix)

  # Estimator of beta and gradients
  beta_hat <- inv_sum_xphix %*% sum_xphiy
  grad_psi1 <- -sum(diag(inv_sum_xphix %*% sum_xphip1phix))
  grad_psi2 <- -sum(diag(inv_sum_xphix %*% sum_xphip2phix))
  grad_rho  <- -sum(diag(inv_sum_xphix %*% sum_xphiprphix))

  mu1 <- c()
  mu2 <- c()
  for (i in 1:n) {
    mu <- t(kronecker(diag(2), X[i,]))%*%beta_hat
    mu1 <- c(mu1, mu[1])
    mu2 <- c(mu2, mu[2])
  }

  A <- sum(1 / (psi1_2 + s1_2))
  B <- sum(1 / (psi2_2 + s2_2))
  C <- sum(rho / sqrt((psi1_2 + s1_2) * (psi2_2 + s2_2)))


  eq1 <- -0.5*(A + grad_psi1 -
                 sum((y1 - mu1)^2 / ((1 - rho^2) * (psi1_2 + s1_2)^2)) +
                 sum(rho * (y1 - mu1) * (y2 - mu2) / ((1 - rho^2) * (psi1_2 + s1_2)^(3/2) * (psi2_2 + s2_2)^(1/2))))

  eq2 <- -0.5*(B + grad_psi2 -
                 sum((y2 - mu2)^2 / ((1 - rho^2) * (psi2_2 + s2_2)^2)) +
                 sum(rho * (y1 - mu1) * (y2 - mu2) / ((1 - rho^2) * (psi1_2 + s1_2)^(1/2) * (psi2_2 + s2_2)^(3/2))))

  eq3 <- -0.5*(-2 * rho * (n) / (1 - rho^2) +
                 grad_rho +
                 2 / (1 - rho^2)^2 * sum(
                   rho * ((y1 - mu1)^2 / (psi1_2 + s1_2) + (y2 - mu2)^2 / (psi2_2 + s2_2)) -
                     (1 + rho^2) * (y1 - mu1) * (y2 - mu2) / sqrt((psi1_2 + s1_2) * (psi2_2 + s2_2))
                 ))
  d_logpsi1 <- psi1_2 * eq1      # *d/d log(psi1^2) by chain rule
  d_logpsi2 <- psi2_2 * eq2      # *d/d log(psi2^2) by chain rule
  d_zrho    <- ((1 - rho^2) / 2) * eq3  # *d/d z_rho by chain rule

  return(c(d_logpsi1, d_logpsi2, d_zrho))
}

reml_hessian <- function(data, reml.root){
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

  compute_matrices.second <- function(i) {
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

    # First derivatives of Phi
    Phi.psi1 <- matrix(c(
      1, 0.5 * rho * sqrt_psi2s / sqrt_psi1s,
      0.5 * rho * sqrt_psi2s / sqrt_psi1s, 0
    ), nrow = 2)

    Phi.psi2 <- matrix(c(
      0, 0.5 * rho * sqrt_psi1s / sqrt_psi2s,
      0.5 * rho * sqrt_psi1s / sqrt_psi2s, 1
    ), nrow = 2)

    Phi.rho <- matrix(c(
      0, sqrt_psi1s * sqrt_psi2s,
      sqrt_psi1s * sqrt_psi2s, 0
    ), nrow = 2)

    # Second derivatives of Phi
    Phi.psi1psi1 <- matrix(c(
      0, -0.25 * rho * sqrt_psi2s / (sqrt_psi1s^3),
      -0.25 * rho * sqrt_psi2s / (sqrt_psi1s^3), 0
    ), nrow = 2)

    Phi.psi1psi2 <- matrix(c(
      0, 0.25 * rho / sqrt_psi2s / sqrt_psi1s,
      0.25 * rho / sqrt_psi2s / sqrt_psi1s, 0
    ), nrow = 2)

    Phi.psi1rho <- matrix(c(
      0, 0.5 * sqrt_psi2s / sqrt_psi1s,
      0.5 * sqrt_psi2s / sqrt_psi1s, 0
    ), nrow = 2)

    Phi.psi2psi2 <- matrix(c(
      0, -0.25 * rho * sqrt_psi1s / (sqrt_psi2s^3),
      -0.25 * rho * sqrt_psi1s / (sqrt_psi2s^3), 0
    ), nrow = 2)

    Phi.psi2rho <- matrix(c(
      0, 0.5 * sqrt_psi1s / sqrt_psi2s,
      0.5 * sqrt_psi1s / sqrt_psi2s, 0
    ), nrow = 2)

    Phi.rhorho <- matrix(c(
      0, 0,
      0, 0
    ), nrow = 2)

    # Kronecker product and related matrix operations
    x <- t(kronecker(diag(2), X[i, ]))
    y <- matrix(c(y1[i], y2[i]), nrow = 2)
    xPhiInv <- t(x) %*% Phi.inv

    qpsi1psi1 = -Phi.inv %*% Phi.psi1 %*% Phi.inv %*% Phi.psi1 %*% Phi.inv - Phi.inv %*% Phi.psi1 %*% Phi.inv %*% Phi.psi1 %*% Phi.inv + Phi.inv %*% Phi.psi1psi1 %*% Phi.inv
    qpsi1psi2 = -Phi.inv %*% Phi.psi1 %*% Phi.inv %*% Phi.psi2 %*% Phi.inv - Phi.inv %*% Phi.psi2 %*% Phi.inv %*% Phi.psi1 %*% Phi.inv + Phi.inv %*% Phi.psi1psi2 %*% Phi.inv
    qpsi1rho = -Phi.inv %*% Phi.psi1 %*% Phi.inv %*% Phi.rho %*% Phi.inv - Phi.inv %*% Phi.rho %*% Phi.inv %*% Phi.psi1 %*% Phi.inv + Phi.inv %*% Phi.psi1rho %*% Phi.inv
    qpsi2psi2 = -Phi.inv %*% Phi.psi2 %*% Phi.inv %*% Phi.psi2 %*% Phi.inv - Phi.inv %*% Phi.psi2 %*% Phi.inv %*% Phi.psi2 %*% Phi.inv + Phi.inv %*% Phi.psi2psi2 %*% Phi.inv
    qpsi2rho = -Phi.inv %*% Phi.psi2 %*% Phi.inv %*% Phi.rho %*% Phi.inv - Phi.inv %*% Phi.rho %*% Phi.inv %*% Phi.psi2 %*% Phi.inv + Phi.inv %*% Phi.psi2rho %*% Phi.inv
    qrhorho = -Phi.inv %*% Phi.rho %*% Phi.inv %*% Phi.rho %*% Phi.inv - Phi.inv %*% Phi.rho %*% Phi.inv %*% Phi.rho %*% Phi.inv + Phi.inv %*% Phi.rhorho %*% Phi.inv

    list(
      xphix       = xPhiInv %*% x,
      xphiy       = xPhiInv %*% y,
      xphip1phix  = xPhiInv %*% Phi.psi1 %*% Phi.inv %*% x,
      xphip2phix  = xPhiInv %*% Phi.psi2 %*% Phi.inv %*% x,
      xphiprphix  = xPhiInv %*% Phi.rho %*% Phi.inv %*% x,
      xqpsi1psi1x = t(x) %*% qpsi1psi1 %*% x,
      xqpsi1psi2x = t(x) %*% qpsi1psi2 %*% x,
      xqpsi1rhox  = t(x) %*% qpsi1rho %*% x,
      xqpsi2psi2x = t(x) %*% qpsi2psi2 %*% x,
      xqpsi2rhox  = t(x) %*% qpsi2rho %*% x,
      xqrhorhox   = t(x) %*% qrhorho %*% x
    )
  }
  # Compute all results
  results <- lapply(1:n, compute_matrices.second)

  # Reduce sums
  sum_xphix      <- Reduce("+", lapply(results, `[[`, "xphix"))
  sum_xphiy      <- Reduce("+", lapply(results, `[[`, "xphiy"))
  sum_xphip1phix <- Reduce("+", lapply(results, `[[`, "xphip1phix"))
  sum_xphip2phix <- Reduce("+", lapply(results, `[[`, "xphip2phix"))
  sum_xphiprphix <- Reduce("+", lapply(results, `[[`, "xphiprphix"))

  sum_xqpsi1psi1x<- Reduce("+", lapply(results, `[[`, "xqpsi1psi1x"))
  sum_xqpsi1psi2x <- Reduce("+", lapply(results, `[[`, "xqpsi1psi2x"))
  sum_xqpsi1rhox <- Reduce("+", lapply(results, `[[`, "xqpsi1rhox"))
  sum_xqpsi2psi2x <- Reduce("+", lapply(results, `[[`, "xqpsi2psi2x"))
  sum_xqpsi2rhox <- Reduce("+", lapply(results, `[[`, "xqpsi2rhox"))
  sum_xqrhorhox <- Reduce("+", lapply(results, `[[`, "xqrhorhox"))

  inv_sum_xphix <- solve(sum_xphix)
  beta_hat <- inv_sum_xphix %*% sum_xphiy
  grad_psi1psi1 <- -sum(diag(inv_sum_xphix %*% sum_xphip1phix %*% inv_sum_xphix %*% sum_xphip1phix)) - sum(diag(inv_sum_xphix %*% sum_xqpsi1psi1x))
  grad_psi1psi2 <- -sum(diag(inv_sum_xphix %*% sum_xphip1phix %*% inv_sum_xphix %*% sum_xphip2phix)) - sum(diag(inv_sum_xphix %*% sum_xqpsi1psi2x))
  grad_psi1rho <- -sum(diag(inv_sum_xphix %*% sum_xphip1phix %*% inv_sum_xphix %*% sum_xphiprphix)) - sum(diag(inv_sum_xphix %*% sum_xqpsi1rhox))
  grad_psi2psi2 <- -sum(diag(inv_sum_xphix %*% sum_xphip2phix %*% inv_sum_xphix %*% sum_xphip2phix)) - sum(diag(inv_sum_xphix %*% sum_xqpsi2psi2x))
  grad_psi2rho <- -sum(diag(inv_sum_xphix %*% sum_xphip2phix %*% inv_sum_xphix %*% sum_xphiprphix)) - sum(diag(inv_sum_xphix %*% sum_xqpsi2rhox))
  grad_rhorho <- -sum(diag(inv_sum_xphix %*% sum_xphiprphix %*% inv_sum_xphix %*% sum_xphiprphix)) - sum(diag(inv_sum_xphix %*% sum_xqrhorhox))

  # calculate mu1 and mu2
  mu1 <- c()
  mu2 <- c()
  for (i in 1:n) {
    mu <- t(kronecker(diag(2), X[i,]))%*%beta_hat
    mu1 <- c(mu1, mu[1])
    mu2 <- c(mu2, mu[2])
  }

  # Compute d2_ell_R_psi1_2_2
  term1_psi1_2_2 <- sum(1 / (psi1_2 + s1_2)^2)
  term4_psi1_2_2 <- sum((y1 - mu1)^2 / ((1 - rho^2) * (psi1_2 + s1_2)^3))
  term5_psi1_2_2 <- (3/2) * sum(rho * (y1 - mu1) * (y2 - mu2) / ((1 - rho^2) * (psi1_2 + s1_2)^(5/2) * (psi2_2 + s2_2)^(1/2)))

  d2_ell_R_psi1_2_2 <- -0.5 * (-term1_psi1_2_2 + grad_psi1psi1 + 2 * term4_psi1_2_2 - term5_psi1_2_2)

  # Compute d2_ell_R_psi2_2_2
  term1_psi2_2_2 <- sum(1 / (psi2_2 + s2_2)^2)
  term4_psi2_2_2 <- sum((y2 - mu2)^2 / ((1 - rho^2) * (psi2_2 + s2_2)^3))
  term5_psi2_2_2 <- (3/2) * sum(rho * (y1 - mu1) * (y2 - mu2) / ((1 - rho^2) * (psi1_2 + s1_2)^(1/2) * (psi2_2 + s2_2)^(5/2)))

  d2_ell_R_psi2_2_2 <- -0.5 * (-term1_psi2_2_2 + grad_psi2psi2 + 2 * term4_psi2_2_2 - term5_psi2_2_2)

  # Compute d2_ell_R_rho_2
  term1_rho_2 <- -2 * n * (1 + rho^2) / (1 - rho^2)^2
  term4_rho_2 <- 8 * rho / (1 - rho^2)^3 * sum(rho * ((y1 - mu1)^2 / (psi1_2 + s1_2) + (y2 - mu2)^2 / (psi2_2 + s2_2)) -
                                                 (1 + rho^2) * (y1 - mu1) * (y2 - mu2) / sqrt((psi1_2 + s1_2) * (psi2_2 + s2_2)))
  term5_rho_2 <- 2 / (1 - rho^2)^2 * sum((y1 - mu1)^2 / (psi1_2 + s1_2) + (y2 - mu2)^2 / (psi2_2 + s2_2) -
                                           2 * rho * (y1 - mu1) * (y2 - mu2) / sqrt((psi1_2 + s1_2) * (psi2_2 + s2_2)))
  d2_ell_R_rho_2 <- -0.5 * (term1_rho_2 + grad_rhorho + term4_rho_2 + term5_rho_2)

  # Compute d2_ell_R_psi1_2_psi2_2
  term3_psi1_2_psi2_2 <- -0.5*sum(rho * (y1 - mu1) * (y2 - mu2) / ((1 - rho^2) * (psi1_2 + s1_2)^(3/2) * (psi2_2 + s2_2)^(3/2)))
  d2_ell_R_psi1_2_psi2_2 <- -0.5 * (grad_psi1psi2 + term3_psi1_2_psi2_2)

  # Compute d2_ell_R_psi1_2_rho
  term3_psi1_2_rho <- (2/(1 - rho^2)^2) * sum(-rho*(y1 - mu1)^2 / (psi1_2 + s1_2)^2
                                              + (1 + rho^2) * (y1 - mu1) * (y2 - mu2) / (2*(psi1_2 + s1_2)^(3/2) * (psi2_2 + s2_2)^(1/2)))
  d2_ell_R_psi1_2_rho <- -0.5 * (grad_psi1rho + term3_psi1_2_rho)

  # Compute d2_ell_R_psi2_2_rho
  term3_psi2_2_rho <- (2/(1 - rho^2)^2) * sum(-rho*(y2 - mu2)^2 / (psi2_2 + s2_2)^2 +
                                                (1 + rho^2) * (y1 - mu1) * (y2 - mu2) / (2 * (psi1_2 + s1_2)^(1/2) * (psi2_2 + s2_2)^(3/2)))
  d2_ell_R_psi2_2_rho <- -0.5 * (grad_psi2rho + term3_psi2_2_rho)

  hessian.reml <- matrix(0, nrow = 3, ncol = 3)
  hessian.reml[1,1] <- d2_ell_R_psi1_2_2
  hessian.reml[2,2] <- d2_ell_R_psi2_2_2
  hessian.reml[3,3] <- d2_ell_R_rho_2
  hessian.reml[1,2] <- hessian.reml[2,1] <- d2_ell_R_psi1_2_psi2_2
  hessian.reml[2,3] <- hessian.reml[3,2] <- d2_ell_R_psi2_2_rho
  hessian.reml[1,3] <- hessian.reml[3,1] <- d2_ell_R_psi1_2_rho

  J <- diag(c(psi1_2, psi2_2, (1-rho^2)/2 ))
  H.trans <- t(J)%*%hessian.reml%*%J
  return(list(H.orig = hessian.reml,
              H.trans = H.trans))
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

  # Check convergence
  if (is.null(res) || res$termcd > 2 || any(abs(res$fvec) > 1e-5)) {  # 1 or 2 mean success
    if (verbose) message("reml_solve(): did not converge. termcd = ", res$termcd)
    return(NULL)
  }

  if (verbose) message("reml_solve(): converged successfully.")
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
