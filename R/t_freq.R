#' @title Frequentist estimation for BRMA model parameters
#' @description Computes point estimates and confidence intervals for the Bivariate Random-Effects Meta-Analysis (BRMA) model using frequentist methods.
#' @importFrom rootSolve multiroot
#' @importFrom MASS mvrnorm

#' @param data A list containing transformed treatment effect estimates and corresponding standard errors for both endpoints.
#'             Required columns: `y1`, `se1`, `y2`, `se2`.
#' @param method Character string specifying the frequentist estimation method. Supported options: \code{c("REML" , "ML")}. Default is \code{"REML"}.
#' @param interval.method Character string specifying the method to construct confidence intervals. Supported options: \code{c("wald", "bootstrap", "auto")}.
#'                        Default is \code{"auto"}, which chooses appropriate method based on sample size.
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
#' data <- list(y1 = c(0.263, -0.007, 0.481, 1.006, 0.734, 0.436, 0.097, -0.191),
#'              se1 = c(0.053, 0.032, 0.047, 0.037, 0.016, 0.072, 0.085, 0.074),
#'              y2 = c(0.661, 0.777, 0.798, 1.061, 1.225, 0.728, 0.036, 0.610),
#'              se2 = c(0.049, 0.023, 0.061, 0.065, 0.012, 0.034, 0.083, 0.080))
#' t_freq(data) # parameters estimated by REML;
#'              # confidence intervals constructed by bootstrap (sample size <10)
#' t_freq(data, method = "ML", interval.method = "Wald", verbose = TRUE)
#' @export

t_freq <- function(data,
                   method = c("REML", "ML"),
                   interval.method = c("auto", "Wald", "bootstrap"),
                   bootstrap.times = 1000,
                   pars.start = NULL,
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

  # Run selected method
  if (method == "ML") {
    if (verbose) message("Running ML estimation...")
    old_warn <- getOption("warn"); options(warn = -1)
    ml.res <- tryCatch(ml.solve(data, pars.start = pars.start, verbose = verbose, alpha = alpha), error = function(e) NULL)
    options(warn = old_warn)
    if (is.null(ml.res)) return(NULL)

    if (interval.method == "bootstrap") {
      if (verbose) message(sprintf("Running ML bootstrap (%d resamples)...", bootstrap.times))
      ML.boot <- ml.param.boot(data, ml.res, bootstrap.times)
      ml.res[[ci_colname]] <- sprintf("(%.3f, %.3f)", ML.boot$bootCI_lower, ML.boot$bootCI_upper)
      ml.res <- ml.res[, c("Parameter", "Estimates", ci_colname)]  # hide se & Wald CI
    } else if (interval.method == "Wald") {
      if (verbose) message("Computing ML Wald CIs...")
      colnames(ml.res) <- c("Parameter", "Estimates", "SE", ci_colname)
    }

    text_summary <- paste0(
      "Point estimates obtained using Maximum Likelihood (ML).\n",
      sprintf("%s%% confidence intervals computed using the %s method%s.",
              round((1 - alpha) * 100),
              interval.method,
              if (interval.method == "bootstrap") sprintf(" with %d bootstrap replicates", bootstrap.times) else "")
    )

    return(list(table = ml.res, summary = text_summary))
  }

  if (method == "REML") {
    if (verbose) message("Running REML estimation...")
    old_warn <- getOption("warn"); options(warn = -1)
    reml.res <- tryCatch(reml.solve(data, pars.start = pars.start, verbose = verbose, alpha = alpha), error = function(e) NULL)
    options(warn = old_warn)
    if (is.null(reml.res)) return(NULL)

    if (interval.method == "bootstrap") {
      if (verbose) message(sprintf("Running REML bootstrap (%d resamples)...", bootstrap.times))
      REML.boot <- reml.param.boot(data, reml.res, bootstrap.times)
      reml.res[[ci_colname]] <- sprintf("(%.3f, %.3f)", REML.boot$bootCI_lower, REML.boot$bootCI_upper)
      reml.res <- reml.res[, c("Parameter", "Estimates", ci_colname)]  # hide se & Wald CI
    } else if (interval.method == "Wald") {
      if (verbose) message("Computing REML Wald CIs...")
      colnames(reml.res) <- c("Parameter", "Estimates", "SE", ci_colname)
    }

    text_summary <- paste0(
      "Point estimates obtained using Restricted Maximum Likelihood (REML). ",
      sprintf("%s%% CI computed using the %s method%s.",
              round((1 - alpha) * 100),
              interval.method,
              if (interval.method == "bootstrap") sprintf(" with %d bootstrap replicates", bootstrap.times) else "")
    )

    return(list(Results = reml.res, Summary = text_summary))
  }
}


### Not export
ml.derivatives.gen <- function(data, x) {
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

  res <- c(score1, score2, score3, score4, score5)
  return(res)
}
ml.hessian.gen <- function(data, ml.est){
  n <- length(data$y1)
  p <- 1

  # Parameters
  beta1.hat <- ml.est[1:p]
  beta2.hat <- ml.est[(p+1):(2*p)]
  psi1_2.hat <- exp(ml.est[2*p + 1])
  psi2_2.hat <- exp(ml.est[2*p + 2])
  rho.hat <- (2 * plogis(ml.est[2*p + 3])) - 1  # inverse logit scaled to (-1, 1)

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
  return(list(est = c(beta1.hat, beta2.hat, psi1_2.hat, psi2_2.hat, rho.hat),
              var = diag(solve(-H))))
}
reml.derivatives.gen <- function(data, x) {
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

  return(c(eq1, eq2, eq3))
}

reml.hessian.gen <- function(data, reml.est){
  n <- length(data$y1)
  p <- 1
  X <- matrix(1, nrow = n, ncol = 1)

  # Parameters
  psi1_2 <- exp(reml.est[1])
  psi2_2 <- exp(reml.est[2])
  rho <- (2 * plogis(reml.est[3])) - 1  # inverse logit scaled to (-1, 1)

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

  return(list(est = c(beta_hat, psi1_2, psi2_2, rho),
              var = c(diag(inv_sum_xphix),
                      diag(solve(-hessian.reml)))))
}
# xx is the output of .hessian.gen
WaldCI.func <- function(xx, alpha) {
  z <- qnorm(1 - alpha/2)

  # beta
  lower.beta <- xx$est[1:2] - z * sqrt(xx$var[1:2])
  upper.beta <- xx$est[1:2] + z * sqrt(xx$var[1:2])

  # psi^2 (log scale)
  log.psi2 <- log(xx$est[3:4])
  var.log.psi2 <- xx$var[3:4] / (xx$est[3:4]^2)
  lower.psi2 <- exp(log.psi2 - z * sqrt(var.log.psi2))
  upper.psi2 <- exp(log.psi2 + z * sqrt(var.log.psi2))

  # rho (logit transform of (rho + 1)/2)
  rho.transformed <- qlogis((xx$est[5] + 1)/2)
  var.rho.transformed <- (2 / (1 - xx$est[5]^2))^2 * xx$var[5]
  lower.rho.trans <- rho.transformed - z * sqrt(var.rho.transformed)
  upper.rho.trans <- rho.transformed + z * sqrt(var.rho.transformed)
  lower.rho <- 2 * plogis(lower.rho.trans) - 1
  upper.rho <- 2 * plogis(upper.rho.trans) - 1

  # combine
  lower <- c(lower.beta, lower.psi2, lower.rho)
  upper <- c(upper.beta, upper.psi2, upper.rho)
  waldci <- sprintf("(%.3f, %.3f)", lower, upper)

  return(waldci)
}
ml.solve <- function(data, pars.start = NULL, maxiter = 50, verbose = TRUE, alpha = 0.05) {
  # Input checks
  stopifnot(is.list(data))
  required <- c("y1", "y2", "se1", "se2")
  if (!all(required %in% names(data))) stop("data must contain y1, y2, se1, se2")
  if (length(data$y1) != length(data$y2)) stop("y1 and y2 must have the same length")

  p <- 1

  # Default initial values
  if (is.null(pars.start)) {
    pars.start <- c(mean(data$y1), mean(data$y2),
                    log(var(data$y1)), log(var(data$y2)), 0)
  }

  # Input check
  if (!is.numeric(pars.start) || length(pars.start) != 5) {
    stop("pars.start must be numeric vector of length 5: (beta1, beta2, psi1, psi2, rho)")
  }

  # Root solving
  res <- tryCatch({
    old_warn <- getOption("warn")
    options(warn = 2)
    on.exit(options(warn = old_warn), add = TRUE)

    ml.root <- suppressMessages(suppressWarnings(
      multiroot(
        f = ml.derivatives.gen,
        start = pars.start,
        maxiter = maxiter,
        atol = 1e-6,
        data = data
      )
    ))

    # Check convergence
    if (is.null(ml.root$estim.precis) || is.na(ml.root$estim.precis) || ml.root$estim.precis >= 1e-6) {
      if (verbose) message("ml.solve(): root finding did not converge.")
      return(NULL)
    }

    # Theoretical variance and CI
    ml.hessian.res <- ml.hessian.gen(data, ml.root$root)
    ml.ci <- WaldCI.func(ml.hessian.res, alpha = alpha)

    res <- data.frame(Parameter = c("beta1", "beta2", "psi1_2", "psi2_2", "rho"),
                      Estimates = ml.hessian.res$est,
                      SE = sqrt(ml.hessian.res$var),
                      CI = ml.ci)

    return(res)

  }, error = function(e) {
    if (verbose) {
      message("ml.solve() error: ", conditionMessage(e))
      tb <- capture.output(traceback())
      message(paste(tb, collapse = "\n"))
    }
    return(NULL)
  })

  # Ensure consistent output
  if (is.null(res)) {
    return(data.frame(
      Parameter = c("beta1", "beta2", "psi1_2", "psi2_2", "rho"),
      Estimates = NA_real_,
      SE = NA_real_,
      CI = NA_character_
    ))
  }
  return(res)
}


reml.solve <- function(data, pars.start = NULL, maxiter = 50, verbose = TRUE, alpha = 0.05) {
  # Input checks
  stopifnot(is.list(data))
  required <- c("y1", "y2", "se1", "se2")
  if (!all(required %in% names(data))) stop("data must contain y1, y2, se1, se2")
  if (length(data$y1) != length(data$y2)) stop("y1 and y2 must have the same length")

  p <- 1

  # Default initial values
  if (is.null(pars.start)) {
    pars.start <- c(log(var(data$y1)), log(var(data$y2)), 0)
  }

  # Input check
  if (!is.numeric(pars.start) || length(pars.start) != 3) {
    stop("pars.start must be numeric vector of length 3: (psi1, psi2, rho)")
  }

  # Root solving
  res <- tryCatch({
    old_warn <- getOption("warn")
    options(warn = 2)
    on.exit(options(warn = old_warn), add = TRUE)

    reml.root <- suppressMessages(suppressWarnings(
      multiroot(
        f = reml.derivatives.gen,
        start = pars.start,
        maxiter = maxiter,
        atol = 1e-6,
        data = data
      )
    ))

    # Check convergence
    if (is.null(reml.root$estim.precis) || is.na(reml.root$estim.precis) || reml.root$estim.precis >= 1e-6) {
      if (verbose) message("reml.solve(): root finding did not converge.")
      return(NULL)
    }

    # Theoretical variance and CI
    reml.hessian.res <- reml.hessian.gen(data, reml.root$root)
    reml.ci <- WaldCI.func(reml.hessian.res, alpha = alpha)

    res <- data.frame(Parameter = c("beta1", "beta2", "psi1_2", "psi2_2", "rho"),
                      Estimates = reml.hessian.res$est,
                      SE = sqrt(reml.hessian.res$var),
                      CI = reml.ci)

    return(res)

  }, error = function(e) {
    if (verbose) {
      message("reml.solve() error: ", conditionMessage(e))
      tb <- capture.output(traceback())
      message(paste(tb, collapse = "\n"))
    }
    return(NULL)
  })

  # Ensure consistent output
  if (is.null(res)) {
    return(data.frame(
      Parameter = c("beta1", "beta2", "psi1_2", "psi2_2", "rho"),
      Estimates = NA_real_,
      SE = NA_real_,
      CI = NA_character_
    ))
  }

  return(res)
}

solve_remlBeta <- function(reml.root, data){
  n <- length(data$y1)
  p <- 1
  X <- matrix(1, nrow = n, ncol = 1)
  # Parameters
  psi1_2 <- reml.root[1]^2
  psi2_2 <- reml.root[2]^2
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
  return(beta_hat)
}

# bootstrap codes
generate.data.general <- function(n, beta1, beta2, rho, psi1_2, psi2_2, se1, se2) {
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
ml.param.boot <- function(data, ml.results, bootstrap.times, alpha = 0.05) {
  n <- length(data$y1)
  p <- 1
  pars.start <- c(mean(data$y1), mean(data$y2), log(var(data$y1)), log(var(data$y2)), 0)

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

    ml.param.boot.samples <- generate.data.general(n = n,
                                                   beta1 = beta1.est, beta2 = beta2.est, rho = rho.est,
                                                   psi1_2 = psi1_2.est, psi2_2 = psi2_2.est,
                                                   se1 = se1, se2 = se2)

    result <- tryCatch({
      suppressWarnings({
        invisible(capture.output({
          out <- multiroot(
            f = ml.derivatives.gen,
            start = pars.start,
            maxiter = 20,
            atol = 1e-6,
            data = ml.param.boot.samples
          )
        }, type = "output"))

        root <- out$root
        if (!is.null(root) &&
            out$estim.precis < 1e-6 &&
            all(is.finite(root)) &&
            all(abs(root) < 10)) {
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

reml.param.boot <- function(data, reml.results, bootstrap.times, alpha = 0.05) {
  pars.start <- c(log(var(data$y1)), log(var(data$y2)), 0)
  p <- 1
  n <- length(data$y1)

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
    boot.data <- generate.data.general(
      n = n,
      beta1 = beta1.est, beta2 = beta2.est, rho = rho.est,
      psi1_2 = psi1_2.est, psi2_2 = psi2_2.est,
      se1 = se1, se2 = se2
    )

    result <- tryCatch({
      suppressWarnings({
        invisible(capture.output({
          out <- multiroot(
            f = reml.derivatives.gen,
            start = pars.start,
            maxiter = 20,
            atol = 1e-6,
            data = boot.data
          )
        }, type = "output"))

        root <- out$root

        if (!is.null(root) &&
            out$estim.precis < 1e-6 &&
            all(is.finite(root)) &&
            all(abs(root) < 10)) {
          reml.beta <- solve_remlBeta(root, boot.data)
          c(reml.beta, root)
        } else NULL
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
