#' @title Bayesian Estimation for Bivariate Random-Effects Meta-Analysis (BRMA) Model
#' @description
#' Performs Bayesian estimation of the Bivariate Random-Effects Meta-Analysis (BRMA) model
#' using \pkg{cmdstanr}. Supports flexible prior choices for the between-study variances
#' (\eqn{\psi_1^2}, \eqn{\psi_2^2}) and the correlation parameter (\eqn{\rho}).
#' Returns posterior summaries and credible intervals.
#' @importFrom cmdstanr cmdstan_model
#' @importFrom posterior as_draws_df subset_draws summarise_draws
#' @importFrom stats quantile plogis qlogis qnorm sd var
#' @importFrom utils capture.output modifyList
#' @param data A list containing transformed treatment effect estimates and corresponding standard errors for both endpoints.
#'             Required columns: `y1`, `se1`, `y2`, `se2`.
#' @param psi.prior.dist Character string specifying the prior for the between-study variances.
#'   Options:
#'   \describe{
#'     \item{\code{"half-t"}}{Half-t prior on the standard deviation (\eqn{\psi}) [default].}
#'     \item{\code{"inverse-gamma"}}{Inverse-gamma prior on the variance (\eqn{\psi^2}).}
#'   }
#' @param ig_shape Numeric. Shape parameter for the inverse-gamma prior if \code{psi.prior.dist = "inverse-gamma"}. Default is 2.1.
#' @param ig_scale Numeric. Scale parameter for the inverse-gamma prior if \code{psi.prior.dist = "inverse-gamma"}. Default is 1.
#' @param half_t_df Numeric. Degrees of freedom for the half-t prior if \code{psi.prior.dist = "half-t"}. Default is 3.
#' @param half_t_scale Numeric. Scale parameter for the half-t prior if \code{psi.prior.dist = "half-t"}. Default is 2.5.
#' @param rho.prior.dist Character string specifying the prior for the correlation parameter \eqn{\rho}.
#'   Options:
#'   \describe{
#'     \item{\code{"fisherz"}}{Fisher z-transform normal prior [default].}
#'     \item{\code{"uniform"}}{Uniform prior on \eqn{\rho} in [-1,1].}
#'   }
#' @param nuts_params A named list of settings controlling the
#'   No-U-Turn Sampler (NUTS) used when fitting the model via cmdstanr.
#'
#'   These settings are passed internally to the
#'   \code{\link[cmdstanr]{sample}} method of a
#'   \code{\link[cmdstanr]{CmdStanModel}} object (created via
#'   \code{\link[cmdstanr]{cmdstan_model}}).
#'
#'   The default is:
#'   \preformatted{
#'   list(chains = 4,
#'     iter_warmup = 1000,
#'     iter_sampling = 1000,
#'     adapt_delta = 0.99
#'   )
#'   }
#'   Typical components include:
#'   \describe{
#'     \item{\code{chains}}{Integer. Number of parallel MCMC chains. Default: 4.}
#'     \item{\code{iter_warmup}}{Integer. Number of warm-up iterations per chain. Default: 1000.}
#'     \item{\code{iter_sampling}}{Integer. Number of post-warm-up sampling iterations per chain. Default: 1000.}
#'     \item{\code{adapt_delta}}{Numeric. Target acceptance probability for the NUTS step size adaptation.
#'     Higher values (closer to 1) reduce the risk of divergent transitions, at the cost of longer runtime. Default: 0.99.}
#'   }
#'
#'   Other arguments supported by \code{\link[cmdstanr]{sample}}, such as
#'   \code{seed}, \code{init}, \code{max_treedepth}, can also be included.
#'   Any arguments not specified in \code{nuts_params} will retain the
#'   default values from \code{\link[cmdstanr]{sample}}.
#' @param verbose Logical. If \code{TRUE}, prints progress messages from Stan. Default is \code{TRUE}.
#' @param alpha Numeric. Significance level for posterior credible intervals. Default is 0.05.
#' @return A list with two elements:
#' \describe{
#'   \item{\code{Results}}{Data frame of posterior summaries for each parameter including:
#'     Posterior Mean, Posterior SD, and Credible Interval.}
#'   \item{\code{prior_summary}}{Character vector describing the priors used for each parameter.}
#' }
#' @examples
#' data <- list(y1 = c(0.198, -0.091, 0.345, 1.193, 0.828, 1.671, 0.152, 1.416),
#'              se1 = c(0.076, 0.081, 0.062, 0.387, 0.497, 0.195, 0.432, 0.254),
#'              y2 = c(-0.234, -0.614, -0.054, -0.670, 0.105, 1.045, -0.910, 0.485),
#'              se2 = c(0.465, 0.408, 0.240, 0.384, 0.239, 0.434, 0.141, 0.159))
#' t_bayes(data)
#' t_bayes(data, rho.prior.dist = "uniform")
#' @export

t_bayes <- function(data,
                    psi.prior.dist = c("half-t", "inverse-gamma"),
                    ig_shape = 2.1,
                    ig_scale = 1,
                    half_t_df = 3,
                    half_t_scale = 2.5,
                    rho.prior.dist = c("fisherz", "uniform"),
                    nuts_params = list(chains = 4,
                                       iter_warmup = 1000,
                                       iter_sampling = 1000,
                                       adapt_delta = 0.99),
                    verbose = TRUE,
                    alpha = 0.05) {
  # Input checks
  stopifnot(is.list(data))
  required <- c("y1", "y2", "se1", "se2")
  if (!all(required %in% names(data))) stop("data must contain y1, y2, se1, se2")
  if (length(data$y1) != length(data$y2)) stop("y1 and y2 must have the same length")
  psi.prior.dist <- match.arg(psi.prior.dist)
  rho.prior.dist <- match.arg(rho.prior.dist)

  # source stan file
  mod_path <- system.file("BayesModel", "brma.stan", package = "tassociation")

  # Prepare Stan data
  stan_data <- create_stan_data(
    data,
    psi.prior.dist = psi.prior.dist,
    rho.prior.dist = rho.prior.dist,
    ig_shape = ig_shape,
    ig_scale = ig_scale,
    half_t_df = half_t_df,
    half_t_scale = half_t_scale
  )

  # Compile model
  mod <- cmdstan_model(mod_path)

  # import nuts parameters from the list
  # first set default NUTS parameters
  default_nuts_params <- list(
    chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.99,
    refresh = 0
  )

  # Merge with user-specified nuts_params
  nuts_params <- modifyList(default_nuts_params, nuts_params)

  # Add required elements
  nuts_params$data <- stan_data
  nuts_params$parallel_chains <- nuts_params$chains
  nuts_params$show_messages <- verbose

  # Fit model
  fit <- tryCatch({
    do.call(mod$sample, nuts_params)
  }, error = function(e) e)

  if (!inherits(fit, "CmdStanMCMC")) {
    warning("MCMC failed: ", fit$message)
    return(NULL)
  }

  # Quick convergence check
  df <- as_draws_df(fit)
  param_names <- c("beta1","beta2","psi1_2","psi2_2","rho")
  df_sub <- subset_draws(df, variable = param_names)
  diag_stats <- summarise_draws(df_sub)
  rhat_max <- max(diag_stats$rhat, na.rm = TRUE)
  ess_min  <- min(diag_stats$ess_bulk, na.rm = TRUE)
  if (rhat_max > 1.01 | ess_min < 200) {
    warning(sprintf("Quick convergence check failed: Rhat_max=%.3f, ESS_min=%.0f",
                    rhat_max, ess_min))
  }

  ci_col <- paste0((1 - alpha) * 100, "% CrI")

  results <- data.frame(
    Parameter = param_names,
    `Posterior Mean` = NA_real_,
    `Posterior SD` = NA_real_,
    check.names = FALSE
  )
  results[[ci_col]] <- NA_character_  # add column for CI as string

  for (p in param_names) {
    var_draws <- df[[p]]
    ci_values <- round(quantile(var_draws, c(alpha/2, 1 - alpha/2)), 3)
    results[results$Parameter == p, "Posterior Mean"] <- round(mean(var_draws), 3)
    results[results$Parameter == p, "Posterior SD"] <- round(sd(var_draws), 3)
    results[results$Parameter == p, ci_col] <- paste0("(", ci_values[1], ", ", ci_values[2], ")")
  }
  # Prior summary
  prior_summary <- c()

  # Psi prior
  if (psi.prior.dist == "inverse-gamma") {
    prior_summary <- c(prior_summary,
                       paste0("Psi prior: inverse-gamma(",ig_shape, ", ",ig_scale,") on variance (psi^2)"))
  } else if (psi.prior.dist == "half-t") {
    prior_summary <- c(prior_summary,
                       paste0("Psi prior: half-t(", half_t_df, ", ",half_t_scale, ") on SD (psi)"))
  }

  # Rho prior
  prior_summary <- c(prior_summary,
                     paste0("Rho prior: ", rho.prior.dist))

  return(list(
    Results = results,
    prior_summary = prior_summary
  ))
}

## not export
create_stan_data <- function(data,
                             psi.prior.dist = c("half-t", "inverse-gamma"),
                             ig_shape = 2.1,
                             ig_scale = 1,
                             half_t_df = 3,
                             half_t_scale = 2.5,
                             rho.prior.dist = c("fisherz", "uniform")) {

  stopifnot(is.list(data))
  required <- c("y1", "y2", "se1", "se2")
  if (!all(required %in% names(data))) stop("data must contain y1, y2, se1, se2")
  if (length(data$y1) != length(data$y2)) stop("y1 and y2 must have the same length.")

  n <- length(data$y1)

  psi.prior.dist <- match.arg(psi.prior.dist)
  rho.prior.dist <- match.arg(rho.prior.dist)

  list(
    N = n,
    y = cbind(data$y1, data$y2),
    se1 = data$se1,
    se2 = data$se2,
    use_uniform = ifelse(rho.prior.dist == "uniform", 1L, 0L),
    use_half_t_psi = ifelse(psi.prior.dist == "half-t", 1L, 0L),
    ig_shape = ig_shape,
    ig_scale = ig_scale,
    half_t_df = half_t_df,
    half_t_scale = half_t_scale
  )
}
