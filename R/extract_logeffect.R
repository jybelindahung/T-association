#' @title Extract log-transformed effect and standard error from reported effect and confidence interval
#'
#' @description
#' Converts reported treatment effect estimates (e.g., odds ratios, hazard ratios)
#' to the log scale and computes their approximate standard errors using the reported confidence
#' interval limits. This method assumes that the log-transformed effect is approximately normally
#' distributed. For other effect measures that do not satisfy this assumption,
#' a different transformation may be required.
#'
#' @param data A data frame where each row represents a study, containing the following columns:
#'   \describe{
#'     \item{effect}{The reported effect estimate (e.g., OR, HR).}
#'     \item{effect.lower}{The lower bound of the reported confidence interval.}
#'     \item{effect.upper}{The upper bound of the reported confidence interval.}
#'     \item{ci.level}{The confidence level (e.g., 95 for a 95\% CI).}
#'   }

#' @details
#' The standard error of the log-transformed effect is approximated as:
#' \deqn{
#' \mathrm{SE}(\log(\hat{\theta})) \approx
#' \frac{\log(U) - \log(L)}{2 z_{1-\alpha/2}},
#' }
#' where \eqn{L} and \eqn{U} are the lower and upper confidence limits, and
#' \eqn{z_{1-\alpha/2}} is the quantile of the standard normal distribution
#' corresponding to the confidence level.
#'
#' This approximation assumes that the confidence interval is symmetric on the log scale,
#' which is standard for ratio-type measures.
#'
#' @return A list with two numeric vectors:
#' \describe{
#'   \item{y}{The log-transformed effect estimates.}
#'   \item{se}{The approximate standard errors of the log-transformed effects.}
#' }
#'
#' @examples
#' df <- data.frame(
#'   effect = c(0.71, 1.76),
#'   effect.lower = c(0.46, 1.08),
#'   effect.upper = c(1.06, 2.86),
#'   ci.level = c(95, 95)
#' )
#' extract_logeffect(df)
#'
#' @export
extract_logeffect <- function(data){
  log.effect <- log(data$effect)
  log.effect.lower <- log(data$effect.lower)
  log.effect.upper <- log(data$effect.upper)
  alpha <- 1 - data$ci.level/100
  z <- qnorm(1 - alpha/2)
  se <- (log.effect.upper - log.effect.lower) / (2 * z)
  return(list(y = log.effect, se = se))
}
