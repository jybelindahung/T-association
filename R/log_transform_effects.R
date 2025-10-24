#' @title Log-transform treatment effect estimates
#' @description Transforms treatment effect estimates (e.g., odds ratios, hazard ratios) and their
#' standard errors to the log scale. This is appropriate for effect estimates that are
#' approximately normally distributed after log transformation. For other types of
#' effect estimates, a different transformation may be required
#'
#' @param effects Numeric vector of treatment effect estimates.
#' @param se.effect Numeric vector of standard errors corresponding to `effects`.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{y}{Log-transformed treatment effect estimates.}
#'     \item{se}{Standard errors on the log scale, calculated using the delta method.}
#'   }
#'
#' @examples
#' eff1 <- c(1.301, 0.993, 1.618, 2.735, 2.083, 1.547, 1.102, 0.826)
#' se1 <- c(0.069, 0.032, 0.076, 0.101, 0.033, 0.111, 0.094, 0.061)
#' log_transform_effects(eff1, se1)
#'
#' @export
log_transform_effects <- function(effects, se.effect) {
  if (length(effects) != length(se.effect)) {
    stop("Length of 'effects' and 'se.effect' must be the same")
  }

  list(
    y = log(effects),
    se = se.effect / effects
  )
}
