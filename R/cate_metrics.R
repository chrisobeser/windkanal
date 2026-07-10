#' Standard metrics for CATE estimates against known truth
#'
#' Computes the fixed metric set for judging individual-effect
#' estimates in one call: ranking quality, magnitude error, bias,
#' dispersion, person-level interval coverage, and an honest error
#' flag. Exists as one shared function so that no experiment script
#' can forget a metric again (the plasmode v1 run recorded no PEHE,
#' which left PAI's trivial r = 1.000 in a linear world without its
#' corrective).
#'
#' Interpretation: `r` measures *ranking* only -- it is scale-free and
#' blind to amplitude errors. In worlds with a linear true effect any
#' affine estimator attains r = 1 regardless of how wrong its
#' magnitudes are. `pehe`, the root mean squared error of the
#' individual effects (the magnitude yardstick of the CATE
#' literature), is the corrective. Read them together: r for ranking,
#' PEHE for magnitude.
#'
#' This is a simulation/plasmode metric set: it requires the true
#' individual effects, which real data never provide.
#'
#' @param tau_hat Numeric vector of estimated individual effects, or
#'   `NULL` for fits that produce none (ATE-only tools, failed fits);
#'   `NULL` returns an all-`NA` row with `error = TRUE`.
#' @param truth Numeric vector of true individual effects (mandatory).
#' @param lo,hi Optional numeric vectors of person-level interval
#'   bounds. If supplied, `covered` is the share of persons whose
#'   interval contains their true effect.
#' @return One-row `data.frame` with columns `n`, `r`, `pehe`, `bias`,
#'   `sd_tau_hat`, `covered`, `error`. `r` is `NA` when `tau_hat` has
#'   zero variance (reported honestly, no warning).
#' @export
cate_metrics <- function(tau_hat, truth, lo = NULL, hi = NULL) {
  stopifnot(is.numeric(truth), length(truth) > 0)
  if (is.null(tau_hat)) {
    return(data.frame(n = length(truth), r = NA_real_, pehe = NA_real_,
                      bias = NA_real_, sd_tau_hat = NA_real_,
                      covered = NA_real_, error = TRUE))
  }
  stopifnot(is.numeric(tau_hat), length(tau_hat) == length(truth))
  if (!is.null(lo) || !is.null(hi)) {
    if (is.null(lo) || is.null(hi)) {
      stop("Supply both `lo` and `hi`, or neither.", call. = FALSE)
    }
    stopifnot(length(lo) == length(truth), length(hi) == length(truth))
  }
  data.frame(
    n = length(truth),
    r = suppressWarnings(stats::cor(tau_hat, truth)),
    pehe = sqrt(mean((tau_hat - truth)^2)),
    bias = mean(tau_hat - truth),
    sd_tau_hat = stats::sd(tau_hat),
    covered = if (is.null(lo)) NA_real_ else mean(lo < truth & truth < hi),
    error = FALSE
  )
}
