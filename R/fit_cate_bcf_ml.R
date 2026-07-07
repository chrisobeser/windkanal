#' Estimate personalized effects (multilevel BCF, stochtree)
#'
#' Like [fit_cate_bcf()], but with the **multilevel BCF** from the
#' stochtree package by the BCF authors (Herren, Hahn et al.): the
#' therapist nesting enters the model as a **random intercept**
#' (`rfx_group_ids`) -- the fix, published 2025/26, for exactly the
#' weakness we measured in standard bcf (coverage collapse under
#' cluster-level treatments).
#'
#' Side finding of the integration: stochtree is ~30x faster than
#' the old bcf package and seed-deterministic (no floating-point
#' noise).
#'
#' **Documented preset:** stochtree additionally runs `num_gfr = 5`
#' grow-from-root warm-start iterations (the authors' package
#' default); we adopt it deliberately and document it here instead
#' of letting it run along silently.
#'
#' @param stream Data stream from [sim_stream()].
#' @param nburn,nsim MCMC burn-in and draws (defaults 500/500).
#' @param pihat Propensity estimate as in [fit_cate_bcf()]:
#'   `"constant"` (default) or `"glm"`.
#' @param rfx Include therapist random intercepts (`TRUE`, default --
#'   the multilevel repair). `FALSE` runs stochtree's BCF *without*
#'   random effects: the fast standard-BCF doppelganger, used to
#'   certify implementation equivalence with the old bcf package.
#' @param seed Random seed (mandatory; MCMC is stochastic).
#' @return `data.frame` (`patient_id`, `x`, `z`, `tau_hat`,
#'   `tau_lo`, `tau_hi`) with attribute `ate` (`estimate`, `se`,
#'   `lo`/`hi` = posterior quantiles).
#' @export
fit_cate_bcf_ml <- function(stream, nburn = 500, nsim = 500,
                            pihat = c("constant", "glm"),
                            rfx = TRUE, seed) {
  pihat <- match.arg(pihat)
  if (!requireNamespace("stochtree", quietly = TRUE)) {
    stop("fit_cate_bcf_ml() requires the 'stochtree' package.",
         call. = FALSE)
  }
  if (missing(seed)) {
    stop("`seed` is mandatory (MCMC is stochastic).", call. = FALSE)
  }
  p <- patients(stream)
  x <- as.data.frame(cate_features(p))
  ph <- if (pihat == "glm") {
    stats::glm(z ~ x, family = stats::binomial, data = p)$fitted.values
  } else {
    rep(mean(p$z), nrow(p))
  }
  args <- list(
    X_train = x, Z_train = p$z, y_train = p$score_mean,
    propensity_train = ph,
    num_burnin = nburn, num_mcmc = nsim,
    general_params = list(random_seed = seed, verbose = FALSE))
  if (isTRUE(rfx)) {
    args$rfx_group_ids_train <- p$therapist_id
    args$rfx_basis_train <- matrix(1, nrow(p), 1)
  }
  fit <- do.call(stochtree::bcf, args)
  draws <- fit$tau_hat_train  # n_persons x nsim

  out <- data.frame(
    patient_id = p$patient_id, x = p$x, z = p$z,
    tau_hat = rowMeans(draws),
    tau_lo  = apply(draws, 1, stats::quantile, 0.025),
    tau_hi  = apply(draws, 1, stats::quantile, 0.975)
  )
  ate_draws <- colMeans(draws)
  attr(out, "ate") <- c(estimate = mean(ate_draws),
                        se = stats::sd(ate_draws),
                        lo = unname(stats::quantile(ate_draws, 0.025)),
                        hi = unname(stats::quantile(ate_draws, 0.975)))
  out
}
