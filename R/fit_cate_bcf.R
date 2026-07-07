#' Estimate personalized effects (Bayesian Causal Forest)
#'
#' Like [fit_cate_grf()], but with the *Bayesian Causal Forest*
#' (Hahn, Murray & Carvalho, 2020) -- the estimator with the
#' **shorter leash**: the effect function is regularized separately
#' and strongly toward "everyone equal" ("shrink to homogeneity"),
#' and every person receives a full posterior distribution of their
#' effect (honest ranges instead of point values).
#'
#' **Honest limitation:** The bcf package knows no clusters -- the
#' therapist nesting is ignored here (unlike grf). Exactly this gap
#' is documented as an open methods point in the project; see
#' [fit_cate_bcf_ml()] for the multilevel repair candidate.
#'
#' @param stream Data stream from [sim_stream()].
#' @param nburn,nsim MCMC burn-in and draws (defaults 500/500).
#' @param n_chains Number of MCMC chains (default 1 -- one suffices
#'   for Monte Carlo studies; the bcf default would be 4 parallel
#'   chains = 4x the cost).
#' @param pihat Propensity estimate: `"constant"` (default --
#'   `mean(z)` for everyone; correct under randomization,
#'   **powerless under confounding**) or `"glm"` (logistic
#'   regression `z ~ x` -- activates bcf's guard against targeted
#'   selection).
#' @param seed Random seed (mandatory; MCMC is stochastic).
#' @return `data.frame` (`patient_id`, `x`, `z`, `tau_hat`,
#'   `tau_lo`, `tau_hi` -- 95% credible interval per person) with
#'   attribute `ate` (posterior mean and SD of the average effect,
#'   plus `lo`/`hi` posterior quantiles).
#' @export
fit_cate_bcf <- function(stream, nburn = 500, nsim = 500,
                         n_chains = 1L, pihat = c("constant", "glm"),
                         seed) {
  pihat <- match.arg(pihat)
  if (!requireNamespace("bcf", quietly = TRUE)) {
    stop("fit_cate_bcf() requires the 'bcf' package.", call. = FALSE)
  }
  if (missing(seed)) {
    stop("`seed` is mandatory (MCMC is stochastic).", call. = FALSE)
  }
  set.seed(seed)
  # bcf dumps tree log files into the working directory -- clean up,
  # but only those it newly created itself
  vorher <- list.files(pattern = "^(con|mod)_trees.*[.]txt$")
  on.exit({
    nachher <- list.files(pattern = "^(con|mod)_trees.*[.]txt$")
    unlink(setdiff(nachher, vorher))
  }, add = TRUE)
  p <- patients(stream)
  x <- cate_features(p)
  pihat <- if (pihat == "glm") {
    stats::glm(z ~ x, family = stats::binomial, data = p)$fitted.values
  } else {
    rep(mean(p$z), nrow(p))
  }

  fit <- suppressWarnings(suppressMessages(
    bcf::bcf(y = p$score_mean, z = p$z,
             x_control = x, x_moderate = x, pihat = pihat,
             nburn = nburn, nsim = nsim, n_chains = n_chains,
             no_output = TRUE, verbose = FALSE)
  ))
  draws <- fit$tau  # nsim x n_persons

  out <- data.frame(
    patient_id = p$patient_id, x = p$x, z = p$z,
    tau_hat = colMeans(draws),
    tau_lo  = apply(draws, 2, stats::quantile, 0.025),
    tau_hi  = apply(draws, 2, stats::quantile, 0.975)
  )
  ate_draws <- rowMeans(draws)
  # lo/hi = true posterior quantiles (the normal approximation
  # estimate +- 1.96*se can deviate for a skewed posterior)
  attr(out, "ate") <- c(estimate = mean(ate_draws),
                        se = stats::sd(ate_draws),
                        lo = unname(stats::quantile(ate_draws, 0.025)),
                        hi = unname(stats::quantile(ate_draws, 0.975)))
  out
}
