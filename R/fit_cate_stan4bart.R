#' Estimate personalized effects (multilevel BART, stan4bart)
#'
#' Pruefling 14: the second multilevel implementation besides
#' [fit_cate_bcf_ml()], to test whether the random-intercept repair is
#' implementation-specific. stan4bart (Dorie et al., 2022, Entropy)
#' combines one BART surface with grouped terms estimated in Stan.
#' Model: `score_mean ~ bart(features + z) + (1 | therapist_id)`, with
#' stan4bart's native `treatment` support: the package generates the
#' counterfactual test frame itself, and individual effects are
#' per-draw contrasts `(2 z - 1) (ev_observed - ev_counterfactual)`.
#'
#' Difference to BCF by design: one response surface, no separate
#' effect forest, no propensity guard -- a genuine second
#' *implementation family*, not a stochtree clone.
#'
#' Note: dbarts is attached if necessary (stan4bart's formula
#' interface resolves `bart()` via dbarts).
#'
#' @param stream Data stream from [sim_stream()].
#' @param chains,warmup,draws MCMC settings (defaults 2/500/500 per
#'   chain, i.e. 1,000 posterior samples total).
#' @param seed Random seed (mandatory; MCMC is stochastic).
#' @return As [fit_cate_bcf_ml()]: `data.frame` with `tau_hat`,
#'   `tau_lo`, `tau_hi` and attribute `ate`.
#' @export
fit_cate_stan4bart <- function(stream, chains = 2, warmup = 500,
                               draws = 500, seed) {
  for (pkg in c("stan4bart", "dbarts")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("fit_cate_stan4bart() requires the '", pkg, "' package.",
           call. = FALSE)
    }
  }
  if (!"dbarts" %in% .packages()) {
    suppressMessages(try(attachNamespace("dbarts"), silent = TRUE))
  }
  if (missing(seed)) {
    stop("`seed` is mandatory (MCMC is stochastic).", call. = FALSE)
  }
  p <- patients(stream)
  d <- as.data.frame(cate_features(p))
  feats <- colnames(d)
  d$z <- p$z
  d$therapist_id <- factor(p$therapist_id)
  d$score_mean <- p$score_mean
  f <- stats::as.formula(paste(
    "score_mean ~ bart(", paste(c(feats, "z"), collapse = " + "),
    ") + (1 | therapist_id)"))
  # stan4bart transformiert den UNAUSGEWERTETEN Call (NSE): eine als
  # Variable uebergebene Formel wuerde bart() als echte Funktion
  # ausgewertet (dbarts::bart -> Fehler). Deshalb wird die Formel per
  # bquote in den Call gespleisst, als stuende sie literal da.
  aufruf <- bquote(stan4bart::stan4bart(
    .(f), data = d, treatment = z, chains = chains, warmup = warmup,
    iter = warmup + draws, seed = seed, verbose = -1))
  fit <- eval(aufruf)
  ev_obs <- dbarts::extract(fit, type = "ev", sample = "train")
  ev_cf  <- dbarts::extract(fit, type = "ev", sample = "test")
  tau_draws <- (2 * d$z - 1) * (ev_obs - ev_cf)  # n_persons x n_samples

  out <- data.frame(
    patient_id = p$patient_id, x = p$x, z = p$z,
    tau_hat = rowMeans(tau_draws),
    tau_lo  = apply(tau_draws, 1, stats::quantile, 0.025),
    tau_hi  = apply(tau_draws, 1, stats::quantile, 0.975)
  )
  ate_draws <- colMeans(tau_draws)
  attr(out, "ate") <- c(estimate = mean(ate_draws),
                        se = stats::sd(ate_draws),
                        lo = unname(stats::quantile(ate_draws, 0.025)),
                        hi = unname(stats::quantile(ate_draws, 0.975)))
  out
}
