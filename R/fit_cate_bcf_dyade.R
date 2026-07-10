#' Dyadischer BCF v1: Multilevel-BCF mit explizitem Passungs-Feature
#'
#' Erweiterung von [fit_cate_bcf_ml()] fuer die Matching-Frage
#' (I30 v1a): Der Feature-Raum erhaelt zusaetzlich das **explizite
#' Produkt** `x * therapist_c` -- die Passungs-Komponente als eigene
#' Spalte. Motivation: Baeume muessen Interaktionen sonst ueber
#' verschachtelte Splits finden; bei Feld-Stichproben (n ~ 200,
#' Caseloads 5-10) ist der vorgerechnete Produkt-Term ein einzelner
#' Split und damit unter Shrinkage deutlich leichter zu entdecken.
#' Therapeuten-Random-Intercepts bleiben aktiv (die Nesting-Reparatur).
#'
#' Einordnung: v1a der Entwicklungslinie I30 (dyadischer CATE);
#' [exploratory], Erwartungen im Experiment-Skript
#' `experiments/dyade/mc_dyade_v1.R`.
#'
#' @param stream Datenstrom aus [sim_stream()]; muss mit
#'   `tau_c`/`tau_xc`-fahigem Design erzeugt sein (braucht
#'   `therapist_c`).
#' @param nburn,nsim MCMC Burn-in und Draws (Defaults 500/500).
#' @param pihat Propensity wie in [fit_cate_bcf_ml()].
#' @param seed Zufalls-Seed (Pflicht).
#' @return Wie [fit_cate_bcf_ml()]: `data.frame` mit `tau_hat`,
#'   `tau_lo`, `tau_hi` und Attribut `ate`.
#' @export
fit_cate_bcf_dyade <- function(stream, nburn = 500, nsim = 500,
                               pihat = c("constant", "glm"), seed) {
  pihat <- match.arg(pihat)
  if (!requireNamespace("stochtree", quietly = TRUE)) {
    stop("fit_cate_bcf_dyade() requires the 'stochtree' package.",
         call. = FALSE)
  }
  if (missing(seed)) {
    stop("`seed` is mandatory (MCMC is stochastic).", call. = FALSE)
  }
  p <- patients(stream)
  if (is.null(p$therapist_c)) {
    stop("Dyadischer BCF braucht `therapist_c` im Stream ",
         "(tau_c/tau_xc-Design).", call. = FALSE)
  }
  x <- as.data.frame(cate_features(p))
  x$passung <- p$x * p$therapist_c  # v1a: explizite Produkt-Spalte
  ph <- if (pihat == "glm") {
    stats::glm(z ~ x, family = stats::binomial, data = p)$fitted.values
  } else {
    rep(mean(p$z), nrow(p))
  }
  args <- list(
    X_train = x, Z_train = p$z, y_train = p$score_mean,
    propensity_train = ph,
    num_burnin = nburn, num_mcmc = nsim,
    general_params = list(random_seed = seed, verbose = FALSE),
    rfx_group_ids_train = p$therapist_id,
    rfx_basis_train = matrix(1, nrow(p), 1))
  fit <- do.call(stochtree::bcf, args)
  draws <- fit$tau_hat_train

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
