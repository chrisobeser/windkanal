#' ATE-Schaetzung mit Bayes-LMM (brms) -- Pruefling 13
#'
#' Bayesianisches Gegenstueck zu [fit_z_satt()]: dasselbe gemischte
#' Modell `score ~ z + (1 | therapist_id) + (1 | patient_id)`, aber
#' mit brms/Stan geschaetzt; berichtet werden Posterior-Mittel, -SD
#' und das zentrale 95-%-Kredibilitaetsintervall des z-Koeffizienten.
#'
#' **Fixierte Spezifikation (vor Erwartungs-Lauf, 2026-07-07):**
#' brms-Default-Priors (flacher Prior auf b_z, student_t(3, ...) auf
#' Intercept und SD-Parameter), 2 Chains, 1000 Warmup + 1000 Draws.
#' Compile-Cache: Das Stan-Modell wird einmal kompiliert und pro Welt
#' via `update()` neu befuellt (Determinismus ueber `seed`).
#'
#' @param stream Datenstrom aus [sim_stream()] (oder Snapshot).
#' @param chains,warmup,draws MCMC-Einstellungen (Default 2/1000/1000).
#' @param seed Zufalls-Seed (Pflicht).
#' @return Benannter Vektor `estimate`, `se` (Posterior-SD), `lo`,
#'   `hi` (2.5-/97.5-%-Quantile) -- kompatibel mit `cover()`.
#' @export
fit_z_brms <- local({
  cache <- new.env(parent = emptyenv())
  function(stream, chains = 2, warmup = 1000, draws = 1000, seed) {
    if (!requireNamespace("brms", quietly = TRUE)) {
      stop("fit_z_brms() requires the 'brms' package.", call. = FALSE)
    }
    if (missing(seed)) stop("`seed` is mandatory.", call. = FALSE)
    d <- as.data.frame(stream)
    if (is.null(cache$modell)) {
      cache$modell <- brms::brm(
        score ~ z + (1 | therapist_id) + (1 | patient_id),
        data = d, chains = chains, warmup = warmup,
        iter = warmup + draws, seed = seed, refresh = 0,
        silent = 2)
      fit <- cache$modell
    } else {
      fit <- stats::update(cache$modell, newdata = d,
                           chains = chains, warmup = warmup,
                           iter = warmup + draws, seed = seed,
                           refresh = 0, silent = 2)
    }
    dr <- as.data.frame(fit)[["b_z"]]
    c(estimate = mean(dr), se = stats::sd(dr),
      lo = unname(stats::quantile(dr, 0.025)),
      hi = unname(stats::quantile(dr, 0.975)))
  }
})
