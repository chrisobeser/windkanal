#' Heterogenitaets-Schaetzer: wem hilft das Treatment?
#'
#' Schaetzt die **Interaktion z:x** -- haengt der Treatment-Effekt vom
#' Personen-Merkmal `x` ab? (Wahrheit: `tau_x` in [sim_stream()].)
#' Das ist die "Personalisierungs-Frage" in ihrer einfachsten Form.
#'
#' * `fit_zx_naive()` -- `lm(score ~ z * x)`: zaehlt jede Sitzung als
#'   unabhaengige Beobachtung. `x` ist aber pro Person konstant --
#'   das naive Modell tut so, als haette es n_sessions-mal so viel
#'   Information ueber jede Person. Der Personalisierungs-
#'   Geister-Generator.
#' * `fit_zx_satt()` -- gemischtes Modell mit Satterthwaite-Inferenz:
#'   respektiert Person und Therapeut:in als Cluster.
#'
#' @param snap Ein Snapshot (`data.frame` mit `score`, `z`, `x`,
#'   `therapist_id`, `patient_id`).
#' @return Benannter Vektor `c(estimate, se)` bzw. `c(estimate, se, p)`.
#' @name fit_zx
NULL

#' @rdname fit_zx
#' @export
fit_zx_naive <- function(snap) {
  m <- stats::lm(score ~ z * x, data = snap)
  co <- summary(m)$coefficients
  c(estimate = co["z:x", "Estimate"],
    se       = co["z:x", "Std. Error"])
}

#' @rdname fit_zx
#' @export
fit_zx_satt <- function(snap) {
  if (!requireNamespace("lmerTest", quietly = TRUE)) {
    stop("fit_zx_satt() braucht das Paket 'lmerTest'.", call. = FALSE)
  }
  m <- lmerTest::lmer(
    score ~ z * x + (1 | therapist_id) + (1 | patient_id),
    data = snap, REML = TRUE
  )
  co <- summary(m)$coefficients
  c(estimate = co["z:x", "Estimate"],
    se       = co["z:x", "Std. Error"],
    p        = co["z:x", "Pr(>|t|)"])
}
