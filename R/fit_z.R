#' Treatment-Effekt-Schaetzer: naiv vs. genestet
#'
#' Zwei Schaetzer fuer den Effekt des Cluster-Ebene-Treatments `z`
#' (Wahrheit: `tau` in [sim_stream()]). Hier beisst die Statistik:
#'
#' * `fit_z_naive()` -- `lm(score ~ z)`: behandelt jede Messung als
#'   unabhaengig. Bei geclusterten Daten sind seine SEs fuer
#'   Cluster-Ebene-Praediktoren **zu klein** (anticonservativ) --
#'   der Falsch-Positiv-Generator.
#' * `fit_z_lmm()` -- `lmer(score ~ z + (1 | therapist_id) +
#'   (1 | patient_id))`: respektiert, dass `z` an der Therapeut:in
#'   haengt und Information auf Cluster-Ebene begrenzt ist.
#' * `fit_z_satt()` -- wie `fit_z_lmm()`, aber mit
#'   **Satterthwaite-Freiheitsgraden** (lmerTest): liefert zusaetzlich
#'   einen echten p-Wert. Bei wenigen Clustern ist der Wald-z-Test
#'   selbst des richtigen Modells noch zu liberal -- diese Variante
#'   repariert das.
#'
#' Interface-Konvention: `fit_*`-Funktionen geben mindestens
#' `c(estimate, se)` zurueck; liefern sie zusaetzlich `p`, nutzen
#' [mc_run()], [run_peek()] und [run_gate()] diesen p-Wert statt der
#' Wald-Naeherung.
#'
#' Gleiches Interface wie alle `fit_*`-Funktionen.
#'
#' @param snap Ein Snapshot (`data.frame` mit `score`, `z`,
#'   `therapist_id`, `patient_id`).
#' @return Benannter Vektor `c(estimate, se)`.
#' @name fit_z
NULL

#' @rdname fit_z
#' @export
fit_z_naive <- function(snap) {
  m <- stats::lm(score ~ z, data = snap)
  co <- summary(m)$coefficients
  c(estimate = co["z", "Estimate"],
    se       = co["z", "Std. Error"])
}

#' @rdname fit_z
#' @export
fit_z_satt <- function(snap) {
  if (!requireNamespace("lmerTest", quietly = TRUE)) {
    stop("fit_z_satt() braucht das Paket 'lmerTest'.", call. = FALSE)
  }
  m <- lmerTest::lmer(
    score ~ z + (1 | therapist_id) + (1 | patient_id),
    data = snap, REML = TRUE
  )
  co <- summary(m)$coefficients
  c(estimate = co["z", "Estimate"],
    se       = co["z", "Std. Error"],
    p        = co["z", "Pr(>|t|)"])
}

#' @rdname fit_z
#' @export
fit_z_lmm <- function(snap) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("fit_z_lmm() braucht das Paket 'lme4'.", call. = FALSE)
  }
  m <- lme4::lmer(
    score ~ z + (1 | therapist_id) + (1 | patient_id),
    data = snap, REML = TRUE
  )
  co <- summary(m)$coefficients
  c(estimate = co["z", "Estimate"],
    se       = co["z", "Std. Error"])
}
