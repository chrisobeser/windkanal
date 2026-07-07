#' Genesteter Slope-Schaetzer (Mixed Model)
#'
#' Schaetzt die mittlere Veraenderung pro Sitzung mit einem linearen
#' gemischten Modell: `score ~ session + (1 | therapist_id) +
#' (1 | patient_id)`. Respektiert damit, was [fit_slope()] bewusst
#' ignoriert -- Patient:innen liefern mehrere Messungen und stecken in
#' Therapeut:innen. Gleiches Interface wie alle `fit_*`-Funktionen:
#' nimmt einen Snapshot, gibt `c(estimate, se)` zurueck.
#'
#' Bewusst schlicht gehalten (nur Random Intercepts, keine Random
#' Slopes) -- Erweiterungen kommen als eigene `fit_*`-Varianten, damit
#' Vergleiche im Replay sauber bleiben.
#'
#' @param snap Ein Snapshot (`data.frame` mit `score`, `session`,
#'   `therapist_id`, `patient_id`).
#' @return Benannter Vektor `c(estimate, se)`.
#' @export
fit_lmm <- function(snap) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("fit_lmm() braucht das Paket 'lme4' (install.packages(\"lme4\")).",
         call. = FALSE)
  }
  m <- lme4::lmer(
    score ~ session + (1 | therapist_id) + (1 | patient_id),
    data = snap, REML = TRUE
  )
  co <- summary(m)$coefficients
  c(estimate = co["session", "Estimate"],
    se       = co["session", "Std. Error"])
}
