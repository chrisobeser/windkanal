#' Replay: den Datenstrom in Zeitreihenfolge nachspielen
#'
#' Fittet dasselbe Modell auf eine Folge von Snapshots und sammelt die
#' Schaetzungen ein -- so laesst sich zusehen, wie sich ein Ergebnis
#' mit wachsendem Datenstrom schaerft (oder wandert). Das ist die
#' Grundfigur des lebenden Systems und spaeter die Basis von
#' Peeking-Simulation und Release-Gates.
#'
#' @param stream Datenstrom aus [sim_stream()].
#' @param at_weeks Vektor von Kalenderwochen, zu denen ein Snapshot
#'   gezogen und gefittet wird.
#' @param fit_fn Funktion, die einen Snapshot (`data.frame`) nimmt und
#'   ein benanntes Vektor-/Listenpaar `c(estimate = ..., se = ...)`
#'   zurueckgibt. Default: [fit_slope()] (naives lm, siehe dort).
#' @param min_rows Snapshots mit weniger Zeilen werden uebersprungen
#'   (Schaetzung = NA), Default 10.
#'
#' @return `data.frame` mit einer Zeile pro Snapshot: `snapshot_week`,
#'   `n_obs`, `n_patients`, `n_therapists`, `estimate`, `se`,
#'   `ci_low`, `ci_high` (95%-Wald-Intervall).
#'
#' @examples
#' s <- sim_stream(n_therapists = 10, patients_per_therapist = 10,
#'                 n_sessions = 8, seed = 42)
#' replay(s, at_weeks = c(10, 20, 30, 40, 50, 60))
#'
#' @export
replay <- function(stream, at_weeks, fit_fn = fit_slope,
                   min_rows = 10) {
  stopifnot(is.numeric(at_weeks), length(at_weeks) >= 1)
  at_weeks <- sort(at_weeks)

  one <- function(w) {
    snap <- snapshot(stream, w)
    base <- data.frame(
      snapshot_week = w,
      n_obs         = nrow(snap),
      n_patients    = length(unique(snap$patient_id)),
      n_therapists  = length(unique(snap$therapist_id))
    )
    if (nrow(snap) < min_rows) {
      return(cbind(base, estimate = NA_real_, se = NA_real_,
                   ci_low = NA_real_, ci_high = NA_real_))
    }
    fit <- fit_fn(snap)
    est <- unname(fit[["estimate"]])
    se  <- unname(fit[["se"]])
    cbind(base, estimate = est, se = se,
          ci_low  = est - 1.96 * se,
          ci_high = est + 1.96 * se)
  }

  out <- do.call(rbind, lapply(at_weeks, one))
  rownames(out) <- NULL
  out
}

#' Naiver Slope-Schaetzer (Referenz-Estimator v0.0.1)
#'
#' Schaetzt die mittlere Veraenderung pro Sitzung per einfachem
#' `lm(score ~ session)`. **Bewusst naiv**: ignoriert das Nesting
#' (Patient:innen in Therapeut:innen) und die Messwiederholung --
#' genau der Referenzpunkt, gegen den spaetere Estimator-Wrapper
#' (lme4, bcf) antreten.
#'
#' @param snap Ein Snapshot (`data.frame` mit `score`, `session`).
#' @return Benannter Vektor `c(estimate, se)`.
#' @export
fit_slope <- function(snap) {
  m <- stats::lm(score ~ session, data = snap)
  co <- summary(m)$coefficients
  c(estimate = co["session", "Estimate"],
    se       = co["session", "Std. Error"])
}
