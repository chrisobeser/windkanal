#' windkanal: Simulation Testbed for Living Clinical Data Systems
#'
#' Ein Windkanal fuer lebende Versorgungs-Datensysteme: simuliert
#' wachsende Routine-Datenstroeme (sitzungsweise Outcomes, Patient:innen
#' genestet in Therapeut:innen, gestaffelter Fall-Eingang), damit
#' Schaetzverfahren und Analyse-Governance-Protokolle an bekannter
#' Wahrheit geprueft werden koennen, bevor echte Klinikdaten fliessen.
#'
#' Geplante Schichten (Roadmap in README):
#' \enumerate{
#'   \item Stream-Generator (`sim_stream()`) -- v0.0.1, minimal
#'   \item Snapshot-/Gate-Mechanik (Peeking vs. Release-Gates)
#'   \item Estimator-Wrapper (lme4, bcf, grf) + Diagnosands
#'   \item Kalibrierte Presets ("ambulanz_de") + Plasmode-Modus
#' }
#'
#' @keywords internal
#' @importFrom stats rnorm runif
"_PACKAGE"
