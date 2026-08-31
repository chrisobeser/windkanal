#' windkanal: Stochastic Simulation Testbed for Living Clinical Data Systems
#'
#' Ein Windkanal fuer lebende Versorgungs-Datensysteme: stochastische
#' (Monte-Carlo-)Simulation wachsender Routine-Datenstroeme
#' (sitzungsweise Outcomes, Patient:innen genestet in Therapeut:innen,
#' gestaffelter Fall-Eingang), damit Schaetzverfahren und
#' Analyse-Governance-Protokolle an bekannter Wahrheit geprueft werden
#' koennen, bevor echte Klinikdaten fliessen.
#'
#' Die Schichten des Pakets (alle ausgeliefert; Details im README):
#' \enumerate{
#'   \item Stream-Generator (`sim_stream()`)
#'   \item Snapshot-/Gate-Mechanik (Peeking vs. Release-Gates)
#'   \item Estimator-Wrapper (lme4, bcf, grf) + Diagnosands
#'   \item Kalibrierte Presets ("ambulanz_de") + Plasmode-Modus
#' }
#'
#' @keywords internal
#' @importFrom stats rnorm runif
"_PACKAGE"
