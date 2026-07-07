#' Standbild des Datenstroms zu einer Kalenderwoche
#'
#' Gibt den Datenstrom so zurueck, wie er in Woche `at_week` sichtbar
#' war: nur Messungen mit `week <= at_week`. Ein Snapshot ist das
#' Grundobjekt der Windkanal-Logik -- Replay, Peeking-Simulation und
#' Release-Gates sind alle nur Folgen von Snapshots.
#'
#' @param stream Ein Datenstrom aus [sim_stream()].
#' @param at_week Kalenderwoche des Standbilds (>= 0).
#'
#' @return Der gefilterte Strom. Traegt sein Aufnahmedatum als
#'   Attribut `snapshot_week` mit sich (Provenienz: jedes Ergebnis
#'   soll wissen, auf welchem Standbild es beruht).
#'
#' @examples
#' s <- sim_stream(n_therapists = 5, patients_per_therapist = 4,
#'                 n_sessions = 6, seed = 42)
#' snap <- snapshot(s, at_week = 20)
#' nrow(snap) < nrow(s)
#' attr(snap, "snapshot_week")
#'
#' @export
snapshot <- function(stream, at_week) {
  stopifnot(
    inherits(stream, "data.frame"),
    "week" %in% names(stream),
    is.numeric(at_week), length(at_week) == 1, at_week >= 0
  )
  out <- stream[stream$week <= at_week, , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "snapshot_week") <- at_week
  out
}
