#' Analyse-Regime: Dauer-Peeking vs. Release-Gate
#'
#' Zwei Arten, denselben wachsenden Datenstrom zu befragen:
#'
#' * `run_peek()` -- der Dauer-Gucker: testet alle `every` Wochen und
#'   "stoppt bei Signifikanz" (meldet den ersten Treffer). So arbeitet
#'   ein lebendes System ohne Governance.
#' * `run_gate()` -- das Release-Gate: **ein** praeregistrierter Test
#'   zu **einem** vorab festgelegten Termin auf dem eingefrorenen
#'   Snapshot. So arbeitet das System mit Governance (v0: Einzel-Gate;
#'   Mehrfach-Gates mit Alpha-Spending folgen).
#'
#' Beide nutzen dieselben Snapshots und denselben Schaetzer -- der
#' Unterschied ist reines Protokoll. Genau das soll messbar werden.
#'
#' @param stream Datenstrom aus [sim_stream()].
#' @param fit_fn Eine `fit_*`-Funktion (z. B. [fit_z_lmm()]).
#' @param every Peek-Intervall in Wochen (Default 4).
#' @param at_week Gate-Termin (Kalenderwoche).
#' @param alpha Signifikanzniveau.
#' @param min_rows Snapshots mit weniger Zeilen werden uebersprungen.
#'
#' @return
#' `run_peek()`: Liste mit `looks` (data.frame aller Blicke:
#' `week`, `estimate`, `se`, `significant`), `hit` (wurde je
#' Signifikanz gemeldet?) und `first_hit_week`.
#' `run_gate()`: Liste mit `week`, `estimate`, `se`, `significant`.
#' @name regime
NULL

#' @rdname regime
#' @export
run_peek <- function(stream, fit_fn, every = 4, alpha = 0.05,
                     min_rows = 30) {
  weeks <- seq(min(stream$week), max(stream$week), by = every)

  looks <- do.call(rbind, lapply(weeks, function(w) {
    snap <- snapshot(stream, w)
    if (nrow(snap) < min_rows) return(NULL)
    f <- suppressWarnings(suppressMessages(fit_fn(snap)))
    data.frame(week = w, estimate = f[["estimate"]], se = f[["se"]],
               significant = is_sig(f, alpha))
  }))
  hit <- !is.null(looks) && any(looks$significant)
  list(
    looks = looks,
    hit = hit,
    first_hit_week = if (hit) looks$week[which(looks$significant)[1]]
                     else NA_real_
  )
}

#' @rdname regime
#' @export
run_gate <- function(stream, fit_fn, at_week, alpha = 0.05,
                     min_rows = 30) {
  snap <- snapshot(stream, at_week)
  if (nrow(snap) < min_rows) {
    return(list(week = at_week, estimate = NA_real_, se = NA_real_,
                significant = FALSE))
  }
  f <- suppressWarnings(suppressMessages(fit_fn(snap)))
  list(week = at_week,
       estimate = f[["estimate"]], se = f[["se"]],
       significant = is_sig(f, alpha))
}
