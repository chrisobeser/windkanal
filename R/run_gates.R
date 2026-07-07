#' Mehrere Release-Gates mit verteiltem Irrtums-Budget
#'
#' Realistischer als ein Einzeltermin ([run_gate()]): das System darf
#' zu mehreren **vorab festgelegten** Terminen pruefen. Damit das
#' Gesamt-Irrtumsrisiko bei `alpha` bleibt, wird das Budget auf die
#' Termine verteilt:
#'
#' * `correction = "bonferroni"`: jeder Termin prueft mit `alpha / k`
#'   (einfach, garantiert gueltig, etwas konservativ).
#' * `correction = "none"`: keine Verteilung -- absichtlich falsch,
#'   um zu zeigen, was unkorrigiertes Mehrfach-Pruefen kostet.
#'
#' Feinere Verteilungen (O'Brien-Fleming/Pocock-Alpha-Spending, wie in
#' klinischen Studien ueblich) sind als spaetere Option vorgesehen.
#'
#' @param stream Datenstrom aus [sim_stream()].
#' @param fit_fn Eine `fit_*`-Funktion.
#' @param at_weeks Vorab festgelegte Gate-Termine (Kalenderwochen).
#' @param alpha Gesamt-Irrtumsbudget.
#' @param correction `"bonferroni"` (Default) oder `"none"`.
#' @param min_rows Snapshots mit weniger Zeilen werden uebersprungen.
#'
#' @return Liste mit `gates` (data.frame: `week`, `estimate`, `se`,
#'   `alpha_local`, `significant`), `hit` und `first_hit_week`.
#' @export
run_gates <- function(stream, fit_fn, at_weeks, alpha = 0.05,
                      correction = c("bonferroni", "none"),
                      min_rows = 30) {
  correction <- match.arg(correction)
  at_weeks <- sort(at_weeks)
  k <- length(at_weeks)
  alpha_local <- if (correction == "bonferroni") alpha / k else alpha

  gates <- do.call(rbind, lapply(at_weeks, function(w) {
    snap <- snapshot(stream, w)
    if (nrow(snap) < min_rows) return(NULL)
    f <- suppressWarnings(suppressMessages(fit_fn(snap)))
    data.frame(week = w, estimate = f[["estimate"]], se = f[["se"]],
               alpha_local = alpha_local,
               significant = is_sig(f, alpha_local))
  }))
  hit <- !is.null(gates) && any(gates$significant)
  list(
    gates = gates,
    hit = hit,
    first_hit_week = if (hit) gates$week[which(gates$significant)[1]]
                     else NA_real_
  )
}
