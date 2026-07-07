#' Monte-Carlo-Experiment deklarieren und ausfuehren
#'
#' Wiederholt eine Simulation ueber viele Welten (Seeds) und wendet
#' auf jede Welt einen oder mehrere Schaetzer an. Das Experiment ist
#' damit vollstaendig durch seine Argumente beschrieben -- Windkanal-
#' Regel: ein Ergebnis, das man nicht aus `mc_run()`-Argumenten
#' reproduzieren kann, existiert nicht.
#'
#' @param n_reps Anzahl simulierter Welten.
#' @param sim_args Benannte Liste von Argumenten fuer [sim_stream()]
#'   (ohne `seed` -- die Seeds vergibt `mc_run`).
#' @param fit_fns Benannte Liste von `fit_*`-Funktionen.
#' @param seed_start Erster Seed; Welt i nutzt `seed_start + i - 1`.
#' @param alpha Signifikanzniveau fuer die `significant`-Spalte.
#'
#' @return `data.frame`: `rep`, `seed`, `estimator`, `estimate`,
#'   `se`, `significant`.
#'
#' @examples
#' res <- mc_run(n_reps = 20,
#'               sim_args = list(n_therapists = 10,
#'                               patients_per_therapist = 5,
#'                               n_sessions = 4, icc = 0.2, tau = 0),
#'               fit_fns = list(naiv = fit_z_naive))
#' mc_summary(res, truth = 0)
#'
#' @export
mc_run <- function(n_reps, sim_args, fit_fns,
                   seed_start = 1, alpha = 0.05) {
  stopifnot(n_reps >= 1, is.list(sim_args), is.list(fit_fns),
            !is.null(names(fit_fns)), !"seed" %in% names(sim_args))
  rows <- lapply(seq_len(n_reps), function(i) {
    seed <- seed_start + i - 1
    s <- do.call(sim_stream, c(sim_args, list(seed = seed)))
    do.call(rbind, lapply(names(fit_fns), function(nm) {
      f <- suppressWarnings(suppressMessages(fit_fns[[nm]](s)))
      data.frame(rep = i, seed = seed, estimator = nm,
                 estimate = f[["estimate"]], se = f[["se"]],
                 significant = is_sig(f, alpha))
    }))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Monte-Carlo-Ergebnisse zusammenfassen
#'
#' @param res Ergebnis von [mc_run()].
#' @param truth Wahrer Parameterwert (fuer Bias und Coverage).
#' @return `data.frame` pro Estimator: `reject_rate` (bei `truth = 0`
#'   die Falsch-Positiv-Rate, sonst Power), `bias`, `coverage`
#'   (95%-Wald-Intervall).
#' @export
mc_summary <- function(res, truth) {
  agg <- function(d) data.frame(
    estimator   = d$estimator[1],
    n_reps      = nrow(d),
    reject_rate = mean(d$significant),
    bias        = mean(d$estimate - truth),
    coverage    = mean(d$estimate - 1.96 * d$se < truth &
                       truth < d$estimate + 1.96 * d$se)
  )
  out <- do.call(rbind, lapply(split(res, res$estimator), agg))
  rownames(out) <- NULL
  out
}
