#' Heterogenitaets-Test des Causal Forest (Kalibrierungstest)
#'
#' Beantwortet die Frage, die rohe `tau_hat`-Streuung NICHT
#' beantworten kann: Gibt es hier wirklich Effekt-Heterogenitaet?
#' Nutzt `grf::test_calibration` -- der Koeffizient der
#' "differential forest prediction" misst, ob die vom Forest
#' behauptete Heterogenitaet aus den Daten reproduzierbar ist
#' (~1 = echt, ~0 = Rauschen); der p-Wert testet sie.
#'
#' Standard-`fit_*`-Interface (`estimate`, `se`, `p`) -- damit laeuft
#' der Test direkt in [mc_run()], [run_peek()] und [run_gates()]:
#' Personalisierungs-Behauptungen werden governance-pruefbar wie
#' Haupteffekte.
#'
#' Verlangt binaeres Treatment: `grf::causal_forest` akzeptiert zwar
#' stetiges W (partieller Effekt), aber der Kalibrierungstest ist im
#' Paket nur fuer den Zwei-Arm-Kontrast validiert -- Dosis-Welten
#' (`z_type = "dose"`) werden mit informativem Fehler abgewiesen.
#'
#' @param snap Snapshot/Datenstrom aus [sim_stream()].
#' @param num_trees Baeume (Default 500 -- MC-tauglich schnell).
#' @param forest_seed Interner Forest-Seed (Default 1, deterministisch;
#'   die Reproduzierbarkeit des Experiments traegt der Stream-Seed).
#' @return Benannter Vektor `c(estimate, se, p)` fuer die
#'   differentielle Forest-Vorhersage.
#' @export
fit_het_grf <- function(snap, num_trees = 500, forest_seed = 1L) {
  if (!requireNamespace("grf", quietly = TRUE)) {
    stop("fit_het_grf() braucht das Paket 'grf'.", call. = FALSE)
  }
  p <- patients(snap)
  z_binaer_pruefen(p, "fit_het_grf()")
  cf <- grf::causal_forest(
    X = matrix(p$x, ncol = 1), Y = p$score_mean, W = p$z,
    clusters = p$therapist_id,
    num.trees = num_trees, seed = forest_seed
  )
  tc <- grf::test_calibration(cf)
  row <- "differential.forest.prediction"
  c(estimate = tc[row, "Estimate"],
    se       = tc[row, "Std. Error"],
    p        = tc[row, "Pr(>t)"])
}
