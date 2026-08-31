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
#' **Der Test sieht nur, wonach er sucht.** Bis einschliesslich
#' Version 0.2.1 bekam der Forest hier ausschliesslich `x`
#' (`matrix(p$x, ncol = 1)`), waehrend [fit_cate_grf()] ueber
#' `cate_features()` zusaetzlich `therapist_c` und die
#' `x_noise*`-Spalten sieht. Eine rein **dyadische** Heterogenitaet
#' (`tau_xc != 0` bei `tau_x = 0`) lag damit ausserhalb der
#' Merkmalsmenge und wurde als sauberes Null gemeldet -- in zehn
#' Welten mit wahrem `SD(tau_i) = 1.53` null Treffer, waehrend
#' dieselbe Staerke entlang `x` in fuenf von fuenf Welten mit
#' Schaetzer ~1.00 erkannt wurde. Default ist deshalb seit 0.3.0
#' `features = "alle"`, also dieselbe Merkmalsmenge wie
#' [fit_cate_grf()]. `features = "nur_x"` haelt die alte, enge
#' Variante als **Referenz-Naivling** verfuegbar: sie beantwortet die
#' Frage nur entlang `x`, eine rein dyadische Heterogenitaet bleibt
#' unentdeckt.
#'
#' @param snap Snapshot/Datenstrom aus [sim_stream()].
#' @param num_trees Baeume (Default 500 -- MC-tauglich schnell).
#' @param forest_seed Interner Forest-Seed (Default 1, deterministisch;
#'   die Reproduzierbarkeit des Experiments traegt der Stream-Seed).
#' @param features Merkmalsmenge des Forests: `"alle"` (Default) =
#'   `cate_features()`, also `x`, `therapist_c` (falls der Strom es
#'   traegt) und alle `x_noise*`; `"nur_x"` = allein `x` (das
#'   Verhalten bis 0.2.1, als deklarierter Naivling).
#' @return Benannter Vektor `c(estimate, se, p)` fuer die
#'   differentielle Forest-Vorhersage, mit den Attributen `features`
#'   (gewaehlte Menge) und `feature_namen`.
#' @examples
#' \donttest{
#' if (requireNamespace("grf", quietly = TRUE)) {
#'   s <- sim_stream(n_therapists = 20, patients_per_therapist = 10,
#'                   n_sessions = 4, z_level = "patient",
#'                   tau = 0.5, tau_x = 0.8, seed = 1)
#'   fit_het_grf(s, num_trees = 200)
#' }
#' }
#' @export
fit_het_grf <- function(snap, num_trees = 500, forest_seed = 1L,
                        features = c("alle", "nur_x")) {
  if (!requireNamespace("grf", quietly = TRUE)) {
    stop("fit_het_grf() braucht das Paket 'grf'.", call. = FALSE)
  }
  features <- match.arg(features)
  p <- patients(snap)
  z_binaer_pruefen(p, "fit_het_grf()")
  X <- if (features == "alle") {
    cate_features(p)
  } else {
    matrix(p$x, ncol = 1, dimnames = list(NULL, "x"))
  }
  cf <- grf::causal_forest(
    X = X, Y = p$score_mean, W = p$z,
    clusters = p$therapist_id,
    num.trees = num_trees, seed = forest_seed
  )
  tc <- grf::test_calibration(cf)
  row <- "differential.forest.prediction"
  out <- c(estimate = tc[row, "Estimate"],
           se       = tc[row, "Std. Error"],
           p        = tc[row, "Pr(>t)"])
  attr(out, "features") <- features
  attr(out, "feature_namen") <- colnames(X)
  out
}
