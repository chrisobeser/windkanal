#' Estimate personalized effects (mixed-effects random forest, MERF)
#'
#' Pruefling 16: the mixed-effects random forest (Hajjem et al., 2014,
#' J Stat Comput Simul; Implementierung: LongituRF, CRAN-Archiv 0.9) --
#' die MERF-Schiene der Li-et-al.-Familie. EM-artiger Wechsel zwischen
#' einem random forest fuer den festen Teil f(features, z) und BLUPs
#' fuer Therapeuten-Random-Intercepts. Individuelle Effekte sind
#' kontrafaktische Differenzen der Forest-Vorhersagen (S-Learner-Stil;
#' der Random Intercept kuerzt sich innerhalb des Therapeuten).
#'
#' `sto = "none"` (kein stochastischer Prozess; wir arbeiten auf
#' Personen-Mitteln, nicht auf Verlaeufen); `time` ist dann nur ein
#' Pflicht-Platzhalter (Innerhalb-Cluster-Index).
#'
#' @param stream Data stream from [sim_stream()].
#' @param ntree Baeume des Forests (default 300).
#' @param iter Maximale EM-Iterationen (default 20, LongituRF-nah).
#' @param B Bootstrap replicates for the ATE interval (default 200;
#'   therapist-cluster bootstrap as for the other learners).
#' @param seed Random seed (mandatory).
#' @return `data.frame` (`patient_id`, `x`, `z`, `tau_hat`) with
#'   attribute `ate` (`estimate`, `se`, `lo`, `hi`; cluster
#'   bootstrap).
#' @export
fit_cate_merf <- function(stream, ntree = 300, iter = 20, B = 200,
                          seed) {
  for (pkg in c("LongituRF", "randomForest")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("fit_cate_merf() requires the '", pkg, "' package.",
           call. = FALSE)
    }
  }
  if (missing(seed)) stop("`seed` is mandatory.", call. = FALSE)
  set.seed(seed)
  p <- patients(stream)
  z_binaer_pruefen(p, "fit_cate_merf()")
  tau_fun <- function(d) {
    X <- as.data.frame(cbind(cate_features(d), z = d$z))
    fit <- NULL
    # LongituRF druckt Konvergenz-Meldungen via print(); stumm schalten
    invisible(utils::capture.output(
      fit <- suppressWarnings(LongituRF::MERF(
        X = X, Y = d$score_mean,
        Z = matrix(1, nrow(d), 1),
        id = d$therapist_id,
        time = stats::ave(seq_len(nrow(d)), d$therapist_id,
                          FUN = seq_along),
        sto = "none", ntree = ntree, iter = iter))))
    X1 <- X; X1$z <- 1
    X0 <- X; X0$z <- 0
    as.numeric(stats::predict(fit$forest, X1) -
               stats::predict(fit$forest, X0))
  }
  tau_hat <- tau_fun(p)
  out <- data.frame(patient_id = p$patient_id, x = p$x, z = p$z,
                    tau_hat = tau_hat)
  boot <- cluster_boot_ate(p, tau_fun, B)
  attr(out, "ate") <- c(estimate = mean(tau_hat), boot)
  out
}
