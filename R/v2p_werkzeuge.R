# Werkzeuge der v2P-Schicht: Doppelzentrierung, Kandidaten-Matrix,
# Benennung post hoc, Entscheidungsregel.

# Doppelzentrierung einer Gittermatrix (Zeilen = x-Rand, Spalten =
# c-Rand): M - Zeilenmittel - Spaltenmittel + Gesamtmittel.
dz_gitter <- function(M) {
  M - outer(rowMeans(M), colMeans(M), "+") + mean(M)
}

# Kandidaten-Matrix am Gitter (Spalten doppelzentriert) mit
# Kollinearitaets-Check. Zwei Kandidaten-Arten:
# * `kandidaten`: benannte Liste von Funktionen function(x, c) --
#   Form-Kandidaten, auf jedem Gitter auswertbar.
# * `kandidaten_prod`: list(patient = Matrix nx x A, therapeut =
#   Matrix nc x B, beide mit colnames) -- Attribut-Produkte a_i * b_j
#   am Kreuz; Zeilen muessen zu den Gitter-Raendern passen
#   (praktisch: gitter = "empirisch").
# Kollinearitaets-Check: Nach Doppelzentrierung sind Kandidaten exakt
# linear abhaengig, sobald ihre h-Komponenten uebereinstimmen --
# dz(-(x-c)^2) = 2 * dz(x*c) auf JEDEM Gitter. lm.fit liefert dann
# stumm NA-Koeffizienten; stattdessen brechen wir laut ab und nennen
# die Spalten.
v2p_kandidaten_matrix <- function(qx, qc, kandidaten = NULL,
                                  kandidaten_prod = NULL) {
  spalten <- list()
  if (!is.null(kandidaten)) {
    stopifnot(is.list(kandidaten), length(names(kandidaten)) > 0,
              all(vapply(kandidaten, is.function, TRUE)))
    for (nm in names(kandidaten)) {
      spalten[[nm]] <- as.vector(dz_gitter(outer(qx, qc, kandidaten[[nm]])))
    }
  }
  if (!is.null(kandidaten_prod)) {
    stopifnot(is.list(kandidaten_prod),
              !is.null(kandidaten_prod$patient),
              !is.null(kandidaten_prod$therapeut))
    pat <- as.matrix(kandidaten_prod$patient)
    ther <- as.matrix(kandidaten_prod$therapeut)
    if (is.null(colnames(pat)) || is.null(colnames(ther))) {
      stop("kandidaten_prod: patient/therapeut brauchen colnames.",
           call. = FALSE)
    }
    if (nrow(pat) != length(qx) || nrow(ther) != length(qc)) {
      stop(sprintf(paste0(
        "kandidaten_prod passt nicht zum Gitter: nrow(patient) = %d ",
        "vs. %d x-Randpunkte, nrow(therapeut) = %d vs. %d c-Randpunkte. ",
        "Attribut-Produkte brauchen gitter = \"empirisch\" (Raender = ",
        "Personen/Therapeuten)."),
        nrow(pat), length(qx), nrow(ther), length(qc)), call. = FALSE)
    }
    for (a in colnames(pat)) for (b in colnames(ther)) {
      spalten[[paste0("prod_", a, "_", b)]] <-
        as.vector(dz_gitter(outer(pat[, a], ther[, b])))
    }
  }
  if (!length(spalten)) return(NULL)
  K <- do.call(cbind, spalten)
  qrK <- qr(K)
  if (qrK$rank < ncol(K)) {
    schuldig <- colnames(K)[qrK$pivot[(qrK$rank + 1):ncol(K)]]
    stop(paste0(
      "Kandidaten-Bibliothek nach Doppelzentrierung kollinear ",
      "(Rang ", qrK$rank, " von ", ncol(K), "); abhaengig: ",
      paste(schuldig, collapse = ", "), ". Kandidaten mit identischer ",
      "h-Komponente (z. B. quadratische Kongruenz UND Produkt: ",
      "dz(-(x-c)^2) = 2*dz(x*c)) sind nicht gemeinsam benennbar -- ",
      "einen davon entfernen."), call. = FALSE)
  }
  K
}

# Draw-weise Benennungs-Regression: OLS von h^(s) auf die
# doppelzentrierten Kandidaten, pro Draw -> volle beta-Posterior.
v2p_beta_draws <- function(K, h_draws) {
  nd <- ncol(h_draws)
  beta <- matrix(NA_real_, ncol(K), nd,
                 dimnames = list(colnames(K), NULL))
  for (s in seq_len(nd)) {
    beta[, s] <- stats::coef(stats::lm.fit(K, h_draws[, s]))
  }
  beta
}

benennung_bauen <- function(qx, qc, h_draws, kandidaten = NULL,
                            kandidaten_prod = NULL) {
  K <- v2p_kandidaten_matrix(qx, qc, kandidaten, kandidaten_prod)
  if (is.null(K)) return(NULL)
  beta <- v2p_beta_draws(K, h_draws)
  list(beta = beta,
       p_gross = rowMeans(abs(beta) > 0.02),
       hinweis = paste("p_gross bei delta = .02 (kalibrierte",
                       "v2P-Schwelle); Entscheidungen ueber",
                       "v2p_entscheidung()"))
}

#' Benennung post hoc auf einem v2P-Fit
#'
#' Rechnet die draw-weise Kandidaten-Regression der Benennungs-Schicht
#' NACH dem Fit -- etwa um mehrere Kandidaten-Bibliotheken (K = 5 vs.
#' K = 25) auf demselben Fit zu vergleichen, ohne den BCF neu zu
#' ziehen. Braucht die h-Draws: [fit_cate_dyade_v2p()] mit
#' `halte_h = TRUE` aufrufen.
#'
#' @param fit Ergebnis von [fit_cate_dyade_v2p()] mit `halte_h = TRUE`.
#' @param kandidaten Benannte Liste von Funktionen `function(x, c)`
#'   (Form-Kandidaten), wie in [fit_cate_dyade_v2p()].
#' @param kandidaten_prod Attribut-Produkt-Bibliothek
#'   `list(patient = , therapeut = )`, wie in [fit_cate_dyade_v2p()].
#' @return Liste `beta` (Draws, Kandidaten x nsim), `p_gross`
#'   (P(|beta| > .02), die auf statischen Produkt-Form-Welten
#'   kalibrierte Benennungs-Schwelle; ausserhalb dieses Regimes eine
#'   Uebertragung, keine Garantie -- die Benennungs-Schicht ist
#'   exploratorisch), `hinweis`.
#' @examples
#' \donttest{
#' if (requireNamespace("stochtree", quietly = TRUE)) {
#'   s <- sim_stream(n_therapists = 8, patients_per_therapist = 6,
#'                   n_sessions = 2, z_level = "patient",
#'                   tau = 0.5, tau_xc = 0.35, seed = 5)
#'   b <- fit_cate_dyade_v2p(s, nburn = 50, nsim = 25, gitter_n = 12,
#'                           halte_h = TRUE, seed = 1)
#'   # zwei Bibliotheken auf DEMSELBEN Fit, ohne den BCF neu zu ziehen
#'   v2p_benennung(b, kandidaten = list(prod = function(x, c) x * c))$p_gross
#'   v2p_benennung(b, kandidaten = list(
#'     prod  = function(x, c) x * c,
#'     rausch = function(x, c) sin(3 * x) * c^2))$p_gross
#' }
#' }
#' @export
v2p_benennung <- function(fit, kandidaten = NULL, kandidaten_prod = NULL) {
  z <- attr(fit, "zerlegung")
  if (is.null(z)) {
    stop("`fit` traegt keine `zerlegung` (kein v2P-Ergebnis?).",
         call. = FALSE)
  }
  if (is.null(z$h_draws)) {
    stop("Keine h-Draws im Fit: fit_cate_dyade_v2p(..., halte_h = TRUE) ",
         "aufrufen.", call. = FALSE)
  }
  b <- benennung_bauen(z$gitter$x, z$gitter$c, z$h_draws,
                       kandidaten, kandidaten_prod)
  if (is.null(b)) {
    stop("Keine Kandidaten uebergeben.", call. = FALSE)
  }
  b
}

#' Entscheidungsregel der v2P-Schicht: Detektion und Benennung
#'
#' Wendet die zweiteilige Regel auf einen v2P-Fit an: **Detektion**
#' amplitudenbasiert ueber die Posterior-Wahrscheinlichkeit
#' P(SD_Gitter(h) > delta_sd) > cut_sd (Richtungs-/Korrelationsmasse
#' sind als Detektor wertlos -- unter der Null ist cor(h_mean, xc)
#' trotzdem ~.88); **Benennung** ueber P(|beta_k| > delta_beta) >
#' cut_beta je Kandidat.
#'
#' Die Defaults sind ein kalibriertes Regelpaar (Kalibrier-Experiment
#' 2026-07-16, 400 Welten): Detektion (.07, .52) mit Power .48 bei
#' Fehlalarm <= .05 in jeder Kalibrier-Zelle; Benennung (.02, .77) mit
#' Power .50 (~ Orakel-Niveau .52) bei Fehlalarm <= .05. Kalibriert
#' auf statischen Produkt-Form-Welten (Caseload 5/10, K5/K25); fuer
#' andere Regime gilt die Regel als Uebertragung, nicht als Garantie.
#' Die Schicht ist **exploratorisch**.
#'
#' @param fit Ergebnis von [fit_cate_dyade_v2p()].
#' @param delta_sd,cut_sd Detektions-Schwelle und Posterior-Cut
#'   (Default: kalibriert .07/.52).
#' @param delta_beta,cut_beta Benennungs-Schwelle und Cut (Default:
#'   kalibriert .02/.77; Benennung nur mit Benennungs-Schicht).
#' @param benennung Optional: Ergebnis von [v2p_benennung()]; Default
#'   ist die Benennung aus dem Fit (`kandidaten`-Argument).
#' @return Liste: `detektiert` (logisch), `p_sd`, und -- falls eine
#'   Benennung vorliegt -- `benannt` (Namen der selektierten
#'   Kandidaten) und `p_beta`.
#' @examples
#' \donttest{
#' if (requireNamespace("stochtree", quietly = TRUE)) {
#'   s <- sim_stream(n_therapists = 8, patients_per_therapist = 6,
#'                   n_sessions = 2, z_level = "patient",
#'                   tau = 0.5, tau_xc = 0.35, seed = 5)
#'   b <- fit_cate_dyade_v2p(s, nburn = 50, nsim = 25, gitter_n = 12,
#'                           kandidaten = list(prod = function(x, c) x * c),
#'                           seed = 1)
#'   e <- v2p_entscheidung(b)      # kalibrierte Defaults .07/.52, .02/.77
#'   e$detektiert                  # Detektion: P(SD_Gitter(h) > .07) > .52
#'   e$p_sd
#'   e$benannt                     # Benennung: P(|beta| > .02) > .77
#' }
#' }
#' @export
v2p_entscheidung <- function(fit, delta_sd = 0.07, cut_sd = 0.52,
                             delta_beta = 0.02, cut_beta = 0.77,
                             benennung = NULL) {
  z <- attr(fit, "zerlegung")
  if (is.null(z)) {
    stop("`fit` traegt keine `zerlegung` (kein v2P-Ergebnis?).",
         call. = FALSE)
  }
  p_sd <- mean(z$sd_h > delta_sd)
  out <- list(detektiert = p_sd > cut_sd, p_sd = p_sd)
  if (is.null(benennung)) benennung <- attr(fit, "benennung")
  if (!is.null(benennung) && !is.null(delta_beta) && !is.null(cut_beta)) {
    p_beta <- rowMeans(abs(benennung$beta) > delta_beta)
    out$p_beta <- p_beta
    out$benannt <- names(which(p_beta > cut_beta))
  }
  out
}
