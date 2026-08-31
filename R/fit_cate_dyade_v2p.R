#' Dyadischer BCF v2P: Flaeche plus draw-weise Zerlegungs-Projektion
#'
#' Fit identisch zu [fit_cate_bcf_dyade()] (explizites
#' Produkt-Feature gegen die BCF-tau-Schrumpfung), danach wird JEDER
#' Posterior-Draw der tau-Flaeche am Gitter der Raender
#' doppelzentriert in tau0 + f(x) + g(c) + h(x, c) zerlegt. Kein
#' Refit, keine Residuen: Die Zerlegung ist eine deterministische
#' Projektion des Draws, die Unsicherheit propagiert exakt. Das
#' Personen-Ranking kommt unveraendert aus der Gesamt-Flaeche:
#' `tau_hat` ist mit [fit_cate_bcf_dyade()] bei gleichem Seed
#' identisch.
#'
#' Detektion ist AMPLITUDEN-basiert: Attribut `sd_h` traegt die
#' Posterior-Verteilung von SD_Gitter(h); Richtungs-/Korrelationsmasse
#' sind als Detektor wertlos (der Wald nutzt das Produkt-Feature auch
#' fuer Rauschen). Die Schwellen sind kalibriert (Detektion .07/.52,
#' Benennung .02/.77; siehe [v2p_entscheidung()]) -- die Schicht
#' selbst bleibt **exploratorisch**, und die Kalibrierung gilt fuer
#' statische Produkt-Form-Welten, in anderen Regimen als
#' Uebertragung. Regel-Anwendung: [v2p_entscheidung()]; Benennung
#' post hoc: [v2p_benennung()].
#'
#' Das Referenzgitter variiert nur `(x, c)`; alle uebrigen
#' Trainings-Features stehen auf 0. Das ist ihr Mittelwert nur, wenn
#' sie zentriert sind -- eine Eigenschaft der Datenquelle, nicht der
#' Funktion. Sind sie es nicht, warnt die Funktion mit der
#' Bedingungsklasse `windkanal_warnung_gitter_nullfeature` und nennt
#' die Spalten; die Zerlegung bleibt gueltig, die Flaeche liegt dann
#' nur an einem anderen Punkt dieser Features als im Stichprobenmittel.
#' Betroffen sind allein die `x_noise*`-Spalten, die
#' [cate_features()] neben `x` und `therapist_c` durchreicht. Aus
#' [sim_stream()] sind die per Konstruktion zentriert (`rnorm`), die
#' Wache zielt also auf Stroeme, deren Rausch-/Kovariatenspalten der
#' Aufrufer selbst gesetzt oder ueberschrieben hat (Plasmode, echte
#' Kovariaten).
#'
#' @param stream Datenstrom aus [sim_stream()] (braucht therapist_c).
#' @param nburn,nsim MCMC burn-in / Draws (Default 500/500).
#' @param pihat `"mean"` (Default) oder `"glm"` -- dieselben zwei
#'   Schaetzer wie in [fit_cate_bcf_dyade()]; die erste Option heisst
#'   dort `"constant"`.
#' @param gitter Referenzgitter der Projektion: `"quantil"` (Default;
#'   Quantilgitter `gitter_n` x `gitter_n` der beiden Raender) oder
#'   `"empirisch"` (Produkt der empirischen Raender: ALLE beobachteten
#'   x_i x ALLE Therapeuten-c_j). `"empirisch"` liefert zusaetzlich
#'   `h_real` (h an den realisierten Paaren, draw-weise) und ist
#'   Voraussetzung fuer `kandidaten_prod`.
#' @param gitter_n Quantilgitter-Punkte je Rand (Default 30; nur
#'   `gitter = "quantil"`). Einzelne endliche **ganze** Zahl >= 2: ein
#'   gebrochener Wert setzt die Randpunkte auf `(1:gitter_n)/(gitter_n
#'   + 1)` und damit auf andere Quantile als gemeint.
#' @param kandidaten Optional: benannte Liste von Funktionen
#'   `function(x, c)` (Form-Kandidaten der Benennung, z. B.
#'   `list(prod = function(x, c) x * c)`). Jeder Kandidat wird am
#'   Gitter doppelzentriert und h draw-weise darauf regressiert;
#'   Ergebnis im Attribut `benennung` (`beta`-Draws, `p_gross`).
#'   Kandidaten mit identischer h-Komponente sind nach der
#'   Doppelzentrierung exakt kollinear (dz(-(x-c)^2) = 2*dz(x*c))
#'   und brechen mit lautem Fehler ab.
#' @param kandidaten_prod Optional: Attribut-Produkt-Bibliothek
#'   `list(patient = Matrix n x A, therapeut = Matrix J x B)` (beide
#'   mit colnames) -- alle A*B Produkte a_i * b_j werden am
#'   empirischen Kreuz doppelzentriert und mitbenannt (Spaltenname
#'   `prod_a_b`). Zeilen-Kontrakt: `patient` in der Reihenfolge von
#'   [patients()], `therapeut` in Erstauftretens-Reihenfolge von
#'   `therapist_id`. Braucht `gitter = "empirisch"`.
#' @param halte_h h-Draws behalten (`zerlegung$h_draws`, Matrix
#'   Gitterpunkte x nsim) -- noetig fuer [v2p_benennung()] post hoc;
#'   Default `FALSE` (Speicher).
#' @param rfx_slope Zusaetzlich zum Therapeuten-Random-Intercept einen
#'   **Random Slope auf z** mitschaetzen (rfx-Basis `cbind(1, z)`):
#'   Therapeuten-Level-Heterogenitaet des Behandlungseffekts, die
#'   nicht ueber das gemessene c laeuft, wird dann vom rfx-Term
#'   aufgesogen statt in die (x, c)-Flaeche zu lecken. Default
#'   `FALSE` = bit-identisch zum Verhalten ohne Slope. Gemessen in
#'   einem Mikro-Experiment (16.07.): senkt die Fehlalarmrate der
#'   Beziehungs-Schicht von .16 auf .10 bei vollem Power-Erhalt --
#'   Teil-Mitigation, die Benennungs-Fehlalarme bleiben (Leckage
#'   sitzt auch auf Patienten-Ebene).
#' @param seed Zufalls-Seed (Pflicht).
#' @return `data.frame` wie [fit_cate_bcf_dyade()] (`tau_hat` etc.),
#'   zusaetzlich Attribute `zerlegung` (Liste: `tau0` Draws, `sd_h`
#'   Draws, `h_mean` Gittermatrix, `f_mean`, `g_mean`, `gitter`,
#'   `typ`; bei `gitter = "empirisch"` zusaetzlich `h_real_mean` und
#'   `therapist_of`; mit `halte_h` zusaetzlich `h_draws`), `ate`,
#'   bei `gitter = "empirisch"` `h_real` (Matrix n x nsim, wahres
#'   h-Pendant: [dyade_h_wahr()] mit `therapist_of`) und ggf.
#'   `benennung`.
#' @examples
#' \donttest{
#' if (requireNamespace("stochtree", quietly = TRUE)) {
#'   s <- sim_stream(n_therapists = 8, patients_per_therapist = 6,
#'                   n_sessions = 2, z_level = "patient",
#'                   tau = 0.5, tau_xc = 0.35, seed = 5)
#'   b <- fit_cate_dyade_v2p(s, nburn = 50, nsim = 25, gitter_n = 12,
#'                           seed = 1)
#'   v2p_entscheidung(b)$p_sd
#' }
#' }
#' @export
fit_cate_dyade_v2p <- function(stream, nburn = 500, nsim = 500,
                               pihat = c("mean", "glm"),
                               gitter = c("quantil", "empirisch"),
                               gitter_n = 30, kandidaten = NULL,
                               kandidaten_prod = NULL,
                               halte_h = FALSE, rfx_slope = FALSE,
                               seed) {
  if (!requireNamespace("stochtree", quietly = TRUE)) {
    stop("fit_cate_dyade_v2p() requires 'stochtree'.", call. = FALSE)
  }
  seed_pruefen(seed, missing(seed), "fit_cate_dyade_v2p()",
               "MCMC is stochastic")
  pihat <- match.arg(pihat)
  gitter <- match.arg(gitter)
  if (gitter == "quantil") {
    if (length(gitter_n) != 1L || !is.numeric(gitter_n) ||
        !is.finite(gitter_n) || gitter_n != round(gitter_n) ||
        gitter_n < 2) {
      stop("`gitter_n` muss eine einzelne endliche ganze Zahl >= 2 ",
           "sein (Doppelzentrierung braucht mindestens 2 Punkte je ",
           "Rand; ein gebrochener Wert erzeugt still ein unsymmetrisch ",
           "gesetztes Referenzgitter, weil (1:gitter_n)/(gitter_n + 1) ",
           "dann nicht mehr die gemeinten Quantile trifft).",
           call. = FALSE)
    }
  }
  if (!is.null(kandidaten_prod) && gitter != "empirisch") {
    stop("`kandidaten_prod` (Attribut-Produkte) braucht ",
         "gitter = \"empirisch\".", call. = FALSE)
  }
  p <- patients(stream)
  z_binaer_pruefen(p, "fit_cate_dyade_v2p()")
  if (is.null(p$therapist_c)) {
    stop("v2P braucht `therapist_c` im Stream (tau_c/tau_xc-Design).",
         call. = FALSE)
  }
  x <- as.data.frame(cate_features(p))
  x$passung <- p$x * p$therapist_c
  rest <- setdiff(names(x), c("x", "therapist_c", "passung"))
  if (length(rest)) {
    m <- colMeans(x[rest], na.rm = TRUE)
    s <- vapply(x[rest], stats::sd, numeric(1), na.rm = TRUE)
    schief <- rest[which(abs(m) > pmax(0.25 * s, 4 * s / sqrt(nrow(x))))]
    if (length(schief)) {
      warning(structure(
        class = c("windkanal_warnung_gitter_nullfeature", "warning",
                  "condition"),
        list(call = NULL, message = paste0(
          "Die Projektion haelt alle Features ausser x, therapist_c ",
          "und passung auf 0 fest -- das ist nur bei zentrierten ",
          "Features ihr Mittelwert. Nicht zentriert (Mittel weiter als ",
          "0.25 SD und als 4 SE von 0 entfernt): ",
          paste0(schief, " (Mittel ",
                 formatC(m[schief], format = "f", digits = 3), ", SD ",
                 formatC(s[schief], format = "f", digits = 3), ")",
                 collapse = ", "),
          ". Die Zerlegung bleibt gueltig, die Flaeche liegt aber an ",
          "einem anderen Punkt dieser Features als im Stichproben",
          "mittel. Wer die Referenz am Mittel will, zentriert die ",
          "Features vor dem Fit."))))
    }
  }
  ph <- if (pihat == "glm") {
    stats::glm(z ~ x, family = stats::binomial, data = p)$fitted.values
  } else {
    rep(mean(p$z), nrow(p))
  }
  rfx_basis <- if (rfx_slope) cbind(1, p$z) else matrix(1, nrow(p), 1)
  fit <- do.call(stochtree::bcf, list(
    X_train = x, Z_train = p$z, y_train = p$score_mean,
    propensity_train = ph,
    num_burnin = nburn, num_mcmc = nsim,
    general_params = list(random_seed = seed, verbose = FALSE),
    rfx_group_ids_train = p$therapist_id,
    rfx_basis_train = rfx_basis))
  draws <- fit$tau_hat_train

  # Gitter = Produkt der Raender; alle uebrigen Trainings-Features
  # auf 0, damit nur (x, c) variiert -- das ist ihr Mittel nur, wenn
  # sie zentriert sind (Wache oben)
  if (gitter == "quantil") {
    qx <- unname(stats::quantile(p$x, (1:gitter_n) / (gitter_n + 1)))
    qc <- unname(stats::quantile(p$therapist_c,
                                 (1:gitter_n) / (gitter_n + 1)))
    therapist_of <- NULL
  } else {
    # empirisch: Zeilen = Personen (Reihenfolge patients()), Spalten =
    # Therapeuten in Erstauftretens-Reihenfolge -- der Kontrakt fuer
    # kandidaten_prod und h_real
    tid_levels <- unique(p$therapist_id)
    therapist_of <- match(p$therapist_id, tid_levels)
    qx <- p$x
    qc <- p$therapist_c[match(tid_levels, p$therapist_id)]
  }
  nx <- length(qx); nc_ <- length(qc)
  gr <- expand.grid(x = qx, therapist_c = qc)   # x variiert am schnellsten
  gx <- as.data.frame(matrix(0, nrow(gr), ncol(x),
                             dimnames = list(NULL, names(x))))
  gx$x <- gr$x; gx$therapist_c <- gr$therapist_c
  gx$passung <- gr$x * gr$therapist_c
  # rfx-Angaben sind fuer predict Pflicht, beruehren tau_hat aber nicht
  # (Random Intercepts wirken auf y_hat); Gruppe beliebig fixiert
  pr <- stats::predict(fit, X = gx, Z = rep(1L, nrow(gx)),
                propensity = rep(mean(p$z), nrow(gx)),
                rfx_group_ids = rep(p$therapist_id[1], nrow(gx)),
                rfx_basis = if (rfx_slope) cbind(1, rep(1, nrow(gx))) else
                  matrix(1, nrow(gx), 1))
  tg <- as.matrix(pr$tau_hat)                   # (nx * nc_) x nsim

  # Projektion: Doppelzentrierung pro Draw
  nd <- ncol(tg)
  tau0 <- numeric(nd); sd_h <- numeric(nd)
  Hd <- matrix(NA_real_, nx * nc_, nd)
  f_sum <- numeric(nx); g_sum <- numeric(nc_)
  h_real <- if (gitter == "empirisch") matrix(NA_real_, nx, nd) else NULL
  for (s in seq_len(nd)) {
    M <- matrix(tg[, s], nx, nc_)               # Zeilen = x, Spalten = c
    m0 <- mean(M); fx <- rowMeans(M) - m0; gc_ <- colMeans(M) - m0
    H <- M - outer(fx, gc_, "+") - m0
    tau0[s] <- m0; sd_h[s] <- stats::sd(H)
    Hd[, s] <- as.vector(H)
    f_sum <- f_sum + fx; g_sum <- g_sum + gc_
    if (!is.null(h_real)) {
      h_real[, s] <- H[cbind(seq_len(nx), therapist_of)]
    }
  }

  out <- data.frame(
    patient_id = p$patient_id, x = p$x, z = p$z,
    tau_hat = rowMeans(draws),
    tau_lo  = apply(draws, 1, stats::quantile, 0.025),
    tau_hi  = apply(draws, 1, stats::quantile, 0.975))
  ate_draws <- colMeans(draws)
  attr(out, "ate") <- c(estimate = mean(ate_draws),
                        se = stats::sd(ate_draws),
                        lo = unname(stats::quantile(ate_draws, 0.025)),
                        hi = unname(stats::quantile(ate_draws, 0.975)))
  zerlegung <- list(
    tau0 = tau0, sd_h = sd_h,
    h_mean = matrix(rowMeans(Hd), nx, nc_),
    f_mean = f_sum / nd, g_mean = g_sum / nd,
    gitter = list(x = qx, c = qc), typ = gitter, engine = "v2p")
  if (gitter == "empirisch") {
    zerlegung$h_real_mean <- rowMeans(h_real)
    zerlegung$therapist_of <- therapist_of
    attr(out, "h_real") <- h_real
  }
  if (halte_h) zerlegung$h_draws <- Hd
  attr(out, "zerlegung") <- zerlegung

  # Benennungs-Schicht: Kandidaten am Gitter DOPPELZENTRIEREN (sonst
  # traefe die Regression auch deren Haupteffekt-Anteile), dann OLS
  # von h^(s) auf die zentrierten Kandidaten PRO DRAW -> volle
  # beta-Posterior. Selektions-Schwellen kalibriert (.02/.77),
  # Schicht [exploratory].
  if (!is.null(kandidaten) || !is.null(kandidaten_prod)) {
    attr(out, "benennung") <- benennung_bauen(qx, qc, Hd, kandidaten,
                                              kandidaten_prod)
  }
  out
}
