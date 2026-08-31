#' Latenten Score auf eine begrenzte Fragebogen-Skala abbilden
#'
#' Bildet den stetigen (latenten) `score` eines Datenstroms auf eine
#' begrenzte, ganzzahlige Summen-Skala ab (BDI/PHQ-Stil): optionale
#' monotone Schiefe-Transformation, affine Verankerung, Zensierung an
#' Boden und Decke, Rundung auf ganze Punkte. Die latente Spalte
#' `score` bleibt unveraendert erhalten; die Schicht ist Opt-in --
#' ohne Aufruf aendert sich keine bestehende Welt.
#'
#' Die beobachtete Skala ist der Default-Massstab fuer Audits auf
#' begrenzten Welten. Diese Schicht DEKLARIERT ihn (Attribut
#' `estimand_default`); die Verdrahtung der Schaetzer-Einstiege
#' (`fit_*`, [mc_run()]) auf `score_obs` ist ein eigener, noch
#' offener Posten -- bis dahin muessen Audits `score_obs` explizit
#' waehlen. Die latente Wahrheit wird stets mit ausgewiesen (siehe
#' [wahrheit_skala()]).
#'
#' Abbildung je Wert: `z = score - mitte`; bei `schiefe != 0` der
#' monotone Tilt `z = expm1(schiefe * z) / schiefe` (Grenzfall
#' `schiefe -> 0` ist die Identitaet; `schiefe > 0` streckt die obere
#' und staucht die untere Flanke, auf der Skala also Rechtsschiefe);
#' dann `roh = (minimum + maximum)/2 + punkte_pro_einheit * z`,
#' zensiert auf `[minimum, maximum]`, gerundet auf ganze Punkte.
#'
#' @param stream Datenstrom aus [sim_stream()] (oder ein
#'   `data.frame` mit Spalte `score`).
#' @param minimum Untere Skalengrenze (ganzzahlig, Default 0).
#' @param maximum Obere Skalengrenze (ganzzahlig, Default 27 --
#'   PHQ-9-Format).
#' @param schiefe Schiefe-Parameter des monotonen Tilts (Default 0 =
#'   keine Schiefe). Geltungsgrenze: fuer `schiefe > 0` ist die
#'   Abbildung nach unten hart beschraenkt
#'   (`roh > (minimum + maximum)/2 - punkte_pro_einheit/schiefe`);
#'   die untere Skalengrenze wird strukturell unerreichbar, sobald
#'   `punkte_pro_einheit <= schiefe * ((maximum - minimum)/2 - 0.5)`
#'   -- spiegelbildlich fuer `schiefe < 0` die obere. [skala_werte()]
#'   warnt in diesem Fall. Fuer Boden-/Deckeneffekte MIT Schiefe die
#'   Verankerung `mitte` verschieben statt `schiefe` zu erhoehen.
#' @param mitte Latenter Wert, der auf die Skalenmitte abgebildet
#'   wird (Default 0 -- die Streams sind um 0 zentriert).
#' @param punkte_pro_einheit Skalenpunkte je latenter Einheit vor
#'   Zensierung (Default `(maximum - minimum) / 6`: +-3 latente
#'   Einheiten spannen die Skala auf).
#' @return Der Datenstrom mit zusaetzlicher ganzzahliger Spalte
#'   `score_obs` in `[minimum, maximum]` und Attribut `skala`
#'   (Abbildungsparameter samt `estimand_default = "beobachtet"`).
#'   Alle uebrigen Spalten, Attribute und die Klasse bleiben
#'   unveraendert.
#' @seealso [skala_werte()] fuer die reine Abbildung,
#'   [wahrheit_skala()] fuer die Wahrheit auf beiden Ebenen.
#' @examples
#' s <- sim_stream(n_therapists = 10, patients_per_therapist = 5,
#'                 n_sessions = 3, tau = 1, seed = 42)
#' b <- skala_begrenzen(s, minimum = 0, maximum = 27)
#' range(b$score_obs)
#' # Deckeneffekt: hohe Verankerung sammelt Masse am Maximum
#' d <- skala_begrenzen(s, minimum = 0, maximum = 27, mitte = -3)
#' mean(d$score_obs == 27)
#' @export
skala_begrenzen <- function(stream, minimum = 0, maximum = 27,
                            schiefe = 0, mitte = 0,
                            punkte_pro_einheit =
                              (maximum - minimum) / 6) {
  if (!is.data.frame(stream) || !("score" %in% names(stream))) {
    stop("skala_begrenzen() braucht einen Datenstrom mit Spalte ",
         "`score` (z. B. aus sim_stream()).", call. = FALSE)
  }
  if ("score_obs" %in% names(stream)) {
    warning("Spalte `score_obs` existiert bereits und wird ",
            "ueberschrieben.", call. = FALSE)
  }
  stream$score_obs <- skala_werte(stream$score, minimum, maximum,
                                  schiefe, mitte, punkte_pro_einheit)
  attr(stream, "skala") <- list(
    minimum = minimum, maximum = maximum, schiefe = schiefe,
    mitte = mitte, punkte_pro_einheit = punkte_pro_einheit,
    estimand_default = "beobachtet")
  stream
}

#' Reine Skalen-Abbildung: latente Werte zu begrenzten Punktwerten
#'
#' Die deterministische Kernabbildung hinter [skala_begrenzen()],
#' als eigene Funktion exportiert, damit Gegencheck-Skripte und
#' [wahrheit_skala()] exakt dieselbe Abbildung verwenden.
#'
#' @param x Numerischer Vektor latenter Werte.
#' @inheritParams skala_begrenzen
#' @return Ganzzahliger Vektor in `[minimum, maximum]`.
#' @examples
#' skala_werte(c(-2, 0, 2))
#' # unerreichbare Grenzen bei punkte_pro_einheit = 1: reine Rundung
#' skala_werte(c(-2.4, 0, 2.4), minimum = -1000, maximum = 1000,
#'             punkte_pro_einheit = 1)
#' @export
skala_werte <- function(x, minimum = 0, maximum = 27, schiefe = 0,
                        mitte = 0,
                        punkte_pro_einheit = (maximum - minimum) / 6) {
  stopifnot(is.numeric(x),
            is.numeric(minimum), length(minimum) == 1L,
            is.finite(minimum), minimum == round(minimum),
            is.numeric(maximum), length(maximum) == 1L,
            is.finite(maximum), maximum == round(maximum),
            minimum < maximum,
            is.numeric(schiefe), length(schiefe) == 1L,
            is.finite(schiefe),
            is.numeric(mitte), length(mitte) == 1L, is.finite(mitte),
            is.numeric(punkte_pro_einheit),
            length(punkte_pro_einheit) == 1L,
            is.finite(punkte_pro_einheit), punkte_pro_einheit > 0)
  if (schiefe > 0 &&
      (minimum + maximum) / 2 - punkte_pro_einheit / schiefe >=
        minimum + 0.5) {
    warning("Untere Skalengrenze strukturell unerreichbar (harte ",
            "Tilt-Schranke bei ",
            (minimum + maximum) / 2 - punkte_pro_einheit / schiefe,
            ") -- Bodeneffekt ist damit ausgeschaltet; ",
            "`mitte` verschieben oder `punkte_pro_einheit` erhoehen.",
            call. = FALSE)
  }
  if (schiefe < 0 &&
      (minimum + maximum) / 2 - punkte_pro_einheit / schiefe <=
        maximum - 0.5) {
    warning("Obere Skalengrenze strukturell unerreichbar (harte ",
            "Tilt-Schranke bei ",
            (minimum + maximum) / 2 - punkte_pro_einheit / schiefe,
            ") -- Deckeneffekt ist damit ausgeschaltet; ",
            "`mitte` verschieben oder `punkte_pro_einheit` erhoehen.",
            call. = FALSE)
  }
  z <- x - mitte
  if (schiefe != 0) z <- expm1(schiefe * z) / schiefe
  roh <- (minimum + maximum) / 2 + punkte_pro_einheit * z
  as.integer(round(pmin(maximum, pmax(minimum, roh))))
}

#' Wahrheit auf beiden Ebenen: latenter und beobachteter Effekt
#'
#' Berechnet fuer eine deklarierte Welt die wahren Behandlungseffekte
#' auf der latenten UND der begrenzten beobachteten Skala, per
#' kontrafaktischem Doppellauf: [sim_stream()] wird zweimal mit
#' identischem Seed gezogen (`z_force = 0` bzw. `1`); die Zuweisung
#' wird dabei normal gezogen (identischer RNG-Verbrauch) und erst
#' danach ueberschrieben, sodass beide Laeufe bis auf `z` identisch
#' sind (CRN) und je Zeile beide Potential Outcomes vorliegen. Drei
#' Wachen pruefen den Gleichlauf (`patient_id`, `session`, `x`
#' identisch), eine vierte den RNG-Endzustand beider Laeufe auf
#' z-abhaengigen Zufallszahlen-Verbrauch; schlaegt eine an, stoppt
#' der Lauf.
#'
#' Die Wahrheit ist als Populations-Estimand der DROPOUT-FREIEN Welt
#' definiert: `dropout = 0` wird intern gesetzt.
#' `dropout_informative` ist unter erzwungenem `dropout = 0`
#' wirkungslos (Warnung). Messfehler (`reliability_score < 1`) bleibt
#' Teil des beobachteten Estimands -- die Zensierung wirkt auf den
#' GEMESSENEN Wert --; er kuerzt sich im Doppellauf nicht heraus und
#' erhoeht die MC-Restunschaerfe der Wahrheitsrechnung. `p_treated`
#' und `confounding` duerfen als Weltparameter durchgereicht werden
#' -- sie formen nur die verworfene Zuweisungs-Ziehung und aendern
#' die Potential Outcomes nicht. Die x-Terzile werden auf dem
#' Strom-`x` gebildet, bei `reliability_x < 1` also auf dem
#' fehlerbehafteten x (dasselbe x, das ein Schaetzer saehe). Nur
#' `z_type = "binary"` wird unterstuetzt.
#'
#' @param ... Weltparameter fuer [sim_stream()] (ohne `dropout`,
#'   `z_force`, `seed` -- diese setzt die Funktion selbst).
#' @inheritParams skala_begrenzen
#' @param session Sitzung, fuer die die Wahrheit berechnet wird
#'   (Default: letzte Sitzung des Stroms).
#' @param seed Zufalls-Seed des kontrafaktischen Doppellaufs
#'   (Pflicht).
#' @return Liste: `session`; `ate_latent` und `ate_beobachtet`
#'   (mittlere Effekte je Skala); `effekte_x_terzile` (`data.frame`
#'   der bedingten Effekte je x-Terzil auf beiden Skalen);
#'   `spann_latent` (Spannweite der Terzil-Effekte in **latenten
#'   Einheiten**) und `spann_beobachtet` (dieselbe Spannweite in
#'   **Skalenpunkten**); `schein_heterogenitaet`; `skala` (Parameter).
#'
#'   ⚠ `spann_latent` und `spann_beobachtet` stehen in
#'   **verschiedenen Einheiten** -- zwischen ihnen liegt der Faktor
#'   `punkte_pro_einheit` (Default `(maximum - minimum) / 6`). Ihre
#'   **Differenz ist keine Kennzahl**, sondern ueberwiegend die
#'   Umrechnungskonstante: in einer gemessenen Welt meldet sie
#'   `+4.33`, waehrend die begrenzte Skala die Heterogenitaet dort in
#'   Wahrheit leicht **gestaucht** hat (`-0.09`) -- Vorzeichen und
#'   Groessenordnung beide falsch. Die Schein-Heterogenitaet der Skala
#'   ist die einheitenbereinigte Differenz
#'   `spann_beobachtet / punkte_pro_einheit - spann_latent`, also der
#'   Vergleich beider Spannweiten auf der **latenten** Skala; sie
#'   faehrt als `schein_heterogenitaet` im Rueckgabewert mit (in
#'   Skalenpunkten: mal `punkte_pro_einheit`). Positiv = die begrenzte
#'   Skala blaeht Heterogenitaet auf, negativ = sie staucht sie.
#'   ⚠ Die Umrechnung `/ punkte_pro_einheit` ist die **lineare**
#'   Ruecknahme der Skalen-Abbildung; unter aktivem `schiefe`-Tilt ist
#'   die Abbildung nicht linear, die Kennzahl ist dann eine
#'   Naeherung. Der zweite, von der Umrechnung unabhaengige Weg ist
#'   der Vergleich von `spann_beobachtet` gegen einen zweiten Lauf mit
#'   unerreichbaren Grenzen bei **gleichem** `punkte_pro_einheit`.
#' @examples
#' w <- wahrheit_skala(n_therapists = 20, patients_per_therapist = 10,
#'                     n_sessions = 3, z_level = "patient",
#'                     tau = 0.5, tau_x = 0.6, seed = 4242)
#' w$ate_latent
#' w$ate_beobachtet
#' w$effekte_x_terzile
#' # NICHT spann_beobachtet - spann_latent (verschiedene Einheiten):
#' w$schein_heterogenitaet
#' @export
wahrheit_skala <- function(..., minimum = 0, maximum = 27,
                           schiefe = 0, mitte = 0,
                           punkte_pro_einheit =
                             (maximum - minimum) / 6,
                           session = NULL, seed) {
  seed_pruefen(seed, missing(seed), "wahrheit_skala()",
               "the counterfactual twin run is seeded")
  args <- list(...)
  gesperrt <- intersect(names(args), c("dropout", "z_force", "seed"))
  if (length(gesperrt)) {
    stop("Diese Argumente setzt wahrheit_skala() selbst: ",
         paste(gesperrt, collapse = ", "), call. = FALSE)
  }
  if (identical(args$z_type, "dose")) {
    stop("wahrheit_skala() unterstuetzt nur z_type = \"binary\".",
         call. = FALSE)
  }
  inert <- intersect(names(args), "dropout_informative")
  inert <- inert[vapply(args[inert], function(v) !identical(v, 0),
                        logical(1))]
  if (length(inert)) {
    warning("Unter intern erzwungenem dropout = 0 wirkungslos: ",
            paste(inert, collapse = ", "), call. = FALSE)
  }
  welt <- function(zf) {
    do.call(sim_stream, c(args, list(dropout = 0, z_force = zf,
                                     seed = seed)))
  }
  s0 <- welt(0)
  rng0 <- get(".Random.seed", envir = globalenv())
  s1 <- welt(1)
  rng1 <- get(".Random.seed", envir = globalenv())
  if (!identical(s0$patient_id, s1$patient_id) ||
      !identical(s0$session, s1$session) ||
      !identical(s0$x, s1$x)) {
    stop("Kontrafaktischer Gleichlauf verletzt (CRN-Wache): ",
         "diese Weltkonfiguration zieht Zufallszahlen abhaengig ",
         "von der Zuweisung.", call. = FALSE)
  }
  if (!identical(rng0, rng1)) {
    stop("Kontrafaktischer Gleichlauf verletzt (RNG-Wache): ",
         "diese Weltkonfiguration verbraucht nach der Zuweisung ",
         "z-abhaengig Zufallszahlen.", call. = FALSE)
  }
  if (is.null(session)) session <- max(s0$session)
  zeilen <- s0$session == session
  if (!any(zeilen)) {
    stop("Sitzung ", session, " kommt im Strom nicht vor.",
         call. = FALSE)
  }
  lat0 <- s0$score[zeilen]
  lat1 <- s1$score[zeilen]
  b0 <- skala_werte(lat0, minimum, maximum, schiefe, mitte,
                    punkte_pro_einheit)
  b1 <- skala_werte(lat1, minimum, maximum, schiefe, mitte,
                    punkte_pro_einheit)
  x <- s0$x[zeilen]
  if (length(unique(x)) < 3) {
    stop("wahrheit_skala() braucht mindestens 3 verschiedene ",
         "x-Werte fuer die Terzil-Bildung.", call. = FALSE)
  }
  brueche <- stats::quantile(x, probs = c(0, 1/3, 2/3, 1))
  if (anyDuplicated(brueche)) {
    stop("x-Terzile nicht bildbar (Quantil-Brueche nicht eindeutig).",
         call. = FALSE)
  }
  terzil <- cut(x, brueche, include.lowest = TRUE,
                labels = c("x_unten", "x_mitte", "x_oben"))
  eff_lat <- tapply(lat1 - lat0, terzil, mean)
  eff_beob <- tapply(b1 - b0, terzil, mean)
  list(
    session = session,
    ate_latent = mean(lat1 - lat0),
    ate_beobachtet = mean(b1 - b0),
    effekte_x_terzile = data.frame(
      terzil = names(eff_lat),
      latent = as.numeric(eff_lat),
      beobachtet = as.numeric(eff_beob)),
    spann_latent = diff(range(eff_lat)),
    spann_beobachtet = diff(range(eff_beob)),
    schein_heterogenitaet = diff(range(eff_beob)) /
      punkte_pro_einheit - diff(range(eff_lat)),
    skala = list(minimum = minimum, maximum = maximum,
                 schiefe = schiefe, mitte = mitte,
                 punkte_pro_einheit = punkte_pro_einheit,
                 estimand_default = "beobachtet"))
}
