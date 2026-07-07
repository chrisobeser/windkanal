#' Item-Antworten aus dem latenten Zustand simulieren
#'
#' Erzeugt fuer jede beobachtete Sitzung ordinale Antworten auf die
#' Items einer Batterie (Struktur z. B. aus [read_items()]). Modell
#' pro Item: standardisierter latenter Score * Ladung + Eigenrauschen,
#' dann in `n_options` gleichwahrscheinliche Kategorien geschnitten;
#' Umkehr-Items (`reversed`) werden gespiegelt.
#'
#' Damit wird Reliabilitaet **emergent**: Cronbachs Alpha der
#' simulierten Skala ist ein Messergebnis (siehe [scale_scores()]),
#' keine Annahme mehr -- und Fragen wie "was kostet Kuerzung von 10
#' auf 4 Items?" werden berechenbar.
#'
#' @param stream Datenstrom aus [sim_stream()].
#' @param items `data.frame` mit `item`, `n_options`, `reversed`
#'   (z. B. aus [read_items()]); optionale Spalten `loading`
#'   (Trennschaerfe je Item, sonst `loading`-Argument) und `latent`
#'   (`"score"` = Symptom-Prozess, Default; `"alliance"` = der
#'   Beziehungs-Prozess -- braucht `alliance = TRUE` im Strom).
#'   So laesst sich die Modul-Batterie abbilden: Symptom-Items und
#'   eingebettete Allianz-Items aus einem Guss.
#' @param loading Standard-Ladung aller Items (Default 0.7 --
#'   typische standardisierte Ladung; 10 Items ergeben damit
#'   Alpha ~ .9).
#' @param careless Ermuedungs-Autopilot: Wahrscheinlichkeit pro Item,
#'   dass eine Person ab dort auf Durchklicken schaltet (Straight-
#'   lining: alle restlichen Antworten = letzte Antwort). Default 0.
#'   Kalibrierung: Meade & Craig (2012, Psych Methods 17, 437-455)
#'   finden ~8-12% Careless-Anteil in Low-Stakes-Stichproben;
#'   Umrechnung auf pro-Item-Onset: p = 1-(1-Anteil)^(1/(k-1)),
#'   z. B. 10% bei 10 Items -> careless ~ 0.012. Klinisch motivierte
#'   Stichproben liegen niedriger (BA-Anker: 4.2%). Achtung,
#'   kontraintuitiv: Straightlining **erhoeht** das scheinbare
#'   Alpha -- kaputte Daten sehen reliabler aus.
#' @param seconds_per_item Zeitkosten pro Item (Default 7 s --
#'   im belegten Band: einfache Likert-Items ~6-10 s/Item laut
#'   Web-Survey-Benchmarks; lange Batterien beschleunigen auf
#'   ~10 s/Frage, "Speeding") -- ergibt das Attribut
#'   `battery_seconds`.
#' @param refusal Verweigerung: Wahrscheinlichkeit pro Sitzung, die
#'   **gesamte Batterie** auszulassen (alle Items `NA`). Default 0.
#'   Laengen-Abhaengigkeit ueber Szenarien kalibrieren: Web-Benchmarks
#'   zeigen ~70-80% Completion bei 5-Minuten- und ~60-70% bei
#'   10-Minuten-Umfragen (graue Quellen).
#' @param seed Zufalls-Seed (Pflicht).
#' @return `data.frame`: `patient_id`, `session`, dann eine Spalte
#'   pro Item (Werte 1..n_options; `NA` bei Verweigerung). Attribut
#'   `items` (die Struktur).
#' @export
sim_items <- function(stream, items, loading = 0.7, careless = 0,
                      seconds_per_item = 7, refusal = 0, seed) {
  if (missing(seed)) {
    stop("`seed` ist Pflicht.", call. = FALSE)
  }
  stopifnot(is.data.frame(items),
            all(c("item", "n_options", "reversed") %in% names(items)),
            !anyNA(items$n_options), all(items$n_options >= 2))
  lam <- if ("loading" %in% names(items)) items$loading else
    rep(loading, nrow(items))
  stopifnot(all(lam > 0 & lam < 1), careless >= 0, careless < 1,
            refusal >= 0, refusal < 1)
  src <- if ("latent" %in% names(items)) items$latent else
    rep("score", nrow(items))
  stopifnot(all(src %in% c("score", "alliance")))
  if (any(src == "alliance") && !"alliance" %in% names(stream)) {
    stop("Items mit latent='alliance' brauchen einen Strom mit ",
         "alliance = TRUE.", call. = FALSE)
  }
  set.seed(seed)

  z <- as.vector(scale(stream$score))
  za <- if (any(src == "alliance"))
    as.vector(scale(stream$alliance)) else NULL
  n <- length(z)
  out <- data.frame(patient_id = stream$patient_id,
                    session = stream$session)
  for (j in seq_len(nrow(items))) {
    quelle <- if (src[j] == "alliance") za else z
    eta <- lam[j] * quelle + sqrt(1 - lam[j]^2) * rnorm(n)
    ko <- items$n_options[j]
    cuts <- stats::qnorm(seq_len(ko - 1) / ko)
    resp <- findInterval(eta, cuts) + 1L
    if (isTRUE(items$reversed[j])) resp <- ko + 1L - resp
    out[[items$item[j]]] <- resp
  }
  # Ermuedungs-Autopilot: ab Onset-Position Straightlining
  k <- nrow(items)
  if (careless > 0 && k > 1) {
    onset <- stats::rgeom(n, careless) + 2L  # fruehestens ab Item 2
    resp_cols <- items$item
    for (r in which(onset <= k)) {
      o <- onset[r]
      lock <- out[[resp_cols[o - 1L]]][r]
      for (j in o:k) {
        # in den Wertebereich des Items geklemmt
        out[[resp_cols[j]]][r] <-
          min(max(lock, 1L), items$n_options[j])
      }
    }
  }
  # Verweigerung: ganze Sitzungs-Batterie ausgelassen (alle NA)
  if (refusal > 0) {
    skip <- stats::runif(n) < refusal
    if (any(skip)) {
      for (j in items$item) out[[j]][skip] <- NA_integer_
    }
  }
  attr(out, "items") <- items
  attr(out, "battery_seconds") <- k * seconds_per_item
  out
}

#' Skalenwerte + emergente Reliabilitaet aus Item-Antworten
#'
#' Polt Umkehr-Items zurueck, summiert je Skala und berechnet
#' Cronbachs Alpha als **Messergebnis** der simulierten Batterie.
#'
#' @param item_data Ergebnis von [sim_items()].
#' @param items Struktur (Default: `items`-Attribut von `item_data`).
#' @return `data.frame` mit `patient_id`, `session` und einer
#'   Summen-Spalte je Skala; Attribut `alpha` (benannter Vektor:
#'   Cronbachs Alpha je Skala, NA bei < 2 Items).
#' @export
scale_scores <- function(item_data, items = attr(item_data, "items")) {
  stopifnot(!is.null(items), "scale" %in% names(items))
  out <- item_data[, c("patient_id", "session")]
  alphas <- c()
  for (sc in unique(items$scale)) {
    it <- items[items$scale == sc, ]
    m <- as.matrix(item_data[, it$item, drop = FALSE])
    # Umkehr-Items zurueckpolen (auf gemeinsame Richtung)
    for (j in seq_len(nrow(it))) {
      if (isTRUE(it$reversed[j])) {
        m[, j] <- it$n_options[j] + 1L - m[, j]
      }
    }
    out[[sc]] <- rowSums(m)  # NA bei Verweigerung (ehrlich fehlend)
    alphas[sc] <- if (ncol(m) < 2) NA_real_ else {
      S <- stats::cov(m, use = "pairwise.complete.obs")
      k <- ncol(m)
      unname(k / (k - 1) * (1 - sum(diag(S)) / sum(S)))
    }
  }
  attr(out, "alpha") <- alphas
  out
}
