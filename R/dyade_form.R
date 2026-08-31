#' Kongruenz-Formen der dyadischen Passung (zentriert, normiert)
#'
#' Die Formen-Familie des dyadischen Terms: Wie wirkt die Paarung
#' (x, c) -- als Produkt, als quadratische Kongruenz, als
#' Betrags-Aehnlichkeit oder als Passungs-Zone (Fenster)? Jede Form
#' wird **zentriert** (Erwartungswert 0 unter unabhaengigen
#' Standardnormal-Raendern, damit `tau` nicht stillschweigend
#' verschoben wird) und **normiert** (SD 1 unter demselben
#' Referenzmass, damit `tau_xc` ueber Formen vergleichbar bleibt).
#'
#' Analytische Konstanten (X, C ~ N(0,1) unabhaengig; D = X - C ~ N(0, 2)).
#' E und SD sind exakt; h-Anteile und Produktnaehe haengen vom
#' Referenzmass ab (Populationswert vs. Quantilgitter) und sind darum
#' als Spannen angegeben:
#' * `misfit_quad` -(x-c)^2: E = -2, SD = sqrt(8); der h-Anteil ist
#'   EXAKT 2xc (Produktform) -- 50 % der Varianz sind Haupteffekte.
#' * `similarity_abs` -|x-c|: E = -2/sqrt(pi), SD = sqrt(2 - 4/pi);
#'   h-Anteil ~.55-.60, Produktnaehe Kor(h, xc) ~ .89-.91.
#' * `fenster` `1{|x-c| <= delta}`: E = p = 2*Phi(delta/sqrt(2)) - 1,
#'   SD = sqrt(p(1-p)); h-Anteil ~.90, Produktnaehe ~ .31-.35 --
#'   der Stresstest der Familie.
#'
#' @param x Patientenmerkmal (numerisch, Standardnormal-Skala).
#' @param c Therapeutenmerkmal derselben Zeilen (numerisch).
#' @param form `"product"` (Default), `"misfit_quad"`,
#'   `"similarity_abs"` oder `"fenster"`.
#' @param fenster_delta Fensterbreite delta (nur `form = "fenster"`;
#'   Default 0.5). Eine **einzelne** endliche Zahl -- ein Vektor wuerde
#'   still gegen `x` und `c` recyceln. Zulaessig ist (0, 3]: jenseits
#'   von delta = 3 faengt
#'   das Fenster fast jedes Paar (p > .96) und die SD-Normierung
#'   teilt durch nahe null -- die Form degeneriert zur Konstanten.
#' @return Numerischer Vektor: zentrierter, normierter Formwert; die
#'   Beitrags-Formel im Generator ist `tau_xc * dyade_form(...)`.
#' @examples
#' set.seed(1)
#' x <- rnorm(1e5); c <- rnorm(1e5)
#' # alle Formen sind zentriert und auf SD 1 normiert (n gross genug,
#' # damit der SD-Schaetzer des Produkts das auch zeigt)
#' sapply(c("product", "misfit_quad", "similarity_abs", "fenster"),
#'        function(f) sd(dyade_form(x, c, f)))
#' @export
dyade_form <- function(x, c,
                       form = c("product", "misfit_quad",
                                "similarity_abs", "fenster"),
                       fenster_delta = 0.5) {
  form <- match.arg(form)
  switch(form,
    product = x * c,                       # E = 0, SD = 1: unveraendert
    misfit_quad = (-(x - c)^2 + 2) / sqrt(8),
    similarity_abs = (-abs(x - c) + 2 / sqrt(pi)) / sqrt(2 - 4 / pi),
    fenster = {
      fenster_delta_pruefen(fenster_delta)
      p <- 2 * stats::pnorm(fenster_delta / sqrt(2)) - 1
      (as.numeric(abs(x - c) <= fenster_delta) - p) / sqrt(p * (1 - p))
    })
}

# internal: the window width is a single number, not a vector. A
# vector passes `stopifnot(fenster_delta > 0, fenster_delta <= 3)`
# element-wise and then recycles silently against x and c -- the run
# would carry two alternating window widths, a world nobody declared.
fenster_delta_pruefen <- function(fenster_delta) {
  if (!is.numeric(fenster_delta) || length(fenster_delta) != 1L ||
      !is.finite(fenster_delta)) {
    stop("fenster_delta muss eine einzelne endliche Zahl sein (ein ",
         "Vektor wuerde in 1{|x - c| <= fenster_delta} still gegen x ",
         "und c recyceln und eine Welt mit wechselnden Fensterbreiten ",
         "erzeugen, die so nicht gemeint war).", call. = FALSE)
  }
  if (fenster_delta <= 0 || fenster_delta > 3) {
    stop("fenster_delta muss in (0, 3] liegen: jenseits von delta = 3 ",
         "faengt das Fenster fast jedes Paar (p > .96) und die ",
         "SD-Normierung teilt durch nahe null.", call. = FALSE)
  }
  invisible(TRUE)
}
