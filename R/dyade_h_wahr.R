#' Wahre h-Komponente einer Kongruenz-Form (Doppelzentrierung)
#'
#' Die Wahrheits-Referenz fuer dyadische Experimente: Die FORM selbst
#' (etwa `-(x-c)^2`) enthaelt Haupteffekt-Anteile -- bei der
#' quadratischen Kongruenz zur Haelfte. Matching-Kennzahlen muessen
#' daher gegen die h-KOMPONENTE gerechnet werden, nie gegen die rohe
#' Form. Diese Funktion berechnet sie numerisch:
#' `tau_xc * dyade_form(...)` am Produktgitter der uebergebenen
#' Raender, doppelzentriert.
#'
#' @param x Patientenwerte (Vektor; die Gitter-Zeilen).
#' @param c Therapeutenwerte (Vektor, ein Wert je Therapeut:in;
#'   die Gitter-Spalten).
#' @param form,fenster_delta Wie in [dyade_form()].
#' @param tau_xc Skalierung (Default 1).
#' @param therapist_of Optional: Index des Therapeuten je Patient
#'   (Laenge wie `x`, Werte in `1:length(c)`); dann enthaelt das
#'   Ergebnis zusaetzlich `realisiert`, das wahre h der realisierten
#'   Paare -- die Referenz fuer personenbezogene Metriken.
#' @return Liste: `h` (Matrix length(x) x length(c)), `f`, `g`,
#'   `tau0`, optional `realisiert`.
#' @examples
#' set.seed(2)
#' x <- rnorm(60); c <- rnorm(8)
#' w <- dyade_h_wahr(x, c, "misfit_quad", tau_xc = 0.35,
#'                   therapist_of = rep(1:8, length.out = 60))
#' max(abs(rowMeans(w$h)))   # doppelzentriert: ~0
#' length(w$realisiert)
#' @export
dyade_h_wahr <- function(x, c, form = c("product", "misfit_quad",
                                        "similarity_abs", "fenster"),
                         fenster_delta = 0.5, tau_xc = 1,
                         therapist_of = NULL) {
  form <- match.arg(form)
  M <- tau_xc * outer(x, c, function(a, b)
    dyade_form(a, b, form, fenster_delta))
  m0 <- mean(M); fx <- rowMeans(M) - m0; gc_ <- colMeans(M) - m0
  H <- M - outer(fx, gc_, "+") - m0
  out <- list(h = H, f = fx, g = gc_, tau0 = m0)
  if (!is.null(therapist_of)) {
    stopifnot(length(therapist_of) == length(x),
              all(therapist_of %in% seq_along(c)))
    out$realisiert <- H[cbind(seq_along(x), therapist_of)]
  }
  out
}
