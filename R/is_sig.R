#' Signifikanz-Entscheidung fuer fit_*-Ergebnisse
#'
#' Nutzt den p-Wert des Schaetzers, falls vorhanden (z. B.
#' Satterthwaite bei [fit_z_satt()]), sonst die Wald-z-Naeherung
#' `|estimate/se| > z_krit`.
#'
#' @param fit Benannter Vektor aus einer `fit_*`-Funktion.
#' @param alpha Signifikanzniveau.
#' @return `TRUE`/`FALSE`.
#' @export
is_sig <- function(fit, alpha = 0.05) {
  if ("p" %in% names(fit)) {
    unname(fit[["p"]] < alpha)
  } else {
    unname(abs(fit[["estimate"]] / fit[["se"]]) >
             stats::qnorm(1 - alpha / 2))
  }
}
