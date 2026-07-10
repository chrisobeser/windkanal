#' Plasmode-Welt: echte Kovariaten und Therapeuten-Struktur,
#' injizierte bekannte Wahrheit
#'
#' Baut einen sim_stream-kompatiblen Datenstrom aus einer ECHTEN
#' Kohorten-Tabelle (z. B. TONI/PsyTOM): echte Patient:innen-Merkmale
#' (x), echte Therapeuten-Zuordnung samt realer Caseload-Verteilung --
#' aber Treatment, Effekt und Outcome werden nach der windkanal-
#' Gleichung mit BEKANNTER Wahrheit erzeugt. So prueft man Schaetzer
#' unter echter Kovariaten-Verteilung (Schiefe!) und echter
#' Cluster-Struktur, ohne die Wahrheit zu verlieren (Kovariaten-
#' Plasmode; Rezept-Prinzip: die Kohorte selbst bleibt lokal und
#' wird nie mitverteilt).
#'
#' @param kohorte data.frame mit mindestens `patient_id`,
#'   `therapist_id` und der in `x_spalte` benannten Kovariate
#'   (z. B. TONI: `t0_ace_summe`). Zeilen mit NA in x fliegen mit
#'   Warnung heraus (dokumentierte Politik E-1: complete case je
#'   Analyse).
#' @param x_spalte Name der Kovariate, die als Moderator x dient;
#'   wird z-standardisiert.
#' @param n_sessions,icc,tau,tau_x,tau_x_form,seed wie in
#'   [sim_stream()]; die Wahrheit ist tau + tau_x * h(x).
#' @param z_level "patient" (Dyaden-Achse, Default) oder "therapist".
#' @return data.frame mit den sim_stream-Spalten (`patient_id`,
#'   `therapist_id`, `session`, `z`, `x`, `score`), Attribut
#'   `wahrheit` (Vektor der wahren Individual-Effekte) und
#'   `plasmode_info`.
#' @export
plasmode_world <- function(kohorte, x_spalte, n_sessions = 4,
                           icc = 0.10, tau = 0.5, tau_x = 0,
                           tau_x_form = c("linear", "step", "quadratic"),
                           z_level = c("patient", "therapist"), seed) {
  tau_x_form <- match.arg(tau_x_form)
  z_level <- match.arg(z_level)
  if (missing(seed)) stop("`seed` is mandatory.", call. = FALSE)
  stopifnot(all(c("patient_id", "therapist_id", x_spalte) %in%
                names(kohorte)))
  k <- kohorte[!is.na(kohorte[[x_spalte]]), ]
  k <- k[order(k$patient_id), ]  # Ordnung = patients()-Ordnung (sort by id)
  if (nrow(k) < nrow(kohorte)) {
    warning(sprintf("%d Zeilen ohne %s entfernt (E-1 complete case).",
                    nrow(kohorte) - nrow(k), x_spalte))
  }
  set.seed(seed)
  x <- as.numeric(scale(k[[x_spalte]]))
  n <- nrow(k)
  ther <- as.integer(factor(k$therapist_id))
  J <- max(ther)
  h <- switch(tau_x_form, linear = x, step = as.numeric(x > 0),
              quadratic = x^2 - 1)
  # Treatment: randomisiert auf gewaehlter Ebene (bekannt, fair)
  z <- if (z_level == "therapist") {
    zt <- stats::rbinom(J, 1, 0.5); zt[ther]
  } else {
    stats::rbinom(n, 1, 0.5)
  }
  # Outcome nach der windkanal-Gleichung auf echter Struktur:
  sd_u <- sqrt(icc / (1 - icc))          # Residual-SD = 1
  u <- stats::rnorm(J, 0, sd_u)[ther]
  b <- stats::rnorm(n, 0, 1)
  tau_i <- tau + tau_x * h
  out <- do.call(rbind, lapply(seq_len(n_sessions), function(s_idx) {
    data.frame(patient_id = k$patient_id, therapist_id = k$therapist_id,
               session = s_idx, z = z, x = x,
               score = b + u + tau_i * z + stats::rnorm(n, 0, 1))
  }))
  attr(out, "wahrheit") <- tau_i
  attr(out, "plasmode_info") <- list(
    n_patienten = n, n_therapeuten = J,
    caseload_median = stats::median(table(ther)),
    x_quelle = x_spalte, x_schiefe = round(mean((x - mean(x))^3), 3))
  out
}
