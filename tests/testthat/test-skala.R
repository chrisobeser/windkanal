# Bounded/Skewed Outcomes (R/skala.R): Begrenzung, Schiefe, Wahrheit
# auf beiden Ebenen (beobachtete Skala als Default-Massstab, latente
# Wahrheit stets mit ausgewiesen)

test_that("score_obs liegt in den Grenzen und ist ganzzahlig", {
  d <- sim_stream(n_therapists = 6, patients_per_therapist = 5,
                  n_sessions = 3, tau = 1, seed = 42)
  b <- skala_begrenzen(d, minimum = 0, maximum = 27)
  expect_true(all(b$score_obs >= 0 & b$score_obs <= 27))
  expect_type(b$score_obs, "integer")
})

test_that("Schicht ist Opt-in: Bestand unveraendert, nur Spalte + Attribut neu", {
  d <- sim_stream(n_therapists = 4, patients_per_therapist = 3,
                  n_sessions = 2, tau = 0.5, seed = 7)
  b <- skala_begrenzen(d)
  expect_identical(b[names(d)], d[names(d)])
  expect_identical(class(b), class(d))
  # jedes vorbestehende Attribut ueberlebt unveraendert
  alt <- attributes(d); alt$names <- NULL
  expect_identical(attributes(b)[names(alt)], alt)
  expect_identical(setdiff(names(b), names(d)), "score_obs")
  expect_identical(attr(b, "skala")$estimand_default, "beobachtet")
})

test_that("Abbildung ist deterministisch und monoton", {
  d <- sim_stream(n_therapists = 5, patients_per_therapist = 4,
                  n_sessions = 2, seed = 3)
  b1 <- skala_begrenzen(d, schiefe = 0.4, punkte_pro_einheit = 9)
  b2 <- skala_begrenzen(d, schiefe = 0.4, punkte_pro_einheit = 9)
  expect_identical(b1$score_obs, b2$score_obs)
  o <- order(d$score)
  expect_true(all(diff(b1$score_obs[o]) >= 0L))
})

test_that("weites Fenster ohne Schiefe = Rundung der latenten Werte", {
  d <- sim_stream(n_therapists = 4, patients_per_therapist = 4,
                  n_sessions = 2, seed = 11)
  b <- skala_begrenzen(d, minimum = -1000, maximum = 1000,
                       punkte_pro_einheit = 1)
  expect_identical(b$score_obs, as.integer(round(d$score)))
})

test_that("schiefe > 0 erzeugt Rechtsschiefe und warnt vor der Tilt-Schranke", {
  d <- sim_stream(n_therapists = 20, patients_per_therapist = 10,
                  n_sessions = 2, seed = 5)
  schief <- function(v) mean((v - mean(v))^3) / stats::sd(v)^3
  b0 <- skala_begrenzen(d, minimum = 0, maximum = 100,
                        punkte_pro_einheit = 5, mitte = 0)
  # ppe = 5 < schiefe * 49.5: untere Grenze strukturell unerreichbar
  expect_warning(
    b1 <- skala_begrenzen(d, minimum = 0, maximum = 100,
                          punkte_pro_einheit = 5, mitte = 0,
                          schiefe = 0.8),
    "unerreichbar")
  expect_gt(schief(as.numeric(b1$score_obs)),
            schief(as.numeric(b0$score_obs)))
  # ppe oberhalb der Schranke: keine Warnung
  expect_silent(skala_begrenzen(d, minimum = 0, maximum = 27,
                                punkte_pro_einheit = 9,
                                schiefe = 0.3))
})

test_that("Deckeneffekt: hohe Verankerung sammelt Masse am Maximum", {
  d <- sim_stream(n_therapists = 10, patients_per_therapist = 10,
                  n_sessions = 2, seed = 13)
  b <- skala_begrenzen(d, minimum = 0, maximum = 27, mitte = -3)
  expect_gt(mean(b$score_obs == 27), 0.2)
  expect_lt(mean(b$score_obs == 0), 0.01)
})

test_that("Fehlerfaelle stoppen sauber", {
  d <- sim_stream(n_therapists = 3, patients_per_therapist = 2,
                  n_sessions = 2, seed = 1)
  expect_error(skala_begrenzen(d, minimum = 5, maximum = 5))
  expect_error(skala_begrenzen(d, minimum = 0, maximum = 10.5))
  expect_error(skala_begrenzen(d, punkte_pro_einheit = 0))
  expect_error(skala_begrenzen(data.frame(a = 1)))
  expect_error(wahrheit_skala(n_therapists = 3, seed = 1,
                              dropout = 0.1),
               "setzt wahrheit_skala")
  expect_error(wahrheit_skala(n_therapists = 3, seed = 1,
                              z_force = 1),
               "setzt wahrheit_skala")
  expect_error(wahrheit_skala(n_therapists = 3, z_type = "dose",
                              seed = 1),
               "binary")
  expect_error(wahrheit_skala(n_therapists = 3), "`seed` is mandatory")
  expect_error(wahrheit_skala(n_therapists = 1,
                              patients_per_therapist = 2,
                              n_sessions = 2, seed = 2),
               "Terzil")
  expect_warning(wahrheit_skala(n_therapists = 6,
                                patients_per_therapist = 4,
                                n_sessions = 2,
                                dropout_informative = 0.5, seed = 4),
                 "wirkungslos")
})

test_that("z_force-Doppellauf ist zeilenweise exakt", {
  a <- list(n_therapists = 8, patients_per_therapist = 6,
            n_sessions = 4, tau = 1.2, tau_x = 0.5, tau_c = 0.3,
            tau_xc = 0.4, x_effect = 0.7, dropout = 0)
  s0 <- do.call(sim_stream, c(a, list(z_force = 0, seed = 21)))
  s1 <- do.call(sim_stream, c(a, list(z_force = 1, seed = 21)))
  te <- 1.2 + 0.5 * s0$x + (0.3 + 0.4 * s0$x) * s0$therapist_c
  expect_equal(s1$score - s0$score, te, tolerance = 1e-12)
})

test_that("z_force ist Opt-in und verlangt binaeres z", {
  a <- list(n_therapists = 5, patients_per_therapist = 4,
            n_sessions = 2, tau = 0.5)
  expect_identical(do.call(sim_stream, c(a, list(seed = 3))),
                   do.call(sim_stream, c(a, list(z_force = NULL,
                                                 seed = 3))))
  expect_error(do.call(sim_stream, c(a, list(z_type = "dose",
                                             z_force = 1, seed = 3))),
               "binary")
  expect_error(do.call(sim_stream, c(a, list(z_force = 2, seed = 3))))
  # der erzwungene Arm steht wirklich ueberall
  s1 <- do.call(sim_stream, c(a, list(z_force = 1, seed = 3)))
  expect_true(all(s1$z == 1L))
})

test_that("Allianz- und Dosis-Welten durch die Schicht", {
  da <- sim_stream(n_therapists = 4, patients_per_therapist = 3,
                   n_sessions = 3, alliance = TRUE, seed = 9)
  ba <- skala_begrenzen(da)
  expect_identical(ba$alliance, da$alliance)
  dd <- sim_stream(n_therapists = 4, patients_per_therapist = 3,
                   n_sessions = 2, z_type = "dose", seed = 9)
  expect_true(all(skala_begrenzen(dd)$score_obs >= 0L))
  expect_identical(skala_werte(c(1, NA, 2)), c(18L, NA, 22L))
  expect_warning(skala_begrenzen(ba), "ueberschrieben")
})

test_that("CRN-Doppellauf: latenter Effekt exakt tau + tau_x * x", {
  w <- wahrheit_skala(n_therapists = 8, patients_per_therapist = 6,
                      n_sessions = 3, tau = 1.2, tau_x = 0.5,
                      minimum = -1000, maximum = 1000,
                      punkte_pro_einheit = 1, seed = 21)
  expect_equal(w$ate_latent, 1.2, tolerance = 0.15)
  expect_gt(w$spann_latent, 0)
})

test_that("ohne aktive Grenze: beobachteter ATE = Punkte * latenter ATE", {
  w <- wahrheit_skala(n_therapists = 30, patients_per_therapist = 10,
                      n_sessions = 2, tau = 1,
                      minimum = -1000, maximum = 1000,
                      punkte_pro_einheit = 3, seed = 8)
  expect_equal(w$ate_beobachtet, 3 * w$ate_latent, tolerance = 0.05)
})

test_that("aktive Decke daempft den beobachteten Effekt", {
  w <- wahrheit_skala(n_therapists = 30, patients_per_therapist = 10,
                      n_sessions = 2, tau = 1,
                      minimum = 0, maximum = 27, mitte = -2,
                      seed = 8)
  linear_erwartet <- w$skala$punkte_pro_einheit * w$ate_latent
  expect_lt(abs(w$ate_beobachtet), abs(linear_erwartet))
})

test_that("Schein-h-Zelle: latent homogen, beobachtet heterogen unter Decke", {
  # wahres tau_x = 0 -- auf der latenten Skala hilft die Behandlung
  # allen gleich; die aktive Decke erzeugt beobachtete Heterogenitaet
  w_decke <- wahrheit_skala(n_therapists = 40,
                            patients_per_therapist = 12,
                            n_sessions = 2, tau = 1.5, tau_x = 0,
                            x_effect = 1.5,
                            minimum = 0, maximum = 27, mitte = -2,
                            seed = 26)
  w_frei <- wahrheit_skala(n_therapists = 40,
                           patients_per_therapist = 12,
                           n_sessions = 2, tau = 1.5, tau_x = 0,
                           x_effect = 1.5,
                           minimum = -1000, maximum = 1000,
                           punkte_pro_einheit = w_decke$skala$punkte_pro_einheit,
                           seed = 26)
  expect_lt(w_decke$spann_latent, 0.05)
  expect_lt(w_frei$spann_beobachtet, 0.6)
  expect_gt(w_decke$spann_beobachtet, 2 * w_frei$spann_beobachtet)
})

test_that("Schein-Heterogenitaet ist einheitenbereinigt", {
  w <- wahrheit_skala(n_therapists = 20, patients_per_therapist = 30,
                      n_sessions = 6, z_level = "patient",
                      tau = 0.5, tau_x = 0.6, seed = 4242)
  # die naive Lesart subtrahiert Skalenpunkte von latenten Einheiten
  expect_equal(w$spann_beobachtet - w$spann_latent, 4.33, tolerance = 0.01)
  # die einheitenbereinigte Kennzahl dreht Vorzeichen UND Groessenordnung
  expect_equal(w$schein_heterogenitaet, -0.09, tolerance = 0.02)
  expect_lt(w$schein_heterogenitaet, 0)
  expect_equal(w$schein_heterogenitaet,
               w$spann_beobachtet / w$skala$punkte_pro_einheit -
                 w$spann_latent)
  # Negativfall: gleiche Einheiten, unerreichbare Grenzen -> nur Rundung
  w2 <- wahrheit_skala(n_therapists = 20, patients_per_therapist = 30,
                       n_sessions = 6, z_level = "patient",
                       tau = 0.5, tau_x = 0.6,
                       minimum = -1000, maximum = 1000,
                       punkte_pro_einheit = 1, seed = 4242)
  expect_lt(abs(w2$schein_heterogenitaet), 0.05)
})

test_that("seed = NULL meldet den Seed, nicht die CRN-Wache", {
  msg <- tryCatch(wahrheit_skala(n_therapists = 6,
                                 patients_per_therapist = 5,
                                 n_sessions = 3, seed = NULL),
                  error = conditionMessage)
  expect_match(msg, "`seed` is mandatory")
  expect_false(grepl("CRN", msg))
  # Negativfall: mit Seed rechnet dieselbe Welt durch
  expect_type(wahrheit_skala(n_therapists = 6, patients_per_therapist = 5,
                             n_sessions = 3, seed = 1), "list")
})
