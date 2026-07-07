test_that("Bit-Identitaet: Defaults (x_effect=0, confounding=0) unveraendert", {
  # Referenzwerte der Version VOR dem Einbau (2026-07-04 gesichert)
  s <- sim_stream(n_therapists = 5, patients_per_therapist = 4,
                  n_sessions = 3, icc = 0.10, z_level = "patient",
                  tau = 0.4, tau_x = 0.3, seed = 42)
  expect_equal(sum(s$score), -2.8151341017, tolerance = 1e-9)
  expect_equal(s$score[1], 2.4850529231, tolerance = 1e-9)
  expect_equal(sum(patients(s)$z), 11)
  expect_equal(patients(s)$x[7], 0.7681787378, tolerance = 1e-9)
})

test_that("confounding koppelt Zuweisung an x", {
  s <- sim_stream(n_therapists = 20, patients_per_therapist = 20,
                  n_sessions = 2, z_level = "patient",
                  confounding = 2, seed = 1)
  p <- patients(s)
  expect_gt(cor(p$z, p$x), 0.3)
  expect_gt(mean(p$x[p$z == 1]), mean(p$x[p$z == 0]))
  # ohne confounding: keine Kopplung
  s0 <- sim_stream(n_therapists = 20, patients_per_therapist = 20,
                   n_sessions = 2, z_level = "patient", seed = 1)
  expect_lt(abs(cor(patients(s0)$z, patients(s0)$x)), 0.15)
})

test_that("x_effect ist ein Haupteffekt auf den Score", {
  s <- sim_stream(n_therapists = 20, patients_per_therapist = 20,
                  n_sessions = 4, z_level = "patient", p_treated = 0,
                  x_effect = 1, seed = 2)
  p <- patients(s)
  b <- unname(coef(lm(score_mean ~ x, data = p))["x"])
  expect_gt(b, 0.8)
  expect_lt(b, 1.2)
})

test_that("Konfundierung by indication verzerrt den naiven Vergleich", {
  # x_effect + confounding: Behandelte haben hoeheres x UND x wirkt
  # auf den Score -> naiver z-Vergleich muss nach oben verzerrt sein
  s <- sim_stream(n_therapists = 20, patients_per_therapist = 20,
                  n_sessions = 4, z_level = "patient", tau = 0,
                  x_effect = 1, confounding = 1.5, seed = 3)
  p <- patients(s)
  naiv <- unname(coef(lm(score_mean ~ z, data = p))["z"])
  expect_gt(naiv, 0.3)  # wahrer Effekt ist 0
})

test_that("Fehlerbedingungen", {
  expect_error(sim_stream(confounding = 1, z_level = "therapist",
                          seed = 1), "patient")
  expect_error(sim_stream(confounding = 1, z_level = "patient",
                          p_treated = 1, seed = 1), "p_treated")
})

test_that("fit_cate_bcf pihat='glm' liefert gueltige Struktur", {
  skip_if_not_installed("bcf")
  s <- sim_stream(n_therapists = 10, patients_per_therapist = 10,
                  n_sessions = 3, z_level = "patient", tau = 0.5,
                  x_effect = 1, confounding = 1, seed = 4)
  ct <- suppressWarnings(fit_cate_bcf(s, nburn = 100, nsim = 100,
                                      pihat = "glm", seed = 1))
  expect_named(ct, c("patient_id", "x", "z", "tau_hat",
                     "tau_lo", "tau_hi"))
  expect_true(is.finite(attr(ct, "ate")[["estimate"]]))
})
