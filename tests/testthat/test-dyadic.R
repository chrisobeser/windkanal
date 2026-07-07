test_that("Default: keine therapist_c-Spalte, keine Extra-Ziehungen", {
  s <- sim_stream(n_therapists = 5, patients_per_therapist = 4,
                  n_sessions = 3, icc = 0.10, z_level = "patient",
                  tau = 0.4, tau_x = 0.3, seed = 42)
  expect_null(s$therapist_c)
  # Referenzwerte der Vorversion (wie test-confounding.R)
  expect_equal(sum(s$score), -2.8151341017, tolerance = 1e-9)
})

test_that("tau_xc erzeugt einen echten dyadischen Effekt", {
  s <- sim_stream(n_therapists = 40, patients_per_therapist = 20,
                  n_sessions = 4, z_level = "patient",
                  tau = 0.5, tau_xc = 1, seed = 5)
  p <- patients(s)
  expect_false(is.null(p$therapist_c))
  # Regression rekonstruiert den Interaktions-Term
  b <- coef(lm(score_mean ~ z * x * therapist_c, data = p))
  expect_gt(unname(b["z:x:therapist_c"]), 0.7)
  expect_lt(unname(b["z:x:therapist_c"]), 1.3)
  # ATE bleibt tau (x, c unabhaengig standardnormal)
  expect_lt(abs(unname(b["z"]) - 0.5), 0.2)
})

test_that("tau_c moderiert auf Therapeuten-Ebene", {
  s <- sim_stream(n_therapists = 40, patients_per_therapist = 20,
                  n_sessions = 4, z_level = "patient",
                  tau = 0.5, tau_c = 1, seed = 6)
  p <- patients(s)
  b <- coef(lm(score_mean ~ z * therapist_c, data = p))
  expect_gt(unname(b["z:therapist_c"]), 0.7)
  expect_lt(unname(b["z:therapist_c"]), 1.3)
})

test_that("tau_shape='ramp': Effekt waechst ueber Sitzungen an", {
  base <- list(n_therapists = 40, patients_per_therapist = 20,
               n_sessions = 4, z_level = "patient", tau = 1,
               tau_shape = "ramp", seed = 7)
  s <- do.call(sim_stream, base)
  eff_at <- function(k) {
    p <- patients(s, outcome = "at_session", session = k)
    unname(coef(lm(score_at ~ z, data = p))["z"])
  }
  # Sitzung 1: tau * 1/4; Sitzung 4: tau * 4/4
  expect_lt(eff_at(1), 0.6)
  expect_gt(eff_at(4), 0.7)
  expect_gt(eff_at(4) - eff_at(1), 0.4)
})

test_that("CATE-Schaetzer sehen therapist_c (Struktur)", {
  skip_if_not_installed("grf")
  s <- sim_stream(n_therapists = 10, patients_per_therapist = 10,
                  n_sessions = 3, z_level = "patient",
                  tau = 0.5, tau_xc = 1, seed = 8)
  ct <- fit_cate_grf(s, num_trees = 200, seed = 1)
  expect_true(all(is.finite(ct$tau_hat)))
})

test_that("bcf-ATE-Attribut traegt Posterior-Quantile lo/hi", {
  skip_if_not_installed("bcf")
  s <- sim_stream(n_therapists = 8, patients_per_therapist = 8,
                  n_sessions = 3, z_level = "patient", tau = 0.5,
                  seed = 9)
  ct <- suppressWarnings(fit_cate_bcf(s, nburn = 100, nsim = 100,
                                      seed = 1))
  a <- attr(ct, "ate")
  expect_true(all(c("estimate", "se", "lo", "hi") %in% names(a)))
  expect_lt(a[["lo"]], a[["estimate"]])
  expect_gt(a[["hi"]], a[["estimate"]])
})
