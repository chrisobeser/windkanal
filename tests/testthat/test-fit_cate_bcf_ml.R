skip_if_not_installed("stochtree")

test_that("fit_cate_bcf_ml: Struktur, Seed-Pflicht, Determinismus", {
  s <- sim_stream(n_therapists = 10, patients_per_therapist = 10,
                  n_sessions = 3, icc = 0.20, z_level = "therapist",
                  tau = 0.5, seed = 30)
  expect_error(fit_cate_bcf_ml(s), "seed")
  a <- fit_cate_bcf_ml(s, nburn = 100, nsim = 100, seed = 1)
  expect_named(a, c("patient_id", "x", "z", "tau_hat",
                    "tau_lo", "tau_hi"))
  expect_true(all(a$tau_lo <= a$tau_hat & a$tau_hat <= a$tau_hi))
  ate <- attr(a, "ate")
  expect_true(all(c("estimate", "se", "lo", "hi") %in% names(ate)))
  # seed-deterministisch (anders als altes bcf: exakt, kein Rauschen)
  b <- fit_cate_bcf_ml(s, nburn = 100, nsim = 100, seed = 1)
  expect_identical(a$tau_hat, b$tau_hat)
})

test_that("rfx=FALSE laeuft (Standard-Doppelgaenger) und wirkt", {
  s <- sim_stream(n_therapists = 10, patients_per_therapist = 10,
                  n_sessions = 3, icc = 0.20, z_level = "therapist",
                  tau = 0.5, seed = 31)
  a <- fit_cate_bcf_ml(s, nburn = 100, nsim = 100, seed = 1)
  b <- fit_cate_bcf_ml(s, nburn = 100, nsim = 100, rfx = FALSE, seed = 1)
  expect_true(all(is.finite(b$tau_hat)))
  expect_false(identical(a$tau_hat, b$tau_hat))
})
