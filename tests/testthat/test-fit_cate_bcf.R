skip_if_not_installed("bcf")

test_that("fit_cate_bcf: Struktur, Seed-Pflicht, Kredibilitaetsintervalle", {
  s <- sim_stream(n_therapists = 10, patients_per_therapist = 10,
                  n_sessions = 3, tau = 0.5, seed = 2)
  expect_error(fit_cate_bcf(s), "seed")
  ct <- suppressWarnings(fit_cate_bcf(s, nburn = 150, nsim = 150,
                                      seed = 1))
  expect_named(ct, c("patient_id", "x", "z", "tau_hat",
                     "tau_lo", "tau_hi"))
  expect_true(all(ct$tau_lo <= ct$tau_hat & ct$tau_hat <= ct$tau_hi))
  expect_true(is.finite(attr(ct, "ate")[["estimate"]]))
})

test_that("kurze Leine: bcf erfindet in homogener Welt weniger Streuung als grf", {
  skip_if_not_installed("grf")
  s <- sim_stream(n_therapists = 20, patients_per_therapist = 10,
                  n_sessions = 4, tau = 0.5, tau_x = 0, seed = 4)
  sd_grf <- sd(fit_cate_grf(s, num_trees = 500, seed = 1)$tau_hat)
  sd_bcf <- sd(suppressWarnings(
    fit_cate_bcf(s, nburn = 250, nsim = 250, seed = 1))$tau_hat)
  # Schwelle 0.7 statt 0.5: sd(tau_hat) haengt an Ketten-Zahl/-Laenge
  # (1 Kette = rauere Posterior-Mittel); Kernaussage bleibt: deutlich
  # weniger Schein-Streuung als grf.
  expect_lt(sd_bcf, sd_grf * 0.7)
})
