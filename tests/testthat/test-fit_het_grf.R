skip_if_not_installed("grf")

test_that("fit_het_grf hat das fit_*-Interface mit p-Wert", {
  s <- sim_stream(n_therapists = 15, patients_per_therapist = 10,
                  n_sessions = 3, seed = 2)
  fit <- fit_het_grf(s, num_trees = 300)
  expect_named(fit, c("estimate", "se", "p"))
  expect_true(fit[["p"]] >= 0 && fit[["p"]] <= 1)
})

test_that("echte Heterogenitaet wird erkannt (gut gepowert)", {
  s <- sim_stream(n_therapists = 30, patients_per_therapist = 15,
                  n_sessions = 4, tau = 0.3, tau_x = 0.5, seed = 3)
  fit <- fit_het_grf(s, num_trees = 500)
  expect_lt(fit[["p"]], 0.05)
})
