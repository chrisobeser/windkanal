test_that("run_gates prueft an allen erreichbaren Terminen", {
  s <- sim_stream(n_therapists = 8, patients_per_therapist = 6,
                  n_sessions = 5, seed = 4)
  g <- run_gates(s, fit_z_naive, at_weeks = c(25, 40, max(s$week)))
  expect_equal(nrow(g$gates), 3)
  expect_equal(g$gates$alpha_local, rep(0.05 / 3, 3))
})

test_that("ohne Korrektur ist das lokale Niveau das volle alpha", {
  s <- sim_stream(n_therapists = 8, patients_per_therapist = 6,
                  n_sessions = 5, seed = 4)
  g <- run_gates(s, fit_z_naive, at_weeks = c(25, 40),
                 correction = "none")
  expect_equal(unique(g$gates$alpha_local), 0.05)
})

test_that("Bonferroni-Treffer impliziert None-Treffer (nie umgekehrt)", {
  for (sd in 1:5) {
    s <- sim_stream(n_therapists = 10, patients_per_therapist = 5,
                    n_sessions = 4, icc = 0.2, seed = sd)
    bon <- run_gates(s, fit_z_naive, at_weeks = c(25, 40, 55))
    non <- run_gates(s, fit_z_naive, at_weeks = c(25, 40, 55),
                     correction = "none")
    if (bon$hit) expect_true(non$hit)
  }
})
