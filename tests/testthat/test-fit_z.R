skip_if_not_installed("lme4")

test_that("z ist konstant innerhalb einer Therapeut:in", {
  s <- sim_stream(n_therapists = 6, patients_per_therapist = 4,
                  n_sessions = 3, seed = 1)
  pro_therapeut <- tapply(s$z, s$therapist_id, function(x) length(unique(x)))
  expect_true(all(pro_therapeut == 1))
})

test_that("p_treated steuert den Anteil behandelter Therapeut:innen", {
  s <- sim_stream(n_therapists = 10, patients_per_therapist = 2,
                  n_sessions = 2, p_treated = 0.3, seed = 2)
  z_therapeut <- tapply(s$z, s$therapist_id, unique)
  expect_equal(sum(z_therapeut), 3)
})

test_that("DIE Kernaussage: naive SE ist bei Cluster-Treatment zu klein", {
  s <- sim_stream(n_therapists = 20, patients_per_therapist = 15,
                  n_sessions = 8, icc = 0.20, tau = 0, seed = 11)
  naiv <- fit_z_naive(s)
  lmm  <- fit_z_lmm(s)
  expect_lt(naiv[["se"]], lmm[["se"]] / 2)  # nicht nur kleiner: DRASTISCH kleiner
})

test_that("LMM-CI enthaelt tau bei eingeschaltetem Effekt", {
  s <- sim_stream(n_therapists = 40, patients_per_therapist = 15,
                  n_sessions = 8, icc = 0.10, tau = 0.5, seed = 9)
  fit <- fit_z_lmm(s)
  expect_true(fit[["estimate"]] - 1.96 * fit[["se"]] < 0.5)
  expect_true(0.5 < fit[["estimate"]] + 1.96 * fit[["se"]])
})

test_that("fit_z_satt liefert estimate, se und p", {
  skip_if_not_installed("lmerTest")
  s <- sim_stream(n_therapists = 10, patients_per_therapist = 6,
                  n_sessions = 4, seed = 3)
  fit <- fit_z_satt(s)
  expect_named(fit, c("estimate", "se", "p"))
  expect_true(fit[["p"]] >= 0 && fit[["p"]] <= 1)
})

test_that("Satterthwaite-p ist konservativer als Wald bei wenigen Clustern", {
  skip_if_not_installed("lmerTest")
  s <- sim_stream(n_therapists = 8, patients_per_therapist = 6,
                  n_sessions = 4, icc = 0.2, seed = 6)
  fit <- fit_z_satt(s)
  wald_p <- 2 * pnorm(-abs(fit[["estimate"]] / fit[["se"]]))
  expect_gt(fit[["p"]], wald_p)
})

test_that("z_level='therapist' (Default) laesst alte Welten bit-identisch", {
  a <- sim_stream(n_therapists = 3, patients_per_therapist = 4,
                  n_sessions = 5, seed = 1)
  b <- sim_stream(n_therapists = 3, patients_per_therapist = 4,
                  n_sessions = 5, z_level = "therapist", seed = 1)
  expect_identical(a, b)
})

test_that("z_level='patient': z variiert innerhalb der Therapeut:innen", {
  s <- sim_stream(n_therapists = 10, patients_per_therapist = 10,
                  n_sessions = 3, z_level = "patient", seed = 2)
  per_ther <- tapply(s$z, s$therapist_id,
                     function(v) length(unique(v)))
  expect_gt(mean(per_ther == 2), 0.7)  # fast ueberall beide Gruppen
  # und z ist pro Patient:in konstant
  per_pat <- tapply(s$z, s$patient_id, function(v) length(unique(v)))
  expect_true(all(per_pat == 1))
})
