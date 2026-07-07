skip_if_not_installed("lme4")

test_that("fit_lmm hat das Standard-Interface (estimate, se)", {
  s <- sim_stream(n_therapists = 8, patients_per_therapist = 8,
                  n_sessions = 6, seed = 3)
  fit <- fit_lmm(s)
  expect_named(fit, c("estimate", "se"))
  expect_true(is.finite(fit[["estimate"]]) && fit[["se"]] > 0)
})

test_that("fit_lmm laeuft im replay()", {
  s <- sim_stream(n_therapists = 8, patients_per_therapist = 8,
                  n_sessions = 6, seed = 3)
  r <- replay(s, at_weeks = c(30, max(s$week)), fit_fn = fit_lmm)
  expect_equal(nrow(r), 2)
  expect_false(any(is.na(r$estimate)))
})

# Lehrstueck (der Windkanal widerlegte hier die erste Intuition):
# Fuer den WITHIN-Person-Slope ist das LMM EFFIZIENTER (kleinere SE),
# weil es Personen-/Therapeuten-Varianz aus dem Residuum zieht.
# Die klassische Anticonservativitaet naiver SEs betrifft Praediktoren
# auf CLUSTER-Ebene (z. B. Treatment je Therapeut:in) -- Test dafuer
# folgt, sobald sim_stream() ein Treatment Z kann.
test_that("within-person Slope: LMM ist effizienter als naives lm", {
  s <- sim_stream(n_therapists = 15, patients_per_therapist = 15,
                  n_sessions = 8, icc = 0.30, seed = 11)
  naiv <- fit_slope(s)
  lmm  <- fit_lmm(s)
  expect_lt(lmm[["se"]], naiv[["se"]])
})

test_that("LMM-CI enthaelt die Wahrheit bei gut gepowerter Simulation", {
  s <- sim_stream(n_therapists = 20, patients_per_therapist = 15,
                  n_sessions = 8, mean_slope = -0.15, seed = 9)
  fit <- fit_lmm(s)
  expect_true(fit[["estimate"]] - 1.96 * fit[["se"]] < -0.15)
  expect_true(-0.15 < fit[["estimate"]] + 1.96 * fit[["se"]])
})
