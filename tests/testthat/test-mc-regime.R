test_that("mc_run: eine Zeile pro Welt x Estimator, reproduzierbar", {
  args <- list(n_therapists = 6, patients_per_therapist = 4,
               n_sessions = 3, tau = 0)
  a <- mc_run(5, args, list(naiv = fit_z_naive), seed_start = 1)
  b <- mc_run(5, args, list(naiv = fit_z_naive), seed_start = 1)
  expect_equal(nrow(a), 5)
  expect_identical(a, b)
})

test_that("mc_summary liefert reject_rate/bias/coverage", {
  args <- list(n_therapists = 6, patients_per_therapist = 4,
               n_sessions = 3, tau = 0)
  res <- mc_run(5, args, list(naiv = fit_z_naive))
  sm <- mc_summary(res, truth = 0)
  expect_named(sm, c("estimator", "n_reps", "reject_rate",
                     "bias", "coverage"))
})

test_that("run_peek liefert Blicke und Stop-Logik", {
  s <- sim_stream(n_therapists = 8, patients_per_therapist = 6,
                  n_sessions = 5, seed = 4)
  p <- run_peek(s, fit_z_naive, every = 8)
  expect_true(is.data.frame(p$looks) && nrow(p$looks) > 1)
  expect_type(p$hit, "logical")
  if (p$hit) expect_true(p$first_hit_week %in% p$looks$week)
})

test_that("run_gate testet genau einmal am Termin", {
  s <- sim_stream(n_therapists = 8, patients_per_therapist = 6,
                  n_sessions = 5, seed = 4)
  g <- run_gate(s, fit_z_naive, at_week = max(s$week))
  expect_true(is.finite(g$estimate))
  expect_type(g$significant, "logical")
})

test_that("zu fruehes Gate wird nicht erzwungen (NA statt Fehler)", {
  s <- sim_stream(n_therapists = 3, patients_per_therapist = 2,
                  n_sessions = 2, seed = 4)
  g <- run_gate(s, fit_z_naive, at_week = 0)
  expect_true(is.na(g$estimate) && !g$significant)
})
