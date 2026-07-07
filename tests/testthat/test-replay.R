test_that("replay liefert eine Zeile pro Snapshot-Woche", {
  s <- sim_stream(n_therapists = 5, patients_per_therapist = 5,
                  n_sessions = 6, seed = 1)
  r <- replay(s, at_weeks = c(20, 40, 60))
  expect_equal(nrow(r), 3)
  expect_equal(r$snapshot_week, c(20, 40, 60))
})

test_that("zu kleine Snapshots ergeben NA statt Fehler", {
  s <- sim_stream(n_therapists = 2, patients_per_therapist = 2,
                  n_sessions = 3, seed = 1)
  r <- replay(s, at_weeks = c(0, max(s$week)))
  expect_true(is.na(r$estimate[1]))
  expect_false(is.na(r$estimate[2]))
})

test_that("SE schrumpft mit wachsendem Datenstrom", {
  s <- sim_stream(n_therapists = 10, patients_per_therapist = 10,
                  n_sessions = 8, seed = 42)
  r <- replay(s, at_weeks = c(15, 60))
  expect_lt(r$se[2], r$se[1])
})

test_that("replay findet den wahren Slope ungefaehr wieder", {
  s <- sim_stream(n_therapists = 20, patients_per_therapist = 20,
                  n_sessions = 10, mean_slope = -0.15, seed = 7)
  r <- replay(s, at_weeks = max(s$week))
  expect_true(r$ci_low < -0.15 && -0.15 < r$ci_high)
})
