test_that("sim_stream liefert eine Zeile pro Sitzung", {
  s <- sim_stream(n_therapists = 3, patients_per_therapist = 4,
                  n_sessions = 5, seed = 1)
  expect_s3_class(s, "data.frame")
  expect_equal(nrow(s), 3 * 4 * 5)
  expect_named(s, c("therapist_id", "z", "patient_id", "x", "entry_week",
                    "session", "week", "score"))
})

test_that("gleicher Seed ergibt identischen Strom", {
  a <- sim_stream(n_therapists = 2, patients_per_therapist = 3,
                  n_sessions = 4, seed = 42)
  b <- sim_stream(n_therapists = 2, patients_per_therapist = 3,
                  n_sessions = 4, seed = 42)
  expect_identical(a, b)
})

test_that("verschiedene Seeds ergeben verschiedene Stroeme", {
  a <- sim_stream(n_therapists = 2, patients_per_therapist = 3,
                  n_sessions = 4, seed = 1)
  b <- sim_stream(n_therapists = 2, patients_per_therapist = 3,
                  n_sessions = 4, seed = 2)
  expect_false(identical(a$score, b$score))
})

test_that("seed ist Pflicht", {
  expect_error(sim_stream(n_therapists = 2), "seed")
})

test_that("unzulaessige ICC wird abgelehnt", {
  expect_error(sim_stream(icc = 1, seed = 1))
  expect_error(sim_stream(icc = -0.1, seed = 1))
})

test_that("week = entry_week + session - 1", {
  s <- sim_stream(n_therapists = 2, patients_per_therapist = 2,
                  n_sessions = 3, seed = 7)
  expect_equal(s$week, s$entry_week + s$session - 1L)
})

test_that("hoehere ICC erzeugt mehr Varianz zwischen Therapeuten", {
  lo <- sim_stream(n_therapists = 50, patients_per_therapist = 20,
                   n_sessions = 1, icc = 0.01, seed = 3)
  hi <- sim_stream(n_therapists = 50, patients_per_therapist = 20,
                   n_sessions = 1, icc = 0.40, seed = 3)
  var_between <- function(d) var(tapply(d$score, d$therapist_id, mean))
  expect_gt(var_between(hi), var_between(lo))
})
