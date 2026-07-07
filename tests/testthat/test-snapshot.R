test_that("snapshot enthaelt nur Messungen bis at_week", {
  s <- sim_stream(n_therapists = 3, patients_per_therapist = 4,
                  n_sessions = 6, seed = 1)
  snap <- snapshot(s, at_week = 20)
  expect_true(all(snap$week <= 20))
  expect_lt(nrow(snap), nrow(s))
})

test_that("snapshot kennt sein Aufnahmedatum (Provenienz)", {
  s <- sim_stream(n_therapists = 2, patients_per_therapist = 2,
                  n_sessions = 3, seed = 1)
  expect_equal(attr(snapshot(s, 15), "snapshot_week"), 15)
})

test_that("Snapshot am Horizont == ganzer Strom", {
  s <- sim_stream(n_therapists = 2, patients_per_therapist = 3,
                  n_sessions = 4, seed = 2)
  snap <- snapshot(s, at_week = max(s$week))
  expect_equal(nrow(snap), nrow(s))
})

test_that("Snapshot in Woche 0 ist leer", {
  s <- sim_stream(n_therapists = 2, patients_per_therapist = 2,
                  n_sessions = 3, seed = 5)
  snap <- snapshot(s, at_week = 0)
  expect_equal(nrow(snap), 0)
  expect_equal(attr(snap, "snapshot_week"), 0)
})
