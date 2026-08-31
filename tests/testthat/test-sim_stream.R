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

test_that("Seed-Wache haelt auch bei explizitem NULL", {
  a <- list(n_therapists = 3, patients_per_therapist = 2, n_sessions = 2)
  expect_error(do.call(sim_stream, c(a, list(seed = NULL))),
               "NOT a missing argument")
  expect_error(do.call(sim_stream, c(a, list(seed = "abc"))),
               "single finite number")
  expect_error(do.call(sim_stream, c(a, list(seed = c(1, 2)))),
               "vector of length 2")
  expect_error(do.call(sim_stream, c(a, list(seed = NA))),
               "type logical")
  expect_error(do.call(sim_stream, c(a, list(seed = NA_real_))),
               "NA, NaN or Inf")
  expect_error(do.call(sim_stream, c(a, list(seed = Inf))),
               "NA, NaN or Inf")
  expect_error(do.call(sim_stream, a), "the argument is missing")
  # Negativfall: ein gueltiger Seed laeuft und liefert denselben Strom
  expect_identical(do.call(sim_stream, c(a, list(seed = 7))),
                   do.call(sim_stream, c(a, list(seed = 7L))))
})

test_that("gebrochene Anzahlen stoppen statt NA-Zeilen zu liefern (S2)", {
  a <- list(n_therapists = 3, patients_per_therapist = 2,
            n_sessions = 2, seed = 5)
  expect_error(do.call(sim_stream, modifyList(a, list(n_therapists = 3.7))),
               "n_therapists. must be a single whole number")
  expect_error(do.call(sim_stream,
                       modifyList(a, list(patients_per_therapist = 2.5))),
               "patients_per_therapist. must be a single whole number")
  expect_error(do.call(sim_stream, modifyList(a, list(n_sessions = 3.5))),
               "n_sessions. must be a single whole number")
  expect_error(do.call(sim_stream, modifyList(a, list(weeks_accrual = 2.5))),
               "weeks_accrual. must be a single whole number")
  expect_error(do.call(sim_stream, modifyList(a, list(n_therapists = c(2, 3)))),
               "n_therapists. must be a single whole number")
  # dieselbe Wache deckt soft_matching: die Zusage "caseloads exact"
  # galt vorher nur so weit, wie die Eingaben ganzzahlig waren
  expect_error(do.call(sim_stream, modifyList(a, list(
                 patients_per_therapist = 2.5, tau_c = 0.3,
                 assignment = "soft_matching", assignment_strength = 2))),
               "patients_per_therapist. must be a single whole number")
  # Negativfall: ganzzahlig laeuft und traegt keine NA
  s <- do.call(sim_stream, a)
  expect_false(anyNA(s$therapist_id))
  expect_false(anyNA(s$score))
  expect_identical(nrow(s), 3L * 2L * 2L)
})

test_that("assignment_strength laeuft nicht ueber", {
  a <- list(n_therapists = 6, patients_per_therapist = 5, n_sessions = 3,
            assignment = "soft_matching", tau_c = 0.3, seed = 11)
  for (lam in c(200, 1e4)) {
    s <- do.call(sim_stream, modifyList(a, list(assignment_strength = lam)))
    expect_false(anyNA(s$therapist_id))
    expect_false(anyNA(s$score))
  }
  # im Grenzfall grosser lambda ist die Paarung die reine Maximum-Wahl:
  # Patient 1 trifft auf volle Kapazitaeten und muss zu argmax(x1 * c_j)
  s <- do.call(sim_stream, modifyList(a, list(assignment_strength = 1e4)))
  p <- patients(s)
  c_j <- tapply(p$therapist_c, p$therapist_id, function(v) v[1])
  erster <- p[p$patient_id == 1L, ][1, ]
  expect_identical(as.integer(erster$therapist_id),
                   as.integer(names(c_j)[which.max(erster$x * c_j)]))
  # Negativfall: unterhalb der Schranke bleibt die Arithmetik unberuehrt.
  # Anker auf der Zuordnung selbst -- sum(score) taugt hier NICHT, weil
  # es bei exakten Caseloads und z auf Therapeutenebene gegen jedes
  # Umsortieren invariant ist.
  klein <- do.call(sim_stream, modifyList(a, list(assignment_strength = 1)))
  expect_identical(
    as.integer(patients(klein)$therapist_id),
    c(2L, 5L, 1L, 5L, 4L, 6L, 4L, 3L, 6L, 1L, 6L, 5L, 6L, 6L, 2L,
      5L, 2L, 4L, 3L, 2L, 3L, 3L, 3L, 2L, 1L, 5L, 1L, 4L, 1L, 4L))
  expect_identical(unname(table(patients(klein)$therapist_id)),
                   unname(table(rep(1:6, each = 5))))
  # und lambda wirkt ueberhaupt: gegen die Zufallszuweisung verschieden
  null <- do.call(sim_stream, modifyList(a, list(assignment_strength = 0,
                                                 assignment = "random")))
  expect_false(identical(klein$therapist_id, null$therapist_id))
})
