# Formations-Schalter (assignment = "soft_matching"): Bit-Identitaet,
# Kapazitaeten, Selektions-Signatur, Guards.

test_that("Default bleibt bit-identisch (explizit vs. Default)", {
  a <- sim_stream(n_therapists = 5, patients_per_therapist = 4,
                  n_sessions = 2, tau = 0.5, tau_xc = 0.3, seed = 42)
  b <- sim_stream(n_therapists = 5, patients_per_therapist = 4,
                  n_sessions = 2, tau = 0.5, tau_xc = 0.3,
                  assignment = "random", seed = 42)
  expect_identical(a, b)
})

test_that("soft_matching: Caseloads exakt, Selektion sichtbar, reproduzierbar", {
  s <- sim_stream(n_therapists = 20, patients_per_therapist = 10,
                  n_sessions = 1, tau = 0.5, tau_xc = 0,
                  assignment = "soft_matching",
                  assignment_strength = 2, seed = 11)
  p <- s[!duplicated(s$patient_id), ]
  expect_true(all(table(p$therapist_id) == 10))       # Kapazitaet exakt
  expect_gt(cor(p$x, p$therapist_c), 0.15)            # Kollider da
  expect_true(!is.null(s$therapist_c))                # c immer gezogen
  s2 <- sim_stream(n_therapists = 20, patients_per_therapist = 10,
                   n_sessions = 1, tau = 0.5, tau_xc = 0,
                   assignment = "soft_matching",
                   assignment_strength = 2, seed = 11)
  expect_identical(s, s2)                              # Seed traegt
})

test_that("soft_matching: Guards (icc_slope, confounding, strength)", {
  expect_error(sim_stream(n_therapists = 4, patients_per_therapist = 3,
                          n_sessions = 1, icc_slope = 0.1,
                          assignment = "soft_matching", seed = 1),
               "icc_slope")
  expect_error(sim_stream(n_therapists = 4, patients_per_therapist = 3,
                          n_sessions = 1, confounding = 0.5,
                          z_level = "patient",
                          assignment = "soft_matching", seed = 1),
               "confounding")
  expect_error(sim_stream(n_therapists = 4, patients_per_therapist = 3,
                          n_sessions = 1, assignment = "soft_matching",
                          assignment_strength = -1, seed = 1))
})

test_that("assignment_strength unter \"random\" stoppt statt still zu wirken", {
  a <- list(n_therapists = 5, patients_per_therapist = 4, n_sessions = 2,
            tau = 0.5, tau_xc = 0.3, seed = 42)
  expect_error(do.call(sim_stream, c(a, list(assignment_strength = 1.2))),
               "assignment_strength wirkt nur")
  expect_error(do.call(sim_stream, c(a, list(assignment = "random",
                                             assignment_strength = 1.2))),
               "soft_matching")
  # Typ-Wache vor der Zweigung: NA, Vektor, NULL, Inf, Text
  for (bad in list(NA_real_, c(1, 2), NULL, Inf, "1")) {
    expect_error(do.call(sim_stream, c(a, list(assignment_strength = bad))),
                 "einzelne endliche Zahl")
  }
  # sie schuetzt auch den Zweig, in dem lambda tatsaechlich wirkt
  expect_error(sim_stream(n_therapists = 5, patients_per_therapist = 4,
                          n_sessions = 2, assignment = "soft_matching",
                          assignment_strength = c(1, 2), seed = 42),
               "einzelne endliche Zahl")
  # Negativfall 1: lambda = 0 ist erlaubt und bit-identisch zum Default
  expect_identical(do.call(sim_stream, c(a, list(assignment_strength = 0))),
                   do.call(sim_stream, a))
  # Negativfall 2: unter soft_matching greift die Wache bei keinem lambda
  b <- list(n_therapists = 6, patients_per_therapist = 5, n_sessions = 2,
            assignment = "soft_matching", tau_xc = 0.3, seed = 42)
  expect_silent(do.call(sim_stream, c(b, list(assignment_strength = 0))))
  expect_silent(do.call(sim_stream, c(b, list(assignment_strength = 1.2))))
})
