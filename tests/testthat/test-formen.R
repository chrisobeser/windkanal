# Kongruenz-Formen-Familie (tau_xc_form): Bit-Identitaet, Normierung,
# Integration.

test_that("Produkt-Default bleibt bit-identisch (explizit vs. Default)", {
  a <- sim_stream(n_therapists = 5, patients_per_therapist = 4,
                  n_sessions = 2, tau = 0.5, tau_c = 0.2, tau_xc = 0.35,
                  seed = 42)
  b <- sim_stream(n_therapists = 5, patients_per_therapist = 4,
                  n_sessions = 2, tau = 0.5, tau_c = 0.2, tau_xc = 0.35,
                  tau_xc_form = "product", seed = 42)
  expect_identical(a, b)
})

test_that("tau_xc = 0: alle Formen liefern identische Streams", {
  streams <- lapply(c("product", "misfit_quad", "similarity_abs", "fenster"),
    function(f) sim_stream(n_therapists = 5, patients_per_therapist = 4,
                           n_sessions = 2, tau = 0.5, tau_c = 0.2,
                           tau_xc = 0, tau_xc_form = f, seed = 7))
  for (i in 2:4) expect_equal(streams[[1]], streams[[i]])
})

test_that("dyade_form: zentriert und auf SD 1 normiert (alle Formen)", {
  set.seed(11); n <- 5e5
  x <- rnorm(n); cc <- rnorm(n)
  for (f in c("product", "misfit_quad", "similarity_abs", "fenster")) {
    v <- dyade_form(x, cc, f)
    expect_lt(abs(mean(v)), 0.01)
    expect_lt(abs(sd(v) - 1), 0.01)
  }
})

test_that("fenster_delta wirkt und wird validiert", {
  set.seed(3); x <- rnorm(100); cc <- rnorm(100)
  expect_false(identical(dyade_form(x, cc, "fenster", 0.5),
                         dyade_form(x, cc, "fenster", 1.0)))
  expect_error(sim_stream(n_therapists = 3, patients_per_therapist = 2,
                          n_sessions = 2, tau_xc = 0.3,
                          tau_xc_form = "fenster", fenster_delta = -1,
                          seed = 1))
  # Obergrenze: jenseits delta = 3
  # degeneriert die Form (p > .96, Normierung teilt durch nahe null)
  expect_error(dyade_form(x, cc, "fenster", 4))
  expect_error(sim_stream(n_therapists = 3, patients_per_therapist = 2,
                          n_sessions = 2, tau_xc = 0.3,
                          tau_xc_form = "fenster", fenster_delta = 4,
                          seed = 1))
})

test_that("fenster_delta als Vektor wird abgelehnt statt recycelt", {
  set.seed(3); x <- rnorm(8); cc <- rnorm(8)
  # ein Vektor besteht die elementweise Schranke und recycelt danach
  # lautlos -- eine Welt mit zwei wechselnden Fensterbreiten
  expect_error(dyade_form(x, cc, "fenster", c(0.1, 2)),
               "einzelne endliche Zahl")
  expect_error(sim_stream(n_therapists = 3, patients_per_therapist = 2,
                          n_sessions = 2, tau_xc = 0.3,
                          tau_xc_form = "fenster",
                          fenster_delta = c(0.5, 2.5), seed = 1),
               "einzelne endliche Zahl")
  expect_error(dyade_form(x, cc, "fenster", "0.5"), "einzelne endliche Zahl")
  expect_error(dyade_form(x, cc, "fenster", NA_real_), "einzelne endliche Zahl")
  expect_error(dyade_form(x, cc, "fenster", Inf), "einzelne endliche Zahl")
  expect_error(dyade_h_wahr(x, cc, "fenster", fenster_delta = c(0.5, 2.5)),
               "einzelne endliche Zahl")
  # Negativfall: der Skalar laeuft unveraendert durch
  expect_length(dyade_form(x, cc, "fenster", 0.5), 8L)
})

test_that("Formen aendern den Stream nur ueber den dyadischen Term", {
  a <- sim_stream(n_therapists = 6, patients_per_therapist = 5,
                  n_sessions = 2, tau = 0.5, tau_xc = 0.35,
                  tau_xc_form = "product", seed = 9)
  b <- sim_stream(n_therapists = 6, patients_per_therapist = 5,
                  n_sessions = 2, tau = 0.5, tau_xc = 0.35,
                  tau_xc_form = "fenster", seed = 9)
  expect_false(identical(a$score, b$score))    # Form wirkt
  expect_identical(a$x, b$x)                   # Draws unveraendert
  expect_identical(a$therapist_c, b$therapist_c)
  expect_identical(a$z, b$z)
})
