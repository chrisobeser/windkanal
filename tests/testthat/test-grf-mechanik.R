skip_if_not_installed("grf")

test_that("Defaults unveraendert: neue Schalter reproduzieren alten Aufruf", {
  s <- sim_stream(n_therapists = 10, patients_per_therapist = 10,
                  n_sessions = 3, icc = 0.20, z_level = "patient",
                  tau = 0.5, seed = 20)
  a <- fit_cate_grf(s, num_trees = 300, seed = 1)
  b <- fit_cate_grf(s, num_trees = 300, honesty = TRUE,
                    clusters = TRUE, seed = 1)
  expect_identical(a$tau_hat, b$tau_hat)
  expect_identical(attr(a, "ate"), attr(b, "ate"))
})

test_that("honesty=FALSE und clusters=FALSE laufen und wirken", {
  s <- sim_stream(n_therapists = 10, patients_per_therapist = 10,
                  n_sessions = 3, icc = 0.20, z_level = "therapist",
                  tau = 0.5, seed = 21)
  an  <- fit_cate_grf(s, num_trees = 300, seed = 1)
  ohne_h <- fit_cate_grf(s, num_trees = 300, honesty = FALSE, seed = 1)
  ohne_c <- fit_cate_grf(s, num_trees = 300, clusters = FALSE, seed = 1)
  expect_false(identical(an$tau_hat, ohne_h$tau_hat))
  # ohne Cluster-Korrektur wirkt die Stichprobe groesser -> kleinere SE
  expect_lt(attr(ohne_c, "ate")[["se"]], attr(an, "ate")[["se"]])
})
