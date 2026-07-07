test_that("Default: linear-Form bit-identisch, keine Noise-Spalten", {
  s <- sim_stream(n_therapists = 5, patients_per_therapist = 4,
                  n_sessions = 3, icc = 0.10, z_level = "patient",
                  tau = 0.4, tau_x = 0.3, seed = 42)
  # Referenzwerte der Vorversion (dritte Verankerung desselben Ankers)
  expect_equal(sum(s$score), -2.8151341017, tolerance = 1e-9)
  expect_length(grep("^x_noise", names(s)), 0)
})

test_that("step: Schwellen-Effekt ohne Wirkung unterhalb", {
  s <- sim_stream(n_therapists = 40, patients_per_therapist = 20,
                  n_sessions = 4, z_level = "patient",
                  tau = 0.2, tau_x = 0.8, tau_x_form = "step", seed = 10)
  p <- patients(s)
  eff <- function(d) unname(coef(lm(score_mean ~ z, data = d))["z"])
  hoch  <- eff(p[p$x > 0.2, ])   # klar oberhalb der Schwelle
  tief  <- eff(p[p$x < -0.2, ])  # klar unterhalb
  expect_gt(hoch, 0.6)           # ~ tau + tau_x = 1.0
  expect_lt(tief, 0.6)           # ~ tau = 0.2
  expect_gt(hoch - tief, 0.4)
})

test_that("quadratic: ATE bleibt tau, Kruemmung vorhanden", {
  s <- sim_stream(n_therapists = 40, patients_per_therapist = 20,
                  n_sessions = 4, z_level = "patient",
                  tau = 0.5, tau_x = 0.6, tau_x_form = "quadratic",
                  seed = 11)
  p <- patients(s)
  b <- coef(lm(score_mean ~ z * I(x^2 - 1), data = p))
  expect_lt(abs(unname(b["z"]) - 0.5), 0.2)
  expect_gt(unname(b["z:I(x^2 - 1)"]), 0.4)
  expect_lt(unname(b["z:I(x^2 - 1)"]), 0.8)
})

test_that("n_noise: Spalten da, aber wirkungslos", {
  s <- sim_stream(n_therapists = 20, patients_per_therapist = 20,
                  n_sessions = 3, z_level = "patient", tau = 0.5,
                  n_noise = 3, seed = 12)
  expect_length(grep("^x_noise", names(s)), 3)
  p <- patients(s)
  expect_true(all(c("x_noise1", "x_noise2", "x_noise3") %in% names(p)))
  b <- coef(lm(score_mean ~ x_noise1 + x_noise2 + x_noise3, data = p))
  expect_lt(max(abs(b[-1])), 0.2)
})

test_that("cate_features sammelt x + therapist_c + noise", {
  s <- sim_stream(n_therapists = 10, patients_per_therapist = 10,
                  n_sessions = 3, z_level = "patient", tau = 0.5,
                  tau_xc = 1, n_noise = 2, seed = 13)
  m <- cate_features(patients(s))
  expect_identical(colnames(m),
                   c("x", "therapist_c", "x_noise1", "x_noise2"))
  expect_true(is.numeric(m))
})

test_that("grf laeuft mit Noise-Features (Struktur)", {
  skip_if_not_installed("grf")
  s <- sim_stream(n_therapists = 10, patients_per_therapist = 10,
                  n_sessions = 3, z_level = "patient", tau = 0.5,
                  n_noise = 2, seed = 14)
  ct <- fit_cate_grf(s, num_trees = 200, seed = 1)
  expect_true(all(is.finite(ct$tau_hat)))
})
