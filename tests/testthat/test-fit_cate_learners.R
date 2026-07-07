test_that("T-Learner: Struktur, Seed-Pflicht, Determinismus", {
  skip_if_not_installed("ranger")
  s <- sim_stream(n_therapists = 10, patients_per_therapist = 10,
                  n_sessions = 3, z_level = "patient", tau = 0.5,
                  seed = 40)
  expect_error(fit_cate_tlearner(s), "seed")
  a <- fit_cate_tlearner(s, num_trees = 100, B = 30, seed = 1)
  expect_named(a, c("patient_id", "x", "z", "tau_hat"))
  ate <- attr(a, "ate")
  expect_true(all(c("estimate", "se", "lo", "hi") %in% names(ate)))
  expect_true(all(is.finite(ate)))
  b <- fit_cate_tlearner(s, num_trees = 100, B = 30, seed = 1)
  expect_identical(a$tau_hat, b$tau_hat)
})

test_that("S-Boost: Struktur, Seed-Pflicht, Determinismus", {
  skip_if_not_installed("xgboost")
  s <- sim_stream(n_therapists = 10, patients_per_therapist = 10,
                  n_sessions = 3, z_level = "patient", tau = 0.5,
                  seed = 41)
  expect_error(fit_cate_sboost(s), "seed")
  a <- fit_cate_sboost(s, nrounds = 100, B = 30, seed = 1)
  expect_named(a, c("patient_id", "x", "z", "tau_hat"))
  expect_true(all(is.finite(attr(a, "ate"))))
  b <- fit_cate_sboost(s, nrounds = 100, B = 30, seed = 1)
  expect_identical(a$tau_hat, b$tau_hat)
})

test_that("beide finden grobe Richtung echter Heterogenitaet", {
  skip_if_not_installed("ranger"); skip_if_not_installed("xgboost")
  s <- sim_stream(n_therapists = 20, patients_per_therapist = 10,
                  n_sessions = 4, z_level = "patient", tau = 0.5,
                  tau_x = 0.8, seed = 42)
  wahr <- 0.5 + 0.8 * patients(s)$x
  # bewusst lockere Schwellen (Lektion: stochastische Groessen)
  t1 <- fit_cate_tlearner(s, num_trees = 300, B = 30, seed = 1)
  expect_gt(cor(t1$tau_hat, wahr), 0.2)
  s1 <- fit_cate_sboost(s, B = 30, seed = 1)
  expect_gt(cor(s1$tau_hat, wahr), 0.2)
})

test_that("X-/DR-/R-Learner: Struktur, Determinismus, grobe Richtung", {
  skip_if_not_installed("ranger")
  s <- sim_stream(n_therapists = 20, patients_per_therapist = 10,
                  n_sessions = 4, z_level = "patient", tau = 0.5,
                  tau_x = 0.8, seed = 42)
  wahr <- 0.5 + 0.8 * patients(s)$x
  for (fn in list(fit_cate_xlearner, fit_cate_drlearner,
                  fit_cate_rlearner)) {
    a <- fn(s, num_trees = 150, B = 20, seed = 1)
    expect_named(a, c("patient_id", "x", "z", "tau_hat"))
    expect_true(all(is.finite(a$tau_hat)))
    b <- fn(s, num_trees = 150, B = 20, seed = 1)
    expect_identical(a$tau_hat, b$tau_hat)
    # bewusst lockere Schwelle (Lektion: stochastische Groessen)
    expect_gt(cor(a$tau_hat, wahr), 0.2)
  }
})

test_that("PAI und MOB: Struktur, Determinismus; B=0 ohne Inferenz", {
  skip_if_not_installed("model4you")
  s <- sim_stream(n_therapists = 15, patients_per_therapist = 8,
                  n_sessions = 3, z_level = "patient", tau = 0.5,
                  tau_x = 0.8, seed = 43)
  wahr <- 0.5 + 0.8 * patients(s)$x
  a <- fit_cate_pai(s, B = 20, seed = 1)
  expect_gt(cor(a$tau_hat, wahr), 0.5)  # linear korrekt spezifiziert
  b <- fit_cate_pai(s, B = 20, seed = 1)
  expect_identical(a$tau_hat, b$tau_hat)
  m1 <- fit_cate_mob(s, num_trees = 50, B = 0, seed = 1)
  expect_true(all(is.finite(m1$tau_hat)))
  expect_true(is.na(attr(m1, "ate")[["se"]]))  # B=0: ehrlich NA
  expect_gt(cor(m1$tau_hat, wahr), 0.2)
  m2 <- fit_cate_mob(s, num_trees = 50, B = 0, seed = 1)
  expect_identical(m1$tau_hat, m2$tau_hat)
})
