# Dosis-Modus (z_type = "dose"): stetiges Treatment auf [0, 1]

test_that("Binaer-Default bleibt bit-identisch (explizit vs. Default)", {
  a <- sim_stream(n_therapists = 4, patients_per_therapist = 3,
                  n_sessions = 2, tau = 0.5, seed = 42)
  b <- sim_stream(n_therapists = 4, patients_per_therapist = 3,
                  n_sessions = 2, tau = 0.5, z_type = "binary", seed = 42)
  expect_identical(a, b)
})

test_that("Dosis-Welt: z stetig in [0,1], Ebenen-Semantik korrekt", {
  d <- sim_stream(n_therapists = 5, patients_per_therapist = 4,
                  n_sessions = 2, z_type = "dose",
                  z_level = "therapist", tau = 0.5, seed = 7)
  expect_true(all(d$z >= 0 & d$z <= 1))
  expect_false(all(d$z %in% c(0, 1)))
  je_th <- tapply(d$z, d$therapist_id, function(v) length(unique(v)))
  expect_true(all(je_th == 1))  # Therapeuten-Ebene: konstant je Cluster

  d2 <- sim_stream(n_therapists = 5, patients_per_therapist = 4,
                   n_sessions = 2, z_type = "dose",
                   z_level = "patient", tau = 0.5, seed = 7)
  je_th2 <- tapply(d2$z, d2$therapist_id, function(v) length(unique(v)))
  expect_true(any(je_th2 > 1))  # Patienten-Ebene: variiert im Cluster
})

test_that("Dosis wirkt per Einheit: Steigung ~ tau", {
  d <- sim_stream(n_therapists = 40, patients_per_therapist = 10,
                  n_sessions = 4, z_type = "dose", z_level = "patient",
                  icc = 0.05, tau = 2, tau_x = 0, seed = 11)
  pm <- aggregate(score ~ patient_id + z, data = d, FUN = mean)
  fit <- lm(score ~ z, data = pm)
  expect_gt(coef(fit)[["z"]], 1.5)
  expect_lt(coef(fit)[["z"]], 2.5)
})

test_that("p_treated in Dosis-Welt warnt; confounding stoppt", {
  expect_warning(
    sim_stream(n_therapists = 3, patients_per_therapist = 2,
               n_sessions = 2, z_type = "dose", p_treated = 0.3, seed = 1),
    "ignored")
  expect_error(
    sim_stream(n_therapists = 3, patients_per_therapist = 2,
               n_sessions = 2, z_type = "dose", z_level = "patient",
               confounding = 0.5, seed = 1),
    "not implemented")
})

test_that("Binaer-Waechter: Arm-Schaetzer verweigern Dosis-Welten", {
  d <- sim_stream(n_therapists = 4, patients_per_therapist = 5,
                  n_sessions = 2, z_type = "dose", z_level = "patient",
                  tau = 0.5, seed = 3)
  expect_error(fit_cate_pai(d, B = 0, seed = 1), "binary treatment")
  expect_error(fit_cate_bcf_ml(d, nburn = 10, nsim = 10, seed = 1),
               "binary treatment")
  expect_error(fit_het_grf(d), "binary treatment")
})

test_that("fit_z_satt() traegt die Dosis-Welt: Steigung ~ tau", {
  d <- sim_stream(n_therapists = 40, patients_per_therapist = 10,
                  n_sessions = 4, z_type = "dose", z_level = "patient",
                  icc = 0.05, tau = 2, tau_x = 0, seed = 11)
  fz <- fit_z_satt(d)
  expect_gt(fz[["estimate"]], 1.5)
  expect_lt(fz[["estimate"]], 2.5)
  expect_lt(fz[["p"]], 0.05)
})
