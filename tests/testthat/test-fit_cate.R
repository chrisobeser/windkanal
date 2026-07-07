skip_if_not_installed("grf")

test_that("patients() verdichtet korrekt auf Personen-Ebene", {
  s <- sim_stream(n_therapists = 4, patients_per_therapist = 3,
                  n_sessions = 5, seed = 1)
  p <- patients(s)
  expect_equal(nrow(p), 12)
  expect_equal(p$n_obs, rep(5, 12))
})

test_that("fit_cate_grf: Struktur, Seed-Pflicht, ATE nahe tau", {
  s <- sim_stream(n_therapists = 20, patients_per_therapist = 15,
                  n_sessions = 4, tau = 0.5, tau_x = 0, seed = 2)
  expect_error(fit_cate_grf(s), "seed")
  ct <- fit_cate_grf(s, num_trees = 500, seed = 1)
  expect_equal(nrow(ct), 300)
  ate <- attr(ct, "ate")
  expect_true(ate[["estimate"]] - 2.5 * ate[["se"]] < 0.5 &&
              0.5 < ate[["estimate"]] + 2.5 * ate[["se"]])
})

test_that("echte Heterogenitaet wird gefunden (tau_hat folgt x)", {
  s <- sim_stream(n_therapists = 30, patients_per_therapist = 15,
                  n_sessions = 4, tau = 0.3, tau_x = 0.5, seed = 3)
  ct <- fit_cate_grf(s, num_trees = 500, seed = 1)
  expect_gt(cor(ct$tau_hat, ct$x), 0.5)
})

test_that("homogene Welt: kaum erfundene Streuung in tau_hat", {
  hom <- sim_stream(n_therapists = 30, patients_per_therapist = 15,
                    n_sessions = 4, tau = 0.5, tau_x = 0, seed = 4)
  het <- sim_stream(n_therapists = 30, patients_per_therapist = 15,
                    n_sessions = 4, tau = 0.5, tau_x = 0.5, seed = 4)
  sd_hom <- sd(fit_cate_grf(hom, num_trees = 500, seed = 1)$tau_hat)
  sd_het <- sd(fit_cate_grf(het, num_trees = 500, seed = 1)$tau_hat)
  # Lehrstueck: auch in der homogenen Welt streuen rohe tau_hat
  # sichtbar (~0.27) -- Streuung allein ist KEIN Heterogenitaets-
  # Beleg; dafuer braucht es Inferenz (Kalibrierungstest, spaeter).
  expect_lt(sd_hom, sd_het)
})

test_that("patients() bietet ehrliche Outcome-Varianten", {
  s <- sim_stream(n_therapists = 10, patients_per_therapist = 8,
                  n_sessions = 8, dropout = 0.15,
                  dropout_informative = 1, seed = 5)
  pm <- patients(s)                       # mean (Default)
  pa <- patients(s, outcome = "at_session", session = 8)
  ps <- patients(s, outcome = "slope")
  expect_true("score_mean" %in% names(pm))
  expect_true(all(pa$n_obs == 8))         # nur Vollendete bei Sitzung 8
  expect_lt(nrow(pa), nrow(pm))           # Dropout reduziert
  expect_true("score_slope" %in% names(ps))
})
