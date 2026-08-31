skip_if_not_installed("grf")

test_that("fit_het_grf hat das fit_*-Interface mit p-Wert", {
  s <- sim_stream(n_therapists = 15, patients_per_therapist = 10,
                  n_sessions = 3, seed = 2)
  fit <- fit_het_grf(s, num_trees = 300)
  expect_named(fit, c("estimate", "se", "p"))
  expect_true(fit[["p"]] >= 0 && fit[["p"]] <= 1)
})

test_that("echte Heterogenitaet wird erkannt (gut gepowert)", {
  s <- sim_stream(n_therapists = 30, patients_per_therapist = 15,
                  n_sessions = 4, tau = 0.3, tau_x = 0.5, seed = 3)
  fit <- fit_het_grf(s, num_trees = 500)
  expect_lt(fit[["p"]], 0.05)
})

test_that("dyadische Heterogenitaet wird gefunden (W2)", {
  # tau_x = 0, tau_xc = 1.5: die Heterogenitaet liegt vollstaendig in
  # der Paarung. Bis 0.2.1 sah der Forest nur x und meldete ein
  # sauberes Null (0 von 10 Welten bei wahrem SD(tau_i) = 1.53).
  s <- sim_stream(n_therapists = 40, patients_per_therapist = 12,
                  n_sessions = 6, z_level = "patient",
                  tau = 0.5, tau_x = 0, tau_xc = 1.5, icc = 0.1,
                  seed = 1001)
  alle <- fit_het_grf(s, num_trees = 500)
  expect_lt(alle[["p"]], 0.05)
  expect_gt(alle[["estimate"]], 0.5)
  expect_true("therapist_c" %in% attr(alle, "feature_namen"))
  # Negativfall: der deklarierte Naivling sieht nur x -- und nichts
  nur_x <- fit_het_grf(s, num_trees = 500, features = "nur_x")
  expect_gt(nur_x[["p"]], 0.05)
  expect_identical(attr(nur_x, "feature_namen"), "x")
})

test_that("Kalibrierung entlang x bleibt bei ~1 (W2, Gegenprobe)", {
  s <- sim_stream(n_therapists = 40, patients_per_therapist = 12,
                  n_sessions = 6, z_level = "patient",
                  tau = 0.5, tau_x = 1.5, tau_xc = 0, icc = 0.1,
                  seed = 1001)
  fit <- fit_het_grf(s, num_trees = 500)
  expect_lt(fit[["p"]], 0.05)
  expect_equal(fit[["estimate"]], 1, tolerance = 0.25)
})

test_that("ohne therapist_c/Rauschmerkmale sind beide Mengen dieselbe (W2)", {
  # Struktur-Aussage, an der die Frage haengt, welche gelaufenen Zahlen
  # die Umstellung beruehrt: zieht eine Welt kein therapist_c (also
  # tau_c = tau_xc = 0, assignment = "random") und keine
  # x_noise-Spalten, dann ist cate_features() genau x -- die Umstellung
  # ist dort folgenlos.
  s <- sim_stream(n_therapists = 20, patients_per_therapist = 10,
                  n_sessions = 4, z_level = "patient",
                  tau = 0.5, tau_x = 0.5, icc = 0.1, seed = 1)
  expect_null(patients(s)$therapist_c)
  expect_identical(colnames(cate_features(patients(s))), "x")
  expect_identical(fit_het_grf(s, num_trees = 300),
                   fit_het_grf(s, num_trees = 300, features = "nur_x"),
                   ignore_attr = TRUE)
  # Negativfall: sobald tau_xc wirkt, traegt der Strom therapist_c und
  # die beiden Mengen fallen auseinander
  s2 <- sim_stream(n_therapists = 20, patients_per_therapist = 10,
                   n_sessions = 4, z_level = "patient",
                   tau = 0.5, tau_xc = 1, icc = 0.1, seed = 1)
  expect_false(is.null(patients(s2)$therapist_c))
  expect_false(identical(fit_het_grf(s2, num_trees = 300)[["estimate"]],
                         fit_het_grf(s2, num_trees = 300,
                                     features = "nur_x")[["estimate"]]))
})
