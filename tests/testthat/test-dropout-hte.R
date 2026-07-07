test_that("dropout = 0: niemand bricht ab (alte Welt bleibt gueltig)", {
  s <- sim_stream(n_therapists = 4, patients_per_therapist = 3,
                  n_sessions = 5, seed = 1)
  expect_equal(nrow(s), 4 * 3 * 5)
})

test_that("dropout > 0: weniger Zeilen, aber Sitzung 1 immer da", {
  s <- sim_stream(n_therapists = 6, patients_per_therapist = 5,
                  n_sessions = 8, dropout = 0.15, seed = 2)
  expect_lt(nrow(s), 6 * 5 * 8)
  expect_equal(length(unique(s$patient_id)), 30)  # alle starten
  # keine Luecken: beobachtete Sitzungen je Person sind 1..k
  ok <- tapply(s$session, s$patient_id,
               function(v) identical(as.integer(v), seq_along(v)))
  expect_true(all(ok))
})

test_that("informativer Dropout: wem es schlechter geht, bleibt kuerzer", {
  s <- sim_stream(n_therapists = 20, patients_per_therapist = 10,
                  n_sessions = 10, dropout = 0.10,
                  dropout_informative = 1.5, seed = 3)
  n_obs  <- tapply(s$session, s$patient_id, max)
  start  <- tapply(s$score, s$patient_id, function(v) v[1])
  # hoher Startscore (= schlechter) -> weniger beobachtete Sitzungen
  expect_lt(cor(start, n_obs), -0.1)
})

test_that("x steht im Output und ist pro Person konstant", {
  s <- sim_stream(n_therapists = 4, patients_per_therapist = 3,
                  n_sessions = 4, seed = 4)
  expect_true("x" %in% names(s))
  expect_true(all(tapply(s$x, s$patient_id,
                         function(v) length(unique(v))) == 1))
})

test_that("tau_x wird vom gemischten Modell wiedergefunden", {
  skip_if_not_installed("lmerTest")
  s <- sim_stream(n_therapists = 30, patients_per_therapist = 12,
                  n_sessions = 6, tau = 0.3, tau_x = 0.4, seed = 5)
  fit <- fit_zx_satt(s)
  expect_true(fit[["estimate"]] - 1.96 * fit[["se"]] < 0.4)
  expect_true(0.4 < fit[["estimate"]] + 1.96 * fit[["se"]])
})

test_that("reliability_score < 1 vergroessert die Score-Varianz passend", {
  a <- sim_stream(n_therapists = 20, patients_per_therapist = 10,
                  n_sessions = 6, seed = 8)
  b <- sim_stream(n_therapists = 20, patients_per_therapist = 10,
                  n_sessions = 6, reliability_score = 0.5, seed = 8)
  expect_equal(var(a$score) / var(b$score), 0.5, tolerance = 0.1)
})

test_that("unreliables x daempft tau_x um den Faktor reliability_x", {
  skip_if_not_installed("lmerTest")
  s <- sim_stream(n_therapists = 40, patients_per_therapist = 15,
                  n_sessions = 6, tau = 0.3, tau_x = 0.4,
                  reliability_x = 0.5, seed = 12)
  fit <- fit_zx_satt(s)
  # erwartete Schaetzung: 0.4 * 0.5 = 0.2 -- die Wahrheit 0.4 liegt
  # AUSSERHALB des CIs, die gedaempfte 0.2 innerhalb
  expect_lt(fit[["estimate"]] + 1.96 * fit[["se"]], 0.4)
  expect_true(fit[["estimate"]] - 1.96 * fit[["se"]] < 0.2 &&
              0.2 < fit[["estimate"]] + 1.96 * fit[["se"]])
})

test_that("alliance = FALSE laesst alte Welten bit-identisch", {
  a <- sim_stream(n_therapists = 3, patients_per_therapist = 4,
                  n_sessions = 5, seed = 1)
  b <- sim_stream(n_therapists = 3, patients_per_therapist = 4,
                  n_sessions = 5, alliance = FALSE, seed = 1)
  expect_identical(a, b)
})

test_that("alliance = TRUE liefert den zweiten Strom", {
  s <- sim_stream(n_therapists = 5, patients_per_therapist = 4,
                  n_sessions = 6, alliance = TRUE, seed = 2)
  expect_true("alliance" %in% names(s))
  expect_true(is.numeric(s$alliance))
})

test_that("coupling ohne alliance wird abgelehnt", {
  expect_error(sim_stream(coupling = -0.2, seed = 1), "alliance")
})

test_that("Kopplung wirkt: gute Allianz -> mehr Besserung danach", {
  s <- sim_stream(n_therapists = 30, patients_per_therapist = 15,
                  n_sessions = 8, alliance = TRUE, coupling = -0.3,
                  seed = 3)
  # Allianz in Sitzung s vs. Score in s+1 (innerhalb Person)
  d <- s[order(s$patient_id, s$session), ]
  lag_alli <- ave(d$alliance, d$patient_id,
                  FUN = function(v) c(NA, v[-length(v)]))
  ok <- !is.na(lag_alli)
  expect_lt(cor(lag_alli[ok], d$score[ok]), -0.05)
})

test_that("Allianz ist traege: AR erzeugt hohe Lag-1-Autokorrelation", {
  lag1 <- function(ar) {
    s <- sim_stream(n_therapists = 20, patients_per_therapist = 10,
                    n_sessions = 10, alliance = TRUE,
                    alliance_ar = ar, seed = 6)
    d <- s[order(s$patient_id, s$session), ]
    la <- ave(d$alliance, d$patient_id,
              FUN = function(v) c(NA, v[-length(v)]))
    cor(la, d$alliance, use = "complete.obs")
  }
  expect_gt(lag1(0.75), 0.5)
  expect_gt(lag1(0.75), lag1(0) + 0.3)
})

test_that("coupling_reverse: hohe Belastung senkt die Folge-Allianz", {
  s <- sim_stream(n_therapists = 30, patients_per_therapist = 15,
                  n_sessions = 8, alliance = TRUE,
                  coupling_reverse = -0.5, seed = 7)
  d <- s[order(s$patient_id, s$session), ]
  lscore <- ave(d$score, d$patient_id,
                FUN = function(v) c(NA, v[-length(v)]))
  ok <- !is.na(lscore)
  expect_lt(cor(lscore[ok], d$alliance[ok]), -0.15)
})

test_that("coupling_reverse ohne alliance wird abgelehnt", {
  expect_error(sim_stream(coupling_reverse = -0.1, seed = 1),
               "alliance")
})

test_that("reliability_alliance verrauscht die Messung passend", {
  a <- sim_stream(n_therapists = 20, patients_per_therapist = 10,
                  n_sessions = 6, alliance = TRUE, seed = 8)
  b <- sim_stream(n_therapists = 20, patients_per_therapist = 10,
                  n_sessions = 6, alliance = TRUE,
                  reliability_alliance = 0.5, seed = 8)
  expect_equal(var(a$alliance) / var(b$alliance), 0.5,
               tolerance = 0.1)
})

test_that("reliability_alliance ohne alliance wird abgelehnt", {
  expect_error(sim_stream(reliability_alliance = 0.8, seed = 1),
               "alliance")
})
