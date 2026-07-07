test_that("preset liefert gueltige sim_stream-Argumente + Quellen", {
  p <- preset("ambulanz_de")
  expect_true(all(c("n_sessions", "icc", "dropout") %in% names(p)))
  src <- attr(p, "sources")
  expect_s3_class(src, "data.frame")
  expect_true(all(c("parameter", "quelle", "status") %in% names(src)))
  s <- do.call(sim_stream, c(p, list(n_therapists = 5,
                                     patients_per_therapist = 4,
                                     seed = 1)))
  expect_s3_class(s, "data.frame")
})

test_that("unbekanntes Preset wird klar abgelehnt", {
  expect_error(preset("mars_kolonie"), "Unbekanntes Preset")
})

test_that("kumulative Abbruchquote liegt nahe der KODAP-Marke (~33%)", {
  p <- preset("ambulanz_de")
  s <- do.call(sim_stream, c(p, list(n_therapists = 40,
                                     patients_per_therapist = 25,
                                     seed = 7)))
  fertig <- tapply(s$session, s$patient_id, max) == 24
  abbruch <- 1 - mean(fertig)
  expect_gt(abbruch, 0.25)
  expect_lt(abbruch, 0.41)
})

test_that("loglinear: fruehe Sitzungen bringen mehr als spaete", {
  s <- sim_stream(n_therapists = 30, patients_per_therapist = 15,
                  n_sessions = 24, shape = "loglinear",
                  mean_slope = -0.40, seed = 10)
  m <- tapply(s$score, s$session, mean)
  frueh <- m[2] - m[1]    # Gewinn Sitzung 1 -> 2
  spaet <- m[24] - m[23]  # Gewinn Sitzung 23 -> 24
  expect_lt(frueh, spaet) # frueh staerker negativ (mehr Besserung)
})

test_that("Preset erreicht die Ziel-Effektstaerke ~d=0.9 prae-post", {
  p <- preset("ambulanz_de")
  p$dropout <- 0  # fuer den reinen Verlaufs-Check ohne Selektion
  s <- do.call(sim_stream, c(p, list(n_therapists = 40,
                                     patients_per_therapist = 25,
                                     seed = 11)))
  prae <- s$score[s$session == 1]
  post <- s$score[s$session == 24]
  d <- (mean(prae) - mean(post)) / sd(prae)
  expect_gt(d, 0.7)
  expect_lt(d, 1.1)
})

test_that("icc_slope = 0 laesst alte Welten bit-identisch", {
  a <- sim_stream(n_therapists = 3, patients_per_therapist = 4,
                  n_sessions = 5, seed = 1)
  b <- sim_stream(n_therapists = 3, patients_per_therapist = 4,
                  n_sessions = 5, icc_slope = 0, seed = 1)
  expect_identical(a, b)
})

test_that("icc_slope > 0: Therapeuten unterscheiden sich im Tempo", {
  tempo_var <- function(icc_s) {
    s <- sim_stream(n_therapists = 40, patients_per_therapist = 15,
                    n_sessions = 6, sd_slope = 0.25,
                    icc_slope = icc_s, seed = 9)
    # beobachtetes Tempo pro Person: (letzte - erste Sitzung) / 5
    tempo <- tapply(s$score, s$patient_id,
                    function(v) (v[length(v)] - v[1]) / 5)
    ther  <- tapply(s$therapist_id, s$patient_id, function(v) v[1])
    var(tapply(tempo, ther, mean))
  }
  expect_gt(tempo_var(0.6), tempo_var(0) * 2)
})

test_that("allianz_beierl2021 laesst sich mit ambulanz_de kombinieren", {
  p <- c(preset("ambulanz_de"), preset("allianz_beierl2021"))
  expect_true(!any(duplicated(names(p))))  # keine Kollisionen
  s <- do.call(sim_stream, c(p, list(n_therapists = 6,
                                     patients_per_therapist = 4,
                                     seed = 1)))
  expect_true("alliance" %in% names(s))
})
