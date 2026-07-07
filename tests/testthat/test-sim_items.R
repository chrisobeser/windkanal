mach_items <- function(k, nopt = 4, rev_at = integer(0)) {
  data.frame(item = paste0("sym_", seq_len(k)), scale = "sym",
             n_options = nopt,
             reversed = seq_len(k) %in% rev_at)
}

test_that("sim_items: Struktur, Wertebereich, Reproduzierbarkeit", {
  s <- sim_stream(n_therapists = 5, patients_per_therapist = 4,
                  n_sessions = 3, seed = 1)
  it <- mach_items(6)
  a <- sim_items(s, it, seed = 1)
  b <- sim_items(s, it, seed = 1)
  expect_identical(a, b)
  expect_equal(nrow(a), nrow(s))
  expect_true(all(a$sym_3 %in% 1:4))
})

test_that("Umkehr-Items korrelieren roh negativ mit normalen", {
  s <- sim_stream(n_therapists = 15, patients_per_therapist = 10,
                  n_sessions = 4, seed = 2)
  a <- sim_items(s, mach_items(4, rev_at = 4), seed = 1)
  expect_lt(cor(a$sym_1, a$sym_4), -0.1)
})

test_that("Reliabilitaet ist emergent und plausibel (10 Items ~ .9)", {
  s <- sim_stream(n_therapists = 20, patients_per_therapist = 10,
                  n_sessions = 4, seed = 3)
  sc10 <- scale_scores(sim_items(s, mach_items(10, rev_at = 9:10),
                                 seed = 1))
  sc4  <- scale_scores(sim_items(s, mach_items(4), seed = 1))
  a10 <- attr(sc10, "alpha")[["sym"]]
  a4  <- attr(sc4, "alpha")[["sym"]]
  expect_gt(a10, 0.75); expect_lt(a10, 0.95)
  expect_gt(a10, a4)   # mehr Items -> hoehere Reliabilitaet
})

test_that("Skalensumme bildet den latenten Score ab", {
  s <- sim_stream(n_therapists = 20, patients_per_therapist = 10,
                  n_sessions = 4, seed = 4)
  sc <- scale_scores(sim_items(s, mach_items(10), seed = 1))
  expect_gt(cor(sc$sym, s$score), 0.85)
})

test_that("funktioniert mit read_items()-Output (Fixture)", {
  s <- sim_stream(n_therapists = 5, patients_per_therapist = 4,
                  n_sessions = 3, seed = 5)
  it <- read_items(testthat::test_path("fixtures", "mini_battery.csv"))
  a <- sim_items(s, it, seed = 1)
  sc <- scale_scores(a)
  expect_true(all(c("sym", "all") %in% names(sc)))
})

test_that("careless = 0 laesst alles unveraendert; Zeit-Attribut da", {
  s <- sim_stream(n_therapists = 5, patients_per_therapist = 4,
                  n_sessions = 3, seed = 6)
  it <- mach_items(8)
  a <- sim_items(s, it, seed = 1)
  b <- sim_items(s, it, careless = 0, seed = 1)
  expect_identical(a, b)
  expect_equal(attr(a, "battery_seconds"), 8 * 7)
})

test_that("Autopilot erzeugt Straightlining-Laeufe", {
  s <- sim_stream(n_therapists = 15, patients_per_therapist = 10,
                  n_sessions = 3, seed = 7)
  it <- mach_items(12)
  longstring <- function(d) {
    m <- as.matrix(d[, it$item]); mean(apply(m, 1, function(v)
      max(rle(v)$lengths)))
  }
  sauber <- sim_items(s, it, seed = 1)
  kaputt <- sim_items(s, it, careless = 0.15, seed = 1)
  expect_gt(longstring(kaputt), longstring(sauber) + 1)
})

test_that("die Alpha-Falle: Straightlining blaeht Reliabilitaet auf", {
  s <- sim_stream(n_therapists = 20, patients_per_therapist = 10,
                  n_sessions = 4, seed = 8)
  it <- mach_items(10)
  a_sauber <- attr(scale_scores(sim_items(s, it, seed = 1)),
                   "alpha")[["sym"]]
  a_kaputt <- attr(scale_scores(sim_items(s, it, careless = 0.10,
                                          seed = 1)),
                   "alpha")[["sym"]]
  expect_gt(a_kaputt, a_sauber)
})

test_that("Allianz-Items messen die Allianz, nicht die Symptome", {
  s <- sim_stream(n_therapists = 20, patients_per_therapist = 10,
                  n_sessions = 6, alliance = TRUE, seed = 9)
  it <- data.frame(item = c(paste0("sym_", 1:4), paste0("all_", 1:4)),
                   scale = rep(c("sym", "alli"), each = 4),
                   n_options = 4, reversed = FALSE,
                   latent = rep(c("score", "alliance"), each = 4))
  sc <- scale_scores(sim_items(s, it, seed = 1))
  expect_gt(cor(sc$alli, s$alliance), cor(sc$alli, s$score) + 0.2)
  expect_gt(cor(sc$sym, s$score), cor(sc$sym, s$alliance) + 0.2)
})

test_that("Allianz-Items ohne Allianz-Strom werden abgelehnt", {
  s <- sim_stream(n_therapists = 4, patients_per_therapist = 3,
                  n_sessions = 3, seed = 9)
  it <- data.frame(item = "all_1", scale = "alli", n_options = 4,
                   reversed = FALSE, latent = "alliance")
  expect_error(sim_items(s, it, seed = 1), "alliance")
})

test_that("refusal = 0 unveraendert; refusal > 0 erzeugt NA-Sitzungen", {
  s <- sim_stream(n_therapists = 15, patients_per_therapist = 10,
                  n_sessions = 4, seed = 11)
  it <- mach_items(6)
  a <- sim_items(s, it, seed = 1)
  b <- sim_items(s, it, refusal = 0, seed = 1)
  expect_identical(a, b)
  c <- sim_items(s, it, refusal = 0.25, seed = 1)
  anteil_na <- mean(is.na(c$sym_1))
  expect_gt(anteil_na, 0.15); expect_lt(anteil_na, 0.35)
  # Verweigerung trifft die ganze Batterie, nicht einzelne Items
  expect_identical(is.na(c$sym_1), is.na(c$sym_6))
  # scale_scores ueberlebt NAs
  sc <- scale_scores(c)
  expect_true(is.finite(attr(sc, "alpha")[["sym"]]))
  expect_true(any(is.na(sc$sym)))
})
