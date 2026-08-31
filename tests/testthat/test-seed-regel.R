# Die Seed-Regel gilt an JEDEM Einstiegspunkt, der Zufallszahlen
# zieht -- nicht nur dort, wo sie zuerst eingebaut wurde. `seed = NULL`
# besteht `missing()` und schickt `set.seed()` auf die Uhr; gemessen
# lief derselbe Aufruf dann zweimal verschieden durch (t-/x-/dr-Learner,
# fit_cate_bcf). Dieser Test haelt die Regel an allen Stellen fest.

welt <- function() {
  sim_stream(n_therapists = 4, patients_per_therapist = 3, n_sessions = 2,
             z_level = "patient", tau = 0.4, tau_c = 0.2, tau_xc = 0.2,
             seed = 1)
}

test_that("seed = NULL wird an jedem stochastischen Einstieg abgelehnt", {
  s <- welt()
  # Name -> Pakete, die der Aufruf VOR der Seed-Wache verlangt
  faelle <- list(
    list(f = function(sd) fit_cate_grf(s, seed = sd), pkg = "grf"),
    list(f = function(sd) fit_cate_bcf(s, seed = sd), pkg = "bcf"),
    list(f = function(sd) fit_cate_bcf_ml(s, seed = sd), pkg = "stochtree"),
    list(f = function(sd) fit_cate_bcf_dyade(s, seed = sd), pkg = "stochtree"),
    list(f = function(sd) fit_cate_dyade_v2p(s, seed = sd), pkg = "stochtree"),
    list(f = function(sd) fit_cate_tlearner(s, seed = sd), pkg = "ranger"),
    list(f = function(sd) fit_cate_xlearner(s, seed = sd), pkg = "ranger"),
    list(f = function(sd) fit_cate_drlearner(s, seed = sd), pkg = "ranger"),
    list(f = function(sd) fit_cate_rlearner(s, seed = sd), pkg = "ranger"),
    list(f = function(sd) fit_cate_sboost(s, seed = sd), pkg = "xgboost"),
    list(f = function(sd) fit_cate_mob(s, seed = sd), pkg = "model4you"),
    list(f = function(sd) fit_cate_pai(s, seed = sd), pkg = NULL),
    list(f = function(sd) fit_cate_stan4bart(s, seed = sd),
         pkg = c("stan4bart", "dbarts")),
    list(f = function(sd) fit_cate_gpboost(s, seed = sd), pkg = "gpboost"),
    list(f = function(sd) fit_cate_merf(s, seed = sd),
         pkg = c("LongituRF", "randomForest")),
    list(f = function(sd) fit_z_brms(s, seed = sd), pkg = "brms"),
    list(f = function(sd) sim_items(s, items = NULL, seed = sd), pkg = NULL),
    list(f = function(sd) plasmode_world(
           kohorte = data.frame(patient_id = 1:4, therapist_id = c(1, 1, 2, 2),
                                x = c(-1, 0, 1, 2)),
           x_spalte = "x", seed = sd), pkg = NULL),
    list(f = function(sd) sim_stream(n_therapists = 2, seed = sd), pkg = NULL),
    list(f = function(sd) wahrheit_skala(n_therapists = 3, seed = sd),
         pkg = NULL)
  )
  gedeckt <- 0L
  for (fall in faelle) {
    if (!is.null(fall$pkg) &&
        !all(vapply(fall$pkg, requireNamespace, logical(1),
                    quietly = TRUE))) next
    expect_error(fall$f(NULL), "NOT a missing argument")
    gedeckt <- gedeckt + 1L
  }
  # ohne installierte Suggests bleiben die vier Kern-Einstiege uebrig
  expect_gte(gedeckt, 4L)
})

test_that("die uebrigen Seed-Fehlklassen greifen auch ausserhalb sim_stream", {
  skip_if_not_installed("ranger")
  s <- welt()
  expect_error(fit_cate_tlearner(s, seed = c(1, 2)), "vector of length 2")
  expect_error(fit_cate_tlearner(s, seed = "abc"), "type character")
  expect_error(fit_cate_tlearner(s, seed = NA_real_), "NA, NaN or Inf")
  expect_error(fit_cate_tlearner(s, seed = Inf), "NA, NaN or Inf")
  expect_error(fit_cate_tlearner(s), "the argument is missing")
  # Gegenprobe: mit gueltigem Seed laeuft derselbe Aufruf zweimal gleich
  a <- fit_cate_tlearner(s, num_trees = 50, B = 20, seed = 7)
  b <- fit_cate_tlearner(s, num_trees = 50, B = 20, seed = 7)
  expect_identical(a, b)
})
