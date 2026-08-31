# v2P: Ranking-Erhalt (N1' per Konstruktion) + Zerlegungs-Attribute
test_that("v2P: tau_hat identisch mit v1a; Zerlegung vorhanden", {
  skip_if_not_installed("stochtree")
  s <- sim_stream(n_therapists = 8, patients_per_therapist = 6,
                  n_sessions = 2, z_level = "patient",
                  tau = 0.5, tau_xc = 0.35, seed = 5)
  a <- fit_cate_bcf_dyade(s, nburn = 50, nsim = 25, seed = 1)
  b <- fit_cate_dyade_v2p(s, nburn = 50, nsim = 25, gitter_n = 12, seed = 1)
  expect_equal(a$tau_hat, b$tau_hat)            # N1': Flaeche unangetastet
  z <- attr(b, "zerlegung")
  expect_length(z$sd_h, 25)
  expect_true(all(z$sd_h >= 0))
  expect_equal(dim(z$h_mean), c(12, 12))
  expect_lt(abs(mean(z$h_mean)), 1e-8)          # doppelzentriert
  expect_lt(abs(mean(z$f_mean)), 1e-8)
})

test_that("dyade_h_wahr: doppelzentriert, realisierte Paare, Falle belegt", {
  set.seed(2); x <- rnorm(60); cc <- rnorm(8)
  w <- dyade_h_wahr(x, cc, "misfit_quad", tau_xc = 0.35,
                    therapist_of = rep(1:8, length.out = 60))
  expect_lt(max(abs(rowMeans(w$h))), 1e-12)
  expect_lt(max(abs(colMeans(w$h))), 1e-12)
  expect_length(w$realisiert, 60)
  roh <- 0.35 * dyade_form(x[1], cc[1], "misfit_quad")
  expect_false(isTRUE(all.equal(roh, w$realisiert[1])))  # Form != h
})

test_that("v2P-Benennung: beta-Draws vorhanden, Produkt fuehrt im Signal", {
  skip_if_not_installed("stochtree")
  s <- sim_stream(n_therapists = 8, patients_per_therapist = 6,
                  n_sessions = 2, z_level = "patient",
                  tau = 0.5, tau_xc = 0.35, seed = 5)
  b <- fit_cate_dyade_v2p(s, nburn = 50, nsim = 25, gitter_n = 12,
                          kandidaten = list(
                            prod = function(x, c) x * c,
                            noise = function(x, c) sin(3 * x) * c^2),
                          seed = 1)
  bn <- attr(b, "benennung")
  expect_equal(dim(bn$beta), c(2L, 25L))
  expect_true(all(is.finite(bn$beta)))
  expect_named(bn$p_gross, c("prod", "noise"))
})

test_that("Benennung: p_gross rechnet auf der dokumentierten .02-Schwelle", {
  # deterministisch: h_draws = beta_s * dz(x*c) -> OLS gibt exakt beta_s
  qx <- c(-1, 0, 1); qc <- c(-1, 0, 1)
  kand <- list(prod = function(x, c) x * c)
  K <- v2p_kandidaten_matrix(qx, qc, kand)
  # die Betas klammern .02 eng ein: jede Schwelle ausserhalb
  # [0.018, 0.022) gaebe einen anderen p_gross-Wert
  b <- c(0.018, 0.022, 0.06, 0.10)
  bn <- benennung_bauen(qx, qc, K[, 1] %o% b, kand)
  expect_equal(as.vector(bn$beta), b)
  expect_equal(unname(bn$p_gross), 0.75)             # P(|beta| > .02)
  for (falsch in c(0.015, 0.025, 0.05)) {
    expect_false(isTRUE(all.equal(unname(bn$p_gross),
                                  mean(abs(b) > falsch))))
  }
  expect_match(bn$hinweis, "\\.02")
})

test_that("v2P warnt, wenn die auf 0 gesetzten Features nicht zentriert sind", {
  skip_if_not_installed("stochtree")
  s <- sim_stream(n_therapists = 8, patients_per_therapist = 6,
                  n_sessions = 2, z_level = "patient",
                  tau = 0.5, tau_xc = 0.35, n_noise = 1, seed = 5)
  ohne <- expect_no_warning(fit_cate_dyade_v2p(s, nburn = 10, nsim = 5,
                                               gitter_n = 4, seed = 1))
  s2 <- s
  s2$x_noise1 <- s2$x_noise1 + 3
  expect_warning(fit_cate_dyade_v2p(s2, nburn = 10, nsim = 5,
                                    gitter_n = 4, seed = 1),
                 class = "windkanal_warnung_gitter_nullfeature")
  expect_warning(mit <- fit_cate_dyade_v2p(s2, nburn = 10, nsim = 5,
                                           gitter_n = 4, seed = 1),
                 "x_noise1")
  # die Wache aendert nichts: derselbe Strom, gewarnt oder nicht, gibt
  # denselben Fit (sie zieht keine Zufallszahlen)
  expect_warning(mit2 <- fit_cate_dyade_v2p(s2, nburn = 10, nsim = 5,
                                            gitter_n = 4, seed = 1))
  expect_identical(mit, mit2)
  expect_false(identical(ohne$tau_hat, mit$tau_hat))   # Strom wirklich anders
  # NA-Robustheit: kein Spaltenname "NA", kein Absturz
  s3 <- s2
  s3$x_noise1[1] <- NA_real_
  expect_warning(fit_cate_dyade_v2p(s3, nburn = 10, nsim = 5,
                                    gitter_n = 4, seed = 1),
                 "x_noise1")
})

test_that("v2P-Guards: gitter_n >= 2, Kollinearitaet bricht laut ab", {
  skip_if_not_installed("stochtree")
  s <- sim_stream(n_therapists = 8, patients_per_therapist = 6,
                  n_sessions = 2, z_level = "patient",
                  tau = 0.5, tau_xc = 0.35, seed = 5)
  expect_error(fit_cate_dyade_v2p(s, nburn = 10, nsim = 5,
                                  gitter_n = 1, seed = 1),
               "gitter_n")
  # gebrochen: (1:2.7)/(2.7+1) traefe .270/.541 statt 1/3, 2/3 --
  # ein schiefes Referenzgitter, auf dem die Doppelzentrierung
  # stillschweigend weiterrechnet
  expect_error(fit_cate_dyade_v2p(s, nburn = 10, nsim = 5,
                                  gitter_n = 2.7, seed = 1),
               "ganze Zahl")
  expect_error(fit_cate_dyade_v2p(s, nburn = 10, nsim = 5,
                                  gitter_n = c(4, 9), seed = 1),
               "ganze Zahl")
  expect_error(fit_cate_dyade_v2p(s, nburn = 10, nsim = 5,
                                  gitter_n = Inf, seed = 1),
               "ganze Zahl")
  # dz(-(x-c)^2) = 2*dz(x*c) exakt ->
  # gemeinsame Benennung muss mit informativem Fehler abbrechen
  expect_error(
    fit_cate_dyade_v2p(s, nburn = 50, nsim = 10, gitter_n = 12,
                       kandidaten = list(
                         prod = function(x, c) x * c,
                         quad = function(x, c) -(x - c)^2),
                       seed = 1),
    "kollinear")
})

test_that("v2P empirisches Gitter: Identitaet, h_real, Attribut-Benennung", {
  skip_if_not_installed("stochtree")
  s <- sim_stream(n_therapists = 8, patients_per_therapist = 6,
                  n_sessions = 2, z_level = "patient",
                  tau = 0.5, tau_xc = 0.35, seed = 5)
  p <- patients(s)
  n <- nrow(p); J <- length(unique(p$therapist_id))
  set.seed(11)
  kp <- list(
    patient = cbind(x = p$x, px1 = rnorm(n)),
    therapeut = cbind(
      c = p$therapist_c[match(unique(p$therapist_id), p$therapist_id)],
      tc1 = rnorm(J)))
  a <- fit_cate_bcf_dyade(s, nburn = 50, nsim = 25, seed = 1)
  b <- fit_cate_dyade_v2p(s, nburn = 50, nsim = 25,
                          gitter = "empirisch", kandidaten_prod = kp,
                          halte_h = TRUE, seed = 1)
  expect_equal(a$tau_hat, b$tau_hat)            # N1' auch am Kreuz
  z <- attr(b, "zerlegung")
  expect_equal(dim(z$h_mean), c(n, J))
  expect_lt(max(abs(rowMeans(z$h_mean))), 1e-8) # doppelzentriert
  expect_lt(max(abs(colMeans(z$h_mean))), 1e-8)
  expect_equal(dim(attr(b, "h_real")), c(n, 25L))
  expect_length(z$h_real_mean, n)
  bn <- attr(b, "benennung")
  expect_equal(rownames(bn$beta),
               c("prod_x_c", "prod_x_tc1", "prod_px1_c", "prod_px1_tc1"))
  # wahres Produkt fuehrt im Signal (Rangordnung, keine Schwelle)
  expect_equal(names(which.max(rowMeans(abs(bn$beta)))), "prod_x_c")
  # Attribut-Produkte ohne empirisches Gitter -> Fehler
  expect_error(fit_cate_dyade_v2p(s, nburn = 10, nsim = 5,
                                  kandidaten_prod = kp, seed = 1),
               "empirisch")
})

test_that("v2p_benennung post hoc == Fit-Benennung; v2p_entscheidung", {
  skip_if_not_installed("stochtree")
  s <- sim_stream(n_therapists = 8, patients_per_therapist = 6,
                  n_sessions = 2, z_level = "patient",
                  tau = 0.5, tau_xc = 0.35, seed = 5)
  kand <- list(prod = function(x, c) x * c)
  b <- fit_cate_dyade_v2p(s, nburn = 50, nsim = 25, gitter_n = 12,
                          kandidaten = kand, halte_h = TRUE, seed = 1)
  nach <- v2p_benennung(b, kandidaten = kand)
  expect_equal(nach$beta, attr(b, "benennung")$beta)  # deterministisch
  e <- v2p_entscheidung(b, delta_sd = 0.02, cut_sd = 0.8,
                        delta_beta = 0.05, cut_beta = 0.8)
  expect_type(e$detektiert, "logical")
  expect_true(e$p_sd >= 0 && e$p_sd <= 1)
  expect_named(e$p_beta, "prod")
  # ohne h_draws bricht die post-hoc-Benennung informativ ab
  b2 <- fit_cate_dyade_v2p(s, nburn = 50, nsim = 10, gitter_n = 12,
                           seed = 1)
  expect_error(v2p_benennung(b2, kandidaten = kand), "halte_h")
})
