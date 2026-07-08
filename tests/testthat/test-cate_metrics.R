test_that("perfekte Schaetzung: r=1, pehe=0, bias=0, error=FALSE", {
  truth <- c(0.1, 0.3, 0.5, 0.7)
  m <- cate_metrics(truth, truth)
  expect_equal(m$r, 1)
  expect_equal(m$pehe, 0)
  expect_equal(m$bias, 0)
  expect_equal(m$n, 4L)
  expect_false(m$error)
})

test_that("PAI-Degeneration sichtbar: affiner Schaetzer hat r=1, aber pehe>0", {
  x <- seq(-2, 2, length.out = 50)
  truth <- 0.3 + 0.2 * x
  tau_hat <- 3 * truth + 1   # affin verzerrt: Ranking perfekt, Betraege falsch
  m <- cate_metrics(tau_hat, truth)
  expect_equal(m$r, 1)
  expect_gt(m$pehe, 1)       # massiver Magnitude-Fehler
  expect_gt(m$bias, 1)
})

test_that("NULL (ATE-only / gescheiterter Fit) ergibt NA-Zeile mit error=TRUE", {
  m <- cate_metrics(NULL, c(0.2, 0.4))
  expect_true(m$error)
  expect_true(is.na(m$r) && is.na(m$pehe) && is.na(m$covered))
  expect_equal(m$n, 2L)
})

test_that("Personen-Intervall-Coverage zaehlt korrekt", {
  truth <- c(0, 0, 0, 0)
  tau_hat <- c(0.1, 0.1, 0.1, 0.1)
  lo <- c(-1, -1, 0.05, 0.05)   # 2 von 4 Intervallen decken die 0
  hi <- c(1, 1, 0.2, 0.2)
  m <- cate_metrics(tau_hat, truth, lo, hi)
  expect_equal(m$covered, 0.5)
})

test_that("nur lo ohne hi wirft Fehler; Laengen-Mismatch wirft Fehler", {
  expect_error(cate_metrics(c(1, 2), c(1, 2), lo = c(0, 0)), "lo")
  expect_error(cate_metrics(c(1, 2, 3), c(1, 2)))
})

test_that("konstantes tau_hat: r ehrlich NA ohne Warnung, pehe berechnet", {
  truth <- c(0.1, 0.2, 0.3)
  expect_no_warning(m <- cate_metrics(c(0.2, 0.2, 0.2), truth))
  expect_true(is.na(m$r))
  expect_equal(m$pehe, sqrt(mean((0.2 - truth)^2)))
  expect_false(m$error)
})
