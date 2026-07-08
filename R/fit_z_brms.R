#' ATE estimation with a Bayesian mixed model (brms)
#'
#' Bayesian counterpart to [fit_z_satt()]: the same mixed model
#' `score ~ z + (1 | therapist_id) + (1 | patient_id)`, estimated with
#' brms/Stan. Reports the posterior mean, posterior SD, and the
#' central 95% credible interval of the treatment coefficient.
#'
#' Specification: brms default priors (flat prior on the treatment
#' coefficient, `student_t(3, ...)` on intercept and SD parameters),
#' 2 chains, 1000 warmup + 1000 draws by default. The Stan model is
#' compiled once per session and refilled per dataset via `update()`
#' (deterministic given `seed`). Note that with the default draw
#' count Stan may warn about low effective sample sizes; increase
#' `draws` for final analyses if those warnings matter for your use.
#'
#' @param stream Data stream from [sim_stream()] (or a snapshot).
#' @param chains,warmup,draws MCMC settings (defaults 2/1000/1000).
#' @param seed Random seed (mandatory).
#' @return Named vector `estimate`, `se` (posterior SD), `lo`, `hi`
#'   (2.5/97.5% quantiles).
#' @export
fit_z_brms <- local({
  cache <- new.env(parent = emptyenv())
  function(stream, chains = 2, warmup = 1000, draws = 1000, seed) {
    if (!requireNamespace("brms", quietly = TRUE)) {
      stop("fit_z_brms() requires the 'brms' package.", call. = FALSE)
    }
    if (missing(seed)) stop("`seed` is mandatory.", call. = FALSE)
    d <- as.data.frame(stream)
    if (is.null(cache$modell)) {
      cache$modell <- brms::brm(
        score ~ z + (1 | therapist_id) + (1 | patient_id),
        data = d, chains = chains, warmup = warmup,
        iter = warmup + draws, seed = seed, refresh = 0,
        silent = 2)
      fit <- cache$modell
    } else {
      fit <- stats::update(cache$modell, newdata = d,
                           chains = chains, warmup = warmup,
                           iter = warmup + draws, seed = seed,
                           refresh = 0, silent = 2)
    }
    dr <- as.data.frame(fit)[["b_z"]]
    c(estimate = mean(dr), se = stats::sd(dr),
      lo = unname(stats::quantile(dr, 0.025)),
      hi = unname(stats::quantile(dr, 0.975)))
  }
})
