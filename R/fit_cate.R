#' Condense the data stream to the person level
#'
#' Many questions live at the person level, not the session level.
#' This function condenses the stream: one row per patient with
#' treatment membership, feature `x`, number of observed sessions,
#' and mean score.
#'
#' **Estimand warning:** Under informative dropout, the mean over
#' *observed* sessions is itself selection-biased (whoever stays
#' longer contributes more and better values). Alternatives:
#' `outcome = "at_session"` (score at a fixed session, only persons
#' observed there) or `outcome = "slope"` (personal trajectory
#' slope).
#'
#' @param stream Data stream from [sim_stream()].
#' @param outcome `"mean"` (default), `"at_session"`, or `"slope"`.
#' @param session Fixed session for `outcome = "at_session"`
#'   (default: last planned session in the stream).
#' @return `data.frame`: `patient_id`, `therapist_id`, `z`, `x`,
#'   `n_obs`, `score_mean` (or `score_at`/`score_slope`).
#' @export
patients <- function(stream, outcome = c("mean", "at_session", "slope"),
                     session = max(stream$session)) {
  outcome <- match.arg(outcome)
  agg <- function(v, f) as.vector(tapply(v, stream$patient_id, f))
  out <- data.frame(
    patient_id   = sort(unique(stream$patient_id)),
    therapist_id = agg(stream$therapist_id, function(v) v[1]),
    z            = agg(stream$z, function(v) v[1]),
    x            = agg(stream$x, function(v) v[1]),
    n_obs        = agg(stream$session, length)
  )
  if (!is.null(stream$therapist_c)) {
    out$therapist_c <- agg(stream$therapist_c, function(v) v[1])
  }
  for (nm in grep("^x_noise", names(stream), value = TRUE)) {
    out[[nm]] <- agg(stream[[nm]], function(v) v[1])
  }
  if (outcome == "mean") {
    out$score_mean <- agg(stream$score, mean)
  } else if (outcome == "at_session") {
    hit <- stream$session == session
    idx <- match(out$patient_id, stream$patient_id[hit])
    out$score_at <- stream$score[hit][idx]
    out <- out[!is.na(out$score_at), , drop = FALSE]
  } else {
    sl <- function(v) if (length(v) < 3) NA_real_ else
      stats::coef(stats::lm(v ~ seq_along(v)))[2]
    out$score_slope <- agg(stream$score, sl)
    out <- out[!is.na(out$score_slope), , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

#' Feature matrix for CATE estimators
#'
#' Collects all candidate features at the person level: `x`, the
#' therapist attribute `therapist_c` (if the stream carries it), and
#' all noise features `x_noise*`. Both CATE estimators use the same
#' selection -- a fair comparison.
#'
#' @param p Person table from [patients()].
#' @return Numeric matrix, one row per person.
#' @keywords internal
cate_features <- function(p) {
  cols <- c("x",
            if (!is.null(p$therapist_c)) "therapist_c",
            grep("^x_noise", names(p), value = TRUE))
  m <- as.matrix(p[cols])
  rownames(m) <- NULL
  m
}

#' Estimate personalized effects (causal forest)
#'
#' Estimates the treatment effect (CATE) for **every person** with
#' `grf::causal_forest` -- the "cautious pattern hunter": it may
#' find heterogeneity, but is regularized against inventing
#' patterns. Nesting is respected: therapists enter the estimation
#' as clusters (cluster-robust inference).
#'
#' The outcome is the mean score per person ([patients()]); the
#' covariate is `x`. Truth in the simulator: individual effect
#' `tau + tau_x * x` (with `reliability_x < 1` operating on *true*
#' x -- the estimate on measured x is attenuated accordingly).
#'
#' @param stream Data stream from [sim_stream()].
#' @param num_trees Number of trees (default 1000).
#' @param honesty Honest splitting (grf's hallmark: one half of the
#'   data grows the trees, the other estimates the values -- buys
#'   honest intervals, costs effective sample size). `FALSE`
#'   disables it. Default `TRUE` (= grf default).
#' @param clusters Respect therapists as clusters (cluster-robust
#'   inference). `FALSE` = ignore nesting -- this is how studies
#'   without cluster correction compute. Default `TRUE`.
#' @param seed Random seed for the forest (mandatory).
#' @return `data.frame` (`patient_id`, `x`, `z`, `tau_hat`) with
#'   attributes `ate` (average effect with SE, cluster-robust) and
#'   `num_trees`.
#' @export
fit_cate_grf <- function(stream, num_trees = 1000, honesty = TRUE,
                         clusters = TRUE, seed) {
  if (!requireNamespace("grf", quietly = TRUE)) {
    stop("fit_cate_grf() requires the 'grf' package.", call. = FALSE)
  }
  if (missing(seed)) {
    stop("`seed` is mandatory (the forest is stochastic).",
         call. = FALSE)
  }
  p <- patients(stream)
  X <- cate_features(p)
  args <- list(X = X, Y = p$score_mean, W = p$z,
               num.trees = num_trees, honesty = honesty, seed = seed)
  if (isTRUE(clusters)) args$clusters <- p$therapist_id
  cf <- do.call(grf::causal_forest, args)
  ate <- grf::average_treatment_effect(cf)
  out <- data.frame(patient_id = p$patient_id, x = p$x, z = p$z,
                    tau_hat = as.vector(cf$predictions))
  attr(out, "ate") <- c(estimate = unname(ate["estimate"]),
                        se = unname(ate["std.err"]))
  attr(out, "num_trees") <- num_trees
  out
}
