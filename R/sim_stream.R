#' Simulate a growing routine-care data stream
#'
#' Generates the wind-tunnel data stream: patients are nested in
#' therapists, enter care staggered over an accrual window, and
#' deliver a session-wise outcome score. The truth is known and
#' configurable -- that is exactly what a wind tunnel is for.
#'
#' Data model:
#' \deqn{score_{ijs} = b_i + u_j + (tau + tau_x x_i) z_j
#'       + slope_i (s - 1) + e_{ijs}}
#' with patient baseline \eqn{b_i ~ N(0,1)}, therapist effect
#' \eqn{u_j ~ N(0, \sigma_u)} (via `icc`), person feature
#' \eqn{x_i ~ N(0,1)}, and residual \eqn{e ~ N(0,1)}.
#'
#' **All "stressors" are off by default** (`tau = 0`, `tau_x = 0`,
#' `dropout = 0`): effects, personalization, or dropout must be
#' switched on explicitly.
#'
#' @param n_therapists Number of therapists.
#' @param patients_per_therapist Cases per therapist.
#' @param n_sessions Maximum sessions per case (weekly).
#' @param icc Therapist intraclass correlation, in `[0, 1)`.
#' @param mean_slope Mean change per session (negative = symptom
#'   reduction).
#' @param sd_slope Total spread of the improvement slopes.
#' @param icc_slope Share of slope variance at the **therapist
#'   level**, in `[0, 1)`: therapists then differ not only in level
#'   but in how fast their patients improve. Empirical anchor: ~17%
#'   of rate variance sits with therapists (Lutz et al., 2007).
#'   Default 0.
#' @param weeks_accrual Accrual window (weeks) of staggered entry.
#' @param shape Shape of the improvement curve: `"linear"` (default,
#'   didactically simple) or `"loglinear"` -- the empirically
#'   supported form (negatively accelerated: largest gains early,
#'   Howard tradition / dose-response research). The trajectory term
#'   is `slope * (session - 1)` or `slope * log(session)`.
#' @param p_treated Share with treatment `z = 1`.
#' @param z_type Type of treatment: `"binary"` (default; `z` in
#'   0/1, the classical two-arm contrast) or `"dose"` (continuous
#'   exposure `z` drawn uniform on `[0, 1]` at the chosen
#'   `z_level`). In dose mode the effect terms `tau`, `tau_x`,
#'   `tau_c`, `tau_xc` act per unit dose: the score contribution is
#'   `[tau + tau_x h(x) + ...] * z` with continuous `z`, so `tau` is
#'   the effect of moving from dose 0 to dose 1. Dose worlds model
#'   questions where amount matters more than assignment (session
#'   dose, alliance exposure, degree of patient-therapist matching).
#'   `p_treated` is ignored (with a warning) and `confounding` is
#'   not implemented in dose mode. Note that the arm-based and
#'   propensity-based CATE wrappers require binary treatment and
#'   refuse dose worlds with an informative error.
#' @param z_level Level of treatment assignment: `"therapist"`
#'   (default -- z constant per therapist, e.g. a feedback system or
#'   training; the hard case for naive methods) or `"patient"` -- z
#'   varies **within** therapists (dyad level, e.g. **matching**:
#'   this pairing fits / does not fit). Outcomes keep clustering in
#'   therapists either way.
#' @param tau Mean treatment effect of `z` on the score level.
#' @param tau_x **Effect heterogeneity**: the z-effect for person i
#'   is `tau + tau_x * x_i`. `tau_x = 0` means: helps everyone
#'   equally.
#' @param tau_x_form Functional form of the heterogeneity:
#'   `"linear"` (default -- contribution `tau_x * x`), `"step"` --
#'   threshold effect `tau_x * 1(x > 0)` (hard for smooth models,
#'   easy for trees: the fairness case for grf/bcf), or
#'   `"quadratic"` -- `tau_x * (x^2 - 1)` (centered so the ATE
#'   stays `tau`).
#' @param n_noise Number of pure **noise features** (`x_noise1`,
#'   ..., standard normal, without any effect) handed to the
#'   estimators as candidate features. Tests whether tools get
#'   distracted by irrelevant variables. Default 0.
#' @param tau_c **Therapist moderation**: the z-effect additionally
#'   depends on a therapist attribute `c` (e.g. a competence/style
#'   dimension, standard normal): contribution `tau_c * c_j`.
#'   Default 0.
#' @param tau_xc **Dyadic matching effect** -- the pairing itself
#'   matters: contribution `tau_xc * x_i * c_j`. Positive =
#'   high-x patients benefit more with high-c therapists (and low-x
#'   with low-c). Only this term makes "who benefits from *whom*"
#'   expressible. With `tau_c != 0` or `tau_xc != 0` the stream
#'   carries the column `therapist_c` (observable to estimators).
#'   Default 0.
#' @param tau_shape Time course of the z-effect: `"constant"`
#'   (default -- full size from session 1) or `"ramp"` -- the effect
#'   grows linearly (`session / n_sessions`), reaching full size
#'   only at the last session. ⚠ The ramp shape is plausible but
#'   not (yet) calibrated against literature.
#' @param x_effect **Prognostic effect** of `x` on the score level
#'   (main effect, independent of treatment). Needed so that
#'   selection on `x` creates true confounding. Default 0.
#' @param confounding **Confounding by indication** (only
#'   `z_level = "patient"`): the treatment probability depends on
#'   `x` -- `P(z=1) = plogis(qlogis(p_treated) + confounding * x)`.
#'   Positive = high-x patients are more likely to receive the
#'   treatment ("targeted selection" in the sense of Hahn et al.).
#'   Confounds the naive comparison only together with
#'   `x_effect != 0` or `tau_x != 0`. Default 0 = coin flip.
#' @param dropout Base probability of dropping out after a session
#'   (at score 0), in `[0, 1)`. `0` = nobody drops out.
#' @param dropout_informative Coupling of dropout to the current
#'   score (logit scale): positive = those doing worse (high score)
#'   are more likely to drop out. Only active with `dropout > 0`.
#' @param reliability_score Reliability of the outcome measurement,
#'   in `(0, 1]`. With `< 1`, classical measurement error is added
#'   to the score (costs precision, does not bias).
#' @param alliance `TRUE` switches on the second trajectory stream:
#'   the **therapeutic alliance** per session, with its own
#'   baseline, its own therapist share (`icc_alliance`), and its own
#'   trend (`alliance_slope`). Default `FALSE` (old worlds
#'   unchanged).
#' @param icc_alliance Therapist share of the alliance variance.
#' @param alliance_slope Mean alliance trend per session
#'   (default 0).
#' @param coupling **The mechanism parameter:** effect of the
#'   previous session's alliance on the current session's symptom
#'   score (negative = good alliance -> stronger improvement).
#'   Default 0 -- null-world principle. Realistic range:
#'   |beta| ~ 0.05-0.15 (Beierl et al., 2021: -0.13
#'   therapist-rated).
#' @param coupling_reverse Reverse direction: effect of the previous
#'   session's symptom score on the current alliance (negative =
#'   high distress -> worse alliance; Beierl et al. 2021: -0.12).
#'   Default 0. The real coupling is **bidirectional** -- simulating
#'   only one direction builds a one-way street that does not exist.
#' @param reliability_alliance Reliability of the alliance
#'   measurement, in `(0, 1]`. With `< 1` the `alliance` column
#'   carries classical measurement error; the **true** alliance
#'   keeps acting undisturbed via `coupling` (error in the
#'   measurement, not in the mechanism).
#' @param alliance_ar Inertia of the alliance (AR(1) share): real
#'   alliance strongly copies itself into the next session. Default
#'   **0.75** (calibrated: Beierl et al. 2021, AR beta .75-.79).
#'   Innovation SD is scaled to stationary variance ~1.
#' @param reliability_x Reliability of the measurement of the person
#'   feature `x`, in `(0, 1]`. With `< 1` the column `x` contains
#'   the **error-laden** measurement while the true effect operates
#'   on true x -- personalization estimates (`tau_x`) are thereby
#'   systematically attenuated toward zero (expected estimate =
#'   `tau_x * reliability_x`).
#' @param seed Random seed. Mandatory: reproducibility is not
#'   optional in the wind tunnel.
#'
#' @return `data.frame` with one row per **observed** session:
#'   `therapist_id`, `z`, `patient_id`, `x`, `entry_week`,
#'   `session`, `week`, `score`.
#'
#' @examples
#' s <- sim_stream(n_therapists = 5, patients_per_therapist = 4,
#'                 n_sessions = 6, icc = 0.05, seed = 42)
#' head(s)
#'
#' @export
sim_stream <- function(n_therapists = 10,
                       patients_per_therapist = 20,
                       n_sessions = 12,
                       icc = 0.05,
                       mean_slope = -0.15,
                       sd_slope = 0.10,
                       icc_slope = 0,
                       weeks_accrual = 52,
                       shape = c("linear", "loglinear"),
                       p_treated = 0.5,
                       z_level = c("therapist", "patient"),
                       z_type = c("binary", "dose"),
                       tau = 0,
                       tau_x = 0,
                       tau_x_form = c("linear", "step", "quadratic"),
                       n_noise = 0,
                       tau_c = 0,
                       tau_xc = 0,
                       tau_shape = c("constant", "ramp"),
                       x_effect = 0,
                       confounding = 0,
                       dropout = 0,
                       dropout_informative = 0,
                       reliability_score = 1,
                       reliability_x = 1,
                       alliance = FALSE,
                       icc_alliance = 0.10,
                       alliance_slope = 0,
                       coupling = 0,
                       coupling_reverse = 0,
                       alliance_ar = 0.75,
                       reliability_alliance = 1,
                       seed) {
  shape <- match.arg(shape)
  tau_shape <- match.arg(tau_shape)
  tau_x_form <- match.arg(tau_x_form)
  z_level <- match.arg(z_level)
  z_type <- match.arg(z_type)
  stopifnot(n_noise >= 0)
  if (z_type == "dose") {
    if (confounding != 0) {
      stop("confounding != 0 is not implemented for z_type = \"dose\".",
           call. = FALSE)
    }
    if (p_treated != 0.5) {
      warning("p_treated is ignored when z_type = \"dose\" ",
              "(dose is drawn uniform on [0, 1]).", call. = FALSE)
    }
  }
  if (missing(seed)) {
    stop("`seed` is mandatory: a wind-tunnel run must be reproducible.",
         call. = FALSE)
  }
  stopifnot(
    n_therapists >= 1, patients_per_therapist >= 1,
    n_sessions >= 1, weeks_accrual >= 1,
    icc >= 0, icc < 1,
    icc_slope >= 0, icc_slope < 1,
    p_treated >= 0, p_treated <= 1,
    dropout >= 0, dropout < 1,
    reliability_score > 0, reliability_score <= 1,
    reliability_x > 0, reliability_x <= 1,
    icc_alliance >= 0, icc_alliance < 1
  )
  if ((coupling != 0 || coupling_reverse != 0) && !alliance) {
    stop("coupling/coupling_reverse require alliance = TRUE.",
         call. = FALSE)
  }
  if (confounding != 0) {
    if (z_level != "patient") {
      stop("confounding != 0 requires z_level = \"patient\" ",
           "(therapist level not implemented yet).",
           call. = FALSE)
    }
    if (p_treated <= 0 || p_treated >= 1) {
      stop("confounding != 0 requires p_treated strictly in (0, 1).",
           call. = FALSE)
    }
  }
  stopifnot(alliance_ar >= 0, alliance_ar < 1,
            reliability_alliance > 0, reliability_alliance <= 1)
  if (reliability_alliance < 1 && !alliance) {
    stop("reliability_alliance requires alliance = TRUE.", call. = FALSE)
  }
  set.seed(seed)

  n_patients <- n_therapists * patients_per_therapist

  # therapist effects: sd_u from ICC at residual SD = 1
  sd_u <- sqrt(icc / (1 - icc))
  u <- rnorm(n_therapists, mean = 0, sd = sd_u)

  # treatment: cluster level (default) or dyad/patient level
  if (z_level == "therapist") {
    if (z_type == "binary") {
      n_treated <- round(p_treated * n_therapists)
      z <- integer(n_therapists)
      z[sample.int(n_therapists, n_treated)] <- 1L
    } else {
      # dose mode: continuous exposure per therapist, uniform on [0, 1]
      # (new branch; binary draw order stays bit-identical)
      z <- stats::runif(n_therapists)
    }
  }

  therapist_id <- rep(seq_len(n_therapists),
                      each = patients_per_therapist)
  patient_id   <- seq_len(n_patients)

  # with confounding = 0, z is drawn here (identical draw order as
  # always -> bit identity); with confounding != 0 only AFTER x,
  # because assignment depends on x
  if (z_level == "patient" && confounding == 0) {
    z_pat <- if (z_type == "binary") {
      stats::rbinom(n_patients, 1L, p_treated)
    } else {
      # dose mode: continuous exposure per patient, uniform on [0, 1]
      stats::runif(n_patients)
    }
  }
  baseline   <- rnorm(n_patients, mean = 0, sd = 1)
  # slope decomposition: patient share + (optional) therapist share;
  # with icc_slope = 0 bit-identical to the old version (no extra draw)
  slope      <- rnorm(n_patients, mean = mean_slope,
                      sd = sd_slope * sqrt(1 - icc_slope))
  if (icc_slope > 0) {
    v_slope <- rnorm(n_therapists, 0, sd_slope * sqrt(icc_slope))
    slope   <- slope + v_slope[therapist_id]
  }
  entry_week <- sample.int(weeks_accrual, n_patients, replace = TRUE)
  x          <- rnorm(n_patients, mean = 0, sd = 1)
  if (z_level == "patient" && confounding != 0) {
    p_i   <- stats::plogis(stats::qlogis(p_treated) + confounding * x)
    z_pat <- stats::rbinom(n_patients, 1L, p_i)
  }

  # long format: one row per (planned) session
  idx     <- rep(seq_len(n_patients), each = n_sessions)
  session <- rep(seq_len(n_sessions), times = n_patients)

  # therapist attribute c: only drawn when it acts (otherwise the
  # draw order of old seeds would not be preserved)
  c_th <- NULL
  if (tau_c != 0 || tau_xc != 0) {
    c_th <- rnorm(n_therapists, mean = 0, sd = 1)
  }
  # noise features: inert, only a distraction for estimators
  x_noise <- NULL
  if (n_noise > 0) {
    x_noise <- matrix(rnorm(n_patients * n_noise), ncol = n_noise,
                      dimnames = list(NULL,
                                      paste0("x_noise", seq_len(n_noise))))
  }

  # heterogeneity contribution: functional form of x ("linear"
  # reproduces the old arithmetic bit-identically)
  h <- switch(tau_x_form,
              linear    = x[idx],
              step      = as.numeric(x[idx] > 0),
              quadratic = x[idx]^2 - 1)
  te <- tau + tau_x * h        # person-specific z-effect
  if (!is.null(c_th)) {        # + therapist/dyad share
    te <- te + (tau_c + tau_xc * x[idx]) * c_th[therapist_id[idx]]
  }
  if (tau_shape == "ramp") {   # effect grows over the sessions
    te <- te * (session / n_sessions)
  }
  z_of_row <- if (z_level == "therapist") z[therapist_id[idx]] else
    z_pat[idx]
  progress <- if (shape == "linear") session - 1 else log(session)

  score <- baseline[idx] +
    u[therapist_id[idx]] +
    te * z_of_row +
    slope[idx] * progress +
    rnorm(n_patients * n_sessions, mean = 0, sd = 1)
  # prognostic main effect of x; only touched when x_effect != 0 so
  # the default path stays untouched
  if (x_effect != 0) score <- score + x_effect * x[idx]

  # alliance stream (optional): own process + coupling onto symptoms
  alli <- NULL
  if (alliance) {
    sd_ua <- sqrt(icc_alliance / (1 - icc_alliance))
    u_a   <- rnorm(n_therapists, 0, sd_ua)
    a_base <- rnorm(n_patients, 0, 1)
    alli <- numeric(n_patients * n_sessions)
    sd_innov <- sqrt(1 - alliance_ar^2)  # stationary variance ~ 1
    for (pt in seq_len(n_patients)) {
      off <- (pt - 1L) * n_sessions
      mu_prev <- 0
      for (s in seq_len(n_sessions)) {
        prog <- if (shape == "linear") s - 1 else log(s)
        mu <- a_base[pt] + u_a[therapist_id[pt]] +
          alliance_slope * prog
        if (s == 1) {
          alli[off + s] <- mu + rnorm(1, 0, 1)
        } else {
          # symptom update first (uses previous session's alliance) ...
          if (coupling != 0) {
            score[off + s] <- score[off + s] +
              coupling * alli[off + s - 1]
          }
          # ... then alliance update: inertia + reverse coupling
          alli[off + s] <- mu +
            alliance_ar * (alli[off + s - 1] - mu_prev) +
            coupling_reverse * score[off + s - 1] +
            rnorm(1, 0, sd_innov)
        }
        mu_prev <- mu
      }
    }
  }

  out <- data.frame(
    therapist_id = therapist_id[idx],
    z            = z_of_row,
    patient_id   = patient_id[idx],
    x            = x[idx],
    entry_week   = entry_week[idx],
    session      = session,
    week         = entry_week[idx] + session - 1L,
    score        = score
  )
  if (!is.null(c_th)) out$therapist_c <- c_th[therapist_id[idx]]
  if (!is.null(x_noise)) {
    for (j in seq_len(n_noise)) {
      out[[colnames(x_noise)[j]]] <- x_noise[idx, j]
    }
  }
  if (alliance) {
    if (reliability_alliance < 1) {
      sd_ma <- stats::sd(alli) *
        sqrt((1 - reliability_alliance) / reliability_alliance)
      alli <- alli + rnorm(length(alli), 0, sd_ma)
    }
    out$alliance <- alli
  }

  # measurement error (after all structural draws -> old seeds stay
  # unchanged at reliability = 1)
  if (reliability_score < 1) {
    sd_m <- stats::sd(score) *
      sqrt((1 - reliability_score) / reliability_score)
    score <- score + rnorm(length(score), 0, sd_m)
    out$score <- score
  }
  x_obs <- x
  if (reliability_x < 1) {
    x_obs <- x + rnorm(n_patients, 0,
                       sqrt((1 - reliability_x) / reliability_x))
    out$x <- x_obs[idx]
  }

  # dropout: after each session a dropout die roll; informative via score
  if (dropout > 0) {
    p_drop <- stats::plogis(stats::qlogis(dropout) +
                              dropout_informative * out$score)
    evt <- stats::runif(nrow(out)) < p_drop
    # first dropout session per patient (Inf = stays to the end)
    first_evt <- stats::ave(as.numeric(evt), out$patient_id,
                            FUN = function(v) {
                              w <- which(v > 0)
                              rep(if (length(w)) w[1] else Inf,
                                  length(v))
                            })
    out <- out[out$session <= first_evt, , drop = FALSE]
    rownames(out) <- NULL
  }

  class(out) <- c("windkanal_stream", "data.frame")
  out
}
