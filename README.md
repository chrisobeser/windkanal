# windkanal <img src="man/figures/logo.png" align="right" width="130" />

[![R-CMD-check](https://github.com/chrisobeser/windkanal/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/chrisobeser/windkanal/actions)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Status: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](#status)

**A wind tunnel for psychotherapy statistics: simulated clinics in
which the true effects are known, so you can test your tools before
they meet real data.**

Psychotherapy data have a shape of their own. Patients are nested in
therapists, outcomes arrive session by session, people drop out, and
questionnaires are filled in by tired humans. Statistical tools that
work beautifully elsewhere can fail quietly under exactly these
conditions, and on real data nobody notices, because the truth is
unknown. windkanal builds artificial outpatient clinics in which the
truth is known by construction. Whatever your estimator claims can
then be checked against what was actually put into the world.

## What you can do with it

**Build a clinic.** `sim_stream()` generates one row per attended
session from an explicit data-generating equation. You control the
therapist share of outcome variance (`icc`, `icc_slope`), treatment
assignment at therapist or patient level, average and person-specific
effects with linear, threshold, or quadratic shapes, confounding by
indication, informative dropout, measurement error, and a
session-wise therapeutic alliance that is coupled with symptoms in
both directions.

**Judge estimators fairly.** All estimators see the same worlds, the
same features, and the same seeds. `mc_run()` repeats a question over
many seeded worlds; `mc_summary()` turns the answers into error
rates: bias, interval coverage, false alarms, ranking quality, and
PEHE.

**Test analysis habits, not just estimators.** `snapshot()` and
`replay()` freeze and replay a growing data stream, so continuous
peeking can be compared with preregistered release gates
(`run_peek()`, `run_gates()`) on identical histories.

**Go down to the item level.** `read_items()` ingests formr
questionnaire tables while keeping only their structure (item texts,
which are often copyrighted, never enter the package). `sim_items()`
produces ordinal responses from latent states, including
careless-responding behavior; `scale_scores()` returns scale scores
with reliability as a measured attribute rather than an assumption.

## The estimator zoo

Seventeen estimator wrappers run behind one uniform interface:

- **Naive OLS** (`fit_z_naive`) — treats every session as
  independent; the field's historical default and the built-in
  negative anchor
- **Mixed model with Satterthwaite inference** (`fit_z_satt`) — the
  classical correct tool for average effects under nesting
- **Bayesian mixed model** (`fit_z_brms`; brms/Stan) — the Bayesian
  counterpart to the Satterthwaite model: same nesting structure,
  posterior mean and credible interval, per-session compile cache
- **Interaction variants** (`fit_zx_naive`, `fit_zx_satt`) — target
  the moderation coefficient directly
- **PAI, per-arm regressions** (`fit_cate_pai`) — the classic
  personalized-advantage approach
- **T-learner** (`fit_cate_tlearner`) — two separate random forests,
  one per treatment arm
- **S-learner** (`fit_cate_sboost`) — one gradient-boosting model
  with treatment as a feature
- **X-learner** (`fit_cate_xlearner`) — imputed individual effects
  with a propensity-weighted blend
- **DR-learner** (`fit_cate_drlearner`) — doubly robust
  pseudo-outcomes, nuisances out-of-bag
- **R-learner** (`fit_cate_rlearner`) — orthogonalized residuals in
  the Nie-Wager style
- **Model-based forest** (`fit_cate_mob`, model4you) — splits where
  the treatment coefficient changes
- **Causal forest** (`fit_het_grf`, `fit_cate_grf`; grf) — honesty
  and cluster-robust inference switches included
- **Bayesian causal forests** (`fit_cate_bcf`, `fit_cate_bcf_ml`;
  bcf and stochtree) — with and without therapist random
  intercepts, so the effect of modeling the nesting is itself
  testable
- **Multilevel BART** (`fit_cate_stan4bart`; stan4bart) — one BART
  response surface plus Stan-sampled therapist random intercepts;
  an independently developed second multilevel implementation
- **GP boosting** (`fit_cate_gpboost`; GPBoost) — tree boosting
  with jointly estimated grouped random effects
- **Mixed-effects random forest** (`fit_cate_merf`; LongituRF) —
  EM-style alternation between forest fit and BLUP update
- **Dyadic BCF** (`fit_cate_bcf_dyade`) — multilevel BCF with an
  explicit patient-by-therapist product feature, for
  matching-effect experiments

Average-effect inference for the learners uses a therapist-cluster
bootstrap; performance measures include dual coverage definitions,
rejection indicators that honor degrees-of-freedom corrections,
ranking correlation, and PEHE.

Judging an estimate against known truth is one call:

``` r
cate_metrics(fit$tau_hat, truth, fit$tau_lo, fit$tau_hi)
#>     n     r  pehe   bias sd_tau_hat covered error
#> 1 200 0.904 0.101 0.0125      0.234       1 FALSE
```

`r` measures ranking only and is blind to amplitude errors (any
affine estimator reaches r = 1 under a linear true effect); `pehe`
is the magnitude corrective. The function always reports both, so
no experiment can silently drop one.

## Calibrated presets, including your own

A preset is a set of `sim_stream()` arguments in which every number
carries a citation and a verification status. The register travels
with the object:

```r
p <- preset("ambulanz_de")   # German outpatient training clinic
attr(p, "sources")           # every value, its source, its status
```

You can, and are encouraged to, define presets for **your own
practice, clinic, hospital, or country**: a preset is simply a named
list of `sim_stream()` arguments with a `sources` data frame attached
as an attribute. The contract is the design rule, not the object
class: values that claim realism carry a source, and values you could
not calibrate are labelled open instead of silently guessed.

```r
meine_klinik <- list(
  n_sessions = 12,      # your typical treatment length
  icc        = 0.08,    # your therapist share of outcome variance
  dropout    = 0.025    # your per-session discontinuation rate
)
attr(meine_klinik, "sources") <- data.frame(
  parameter = c("n_sessions", "icc", "dropout"),
  quelle    = c("clinic records 2024", "own MLM estimate", "open"),
  status    = c("belegt", "belegt", "offen")
)

s <- do.call(sim_stream, c(meine_klinik,
       list(n_therapists = 15, patients_per_therapist = 8,
            tau = 0.3, seed = 1)))
```

Presets are combinable, and sample sizes deliberately stay out of
them: how large your simulated clinic is remains a design choice, not
a reality claim.

## Quick start

```r
# install.packages("devtools")
devtools::install_github("chrisobeser/windkanal")
library(windkanal)

s <- sim_stream(n_therapists = 20, patients_per_therapist = 10,
                n_sessions = 4, z_level = "patient",
                tau = 0.5, tau_x = 0.5, icc = 0.10, seed = 1)

fit_z_naive(s)  # pretends sessions are independent
fit_z_satt(s)   # knows about therapists and patients
```

Run both and compare the standard errors: that difference is the
package's founding observation.

## Design principles

- `seed` is a required argument everywhere. Reproducibility is not
  optional.
- Every preset parameter carries a source; values that could not be
  calibrated are labelled open instead of silently guessed.
- New features default to off and leave old worlds bit-identical,
  which the test suite enforces across releases.
- Licensed questionnaire texts never enter the repository.

## Status

**Work in progress.** windkanal is an early pre-release (version
0.1.0) under active development: interfaces may still change, and the
package is currently in its testing phase. That said, development is
test-driven from the start: more than 220 unit tests, continuous
integration on GitHub Actions, `R CMD check` clean, and the
data-generating equation has been independently reimplemented from its
written specification and validated against the package. The first validation study is
complete — sixteen estimators across fifteen pre-specified cells,
62,100 estimator-world fits, every expectation version-controlled
before execution — and its preprint is in preparation. Suitable
today for
methodological experiments and teaching; not yet for unsupervised
production use.

## Roadmap

Planned, in rough order. Suggestions and use cases are welcome via
the issue tracker.

- [x] **Bayesian mixed model wrapper** (`fit_z_brms()`): the Bayesian
      counterpart to the Satterthwaite mixed model, behind the same
      uniform interface (shipped)
- [x] **More estimator wrappers** (shipped: stan4bart, GPBoost,
      MERF — see the estimator zoo)
- [x] **Plasmode mode** (shipped: `plasmode_world()` — real
      covariate tables and cluster structures, injected truth known)
- [ ] **Continuous treatment**: dose as a treatment axis, for
      questions where amount matters more than assignment
- [ ] **Questionnaire layer, next stage**: per-item loadings and
      IRT-style difficulties, additional careless-responding styles
      beyond straightlining, item-level missingness with prorating,
      and response shift (participants recalibrating their answer
      scale over the course of therapy)
- [ ] **Bounded and skewed outcomes**: floor and ceiling effects for
      questionnaire-realistic score distributions
- [ ] **Community presets**: calibrated settings for other services
      and countries, contributed with mandatory sources
- [ ] **Shiny design explorer**: interactive what-if for clustered
      longitudinal designs
- [ ] **CRAN submission** once the interfaces have settled

## Citation

``` r
citation("windkanal")
```

Obeser, C. (2026). *windkanal: A simulation testbed for
treatment-effect estimators under therapist nesting* (R package
version 0.1.0). https://github.com/chrisobeser/windkanal

## License

MIT. Contributions are welcome, see `CONTRIBUTING.md`.
