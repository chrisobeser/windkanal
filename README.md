# windkanal <img src="man/figures/logo.png" align="right" width="130" />

**A wind tunnel for psychotherapy statistics: simulated clinics with
known ground truth, so you can test your tools before they meet real
patients.**

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

**Start from calibrated settings.** `preset("ambulanz_de")` encodes a
German outpatient training clinic. Every number in a preset carries a
citation and a verification status, and the register travels with the
object: `attr(preset("ambulanz_de"), "sources")`.

**Let estimators compete.** Twelve estimator classes run behind one
interface and return comparable output: naive OLS, mixed models with
Satterthwaite inference, per-arm regressions (PAI), T-, S-, X-, DR-
and R-learners, model-based forests, causal forests (grf), and
Bayesian causal forests (stochtree), with and without therapist
random effects. Average-effect inference for the learners uses a
therapist-cluster bootstrap.

**Judge them fairly.** `mc_run()` repeats a question over seeded
worlds and `mc_summary()` turns the answers into error rates: bias,
interval coverage, false alarms, ranking quality, and PEHE.

**Test analysis habits, not just estimators.** `snapshot()` and
`replay()` freeze and replay a growing data stream, so that
continuous peeking can be compared with preregistered release gates
(`run_peek()`, `run_gates()`) on identical histories.

**Go down to the item level.** `read_items()` ingests formr
questionnaire tables while keeping only their structure; item texts,
which are often copyrighted, never enter the package. `sim_items()`
produces ordinal responses from latent states, including
careless-responding behavior, and `scale_scores()` returns scale
scores with reliability as a measured attribute rather than an
assumption.

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

Version 0.1.0. More than 220 unit tests, continuous integration on
GitHub Actions, `R CMD check` clean. The data-generating equation has
been independently reimplemented from its written specification and
validated against the package. Validation studies built on windkanal
are in preparation; the package itself is stable enough for
methodological experiments and teaching.

## License

MIT. Contributions are welcome, see `CONTRIBUTING.md`.
