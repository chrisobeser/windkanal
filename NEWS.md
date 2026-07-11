# windkanal 0.2.0 (2026-07-11)

- Continuous treatment axis: `sim_stream(z_type = "dose")` draws a
  continuous exposure `z` uniform on [0, 1] at the chosen `z_level`
  (therapist or patient/dyad). All effect terms (`tau`, `tau_x`,
  `tau_c`, `tau_xc`) act per unit dose, so `tau` is the effect of
  moving from dose 0 to dose 1. Dose worlds model questions where
  amount matters more than assignment: session dose, alliance
  exposure, degree of patient-therapist matching. Binary worlds
  (`z_type = "binary"`, the default) remain bit-identical to
  earlier versions.
- Guard rails: the arm- and propensity-based estimator wrappers
  (all `fit_cate_*` functions and `fit_het_grf()`) refuse dose
  worlds with an informative error instead of silently
  reinterpreting a continuous `z`. The classical mixed-model
  wrappers (`fit_z_lmm()`, `fit_z_satt()`, `fit_zx_*()`) carry the
  per-unit dose effect; continuous-treatment CATE estimators
  remain on the roadmap.
- In dose mode `p_treated` is ignored (with a warning) and
  `confounding` is not implemented (error); both are documented in
  `?sim_stream`.

# windkanal 0.1.0 (2026-07-10)

- Four new estimator wrappers behind the same uniform interface:
  multilevel BART (`fit_cate_stan4bart`; stan4bart), GP boosting
  (`fit_cate_gpboost`; GPBoost), mixed-effects random forest
  (`fit_cate_merf`; LongituRF), and a dyadic BCF variant with an
  explicit patient-by-therapist product feature (`fit_cate_bcf_dyade`).
- `plasmode_world()`: build worlds from a real covariate table and a
  real cluster structure while keeping the injected treatment effect
  known (generate-treatment framework; ships as recipe code, no data).
- `inst/CITATION` added; `citation("windkanal")` now works.

Shipped in the same release (previously listed under "development
version"):

* New `cate_metrics()`: one shared metric set for CATE estimates against
  known truth (ranking `r`, magnitude `pehe`, bias, dispersion, person-level
  interval coverage, honest error flag). PEHE is always reported: `r` alone
  is blind to amplitude errors (any affine estimator reaches r = 1 under a
  linear true effect).

* New `fit_z_brms()`: Bayesian mixed-model ATE estimator (brms/Stan),
  the Bayesian counterpart to `fit_z_satt()` behind the same interface,
  with a per-session compile cache. Requires the suggested `brms`
  package.

* Added `CITATION.cff` (GitHub citation support).
- First full validation program complete: sixteen estimators across
  fifteen pre-specified cells, 62,100 estimator-world fits, all
  expectations version-controlled before execution. Preprint in
  preparation.

# windkanal 0.1.0 (2026-07-06)

First feature-complete pre-release.

## Simulation engine
* `sim_stream()`: patients nested in therapists, staggered entry, session-wise
  outcomes; treatment at therapist or patient (dyad) level; effect heterogeneity
  with linear, step, or quadratic moderator shapes (`tau_x_form`); therapist
  moderation (`tau_c`) and dyadic matching effects (`tau_xc`); effect ramps
  (`tau_shape`); prognostic covariate effects and confounding by indication;
  inert noise features; informative dropout; measurement error on outcome and
  moderator; optional alliance process with bidirectional coupling.
* Calibrated presets with per-parameter source registers (`preset()`):
  `"ambulanz_de"`, `"allianz_beierl2021"`.
* Time machine: `snapshot()`, `replay()`; analysis regimes `run_peek()`,
  `run_gate()`, `run_gates()`.

## Estimators (uniform `fit_*` interface, shared feature set)
* Classical: naive OLS, mixed models with Satterthwaite inference (main effect,
  interaction, slope).
* CATE: PAI (per-arm linear), T-/S-/X-/DR-/R-learners, model-based forest
  (model4you), causal forest (grf; honesty and cluster switches), legacy BCF
  (bcf), multilevel BCF with therapist random intercepts (stochtree),
  therapist-cluster bootstrap for learner inference.

## Infrastructure
* Item layer: `read_items()` (formr structure only), `sim_items()`,
  `scale_scores()`.
* Monte Carlo driver (`mc_run()`/`mc_summary()`), experiment scripts with
  fixed expectations, per-world checkpointing, and per-row provenance
  (commit hash, package versions).
* `experiments/reproduce.R`: one-command re-execution and byte-level
  verification against committed results; CI via GitHub Actions;
  independent reimplementation check (`verify_independent.R`).
