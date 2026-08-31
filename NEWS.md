# windkanal 0.3.0 (2026-08-31)

Two new capability blocks: a **dyad toolkit** for patient-therapist
matching questions, and **bounded, skewed outcomes** for
questionnaire-realistic score distributions. Both default to off, and
worlds from earlier versions are unchanged: the test suite pins a
reference world against its 0.2.1 values and checks that each new
switch leaves the default path bit-identical.

## Dyad toolkit

- Congruence forms of the dyadic term: `sim_stream(tau_xc_form =)`
  offers `"product"` (the previous arithmetic, unchanged),
  `"misfit_quad"` (`-(x - c)^2`), `"similarity_abs"` (`-|x - c|`)
  and `"fenster"` (matching zone `1{|x - c| <= fenster_delta}`).
  All forms are centred and scaled to SD 1 under independent
  standard-normal margins, so `tau_xc` means the same amplitude
  across forms; `dyade_form()` exposes the mapping directly.
  `fenster_delta` is restricted to (0, 3] -- beyond that the window
  catches almost every pair and the normalization divides by nearly
  zero.
- `dyade_h_wahr()`: the truth reference for these forms. A congruence
  form is not the same thing as its interaction component -- the
  quadratic form is half main effects -- so matching claims must be
  scored against the double-centred `h` component, not against the
  raw form. The function computes it, optionally for the realized
  pairs.
- Formation switch: `sim_stream(assignment = "soft_matching")`
  replaces random allocation with the informal matching of routine
  care -- patient i lands with therapist j with probability
  proportional to `exp(assignment_strength * x_i * c_j)` among
  therapists with free capacity, caseloads exact. This is the
  formation collider and the range restriction on the fit surface
  that any routine-data application has to survive. Requires
  `icc_slope = 0` and `confounding = 0` in this version. A non-zero
  `assignment_strength` under `assignment = "random"` is an error
  rather than a silently ignored argument, and extreme lambda no
  longer overflows `exp()`.
- `fit_cate_dyade_v2p()`: dyadic BCF with a draw-wise decomposition.
  The fit is the one from `fit_cate_bcf_dyade()`; each posterior draw
  of the tau surface is then double-centred on a reference grid into
  `tau0 + f(x) + g(c) + h(x, c)`. No refit, no residuals -- the
  decomposition is a deterministic projection, so uncertainty
  propagates exactly and `tau_hat` stays identical to
  `fit_cate_bcf_dyade()` at the same seed. Detection is
  amplitude-based (`sd_h` draws); directional and correlation
  measures are useless as detectors here, since the forest uses the
  product feature for noise as well.
- Candidate naming on top of it: `v2p_benennung()` regresses the
  h-draws on a named library of candidate forms (or attribute
  products), `v2p_entscheidung()` applies the two-part decision rule.
  The default thresholds are calibrated (detection .07/.52, naming
  .02/.77) on static product-form worlds; outside that regime they
  are a transfer, not a guarantee, and the layer stays exploratory.
  Candidates whose h components coincide are exactly collinear after
  double-centring (`dz(-(x-c)^2) = 2*dz(x*c)` on every grid) and are
  refused with a named error instead of silent `NA` coefficients.

## Bounded and skewed outcomes

- `skala_begrenzen()` maps the latent `score` onto a bounded integer
  questionnaire scale (BDI/PHQ style): optional monotone skew tilt,
  affine anchoring, censoring at floor and ceiling, rounding to whole
  points. The latent column is kept; the layer is opt-in and adds
  `score_obs` plus a `skala` attribute. `skala_werte()` exposes the
  pure mapping.
- `wahrheit_skala()` computes the true treatment effect on **both**
  scales via a counterfactual twin run, so a bounded world can be
  audited without losing the latent truth. The floor/ceiling of a
  real questionnaire manufactures apparent effect heterogeneity where
  the latent effect is homogeneous -- the function reports it as
  `schein_heterogenitaet`, unit-corrected, and documents why the raw
  difference of the two ranges is not a measure of anything.
- The twin run rests on `sim_stream(z_force = )`: the assignment is
  drawn exactly as usual (identical RNG consumption) and only then
  overwritten, so both runs are common-random-number aligned per row.
  `wahrheit_skala()` guards the alignment four ways, including an
  RNG end-state check.

## Fixes

- `fit_het_grf()` now sees the same feature set as `fit_cate_grf()`
  (`features = "alle"`, the new default). Until 0.2.1 it received
  only `x`, so purely dyadic heterogeneity (`tau_xc != 0` at
  `tau_x = 0`) lay outside the feature set and was reported as a
  clean null -- in ten worlds with true `SD(tau_i) = 1.53`, zero
  hits, while the same strength along `x` was found in five of five.
  `features = "nur_x"` keeps the old, narrow variant available as a
  declared naive reference. **This changes results of earlier runs
  of `fit_het_grf()` on worlds carrying `therapist_c` or noise
  features.**
- The seed rule is now enforced uniformly at **all twenty** entry
  points that take a mandatory `seed` -- `sim_stream()`,
  `wahrheit_skala()`, `sim_items()`, `plasmode_world()` and the
  sixteen `fit_*` wrappers. `seed = NULL` passes `missing()`
  unnoticed and would send `set.seed()` to the clock; measured, the
  same call then ran through twice with different results in
  `fit_cate_tlearner()`, `fit_cate_xlearner()`,
  `fit_cate_drlearner()` and `fit_cate_bcf()`. It is now an error
  with an explicit message, as are vectors, `NA` and `Inf`. Until
  0.2.1 the rule held at four of these twenty. (`fit_het_grf()` is
  the one stochastic entry point without a mandatory `seed`: its
  internal `forest_seed` is defaulted and the fit is deterministic
  by construction.)
- Argument guards against silent vector recycling and fractional
  counts: `fenster_delta` must be a single finite number (a vector
  passed the element-wise bound and then recycled into a world with
  alternating window widths), `gitter_n` in
  `fit_cate_dyade_v2p()` must be a single whole number >= 2 (a
  fractional value produced a skewed reference grid), and
  `n_therapists`, `patients_per_therapist`, `n_sessions`,
  `weeks_accrual` must be single whole numbers. The last one is a
  long-standing generator bug: `n_therapists = 3.7` used to run and
  deliver rows carrying `NA` in `therapist_id`, `z` and `score`
  without any warning -- and it silently broke the "caseloads exact"
  promise of `assignment = "soft_matching"`.
- The five arguments new in this version (`tau_xc_form`,
  `fenster_delta`, `assignment`, `assignment_strength`, `z_force`)
  are appended **after** `seed` in the signature of `sim_stream()`,
  so arguments 1 to 33 keep the positions they had in 0.2.1 and
  positional calls written against earlier versions keep their
  meaning.

## Docs

- Positioning made explicit: windkanal is a stochastic (Monte Carlo)
  simulation testbed with controlled randomness -- now stated in the
  DESCRIPTION title and description, the README lead, and
  CITATION.cff (new keyword `stochastic-simulation`).

# windkanal 0.2.1 (2026-07-16)

- New vignette "Dose worlds: when the amount matters more than the
  assignment": walkthrough for the continuous treatment axis shipped
  in 0.2.0 (`z_type = "dose"`), including the nesting trap on the
  dose slope and the guard rails of the arm-based CATE wrappers.
- Documentation release: no changes to package code. Citation and
  README version strings brought up to date.

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
