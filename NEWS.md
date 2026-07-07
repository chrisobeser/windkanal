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
