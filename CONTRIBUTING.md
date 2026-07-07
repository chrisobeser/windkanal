# Contributing to windkanal

Thank you for your interest. The package is young; contributions are welcome.

* **Bugs / questions:** open a GitHub issue with a minimal reproducible
  example (a `sim_stream()` call with its seed is usually enough).
* **New estimators:** implement the `fit_*` interface (take a stream or
  snapshot, return `c(estimate, se)` or a CATE data frame with an `ate`
  attribute; mandatory `seed` argument). One file per estimator; add tests
  mirroring `tests/testthat/test-fit_cate_learners.R`.
* **New presets:** every parameter needs a source (see `attr(preset(...),
  "sources")`). Uncited values are not accepted.
* **Style:** match the existing code; comments state constraints, not
  narration. Run `devtools::test()` before submitting; CI must stay green.
* **Data discipline:** never commit licensed item texts or third-party
  datasets to the repository.
