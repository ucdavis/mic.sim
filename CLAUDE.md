# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this package does

`mic.sim` fits a latent Gaussian mixture model to minimum inhibitory concentration (MIC)
data on the log2 scale. MIC readings are interval-censored (a plate reports `<=2`, `4`,
`>8`), so every model fit is an interval-censored regression. The mixture separates a
wild-type from a non-wild-type subpopulation, and the point of the package is to model
*trends over time* in both the component means (mu) and the mixture weights (pi).

## Development commands

```r
devtools::load_all()      # attach the package for interactive work
devtools::document()      # regenerate man/ and NAMESPACE from roxygen (required after any doc change)
devtools::test()          # run the test suite
devtools::check()         # full R CMD check, as CI runs it
devtools::install()
pkgdown::build_site()     # site is published to https://ajmichaelucd.github.io/mic.sim/
```

Run a single test file:

```r
devtools::load_all(); testthat::test_file("tests/testthat/test-import_mics.R")
```

CI runs `R-CMD-check` on five OS/R combinations with `error-on: '"error"'`, plus
pkgdown, a link check, and `check-non-standard-chars`.

## Architecture

### The pipeline

```
import_mics_with_metadata()  |  simulate_mics()      ->  visible_data
                                        fit_EM()     ->  output (a list)
   plot_fm() / plot_likelihood() / plot_mean() / plot_pi()
```

`fit_EM()` is a dispatcher, not the algorithm: it optionally runs cross-validation to
pick smoothness, builds the mu formulas via `write_all_formulas()`, and then calls
`EM_algorithm()` (full) or `EM_algorithm_reduced()` (reduced).

### Attributes are the data contract

Data frames carry `scale` (`"log"` or `"MIC"`), `source`, `mic_class`, `metadata`,
`lr_col`, and later `model`. Downstream functions read these attributes rather than
taking equivalent arguments (`add_scale()`, `add_attribute_data()`,
`check_scale_and_var()`, every backend branch). A transformation that drops attributes,
or a hand-built tibble, fails far from the edit site with a confusing error. The
`import_mics()` test asserts the full attribute set for this reason.

### `possible_data` is long: one row per observation x component

`visible_data` is expanded with `reframe(.by = everything(), c = as.character(1:2))`
plus `add_obs_id()`. The probability columns are backtick-named and mirror the EM
derivation: `` `P(C=c|t)` ``, `` `P(Y|t,c)` ``, `` `P(c,y|t)` ``, `` `P(Y=y|t)` ``, and
`` `P(C=c|y,t)` `` -- the last is the observation weight passed to every weighted fit.
The names contain `|`, `(`, and `=`, so they must always be backticked.

### Two model families per iteration

- **mu** (component means): one interval-censored regression per estimated component.
  Backend `"surv"` uses `survival::survreg()` with `pspline()`; `"polynomial"` uses
  `survreg()` with a polynomial term; `"mgcv"` uses `mgcv::gam(family = cnorm())`.
- **pi** (mixture weights): `mgcv::gam()` on `c == "2"`, binomial with a `logit` or
  `identity` link, weighted by `` `P(C=c|y,t)` ``.

Convergence is checked on both the log-likelihood delta (`tol_ll`) and the model
coefficients (`model_coefficient_tolerance`); `likelihood_documentation` is a matrix
with fixed `dimnames` that accumulates per-iteration diagnostics and is converted to a
tibble in the output.

### Backend dispatch is by attribute, not S3

Functions named `fit_mu_model.mgcv`, `fit_single_component_model.surv`,
`modify_bounds.mgcv` *look* like S3 methods, but there is no `UseMethod` anywhere in
the package. The bare wrapper branches on `attr(x, "model")` with `if`/`else`. Adding
or changing a backend means editing every branching wrapper
(`fit_single_component_model()`, `fit_all_mu_models()`, `modify_bounds()`,
`add_attribute_data()`, and the checks inside `EM_algorithm()`), not registering a
method.

Note the name remap: the user-facing `model = "pspline"` is rewritten to the internal
`"surv"` at the top of `fit_EM()`. Internal values are `"surv"`, `"polynomial"`,
`"mgcv"`.

### Full vs. reduced approach

`approach = "reduced"` exists for data where one component lies (nearly) entirely
outside the range of tested concentrations, so its mean is not estimable. It skips the
mu model for the side named by `fixed_side` (`"RC"` upper, `"LC"` lower) and uses
`ecoff` and `extra_row` to decide which observations may belong to that fixed
component. `ncomp` still counts the fixed component.

### Cross-validation

When `pre_set_degrees` is `NULL`: `full_cv()` -> `single_cv_all()` ->
`get_fold_likelihood_all()` -> `calculate_fold_likelihood_all()`. Degree combinations
come from `create_degree_combinations_all()` (`degree_sets` is `"matched"` or
`"independent"`), and the winner is the set with the highest summed held-out
log-likelihood. Folds that fail are retried up to `reruns_allowed`.

### The `fit_EM()` output list

`likelihood`, `possible_data`, `mu_model`, `pi_model`, `steps`, `converge`, `ncomp`,
`prior_step_models`, `cv_results`, `mu_formula`, `model`. Every plotting function takes
this list as `output`. `plot_fm(use_prior_step = TRUE)` falls back to
`prior_step_models` when a mu model failed to converge on the final iteration -- that is
why the previous iteration is retained in the output at all.

## Conventions and gotchas

- **`man/` and `NAMESPACE` are generated.** Edit the roxygen block and run
  `devtools::document()`. A hand-edited `.Rd` drifting from its roxygen source has
  already caused one broken example here.
- **Only 15 of ~150 top-level functions are exported.** Most files in `R/` are internal helpers
  with no `@export`. `import_mics()` is internal; `import_mics_with_metadata()` is the
  public entry point (tests call the internal one directly, which works because
  testthat loads the namespace).
- **Roxygen examples are a recurring source of R CMD check failures.** Examples run
  with only `mic.sim` attached, so qualify calls into Imports packages
  (`tibble::tibble()`, `dplyr::rename()`). Wrap examples for non-exported functions in
  `\dontrun{}`, and slow or printing-fragile fits (`fit_EM`) in `\donttest{}`.
- **`≤` (U+2264) in MIC strings is data syntax, not a typo.** `import_mics()` parses
  `≤`, `<=`, `=<`, and `>`. The `check-non-standard-chars` CI job only bans curly quotes
  and en/em dashes, so leave `≤` alone.
- **magrittr pipes** (`%>%`, `%<>%`) throughout, not the native pipe; the `.Rproj` sets
  `UseNativePipeOperator: No`.
- **`misc/` and `data-raw/` are tracked scratch work**, both `.Rbuildignore`d. They hold
  abandoned drafts and exploration scripts, not package code -- do not treat them as
  reference implementations. (`R/scratch_2C_mean_plots.R` is scratch that ended up in
  `R/`.)
