# Fit Model Using EM Algorithm

Function that uses the EM Algorithm to fit a model to output from either
import_mics_with_metadata() or simulate_mics(). The algorithm fits a
regression model to the mean (mu) of each component in the data,
assuming a mixture of gaussians on the log2 scale, and a generalized
additive model to the component weights (pi). The function includes an
option to avoid estimating mu for a component that is entirely or nearly
entirely outside the range of tested concentrations by using a reduced
approach and selecting which side should be fixed outside the tested
range. Mu model fits are done on the log2(MIC) scale. Output is a list
containing the data and observation weights, fitted models, and many of
the settings used to run the function.

## Usage

``` r
fit_EM(
  model = "pspline",
  approach = "full",
  pre_set_degrees = NULL,
  max_degree = 8,
  degree_sets = "matched",
  visible_data,
  nfolds = 10,
  non_linear_term = "t",
  covariates = NULL,
  pi_formula = c == "2" ~ s(t),
  fixed_side = NULL,
  extra_row = FALSE,
  ecoff = NA,
  max_it = 3000,
  ncomp = 2,
  tol_ll = 1e-06,
  pi_link = "logit",
  verbose = 3,
  model_coefficient_tolerance = 1e-05,
  maxiter_survreg = 30,
  initial_weighting = 3,
  sd_initial = 0.2,
  scale = NULL,
  reruns_allowed = 3,
  max_out_break = FALSE
)
```

## Arguments

- model:

  String, "pspline" or "polynomial". Which non-linear term should be
  used in model

- approach:

  String, either "full" or "reduced". Full model attempts to fit mu
  models to each component, reduced models does not attempt to fit a mu
  model for one component. Suitable if one component is entirely or
  nearly entirely outside the range of tested concentrations

- pre_set_degrees:

  Vector of numeric values with length equal to number of components
  where mu is being estimated. NULL if cross validation is being used to
  select degrees (polynomial) or degrees of freedom (pspline) for mu
  model(s).

- max_degree:

  Numeric, must be an integer \> 0. Highest value to test in cross
  validation for degrees (polynomial) or degrees of freedom (pspline)
  for mu model(s).

- degree_sets:

  String: "independent" or "matched", if using cross validation to
  select appropriate degrees (polynomial) or degrees of freedom
  (pspline) for mu model(s), whether to test a grid of all possible
  combinations("independent") or fix the degrees of the components to
  match ("matched")

- visible_data:

  Data frame, data including the left and right bound of the MICs (use
  import_mics_with_metadata to format correctly) and any covariates
  (including the non-linear term)

- nfolds:

  Numeric, number of folds for the cross validation to select
  appropriate degrees (polynomial) or degrees of freedom (pspline) for
  mu model(s)

- non_linear_term:

  String, non-linear term to be included in the model. Variable in the
  pspline term in the pspline model or in the polynomial term in the
  polynomial model.

- covariates:

  String, covariates to be included in mu model aside from the
  non-linear term.

- pi_formula:

  Formula, for the component weight model. Model is fit using mgcv's gam
  function. Nonlinear terms include s() and lo(). Basis of s() function
  can be changed using bs argument to s(). Use c == "2" for the left
  side of the formula.

- fixed_side:

  String, if using a reduced model, specify which component the
  algorithm won't estimate mu for. "RC" corresponds to the upper
  component, "LC" corresponds to the lower. NULL if using the full
  model.

- extra_row:

  Logical, if using a reduced model, the highest ("RC") or lowest ("LC")
  MIC value may be included in the set observations weighted as possibly
  in the fixed component

- ecoff:

  String or numeric, represents the largest MIC value of the WT
  distribution on MIC scale. Model will use this in combination with the
  fixed side to assume the observations on the non-fixed side of the
  ecoff are not part of the component on the fixed side of the ecoff.
  Can be a number on the MIC scale: 32 or a string representation of WT
  classification "\<=32", should be inclusive

- max_it:

  Numeric, maximum number of iterations for the EM algorithm for any
  given model fitting.

- ncomp:

  Numeric, number of components to be fitted. When fitting a reduced
  model where one component is not estimated, that component should
  still contribute to the value in ncomp. E.g. a reduced model where the
  upper component is fixed and mu for the lower component is being
  estimated has a value of ncomp = 2.

- tol_ll:

  Numeric, maximum tolerance for change in likelihood between steps of
  the algorithm for model convergence to be achieved.

- pi_link:

  String: "logit" or "identity", link function for the generalized
  linear model fit for the pi model (component weights).

- verbose:

  Numeric, controls amount of information printed during model fitting

- model_coefficient_tolerance:

  Numeric, maximum tolerance for change in model coefficients (insluding
  spline terms) between steps of the algorithm for model convergence to
  be achieved.

- maxiter_survreg:

  Maximum iterations used in survreg model fitting, default is 30.

- initial_weighting:

  Numeric, For the "full" model fitting: if 1: initial observation
  weights are estimated using linear regression at the highest and
  lowest tested concentrations. If 2, used a randomized start suitable
  for simulation studies on model validity but otherwise not
  recommended. If 3 or greater, fits a linear models to the components
  and estimates intial weights based on this model fit. For the reduced
  model fitting: 1 sets initial weights corresponding to fixed side (and
  extra_row) where observations not outside the range on the side
  corresponding to the fixed side are forced to be in the component
  where mu is being estimated. For initial weighting two a linear model
  is fit for the component still being estimated to provide initial
  observation weights.

- sd_initial:

  Numeric, value greater than 0 and less than 1. Proportion of the range
  from the highest concentration to lowest concentration that is used as
  the initial estimate of sigma for the estimated components. Default is
  0.2

- reruns_allowed:

  Numeric, if the cross-validation for a particular combination of
  degrees (polynomial) or degrees of freedom (pspline) fails, how many
  repeat attempts should be allowed?

- max_out_break:

  Logical, if TRUE when a CV fold reaches maximum iterations it breaks
  the loop and moves to next rerun if applicable

## Examples

``` r
# \donttest{
data = simulate_mics()
result = fit_EM(model = "pspline",
approach = "full",
pre_set_degrees = c(4,4),
visible_data = data,
non_linear_term = "t",
covariates = NULL,
pi_formula = c == "2" ~ s(t),
max_it = 300,
ncomp = 2,
tol_ll = 1e-6,
pi_link = "logit",
verbose = 1,
model_coefficient_tolerance = 0.00001,
initial_weighting = 3,
sd_initial = 0.2
)
#> Stopped on combined LL and parameters
# }
```
