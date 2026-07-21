#' Fit Model Using EM Algorithm
#'
#' Function that uses the EM Algorithm to fit a model to output from either import_mics_with_metadata() or simulate_mics().
#' The algorithm fits a regression model to the mean (mu) of each component in the data, assuming a mixture of gaussians on the log2 scale,
#' and a generalized additive model to the component weights (pi). The function includes an option to avoid estimating mu for a component that is entirely or nearly
#' entirely outside the range of tested concentrations by using a reduced approach and selecting which side should be fixed outside the tested range.
#' Mu model fits are done on the log2(MIC) scale. Output is a list containing the data and observation weights, fitted models, and many of the settings used to run the function.
#'
#'
#' @param model String, "pspline" or "polynomial". Which non-linear term should be used in model
#' @param approach String, either "full" or "reduced". Full model attempts to fit mu models to each component, reduced models does not attempt to fit a mu model for one component. Suitable if one component is entirely or nearly entirely outside the range of tested concentrations
#' @param pre_set_degrees Vector of numeric values with length equal to number of components where mu is being estimated. NULL if cross validation is being used to select degrees (polynomial) or degrees of freedom (pspline) for mu model(s).
#' @param max_degree Numeric, must be an integer > 0. Highest value to test in cross validation for degrees (polynomial) or degrees of freedom (pspline) for mu model(s).
#' @param degree_sets String: "independent" or "matched", if using cross validation to select appropriate degrees (polynomial) or degrees of freedom (pspline) for mu model(s), whether to test a grid of all possible combinations("independent") or fix the degrees of the components to match ("matched")
#' @param visible_data Data frame, data including the left and right bound of the MICs (use import_mics_with_metadata to format correctly) and any covariates (including the non-linear term)
#' @param nfolds Numeric, number of folds for the cross validation to select appropriate degrees (polynomial) or degrees of freedom (pspline) for mu model(s)
#' @param non_linear_term String, non-linear term to be included in the model. Variable in the pspline term in the pspline model or in the polynomial term in the polynomial model.
#' @param covariates String, covariates to be included in mu model aside from the non-linear term.
#' @param pi_formula Formula, for the component weight model. Model is fit using mgcv's gam function. Nonlinear terms include s() and lo(). Basis of s() function can be changed using bs argument to s(). Use c == "2" for the left side of the formula.
#' @param fixed_side String, if using a reduced model, specify which component the algorithm won't estimate mu for. "RC" corresponds to the upper component, "LC" corresponds to the lower. NULL if using the full model.
#' @param extra_row Logical, if using a reduced model, the highest ("RC") or lowest ("LC") MIC value may be included in the set observations weighted as possibly in the fixed component
#' @param ecoff String or numeric, represents the largest MIC value of the WT distribution on MIC scale. Model will use this in combination with the fixed side to assume the observations on the non-fixed side of the ecoff are not part of the component on the fixed side of the ecoff. Can be a number on the MIC scale: 32 or a string representation of WT classification "<=32", should be inclusive
#' @param max_it Numeric, maximum number of iterations for the EM algorithm for any given model fitting.
#' @param ncomp Numeric, number of components to be fitted. When fitting a reduced model where one component is not estimated, that component should still contribute to the value in ncomp. E.g. a reduced model where the upper component is fixed and mu for the lower component is being estimated has a value of ncomp = 2.
#' @param tol_ll Numeric, maximum tolerance for change in likelihood between steps of the algorithm for model convergence to be achieved.
#' @param pi_link String: "logit" or "identity", link function for the generalized linear model fit for the pi model (component weights).
#' @param verbose Numeric, controls amount of information printed during model fitting
#' @param model_coefficient_tolerance Numeric, maximum tolerance for change in model coefficients (insluding spline terms) between steps of the algorithm for model convergence to be achieved.
#' @param maxiter_survreg Maximum iterations used in survreg model fitting, default is 30.
#' @param initial_weighting Numeric, For the "full" model fitting: if 1: initial observation weights are estimated using linear regression at the highest and lowest tested concentrations. If 2, used a randomized start suitable for simulation studies on model validity but otherwise not recommended. If 3 or greater, fits a linear models to the components and estimates intial weights based on this model fit. For the reduced model fitting: 1 sets initial weights corresponding to fixed side (and extra_row) where observations not outside the range on the side corresponding to the fixed side are forced to be in the component where mu is being estimated. For initial weighting two a linear model is fit for the component still being estimated to provide initial observation weights.
#' @param sd_initial Numeric, value greater than 0 and less than 1. Proportion of the range from the highest concentration to lowest concentration that is used as the initial estimate of sigma for the estimated components. Default is 0.2
#' @param reruns_allowed Numeric, if the cross-validation for a particular combination of degrees (polynomial) or degrees of freedom (pspline) fails, how many repeat attempts should be allowed?
#' @param max_out_break Logical, if TRUE when a CV fold reaches maximum iterations it breaks the loop and moves to next rerun if applicable
#'
#' @importFrom purrr map_dfc
#' @importFrom survival survreg coxph.wtest
#' @import purrr
#' @import survival
#'
#' @return
#' @export
#'
#' @examples
#' \donttest{
#' data = simulate_mics()
#' result = fit_EM(model = "pspline",
#' approach = "full",
#' pre_set_degrees = c(4,4),
#' visible_data = data,
#' non_linear_term = "t",
#' covariates = NULL,
#' pi_formula = c == "2" ~ s(t),
#' max_it = 300,
#' ncomp = 2,
#' tol_ll = 1e-6,
#' pi_link = "logit",
#' verbose = 1,
#' model_coefficient_tolerance = 0.00001,
#' initial_weighting = 3,
#' sd_initial = 0.2
#' )
#' }
#'
#'
fit_EM = function(model = "pspline", #"polynomial",
                  approach = "full", #"reduced"
                  pre_set_degrees = NULL, #c(7,7)
                  max_degree = 8,
                  degree_sets = "matched",
                  #"independent"
                  visible_data,
                  nfolds = 10,
                  non_linear_term = "t",
                  covariates = NULL,
                  pi_formula = c == "2" ~ s(t),
                  fixed_side = NULL,
                  extra_row = FALSE,
                  ecoff = NA,
                  max_it = 3000,
                  ncomp = 2, #relevant
                  tol_ll = 1e-6,
                  pi_link = "logit",
                  verbose = 3,
                  model_coefficient_tolerance = 0.00001,
                  maxiter_survreg = 30,
                  initial_weighting = 3,
                  sd_initial = 0.2,
                  scale = NULL,
                  reruns_allowed = 3,
                  max_out_break = FALSE) {
  ##check here if approach is reduced but fixed side is null then we have a problem
 if(model == "pspline"){
   model = "surv"
 }

   if (is.null(pre_set_degrees)) {
    cv_results_intermediate = full_cv(
      model = model,
      approach = approach,
      max_degree = max_degree,
      degree_sets = degree_sets,
      visible_data = visible_data,
      nfolds = nfolds,
      non_linear_term = non_linear_term,
      covariates = covariates,
      pi_formula = pi_formula,
      fixed_side = fixed_side,
      extra_row = extra_row,
      ecoff = ecoff,
      max_it = max_it,
      ncomp = ncomp,
      tol_ll = tol_ll,
      pi_link = pi_link,
      verbose = verbose,
      model_coefficient_tolerance = model_coefficient_tolerance,
      maxiter_survreg = maxiter_survreg,
      initial_weighting = initial_weighting,
      sd_initial = sd_initial,
      scale = scale,
      reruns_allowed = reruns_allowed,
      max_out_break = max_out_break
    )


    cv_results =
      cv_results_intermediate %>%
      summarize(
        .by = starts_with("degree"),
        log_likelihood = sum(fold_likelihood),
        total_repeats = sum(repeats)
      ) %>%
      arrange(desc(log_likelihood))



    mu_formula = write_all_formulas(non_linear_term,
                                    pull_top_degree_set(cv_results),
                                    covariates, model)
  } else{
    if(model == "surv" & any(pre_set_degrees == 1)){
      errorCondition("degrees of freedom for pspline in surv package must be at least 2")
    }
    mu_formula = write_all_formulas(non_linear_term, pre_set_degrees, covariates, model)
    cv_results = NULL
  }
  ##add a choice here between full and reduced models FIX

  if(approach == "full"){
    output = EM_algorithm(
      visible_data = visible_data,
      model = model,
      mu_formula = mu_formula,
      pi_formula = pi_formula,
      max_it = max_it,
      ncomp = ncomp,
      tol_ll = tol_ll,
      browse_at_end = FALSE,
      browse_each_step = FALSE,
      plot_visuals = FALSE,
      prior_step_plot = FALSE,
      pause_on_likelihood_drop = FALSE,
      pi_link = pi_link,
      verbose = verbose,
      model_coefficient_tolerance = model_coefficient_tolerance,
      maxiter_survreg = maxiter_survreg,
      initial_weighting = initial_weighting,
      sd_initial = sd_initial,
      stop_on_likelihood_drop = FALSE,
      n_models = 100,
      seed = NULL,
      randomize = "all",
      scale = scale
    )
  }else if(approach == "reduced" & !is.null(fixed_side)){
    output = EM_algorithm_reduced(fixed_side = fixed_side,
                                  extra_row = extra_row,
                                  ecoff = ecoff,
                                  visible_data = visible_data,
                                  model = model,
                                  mu_formula = mu_formula,
                                  pi_formula = pi_formula,
                                  max_it = max_it,
                                  ncomp = ncomp,
                                  tol_ll = tol_ll,
                                  browse_at_end = FALSE,
                                  browse_each_step = FALSE,
                                  plot_visuals = FALSE,
                                  prior_step_plot = FALSE,
                                  pause_on_likelihood_drop = FALSE,
                                  pi_link = pi_link,
                                  verbose = verbose,
                                  model_coefficient_tolerance = model_coefficient_tolerance,
                                  maxiter_survreg = maxiter_survreg,
                                  initial_weighting = initial_weighting,
                                  sd_initial = sd_initial,
                                  stop_on_likelihood_drop = FALSE,
                                  non_linear_term = non_linear_term,
                                  covariates = covariates,
                                  scale = scale
    )
  }else{
    errorCondition("Values for approach are 'full' and 'reduced', if using reduced model, supply a value for fixed_side (RC or LC) and consider extra_row")
  }

  output$cv_results = cv_results

  output %>% return()

}




pull_top_degree_set = function(cv_results) {
  cv_results %>% select(-c(log_likelihood,total_repeats)) %>% purrr::map_dbl(., head(1)) %>% unname %>% return()
}
