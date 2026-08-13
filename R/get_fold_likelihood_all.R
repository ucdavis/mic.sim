#' Title
#'
#' @param model
#' @param approach
#' @param i
#' @param visible_data
#' @param degrees
#' @param non_linear_term
#' @param covariates
#' @param pi_formula
#' @param fixed_side
#' @param extra_row
#' @param ecoff
#' @param max_it
#' @param ncomp
#' @param tol_ll
#' @param pi_link
#' @param verbose
#' @param model_coefficient_tolerance
#' @param maxiter_survreg
#' @param initial_weighting
#' @param sd_initial
#' @param scale
#'
#' @return
#' @keywords internal
#'
#' @examples
get_fold_likelihood_all = function(model = "surv",
                                   approach = "full",
                                   i,
                                   visible_data,
                                   degrees,
                                   non_linear_term = "t",
                                   covariates = NULL,
                                   pi_formula = c == "2" ~ s(t),
                                   fixed_side = NULL,
                                   extra_row = FALSE,
                                   ecoff = NA,
                                   max_it = 300,
                                   ncomp = 2,
                                   tol_ll = 1e-6,
                                   pi_link = "logit",
                                   verbose = 3,
                                   model_coefficient_tolerance = 0.00001,
                                   maxiter_survreg = 30,
                                   initial_weighting = 3,
                                   sd_initial = 0.2,
                                   scale = NULL,
                                   max_out_break = FALSE) {
  test = i

  if(verbose >= 1){
    message(paste0("fold ", i))
  }

  mu_formula = write_all_formulas(non_linear_term, degrees, covariates, model)

  ##add check for if fold column exists
  training_set = visible_data %>% filter(fold != test)
  testing_set = visible_data %>% filter(fold == test)


  ##OPTION FOR FULL OR REDUCED MODEL HERE
  if(approach == "full"){
    trained_model = EM_algorithm(
      visible_data = training_set,
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
      non_linear_term = non_linear_term,
      covariates = covariates,
      scale = scale
    )
  }else if(approach == "reduced"){
    trained_model = EM_algorithm_reduced(fixed_side = fixed_side,
                                         extra_row = extra_row,
                                         ecoff = ecoff,
                                         visible_data = training_set,
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
    errorCondition("Use 'full' or 'reduced' for approach")
  }

  calculate_fold_likelihood_all(testing_set, trained_mu_model = trained_model$mu_model, trained_pi_model = trained_model$pi_model, approach = approach, fixed_side = fixed_side, extra_row = extra_row, ecoff = ecoff, ncomp = ncomp, converge = trained_model$converge, max_out_break = max_out_break) %>%
    return()  ###MAY NEED TO FIX THIS TOO

}



get_fold_likelihood_all_safe = purrr::possibly(get_fold_likelihood_all, otherwise = NaN)


get_fold_likelihood_all_safe_single_output =
  function(model = "surv",
           approach = "full",
           i,
           visible_data,
           degrees,
           non_linear_term = "t",
           covariates = NULL,
           pi_formula = c == "2" ~ s(t),
           fixed_side = NULL,
           extra_row = FALSE,
           max_it = 3000,
           ncomp = 2,
           tol_ll = 1e-6,
           pi_link = "logit",
           verbose = 3,
           model_coefficient_tolerance = 0.00001,
           maxiter_survreg = 30,
           initial_weighting = 1,
           sd_initial = 0.2,
           scale = NULL){

    out = get_fold_likelihood_all_safe(
      model = model,
      approach = approach,
      i = i,
      visible_data = assign_folds(visible_data, nfolds),
      degrees,
      non_linear_term = non_linear_term,
      covariates = covariates,
      pi_formula = pi_formula,
      fixed_side = fixed_side,
      extra_row = extra_row,
      max_it = max_it,
      ncomp = ncomp,
      tol_ll = tol_ll,
      pi_link = pi_link,
      verbose = verbose,
      model_coefficient_tolerance = model_coefficient_tolerance,
      maxiter_survreg = maxiter_survreg,
      initial_weighting = initial_weighting,
      sd_initial = sd_initial,
      scale = scale
    )

    out$result %>% return()
  }


