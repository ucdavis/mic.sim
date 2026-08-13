#' Title
#'
#'
#' @inheritParams fit_EM
#' @param visible_data
#' @param model
#' @param mu_formula A formula for a survreg object from the survival package, left side of equation should be a surv object using "interval2" format, right side should be the non-linear term (polynomial or pspline) and any covariates. Can be a single formula or a list of formulas where length is equal to the number of components where the trend in the mean is being estimated.
#' @param pi_formula
#' @param max_it
#' @param ncomp
#' @param tol_ll
#' @param browse_at_end For internal model testing
#' @param browse_each_step For internal model testing
#' @param plot_visuals For internal model testing
#' @param prior_step_plot For internal model testing
#' @param pause_on_likelihood_drop For internal model testing
#' @param pi_link
#' @param verbose
#' @param model_coefficient_tolerance
#' @param maxiter_survreg
#' @param initial_weighting Numeric, if 1: initial observation weights are estimated using linear regression at the highest and lowest tested concentrations. If 2, used a randomized start suitable for simulation studies on model validity but otherwise not recommended. If 3 or greater, fits a linear models to the components and estimates intial weights based on this model fit.
#' @param sd_initial
#' @param stop_on_likelihood_drop For internal model testing
#' @param n_models Currently deprecated, will be used in future versions as part of model validation
#' @param seed Currently deprecated, will be used in future versions as part of model validation
#' @param randomize Currently deprecated, will be used in future versions as part of model validation
#'
#' @return
#' @export
#'
#' @importFrom magrittr %<>%
#' @importFrom survival survreg.control survreg Surv pspline coxph.wtest
#'
#' @examples
EM_algorithm = function(
    visible_data,
    model = "surv", #"mgcv", "polynomial"
    mu_formula = Surv(time = left_bound,
                      time2 = right_bound,
                      type = "interval2") ~ pspline(t, df = 4),
    #mu_formula = yi ~ s(t),
    pi_formula = c == "2" ~ s(t), #or: c == "2" ~ lo(t)
    max_it = 3000,
    ncomp = 2,
    tol_ll = 1e-6,
    browse_at_end = FALSE,
    browse_each_step = FALSE,
    plot_visuals = FALSE,
    prior_step_plot = FALSE,
    pause_on_likelihood_drop = FALSE,
    pi_link = "logit",
    verbose = 3,
    model_coefficient_tolerance = 0.00001,
    maxiter_survreg = 30,
    initial_weighting = 1,
    sd_initial = 0.2,
    stop_on_likelihood_drop = FALSE,
    n_models = 100,
    seed = NULL,
    randomize = "all",
    non_linear_term = "t",
    covariates = NULL,
    scale = NULL
){
  if(model == "pspline"){
    model = "surv"
  }

  #add attribute model to visible data

  visible_data = modify_visible_data(visible_data, model)

  visible_data = set_scale_log(visible_data, scale)

  set_algorithm_seed(seed)


  converge = NA_character_

  if(ncomp == 1){
    out = fit_single_component_model(visible_data, mu_formula, maxiter_survreg)
  attr(out$mu_model, "model") = model
  return(out)
    }else{

    #first E step-----
    possible_data = first_E_step(initial_weighting, visible_data, plot_visuals, sd_initial, ncomp, randomize, n_models, model, non_linear_term, covariates, pi_formula, max_it,tol_ll, pi_link, model_coefficient_tolerance)

    #wrapper function
    plot_initial_weighting_regression(possible_data = possible_data)
    #browser?
    likelihood_documentation = set_up_likelihood_matrix(max_it)
    possible_data = modify_bounds(possible_data)

    possible_data = add_scale(possible_data, attr(visible_data, "scale"))


    for(i in 1:max_it){
      message_iteration_count(i, verbose)
      #first M step--------
      #MLE of all parameters

      if(i != 1){
        prior_iteration = save_previous_iteration(mu_models_new, pi_model_new, log_likelihood_new, possible_data)
      }

      mu_models_new = fit_all_mu_models(possible_data, ncomp, mu_formula, approach = "full", fixed_side = NULL, maxiter_survreg)
      if(check_mu_models_convergence(mu_models_new, ncomp)){
        converge = "NO"
        break #if wrapping into an M-step function, add list to list of outputs, then check between steps
      }

      likelihood_documentation = M_step_likelihood_matrix_updates(likelihood_documentation, i, mu_models_new, ncomp, maxiter_survreg) ###use dimnames in matrix


      pi_model_new = fit_mgcv_pi_model_safe_reformat(pi_formula = pi_formula, pi_link = pi_link, possible_data = possible_data)

      if(check_error(pi_model_new)){
        converge = "NO"
        break
      }

      #pi_model_new = reformat_safe(pi_model_new_list)


      if(i > 1){
        model_coefficient_checks_results = model_coefficient_checks(mu_models_new, pi_model_new, prior_iteration$mu_models_old, prior_iteration$pi_model_old, model_coefficient_tolerance, ncomp)
      }




      #Next E step-------------

      likelihood_documentation[i, "m_step_check_new"] <- m_step_check_maximizing(possible_data, mu_models_new, pi_model_new) ###use dimnames in matrix
      if(i > 1){
        likelihood_documentation[i, "m_step_check_old"] <- m_step_check_maximizing(possible_data, prior_iteration$mu_models_old, prior_iteration$pi_model_old) ###use dimnames in matrix
      }else{
        likelihood_documentation[i, "m_step_check_old"] <- NaN
      }



      possible_data = E_step(possible_data, mu_models_new, pi_model_new)

      if(verbose > 2){
        print(pi_model_new)
        print(mu_models_new)
      }



      log_likelihood_new = calculate_log_likelihood(possible_data, verbose)
      likelihood_documentation[i, "loglikelihood"] <- log_likelihood_new ###use dimnames

      if(i > 1 && likelihood_documentation[i, "loglikelihood"] < likelihood_documentation[i - 1, "loglikelihood"] & stop_on_likelihood_drop){
        converge = "likelihood decreased"
        break
      }

      if(i > 1 & pause_on_likelihood_drop && likelihood_documentation[i, "loglikelihood"] < likelihood_documentation[i - 1, "loglikelihood"]){
        browser("likelihood decreased")
      }

      if(browse_each_step & plot_visuals){

        browser(message("End of step ", i))
        plot_fm_step(pi_model_new, mu_models_new, ncomp, possible_data, prior_step_plot, i)

        if(i > 1){
          plot_likelihood(likelihood_documentation, format = "matrix")
        }
      }

      if(browse_each_step & !plot_visuals){browser(message("End of step ", i))}

      if(i != 1){
        ##too much overlap with 103
        check_ll = (log_likelihood_new - prior_iteration$log_likelihood_old)

        if(check_ll < 0 & verbose >= 2){
          warning("Log Likelihood decreased")  ###  HAS BEEN GETTING USED A LOT, WHY IS THE LOG LIKELIHOOD GOING DOWN????
        }



        if(check_ll < tol_ll & model_coefficient_checks_results)
        {
          if(verbose > 0){
            message("Stopped on combined LL and parameters")}
          converge = "YES"
          break
        }

      }



    }


    ####group 149 to 173
    if(i > 1){
      if(i == max_it & !(check_ll < tol_ll & model_coefficient_checks_results)){
        converge = "iterations"
      }
    }

    if(i == 1 & i == max_it & is.na(converge)){
      converge = "iterations"
    }


    if(i == 1){
      prior_iteration = list(
        mu_models_old = NA,
        pi_model_old = NA,
        log_likelihood_old = NA,
        possible_data_old = NA
      )
    }

    if(initial_weighting == 2){
      random_start_set = randomize
    }else{
      random_start_set = NA_character_
    }

    return(
      list(
        likelihood = tibble_like(likelihood_documentation),
        model = model,
        possible_data = possible_data,
        pi_model = pi_model_new,
        mu_model = mu_models_new,
        steps = i,
        converge = converge,
        ncomp = ncomp,
        prior_step_models = prior_iteration,
        seed = seed,
        random_start_set = random_start_set,
        sd_initial = ifelse(initial_weighting >= 7, sd_initial, NaN),
        mu_formula = mu_formula
      )
    )
  }
}

fit_mu_model.mgcv = function(possible_data, pred_comp, mu_formula){
  df = possible_data %>% filter(`P(C=c|y,t)` > 0 & c == pred_comp)
  df$yi = cbind(df$left_bound_mgcv, df$right_bound_mgcv)
  mgcv::gam(mu_formula, family= mgcv::cnorm(link = "identity"), weights = `P(C=c|y,t)`, data=df, method = "REML") %>% return()
}

fit_single_component_model.mgcv = function(visible_data, mu_formula){

  possible_data = visible_data %>%
    modify_bounds.mgcv()
  mu_model = mgcv::gam(mu_formula, family= mgcv::cnorm(link = "identity"), data = possible_data, method = "REML")

  return(list(possible_data = possible_data,
              mu_model = mu_model,
              converge = ifelse(length(mu_model) > 1, "YES", "NO"),
              ncomp = 1))
}

fit_single_component_model.surv = function(visible_data, mu_formula, maxiter_survreg){
  mu_model = visible_data %>% survival::survreg(
    mu_formula,  ##Make this chunk into an argument of the function
    data = .,
    dist = "gaussian",
    control = survreg.control(maxiter = ifelse(!is.null(maxiter_survreg), maxiter_survreg, 30))
    )

  return(list(possible_data = visible_data,
              mu_model = mu_model,
              converge = "YES",
              ncomp = 1,
              likelihood = tibble(step = 1, likelihood = mu_model$loglik[2]))
  )
}

modify_bounds.mgcv = function(data){
  data %>% mutate(
    left_bound_mgcv =
      case_when(
        left_bound == -Inf ~ right_bound,
        TRUE ~ left_bound
      ),
    right_bound_mgcv =
      case_when(
        left_bound == -Inf ~ -Inf,
        TRUE ~ right_bound
      ),
    yi = cbind(left_bound_mgcv, right_bound_mgcv)
  ) %>% return()
}

modify_bounds = function(data){
  method = attr(data, "model")
  if(method == "mgcv"){
    modify_bounds.mgcv(data) %>% return()
  }else{
    return(data)
  }
}

add_attribute_data = function(data, model){
  if(model %in% c("surv", "mgcv", "polynomial")){
    attr(data, "model") <- model
    return(data)
  }else{
    stop("Invalid value for model use 'surv', 'mgcv', or 'polynomial'")
  }
}


fit_single_component_model = function(visible_data, mu_formula, maxiter_survreg){

  if(!is.list(mu_formula)){
   mu_formula = list(mu_formula)
  }


  if(attr(visible_data, "model") == "mgcv"){
   mu_fit = purrr::map(mu_formula, ~fit_single_component_model.mgcv(visible_data, .x)) %>% unlist(., recursive = FALSE)
  }
  if(attr(visible_data, "model") %in% c("surv", "polynomial")){
   mu_fit = purrr::map(mu_formula, ~fit_single_component_model.surv(visible_data, .x, maxiter_survreg)) %>% unlist(., recursive = FALSE)
  }

  attr(mu_fit, "model") = attr(visible_data, "model")
  return(mu_fit)

}


set_up_likelihood_matrix = function(max_it){
  likelihood_documentation <- matrix(data = NA, nrow = max_it, ncol = 7, dimnames = list(NULL, c("step", "loglikelihood", "survreg_maxout", "m_step_check_new", "m_step_check_old", "scale_comp_1", "scale_comp_2")))

  likelihood_documentation [,"step"] <- 1:max_it
  return(likelihood_documentation)
}

message_iteration_count = function(i, verbose = 3){
  if(verbose > 1){
    message("starting iteration number ", i)}
}

save_previous_iteration = function(mu_models_new, pi_model_new, log_likelihood_new, possible_data){
  list(
    mu_models_old = mu_models_new,
    pi_model_old = pi_model_new,
    log_likelihood_old = log_likelihood_new,
    possible_data_old = possible_data
  ) %>% return()
}

M_step_likelihood_matrix_updates.mgcv = function(likelihood_documentation, i, mu_models_new){
  likelihood_documentation[i, "scale_comp_1"] = mu_models_new[[1]]$family$getTheta(TRUE)
  likelihood_documentation[i, "scale_comp_2"] = mu_models_new[[2]]$family$getTheta(TRUE)
  return(likelihood_documentation)
}

M_step_likelihood_matrix_updates.surv = function(likelihood_documentation, i, mu_models_new, ncomp, maxiter_survreg){
  likelihood_documentation[i, "scale_comp_1"] = mu_models_new[[1]]$scale
  likelihood_documentation[i, "scale_comp_2"] = mu_models_new[[2]]$scale
  likelihood_documentation[i, "survreg_maxout"] = check_survreg_iteration_maxout(mu_models_new, ncomp, maxiter_survreg)
  return(likelihood_documentation)
}

M_step_likelihood_matrix_updates = function(likelihood_documentation, i, mu_models_new, ncomp, maxiter_survreg = 30){
  if(attr(mu_models_new, "model") == "mgcv"){
    M_step_likelihood_matrix_updates.mgcv(likelihood_documentation, i, mu_models_new) %>% return()
  }else{
    M_step_likelihood_matrix_updates.surv(likelihood_documentation, i, mu_models_new, ncomp, maxiter_survreg) %>% return()
  }
}

check_error = function(model){
  !is.null(attr(model, "error")) %>% return()
}

check_mu_models_convergence.mgcv = function(mu_models, ncomp){
  purrr::map(1:ncomp, ~any(is.na(mu_models[[.x]]$family$getTheta(TRUE)), is.na(mu_models[[.x]]$coefficients))) %>% return()
}

check_mu_models_convergence.surv = function(mu_models, ncomp){
  purrr::map(1:ncomp, ~any(is.na(mu_models[[.x]]$scale), is.na(mu_models[[.x]]$coefficients)))
}

check_mu_models_convergence = function(mu_models, ncomp){
  if(any(purrr::map_lgl(1:ncomp, ~check_error(model = mu_models[[.x]])))){
   return(TRUE)
  }

  if(attr(mu_models, "model") == "mgcv"){
    check_mu_models_convergence.mgcv(mu_models, ncomp) %>% unlist %>% any %>% return()
  }else{
    check_mu_models_convergence.surv(mu_models, ncomp) %>% unlist %>% any %>% return()
  }
}

fit_mgcv_pi_model = function(pi_formula, pi_link, possible_data){

  if(pi_link == "logit"){
    pi_model = gam(pi_formula, family = binomial(link = "logit"), data = possible_data, weights = `P(C=c|y,t)`, method = "ML") %>% suppressWarnings()
  } else if(pi_link == "identity"){
    pi_model = gam(pi_formula, family = binomial(link = "identity"), data = possible_data, weights = `P(C=c|y,t)`, method = "ML") %>% suppressWarnings()
  }
  if(pi_link == "logit_simple"){
    pi_model = glm(c == "2" ~ t, family = binomial(link = "logit"), data = possible_data, weights = `P(C=c|y,t)`) %>% suppressWarnings()
    }else{ errorCondition("pick logit or identity link function")}

  return(pi_model)
}

fit_mgcv_pi_model_safe = safely(fit_mgcv_pi_model, otherwise = "Error")

fit_mgcv_pi_model_safe_reformat = function(pi_formula, pi_link, possible_data){
pi_model_new_list = fit_mgcv_pi_model_safe(pi_formula = pi_formula, pi_link = pi_link, possible_data = possible_data)

pi_model_new = reformat_safe(pi_model_new_list)
return(pi_model_new)
}


reformat_safe = function(model){
  out = model$result
  attr(out, "error") = model$error
  return(out)
}

#gam_safe = safely(mgcv::gam, otherwise = "Error")
#glm_safe = safely(glm, otherwise = "Error")

model_coefficient_checks.mgcv = function(mu_models_new, pi_model_new, mu_models_old, pi_model_old, model_coefficient_tolerance, ncomp){
  #do the weird coefficients gam returns chaange (only one for s(t) for some reason)
  pi_parametric_coef_check = max((pi_model_new %>% coefficients()) - (pi_model_old %>% coefficients())) < model_coefficient_tolerance

  #these are the actual smoothing things for every observation I believe
  #pi_nonparametric_check = max(pi_model_new$smooth - pi_model_old$smooth) < model_coefficient_tolerance

  #check if the number of coefficients in the mu models changes
  mu_number_of_coef_check = purrr::map(1:ncomp, ~(length(na.omit(mu_models_new[[.x]]$coefficients)) == length(na.omit(mu_models_old[[.x]]$coefficients)))) %>%
    unlist %>% all

  mu_coefficient_check = purrr::map(1:ncomp, ~((mu_models_new[[.x]]$coefficients - mu_models_old[[.x]]$coefficients) %>% abs %>% sum)) %>% unlist %>% max < model_coefficient_tolerance

  mu_sigma_check = purrr::map(1:ncomp, ~((mu_models_new[[.x]]$family$getTheta(TRUE) - mu_models_old[[.x]]$family$getTheta(TRUE)) %>% abs)) %>% unlist %>% max < model_coefficient_tolerance

  #check likelihood at end of E step

  return(all(pi_parametric_coef_check,
             #pi_nonparametric_check,
             mu_number_of_coef_check,
             mu_coefficient_check,
             mu_sigma_check))
}
 #break into sub-functions
model_coefficient_checks.surv = function(mu_models_new, pi_model_new, mu_models_old, pi_model_old, model_coefficient_tolerance, ncomp){
  #do the weird coefficients gam returns chaange (only one for s(t) for some reason)
  pi_parametric_coef_check = max((pi_model_new %>% coefficients()) - (pi_model_old %>% coefficients())) < model_coefficient_tolerance

  #check if the number of coefficients in the mu models changes
  mu_number_of_coef_check = purrr::map(1:ncomp, ~(length(na.omit(mu_models_new[[.x]]$coefficients)) == length(na.omit(mu_models_old[[.x]]$coefficients)))) %>%
    unlist %>% all

  mu_coefficient_check = purrr::map(1:ncomp, ~((mu_models_new[[.x]]$coefficients - mu_models_old[[.x]]$coefficients) %>% abs %>% sum)) %>% unlist %>% max < model_coefficient_tolerance

  mu_sigma_check = purrr::map(1:ncomp, ~((mu_models_new[[.x]]$scale - mu_models_old[[.x]]$scale) %>% abs)) %>% unlist %>% max < model_coefficient_tolerance

  #check likelihood at end of E step

  return(all(pi_parametric_coef_check, mu_number_of_coef_check, mu_coefficient_check, mu_sigma_check))
}
 #break into sub-functions

model_coefficient_checks = function(mu_models_new, pi_model_new, mu_models_old, pi_model_old, model_coefficient_tolerance, ncomp, i){
    if (attr(mu_models_new, "model") == "mgcv") {
      model_coefficient_checks.mgcv(
        mu_models_new,
        pi_model_new,
        mu_models_old,
        pi_model_old,
        model_coefficient_tolerance,
        ncomp
      ) %>% return()
    } else{
      model_coefficient_checks.surv(
        mu_models_new,
        pi_model_new,
        mu_models_old,
        pi_model_old,
        model_coefficient_tolerance,
        ncomp
      ) %>% return()
    }
}

get_scale = function(mu_model){
  if(attr(mu_model, "model") == "mgcv"){
    mu_model$family$getTheta(TRUE) %>% return()
  }else{
    mu_model$scale %>% return()
  }
}

calculate_density_obs =  function(possible_data, mu_models){
  possible_data %>%
    mutate(
      `E[Y|t,c]` = case_when(c == "1" ~ predict(mu_models[[1]], newdata = possible_data),
                             c == "2" ~ predict(mu_models[[2]], newdata = possible_data),
                             TRUE ~ NaN),
      `sd[Y|t,c]` = case_when(c == "1" ~ get_scale(mu_models[[1]]),
                              c == "2" ~ get_scale(mu_models[[2]]), #1,
                              TRUE ~ NaN),
      `P(Y|t,c)` = case_when(
        left_bound == right_bound ~ dnorm(x = left_bound, mean = `E[Y|t,c]`, sd =  `sd[Y|t,c]`),
        left_bound <= `E[Y|t,c]` ~ pnorm(right_bound, mean = `E[Y|t,c]`, sd =  `sd[Y|t,c]`) -
          pnorm(left_bound, mean = `E[Y|t,c]`, sd =  `sd[Y|t,c]`),
        TRUE ~ pnorm(left_bound, mean = `E[Y|t,c]`, sd =  `sd[Y|t,c]`, lower.tail = FALSE) -
          pnorm(right_bound, mean = `E[Y|t,c]`, sd =  `sd[Y|t,c]`, lower.tail = FALSE)
      )
    ) %>% return()
}

pi_model_predictions = function(possible_data, pi_model_new){
  possible_data %>% mutate(
    `P(C=c|t)` = case_when(
      c == "2" ~ predict(pi_model_new, newdata = possible_data, type = "response"),
      c == "1" ~ 1 - predict(pi_model_new, newdata = possible_data, type = "response")
    )
  ) %>% return()
}

calculate_new_weights = function(possible_data) {
  possible_data %>% mutate(`P(c,y|t)` = `P(C=c|t)` * `P(Y|t,c)`) %>%
    mutate(.by = obs_id,
           `P(Y=y|t)` = sum(`P(c,y|t)`)) %>%
    mutate(`P(C=c|y,t)` = `P(c,y|t)` / `P(Y=y|t)`)
}


E_step = function(possible_data, mu_models, pi_model){
  possible_data %<>% calculate_density_obs(., mu_models) %>%
    pi_model_predictions(., pi_model) %>%
    calculate_new_weights() %>% return()
}

tibble_like <- function(likelihood_documentation, model = "surv"){
  likelihood_documentation %>% as_tibble %>% filter(!is.na(.data$loglikelihood)) %>% return()
}

modify_visible_data = function(visible_data, model, scale = NULL){
  visible_data = add_attribute_data(visible_data, model)
  visible_data = add_obs_id(visible_data)
  visible_data = add_scale(visible_data, scale)
  return(visible_data)
}

add_scale = function(data, scale){
  if(is.null(attr(data, "scale"))){
    if(is.null(scale)){
      errorCondition("This data set does not have a scale attribute already attached, please use 'scale' argument to provide one: acceptable values are log, fold, log2, MIC")
    }else{
      scale_modified = case_when(
        tolower(scale) %in% c("log", "fold", "log2") ~ "log",
        tolower(scale) %in% c("mic", "concentration") ~ "MIC",
        TRUE ~ "error"
      )

      if(scale_modified == "error"){errorCondition("Invalid value of scale, please use 'log', 'fold', 'log2', or 'MIC'")}

      attr(data, "scale") = scale
    }
  }
  return(data)
}

set_algorithm_seed = function(seed){
  if(!is.null(seed)){
    set.seed(seed)
  }
}




































m_step_check_maximizing = function(possible_data, mu_models, pi_model){
  possible_data %>%
    mutate(
      `E[Y|t,c]` = case_when(c == "1" ~ predict(mu_models[[1]], newdata = possible_data),
                             c == "2" ~ predict(mu_models[[2]], newdata = possible_data),
                             TRUE ~ NaN),
      `sd[Y|t,c]` = case_when(c == "1" ~ get_scale(mu_models[[1]]),
                              c == "2" ~ get_scale(mu_models[[2]]), #1,
                              TRUE ~ NaN),
      `P(Y|t,c)` = case_when(
        left_bound == right_bound ~ dnorm(x = left_bound, mean = `E[Y|t,c]`, sd =  `sd[Y|t,c]`),
        left_bound <= `E[Y|t,c]` ~ pnorm(right_bound, mean = `E[Y|t,c]`, sd =  `sd[Y|t,c]`) -
          pnorm(left_bound, mean = `E[Y|t,c]`, sd =  `sd[Y|t,c]`),
        TRUE ~ pnorm(left_bound, mean = `E[Y|t,c]`, sd =  `sd[Y|t,c]`, lower.tail = FALSE) -
          pnorm(right_bound, mean = `E[Y|t,c]`, sd =  `sd[Y|t,c]`, lower.tail = FALSE)
      ),
      `P(C=c|t)` = case_when(
        c == "2" ~ predict(pi_model, newdata = tibble(t = t), type = "response"),
        c == "1" ~ 1 - predict(pi_model, newdata = tibble(t = t), type = "response")
      ),
      `P(c,y|t)` = `P(C=c|t)` * `P(Y|t,c)`
    ) %>% select(obs_id, c, `P(c,y|t)`, `P(C=c|y,t)`) %>%
    mutate(m_step_check = `P(C=c|y,t)` * log(`P(c,y|t)`)) %>% pull(m_step_check) %>% sum %>% return()



}





# fit_pi_model = function(pi_formula, pi_link, possible_data){
#
#   if(pi_link == "logit"){
#     pi_model = gam(pi_formula, family = binomial(link = "logit"), data = possible_data, weights = `P(C=c|y,t)`)
#   } else if(pi_link == "identity"){
#
#     pi_model = gam(pi_formula, family = binomial(link = "identity"), data = possible_data, weights = `P(C=c|y,t)`)
#   }else{ errorCondition("pick logit or identity link function")}
#
#
#   return(pi_model)
# }

check_mu_models_convergence_surv = function(mu_models, ncomp){
  purrr::map(1:ncomp, ~any(is.na(mu_models[[.x]]$scale), is.na(mu_models[[.x]]$coefficients)))
}








m_step_check_maximizing_mgcv = function(possible_data, mu_models, pi_model){
  possible_data %>%
    mutate(
      `E[Y|t,c]` = case_when(c == "1" ~ predict(mu_models[[1]], newdata = possible_data),
                             c == "2" ~ predict(mu_models[[2]], newdata = possible_data),
                             TRUE ~ NaN),
      #predict(model, newdata = possible_data),
      `sd[Y|t,c]` = case_when(c == "1" ~ mu_models[[1]]$family$getTheta(TRUE),
                              c == "2" ~ mu_models[[2]]$family$getTheta(TRUE), #1,
                              TRUE ~ NaN),
      #model$scale[c], #####QUESTION HERE????????????????????????????
      # `Var[Y|t,c]` = `sd[Y|t,c]`^2,

      `P(Y|t,c)` = case_when(
        left_bound == right_bound ~ dnorm(x = left_bound, mean = `E[Y|t,c]`, sd =  `sd[Y|t,c]`),
        left_bound <= `E[Y|t,c]` ~ pnorm(right_bound, mean = `E[Y|t,c]`, sd =  `sd[Y|t,c]`) -
          pnorm(left_bound, mean = `E[Y|t,c]`, sd =  `sd[Y|t,c]`),
        TRUE ~ pnorm(left_bound, mean = `E[Y|t,c]`, sd =  `sd[Y|t,c]`, lower.tail = FALSE) -
          pnorm(right_bound, mean = `E[Y|t,c]`, sd =  `sd[Y|t,c]`, lower.tail = FALSE)
      ),
      `P(C=c|t)` = case_when(
        c == "2" ~ predict(pi_model, newdata = tibble(t = t), type = "response"),
        c == "1" ~ 1 - predict(pi_model, newdata = tibble(t = t), type = "response")
      ),
      `P(c,y|t)` = `P(C=c|t)` * `P(Y|t,c)`
    ) %>% select(obs_id, c, `P(c,y|t)`, `P(C=c|y,t)`) %>%
    mutate(m_step_check = `P(C=c|y,t)` * log(`P(c,y|t)`)) %>% pull(m_step_check) %>% sum %>% return()
}

set_scale_log = function(visible_data, scale){
  visible_data = add_scale(visible_data, scale)
  if(attr(visible_data, "scale") == "MIC"){
    visible_data = visible_data %>% mutate(left_bound = log2(left_bound),
                                           right_bound = log2(right_bound))
  }
  return(visible_data)

}

