#' Plot Output of Model Fitting
#'
#' Takes output of a fit_EM() run and plots it.
#' Logistic regression fitting is in progress (currently does not accept covariates and non-linear term must be "t")
#'
#' @param output List, output of a fit_EM()
#' @param title String, title of plot
#' @param add_log_reg Logical, add a curve to the pi (component weight) plot showing proportion resistant and susceptible
#' @param ecoff Numeric or String, represents an ECOFF on the MIC Scale to define the upper limit of the WT component, A number is interpreted as WT being <= ECOFF.
#' @param s_breakpoint String, represents S breakpoint (e.g. <= 2) on MIC scale (not log2 scale)
#' @param r_breakpoint String, represents R breakpoint (e.g. >= 32) on MIC scale (not log2 scale)
#' @param visual_split Numeric or String, represents a visual split point on the MIC scale to define the upper limit of the WT component. A number is interpreted as WT being <= visual_split
#' @param use_prior_step Logical, if one mu model did not converge, can try plotting mu models from previous step by setting this to TRUE
#' @param range_zoom Logical, zoom y axis to range of tested concentrations
#' @param plot_range Vector of length 2, minimum and maximum values of y axis of plot
#' @param start_date Integer, value at which x axis should start (year).
#' @param  x_axis_t_breaks Numerical vector, vector of values on the scale of t, the time variable in years from the start of the study. Helpful to use seq(0,t_max, by = spacing) where t_max is the length of study period and spacing is how many years to separate major ticks by
#' @param skip Vector, vector of either "ecoff", "bkpts", or c("ecoff", "bkpts"), to describe any splits for which logistic regression should not be plotted if another logistic regression is being plotted. If only one divider is used, just turn off add_log_reg
#'
#' @import ggplot2
#' @import ggnewscale
#' @importFrom patchwork wrap_plots
#' @importFrom latex2exp TeX
#'
#' @return
#' @export
#'
#' @examples
#' data = simulate_mics()
#' output = fit_EM(model = "pspline",
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
#' plot_fm(output = output, title = "Example", add_log_reg = TRUE, s_breakpoint = "<=1", r_breakpoint = ">=4")
#'
#'
plot_fm <- function(output, title ="", add_log_reg = FALSE, ecoff = NA, s_breakpoint = NA, r_breakpoint = NA, visual_split = NA, use_prior_step = FALSE, range_zoom = FALSE, plot_range = NULL, start_date = 0, x_axis_t_breaks = NULL, skip = NULL){
  #assumed_components is how many components were intended to be in the distribution, even if not estimated (e.g. for a reduced model where the NWT component is unestimated, assumed_components would still be 2)
  assumed_components = output$ncomp

  if(!is.null(output$prior_step_models) & use_prior_step){
    output$mu_model = output$prior_step_models$mu_models
  }

  if(!is.null(output$fixed_side)){
    if(output$fixed_side == "LC"){
      check = check_comp_conv(output$mu_model[[1]])
      dnc = c(TRUE, check)
      n_fitted_components = 1 - check
      results = tibble(c = 1:2, dnc)
      fitted_comp = output$mu_model[[1]]
    }else if(output$fixed_side == "RC"){
      check = check_comp_conv(output$mu_model[[1]])
      dnc = c(check, TRUE)
      n_fitted_components = 1 - check
      results = tibble(c = 1:2, dnc)
      fitted_comp = output$mu_model[[1]]
    }else{
      errorCondition("Invalid value for fixed_side")
    }
  }else{

    if(assumed_components == 2){
      results <- tibble(c = 1:2, dnc = purrr::map_lgl(output$mu_model, ~check_comp_conv(.x)))
      if(nrow(results %>% filter(dnc)) > 0){
        fitted_comp = output$mu_model[[results %>% filter(!dnc) %>% pull(c)]]
        n_fitted_components = 1
      } else{
        n_fitted_components = 2
        fitted_comp = NULL
      }
    }else{
      fitted_comp = output$mu_model
      n_fitted_components = 1
    }
  }

  df = output$possible_data %>% mutate(cens =
                                         case_when(
                                           left_bound == -Inf ~ "lc",
                                           right_bound == Inf ~ "rc",
                                           TRUE ~ "int"
                                         ),
                                       mid =
                                         case_when(
                                           left_bound == -Inf ~ right_bound - 0.5,
                                           right_bound == Inf ~ left_bound + 0.5,
                                           TRUE ~ (left_bound + right_bound) / 2
                                         ))

  attr(df, "model") <- attr(output$possible_data, "model")



  # Set the range of the plot. If not specified, use the plot_bounds function, otherwise just use the specified values
  if(is.null(plot_range)){
    plot_min <- plot_bounds(output$possible_data, "min", n_fitted_components, range_zoom, output, fitted_comp)
    plot_max <- plot_bounds(output$possible_data, "max", n_fitted_components, range_zoom, output, fitted_comp)
  }else{
    plot_min = plot_range[1]
    plot_max = plot_range[2]
  }

  mean = plot_mean(output = output, df = df, results = results, start_date = start_date, fitted_comp = fitted_comp, title = title, plot_min = plot_min, plot_max = plot_max, ecoff = ecoff, s_breakpoint = s_breakpoint, r_breakpoint = r_breakpoint, visual_split = visual_split, x_axis_t_breaks = x_axis_t_breaks, assumed_components = assumed_components, n_fitted_components = n_fitted_components)

  if(assumed_components > 1){
    pi = plot_pi(output = output, df = df, start_date = start_date, add_log_reg = add_log_reg, ecoff = ecoff, s_breakpoint = s_breakpoint, r_breakpoint = r_breakpoint, visual_split = visual_split, skip = skip, x_axis_t_breaks = x_axis_t_breaks)
    return(patchwork::wrap_plots(mean,pi, ncol = 1))
  }else{
    return(patchwork::wrap_plots(mean, ncol = 1))
  }

}


plot_bounds = function(df, side, n_fitted_components, range_zoom = FALSE, output, fitted_comp = NULL){
  if(side == "min"){
    if(nrow(df %>% filter(left_bound == -Inf)) > 0){
      plot_min_1 <- (df %>% filter(left_bound == -Inf) %>% pull(right_bound) %>% min(., na.rm = TRUE)) - 2.5
    }else{
      plot_min_1 <- (df %>% pull(left_bound) %>% min(., na.rm = TRUE)) - 2.5
    }

if(attr(df, "model") != "mgcv"){
    if(n_fitted_components == 2){
      plot_min_2 <- min(sim_pi_survreg_boot(df, fit = output$mu_model[[1]], alpha = 0.05, nSims = 10000) %>% pull(lwr) %>% min(., na.rm = TRUE) - 0.2,
                        sim_pi_survreg_boot(df, fit = output$mu_model[[2]], alpha = 0.05, nSims = 10000) %>% pull(lwr) %>% min(., na.rm = TRUE) - 0.2)
    } else if(n_fitted_components == 1){
      plot_min_2 <- sim_pi_survreg_boot(df, fit = fitted_comp, alpha = 0.05, nSims = 10000) %>% pull(lwr) %>% min(., na.rm = TRUE) - 0.2
    }else{
      plot_min_2 = plot_min_1
    }
}else{
  plot_min_2 = plot_min_1
}

    plot_min = min(plot_min_1, plot_min_2, na.rm = TRUE)
    if(range_zoom){
      return(plot_min_1)
    }else{
    return(plot_min)
    }
  } else if(side == "max"){
    if(nrow(df %>% filter(right_bound == Inf)) > 0){
      plot_max_1 <- (df %>% filter(right_bound == Inf) %>% pull(left_bound) %>% max(., na.rm = TRUE)) + 2.5
    }else{
      plot_max_1 <- (df %>% pull(right_bound) %>% max(., na.rm = TRUE)) + 2.5
    }

    if(attr(df, "model") != "mgcv"){

    if(n_fitted_components == 1){
      plot_max_2 <- sim_pi_survreg_boot(df, fit = fitted_comp, alpha = 0.05, nSims = 10000) %>% pull(upr) %>% max(., na.rm = TRUE) + 0.2
    } else if(n_fitted_components == 2){
      plot_max_2 <- max(sim_pi_survreg_boot(df, fit = output$mu_model[[1]], alpha = 0.05, nSims = 10000) %>% pull(upr) %>% max(., na.rm = TRUE) + 0.2,
                        sim_pi_survreg_boot(df, fit = output$mu_model[[2]], alpha = 0.05, nSims = 10000) %>% pull(upr) %>% max(., na.rm = TRUE) + 0.2)
    }else{
      plot_max_2 = plot_max_1
    }

    }else{
      plot_max_2 = plot_max_1
    }

    plot_max = max(plot_max_1, plot_max_2, na.rm = TRUE)
    if(range_zoom){
      return(plot_max_1)
    }else{
      return(plot_max)
    }
  }else{
    errorCondition("choose 'min' or 'max'")
  }
}

check_comp_conv = function(models){
  if((length(models) == 1 && models == "Error")){
    return(TRUE)
  }else{is.na(models$scale) | (tibble(a = models$coefficients) %>% filter(is.na(a)) %>% nrow) > 0}
}

get_two_comp_ci = function(output){
  tibble(t = rep(seq(0, max(output$possible_data$t), len = 300), 2)) %>%
    mutate(
      c1pred = predict(output$mu_model[[1]], tibble(t), se = T)$fit,
      c1pred_se = predict(output$mu_model[[1]], tibble(t), se = T)$se.fit,
      c1pred_lb = c1pred - 1.96 * c1pred_se,
      c1pred_ub = c1pred + 1.96 * c1pred_se,
      c2pred = predict(output$mu_model[[2]], tibble(t), se = T)$fit,
      c2pred_se = predict(output$mu_model[[2]], tibble(t), se = T)$se.fit,
      c2pred_lb = c2pred - 1.96 * c2pred_se,
      c2pred_ub = c2pred + 1.96 * c2pred_se,
    ) %>% return()}

offset_time_as_date_in_df = function(df, start_date){
  df %>% mutate(t = offset_time_as_date(t, start_date)) %>% return()
}


inverse_logit = function(x){1 / (1 + (exp(-x)))}

logit = function(p){log(p / (1 - p))}



set_y_labels = function(value){
  case_when(2^value >= 0.12 ~ 2^value,
    0.03 > 2^value  & 2^value >= 0.015 ~ signif(2^value, 2),
    TRUE ~ signif(2^value, 1)
    )
}

mu.se.brd <- function(t, c, z){predict(output$mu_model[[c]], data.frame(t = t)) + (z * predict(output$mu_model[[c]], data.frame(t = t), se = TRUE)$se.fit)}
mu.se.brd.fms <- function(t, z){predict(fitted_comp, data.frame(t = t)) + (z * predict(fitted_comp, data.frame(t = t), se = TRUE)$se.fit)}
