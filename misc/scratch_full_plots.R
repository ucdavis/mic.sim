plot_out = function(output, title ="", add_log_reg = FALSE, ecoff = NA, s_breakpoint = NA, r_breakpoint = NA, visual_split = NA, use_prior_step = FALSE, range_zoom = FALSE, plot_range = NULL, start_date = 0, x_axis_t_breaks = NULL, skip = NULL){
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

mean = plot_mean(output = output, df = df, results = results, start_date = start_date, fitted_comp = fitted_comp, title = title, plot_min = plot_min, plot_max = plot_max, ecoff = ecoff, s_breakpoint = s_breakpoint, r_breakpoint = r_breakpoint, visual_split = visual_split, x_axis_t_breaks = x_axis_t_breaks)

 if(assumed_components > 1){
   pi = plot_pi(output = output, df = df, start_date = start_date, add_log_reg = add_log_reg, ecoff = ecoff, s_breakpoint = s_breakpoint, r_breakpoint = r_breakpoint, visual_split = visual_split, skip = skip)
   return(patchwork::wrap_plots(mean,pi, ncol = 1))
   }else{
     return(patchwork::wrap_plots(mean, ncol = 1))
 }

}
