#' Title
#'
#' @param output
#' @param df
#' @param start_date
#' @param add_log_reg
#' @param ecoff
#' @param s_breakpoint
#' @param r_breakpoint
#' @param visual_split
#' @param skip
#' @param x_axis_t_breaks
#'
#' @keywords internal
#'
#' @return
#' @export
#'
#' @examples
plot_pi = function(output, df, start_date, add_log_reg, ecoff, s_breakpoint, r_breakpoint, visual_split, skip, x_axis_t_breaks){

  pi_bounds = tibble(t = seq(0, max(output$possible_data$t), len = 300),
                     pi2 = predict(output$pi_model, newdata = data.frame(t = t), type = "response"),
                     pi1 = 1 - pi2,
                     pi_se = predict(output$pi_model, newdata = data.frame(t = t), se.fit = TRUE)$se.fit,
                     pi_1_lp = logit(1 - predict(output$pi_model, newdata = data.frame(t = t), type = "response")),
                     pi_2_lp = predict(output$pi_model, newdata = data.frame(t = t)),
                     pi_1_lb = inverse_logit(pi_1_lp - (1.96 * pi_se)),
                     pi_1_ub = inverse_logit(pi_1_lp + (1.96 * pi_se)),
                     pi_2_lb = inverse_logit(pi_2_lp - (1.96 * pi_se)),
                     pi_2_ub = inverse_logit(pi_2_lp + (1.96 * pi_se))
  )


  pi = df %>%
    offset_time_as_date_in_df(., start_date) %>%
    ggplot(aes(x = t)) +
    geom_line(aes(x = offset_time_as_date(t, start_date), y = pi1, color = "Component 1 Proportion", linetype = "Component 1 Proportion"), data = pi_bounds) +
    geom_line(aes(x = offset_time_as_date(t, start_date), y = pi2, color = "Component 2 Proportion", linetype = "Component 2 Proportion"), data = pi_bounds) +
    #geom_function(fun = function(t){(1 - predict(output$pi_model, newdata = data.frame(t = as_offset_time(x = t, start_date)), type = "response"))}, aes(color = "Component 1 Proportion", linetype = "Fitted Model")) +
    #geom_function(fun = function(t){predict(output$pi_model, newdata = data.frame(t = as_offset_time(x = t, start_date)), type = "response")}, aes(color = "Component 2 Proportion", linetype = "Fitted Model")) +
    #scale_color_manual(breaks = c("Component 1 Proportion", "Component 2 Proportion"), values = c("#e4190b", "#00BFC4"), name = "Component Prevalence") +
    ylim(0,1)  +
    xlab("Time") + ylab("Proportion of Isolates") + theme_minimal() +
    scale_color_manual(breaks = c("Component 1 Proportion", "Component 2 Proportion"), values = c("#e4190b", "#00999d"), labels = c(TeX(r'(Component 1 Prevalence: $\hat{\pi}_{1,t}$)'), TeX(r'(Component 2 Prevalence: $\hat{\pi}_{2,t}$)')), name = "Component Prevalence") +
    scale_linetype_manual(breaks = c("Component 1 Proportion", "Component 2 Proportion"), values = c(1, 1), labels = c(TeX(r'(Component 1 Prevalence: $\hat{\pi}_{1,t}$)'), TeX(r'(Component 2 Prevalence: $\hat{\pi}_{2,t}$)')), name = "Component Prevalence")


  ## now add for the extra LRs
  ##maybe an if here? tbd

  if(add_log_reg && (!is.na(ecoff) | (!is.na(s_breakpoint) & !is.na(r_breakpoint)) | !is.na(visual_split))){
    if(!is.na(s_breakpoint) & !is.na(r_breakpoint)){
      lr_output_bkpt = log_reg(output$possible_data, split_by = "r_breakpoint", data_type = "possible_data", drug = NULL, date_col = "t", date_type = "decimal", first_year = NULL, s_breakpoint = s_breakpoint, r_breakpoint = r_breakpoint)
    }
    if(!is.na(ecoff)){
      lr_output_ecoff = log_reg(output$possible_data, split_by = "ecoff", data_type = "possible_data", drug = NULL, date_col = "t", date_type = "decimal", first_year = NULL, ecoff = ecoff)
    }
    if(!is.na(visual_split)){
      lr_output_visual_split = log_reg(output$possible_data, split_by = "visual_split", data_type = "possible_data", drug = NULL, date_col = "t", date_type = "decimal", first_year = NULL, visual_split = visual_split)
    }

    pi = pi +
      ggnewscale::new_scale_color() +
      ggnewscale::new_scale("linetype")

    breaks_list = c()
    color_values_list = c()
    linetypes_list = c()




    if(!is.na(r_breakpoint) & (is.null(skip) | (!is.null(skip) && !("bkpts" %in% skip) ) ) ){
      pi_bounds = pi_bounds %>%
        mutate(.by = t,
               susceptible = 1 - predict(lr_output_bkpt, newdata = tibble(t = t), type = "response"),
               resistant = predict(lr_output_bkpt, newdata = tibble(t = t), type = "response") )

      pi = pi +
        geom_line(aes(x = offset_time_as_date(t, start_date), y = susceptible, color = "Susceptible", linetype = "Susceptible"), data = pi_bounds) +
        geom_line(aes(x = offset_time_as_date(t, start_date), y = resistant, color = "Resistant", linetype = "Resistant"), data = pi_bounds)

      breaks_list = breaks_list %>% append(c("Susceptible", "Resistant"))
      color_values_list = color_values_list %>% append(c("#7CAE00", "#C77CFF"))
      linetypes_list = linetypes_list %>% append(c(2,2))
    }

    if(!is.na(ecoff) & (is.null(skip) | (!is.null(skip) && !("ecoff" %in% skip) ) ) ){
      pi_bounds = pi_bounds %>%
        mutate(.by = t,
               wt = 1 - predict(lr_output_ecoff, newdata = tibble(t = t), type = "response"),
               nwt = predict(lr_output_ecoff, newdata = tibble(t = t), type = "response") )

      pi = pi +
        geom_line(aes(x = offset_time_as_date(t, start_date), y = wt, color = "WT (ECOFF)", linetype = "WT (ECOFF)"), data = pi_bounds) +
        geom_line(aes(x = offset_time_as_date(t, start_date), y = nwt, color = "NWT (ECOFF)", linetype = "NWT (ECOFF)"), data = pi_bounds)

      breaks_list = breaks_list %>% append(c("WT (ECOFF)", "NWT (ECOFF)"))
      color_values_list = color_values_list %>% append(c("#fcbf07", "#0211a3"))
      linetypes_list = linetypes_list %>% append(c(2,2))
    }

    if(!is.na(visual_split) ){
      pi_bounds = pi_bounds %>%
        mutate(.by = t,
               c1vs = 1 - predict(lr_output_visual_split, newdata = tibble(t = t), type = "response"),
               c2vs = predict(lr_output_visual_split, newdata = tibble(t = t), type = "response") )

      pi = pi +
        geom_line(aes(x = offset_time_as_date(t, start_date), y = wt, color = "Below Split", linetype = "Below Split"), data = pi_bounds) +
        geom_line(aes(x = offset_time_as_date(t, start_date), y = nwt, color = "Above Split", linetype = "Above Split"), data = pi_bounds)

      breaks_list = breaks_list %>% append(c("Below Split", "Above Split"))
      color_values_list = color_values_list %>% append(c("#DF4601", "#000000"))
      linetypes_list = linetypes_list %>% append(c(2,2))
    }



    pi = pi + scale_color_manual(breaks = breaks_list,
                                 values = color_values_list,
                                 name = "Logistic Regression Models") +  #+ guides(linetype = "none")
      scale_linetype_manual(
        values = c(linetypes_list),
        breaks = c(breaks_list),
        name = "Logistic Regression Models"
      )


  }

  pi = pi +
    geom_ribbon(aes(ymin = pi_1_lb, ymax = pi_1_ub, x = offset_time_as_date(t, start_date), fill = "Component 1 Proportion"), data = pi_bounds, alpha = 0.2) +
    geom_ribbon(aes(ymin = pi_2_lb, ymax = pi_2_ub, x = offset_time_as_date(t, start_date), fill = "Component 2 Proportion"), data = pi_bounds, alpha = 0.2) +
    scale_fill_manual(breaks = c("Component 1 Proportion", "Component 2 Proportion"), values = c("#e4190b", "#00999d"), labels = c(TeX(r'(Component 1 Prevalence: $\hat{\pi}_{1,t}$)'), TeX(r'(Component 2 Prevalence: $\hat{\pi}_{2,t}$)')), name = "Component Prevalence")

  if(!is.null(x_axis_t_breaks)){
    pi = pi + scale_x_continuous(breaks = x_axis_t_breaks %>% offset_time_as_date(., start_date = start_date), labels = x_axis_t_breaks %>% offset_time_as_date(., start_date = start_date) %>% year())

  }

  return(pi)

}
