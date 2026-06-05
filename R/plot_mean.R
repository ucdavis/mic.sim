#' Title
#'
#' @param output
#' @param df
#' @param results
#' @param start_date
#' @param fitted_comp
#' @param title
#' @param plot_min
#' @param plot_max
#' @param ecoff
#' @param s_breakpoint
#' @param r_breakpoint
#' @param visual_split
#' @param x_axis_t_breaks
#' @param assumed_components
#' @param n_fitted_components
#'
#' @keywords internal
#'
#' @return
#' @export
#'
#' @examples
plot_mean = function(output, df, results, start_date, fitted_comp, title, plot_min, plot_max, ecoff, s_breakpoint, r_breakpoint, visual_split, x_axis_t_breaks, assumed_components, n_fitted_components){
  if(n_fitted_components == 2){
    mean = make_2C_mean_plot(output = output, df = df, start_date = start_date, title = title, plot_min = plot_min, plot_max = plot_max)
    mean = mean_plot_add_splits(mean = mean, ecoff = ecoff, s_breakpoint = s_breakpoint, r_breakpoint = r_breakpoint, visual_split = visual_split)
    mean = readjust_scales_mean_plot(mean = mean, x_axis_t_breaks = x_axis_t_breaks, start_date = start_date)
  }else if(assumed_components == 2 & n_fitted_components == 1){
    mean = set_up_reduced_mean_plot(output = output, results = results, df = df, start_date = start_date, fitted_comp = fitted_comp, title = title, plot_min = plot_min, plot_max = plot_max)
    mean = mean_plot_add_splits(mean = mean, ecoff = ecoff, s_breakpoint = s_breakpoint, r_breakpoint = r_breakpoint, visual_split = visual_split)
    mean = readjust_scales_mean_plot(mean = mean, x_axis_t_breaks = x_axis_t_breaks, start_date = start_date)
  }else if(assumed_components == 1){
    mean = set_up_1C_mean_plot(output = output, df = df, start_date = start_date, fitted_comp = fitted_comp, title = title, plot_min = plot_min, plot_max = plot_max)
    mean = mean_plot_add_splits(mean = mean, ecoff = ecoff, s_breakpoint = s_breakpoint, r_breakpoint = r_breakpoint, visual_split = visual_split)
    mean = readjust_scales_mean_plot(mean = mean, x_axis_t_breaks = x_axis_t_breaks, start_date = start_date)
  }else{
    warningCondition("No components converged")
    mean = df %>% ggplot() +
      geom_segment(aes(x = t, xend = t, y = left_bound, yend = right_bound, color = cens), data = (df %>% filter(cens == "int")), alpha = 0.3) +
      geom_segment(aes(x = t, xend = t, y = right_bound, yend = left_bound, color = cens), data = (df %>% filter(cens == "lc") %>% mutate(left_bound = plot_min)), arrow = arrow(length = unit(0.03, "npc")), alpha = 0.3) +
      geom_segment(aes(x = t, xend = t, y = left_bound, yend = right_bound, color = cens), data = (df %>% filter(cens == "rc") %>% mutate(right_bound = plot_max)), arrow = arrow(length = unit(0.03, "npc")), alpha = 0.3) +
      geom_point(aes(x = t, y = left_bound,  color = cens), data = df %>% filter(left_bound != -Inf), alpha = 0.3) +
      geom_point(aes(x = t, y = right_bound,  color = cens), data = df %>% filter(right_bound != Inf), alpha = 0.3) +
      theme_minimal()
  }
  return(mean)
}


mean_plot_add_splits = function(mean, ecoff, s_breakpoint, r_breakpoint, visual_split){
  if((!is.na(ecoff) | (!is.na(s_breakpoint) & !is.na(r_breakpoint)) | !is.na(visual_split))){

    mean = mean +
      ggnewscale::new_scale_color()

    if(!is.na(s_breakpoint)){
      s_line = case_when(
        grepl("≤",s_breakpoint) ~ s_breakpoint %>% as.character() %>% parse_number() %>% log2,
        grepl("=",s_breakpoint) ~ s_breakpoint %>% as.character() %>% parse_number() %>% log2,
        grepl("<",s_breakpoint) ~ s_breakpoint %>% as.character() %>% parse_number() %>% log2 - 1,
        TRUE ~ s_breakpoint %>% as.character() %>% parse_number() %>% log2
      )

    }

    if(!is.na(r_breakpoint)){
      r_line = case_when(
        grepl("≥",r_breakpoint) ~ r_breakpoint %>% as.character() %>% parse_number() %>% log2 - 1,
        grepl("=",r_breakpoint) ~ r_breakpoint %>% as.character() %>% parse_number() %>% log2 - 1,
        grepl(">",r_breakpoint) ~ r_breakpoint %>% as.character() %>% parse_number() %>% log2,
        TRUE ~ r_breakpoint %>% as.character() %>% parse_number() %>% log2 - 1
      )
    }

    if(!is.na(ecoff)){
      ecoff_line = case_when(
        grepl("≤", ecoff) ~ ecoff %>% as.character() %>% parse_number() %>% log2,
        grepl("=", ecoff) ~ ecoff %>% as.character() %>% parse_number() %>% log2,
        grepl("<", ecoff) ~ ecoff %>% as.character() %>% parse_number() %>% log2 - 1,
        TRUE ~ ecoff %>% as.character() %>% parse_number() %>% log2
      )

      if(!is.na(s_breakpoint)){
        if(s_line == ecoff_line){
          s_line = s_line + 0.05
        }
      }
    }

    if(!is.na(visual_split)){
      visual_split_line = case_when(
        grepl("≤", visual_split) ~ visual_split %>% as.character() %>% parse_number() %>% log2,
        grepl("=", visual_split) ~ visual_split %>% as.character() %>% parse_number() %>% log2,
        grepl("<", visual_split) ~ visual_split %>% as.character() %>% parse_number() %>% log2 - 1,
        TRUE ~ visual_split %>% as.character() %>% parse_number() %>% log2
      )
      if(!is.na(ecoff)){
        if(visual_split_line == ecoff_line){
          visual_split_line = visual_split_line + 0.05
        }
      }
      if(!is.na(s_breakpoint)){
        if(s_line == visual_split_line){
          s_line = s_line + 0.05
        }
      }
    }


    included = tibble(
      breaks = c("Susceptible Breakpoint", "Resistant Breakpoint", "ECOFF", "Visual Split"),
      labels = c(TeX(paste0("Susceptible Breakpoint: ", s_breakpoint,r'($\mu$)',"g/mL")), TeX(paste0("Resistant Breakpoint: ", r_breakpoint,r'($\mu$)',"g/mL")), TeX(paste0("ECOFF: ", ecoff,r'($\mu$)',"g/mL")), TeX(paste0("Visual Split: ", visual_split,r'($\mu$)',"g/mL"))),
      values = c("#7CAE00", "#C77CFF", "darkred", "#DF4601"),
      linetypes = c(5,5,2,2),
      entries = c(s_breakpoint, r_breakpoint, ecoff, visual_split)

    ) %>% filter(!is.na(entries))

    if(!is.na(s_breakpoint)){
      mean = mean + geom_hline(aes(yintercept = s_line, color = "Susceptible Breakpoint", linetype =  "Susceptible Breakpoint"))
    }
    if(!is.na(r_breakpoint)){
      mean = mean + geom_hline(aes(yintercept = r_line, color = "Resistant Breakpoint", linetype =  "Resistant Breakpoint"))
    }
    if(!is.na(ecoff)){
      mean = mean + geom_hline(aes(yintercept = ecoff_line, color = "ECOFF", linetype =  "ECOFF"))
    }
    if(!is.na(visual_split)){
      mean = mean + geom_hline(aes(yintercept = visual_split_line, color = "Visual Split", linetype =  "Visual Split"))
    }

    mean = mean + scale_color_manual(
      breaks = included %>% pull(breaks),
      labels = included %>% pull(labels),
      values = included %>% pull(values),
      name = "Breakpoints and Cutoffs") +
      scale_linetype_manual(
        breaks = included %>% pull(breaks),
        labels = included %>% pull(labels),
        values = included %>% pull(linetypes),
        name = "Breakpoints and Cutoffs")  #+ guides(linetype = "none", color = "none")


  }else{
    mean = mean + guides(linetype = "none")
  }
  return(mean)
}

readjust_scales_mean_plot = function(mean, x_axis_t_breaks, start_date){

  if(!is.null(x_axis_t_breaks)){

    mean = mean + scale_x_continuous(breaks = x_axis_t_breaks %>% offset_time_as_date(., start_date = start_date), labels = x_axis_t_breaks %>% offset_time_as_date(., start_date = start_date) %>% year())

  }

  mean = mean + scale_y_continuous( breaks = function(limits) seq(floor(limits[1]), ceiling(limits[2]), by = 1), minor_breaks = NULL, name = latex2exp::TeX(r'($\log_2$ MIC (Fold Scale))'),
                                    sec.axis = sec_axis(
                                      trans = ~ ., # Transformation: multiply by 5
                                      labels = set_y_labels,
                                      name = TeX(r'(MIC (Logarithmic Spacing) [$\mu$g/mL])'),
                                      breaks = function(limits) seq(floor(limits[1]), ceiling(limits[2]), by = 1)
                                    ))
  return(mean)
}

set_up_1C_mean_plot = function(output, df, start_date, fitted_comp, title, plot_min, plot_max){
  ci_data <- tibble(t = rep(seq(0, max(output$possible_data$t), len = 300), 2)) %>%
    mutate(
      c1pred = predict(fitted_comp, tibble(t), se = T)$fit,
      c1pred_se = predict(fitted_comp, tibble(t), se = T)$se.fit,
      c1pred_lb = c1pred - 1.96 * c1pred_se,
      c1pred_ub = c1pred + 1.96 * c1pred_se
    )

  #fitted_comp$scale %>% print()
  mean <-
    df %>%
    offset_time_as_date_in_df(., start_date) %>%
    ggplot(aes(x = t)) +
    geom_segment(aes(x = t, xend = t, y = left_bound, yend = right_bound, color = "Observations"), data = (df %>% filter(cens == "int") %>% offset_time_as_date_in_df(., start_date)), alpha = 0.3) +
    geom_segment(aes(x = t, xend = t, y = right_bound, yend = left_bound, color = "Observations"), data = (df %>% filter(cens == "lc") %>% mutate(left_bound = plot_min) %>% offset_time_as_date_in_df(., start_date)), arrow = arrow(length = unit(0.03, "npc")), alpha = 0.3) +
    geom_segment(aes(x = t, xend = t, y = left_bound, yend = right_bound, color = "Observations"), data = (df %>% filter(cens == "rc") %>% mutate(right_bound = plot_max) %>% offset_time_as_date_in_df(., start_date)), arrow = arrow(length = unit(0.03, "npc")), alpha = 0.3) +
    geom_point(aes(x = t, y = left_bound,  color = "Observations"), data = df %>% filter(left_bound != -Inf) %>% offset_time_as_date_in_df(., start_date), alpha = 0.3) +
    geom_point(aes(x = t, y = right_bound,  color = "Observations"), data = df %>% filter(right_bound != Inf) %>% offset_time_as_date_in_df(., start_date), alpha = 0.3) +
    scale_colour_manual(values = c("Observations" = "#F8766D"), guide = "none") +
    ggnewscale::new_scale_color() +
    #        geom_function(fun = function(t){predict(fitted_comp, newdata = data.frame(t = as_offset_time(x = t, start_date)))}, aes(color = "Component 1 Mu", linetype = "Fitted Model")) +
    #        geom_function(fun = function(t){mu.se.brd.fms(t, z = 1.96)}, aes(color = "Component Mu", linetype = "Fitted Model SE"), size = 0.6, alpha = 0.6) +
    #        geom_function(fun = function(t){mu.se.brd.fms(t, z = -1.96)}, aes(color = "Component Mu", linetype = "Fitted Model SE"), size = 0.6, alpha = 0.6) +
    geom_line(aes(x = offset_time_as_date(t, start_date), y = c1pred, color = "Component 1 Mu"), data = ci_data) +
    geom_ribbon(aes(ymin = c1pred_lb, ymax = c1pred_ub, x = offset_time_as_date(t, start_date), fill = "Component 1 Mu"), data = ci_data, alpha = 0.2) +
    geom_ribbon(aes(ymin = lwr, ymax = upr, x = t, fill = "Component 1 Mu"), data = sim_pi_survreg_boot(df, fit = fitted_comp, alpha = 0.05, nSims = 10000) %>% offset_time_as_date_in_df(., start_date), alpha = 0.15) +
    scale_color_manual(breaks = c("Component 1 Mu"), values = c("#e4190b"), labels = c(TeX(r'(Component 1 Mean: $\hat{\mu}_{1,t}$)')), name = "Component Mean") +
    scale_fill_manual(breaks = c("Component 1 Mu"), values = c("#e4190b"), labels = c(TeX(r'(Component 1 Mean: $\hat{\mu}_{1,t}$)')), name = "Component Mean") +
    #scale_linetype_manual(values = c("Fitted Model" = 1), guide = "none") +
    ggtitle(title) +
    xlab("Time") +
    # ylab(TeX(r'(MIC ($\mu$g/mL))')) +
    ylim(plot_min - 1, plot_max + 1) +
    scale_y_continuous(breaks = scales::breaks_extended((plot_max - plot_min)/1.5)) +
    #scale_x_continuous(breaks = scales::breaks_extended(6)) +
    theme_minimal() +
    theme(legend.position = "bottom")
  return(mean)
}

set_up_reduced_mean_plot = function(output, results, df, start_date, fitted_comp, title, plot_min, plot_max){
  ci_data <- tibble(t = rep(seq(0, max(output$possible_data$t), len = 300), 2)) %>%
    mutate(
      c1pred = predict(fitted_comp, tibble(t), se = T)$fit,
      c1pred_se = predict(fitted_comp, tibble(t), se = T)$se.fit,
      c1pred_lb = c1pred - 1.96 * c1pred_se,
      c1pred_ub = c1pred + 1.96 * c1pred_se
    )
  if((results %>% filter(!dnc) %>% pull(c)) == 1){
    corresponding_color = "#e4190b"
    corresponding_label = TeX(r'(Component 1 Mean: $\hat{\mu}_{1,t}$)')
  }else{ #comp1
    corresponding_color = "#00999d"
    corresponding_label = TeX(r'(Component 2 Mean: $\hat{\mu}_{2,t}$)')
  }#comp2


  mean <- df %>%
    offset_time_as_date_in_df(., start_date) %>%
    ggplot(aes(x = t)) +
    scale_colour_gradient2(high = "#00BFC4", low = "#F8766D", mid = "green", midpoint = 0.5, name = "P(C=2|y,t)") +
    geom_segment(aes(x = t, xend = t, y = left_bound, yend = right_bound, color = `P(C=c|y,t)`), data = (df %>% filter(cens == "int" & c == "2") %>% offset_time_as_date_in_df(., start_date)), alpha = 0.2) +
    geom_segment(aes(x = t, xend = t, y = right_bound, yend = left_bound, color = `P(C=c|y,t)`), data = (df %>% filter(cens == "lc" & c == "2") %>% mutate(left_bound = plot_min) %>% offset_time_as_date_in_df(., start_date)), arrow = arrow(length = unit(0.03, "npc")), alpha = 0.2) +
    geom_segment(aes(x = t, xend = t, y = left_bound, yend = right_bound, color = `P(C=c|y,t)`), data = (df %>% filter(cens == "rc"& c == "2") %>% mutate(right_bound = plot_max) %>% offset_time_as_date_in_df(., start_date)), arrow = arrow(length = unit(0.03, "npc")), alpha = 0.2) +
    geom_point(aes(x = t, y = left_bound,  color = `P(C=c|y,t)`), data = (df %>% filter(left_bound != -Inf & c == "2") %>% offset_time_as_date_in_df(., start_date)), alpha = 0.2) +
    geom_point(aes(x = t, y = right_bound,  color = `P(C=c|y,t)`), data = (df %>% filter(right_bound != Inf & c == "2") %>% offset_time_as_date_in_df(., start_date)), alpha = 0.2) +
    ggnewscale::new_scale_color() +
    #geom_function(fun = function(t){predict(fitted_comp, newdata = data.frame(t = as_offset_time(x = t, start_date)))}, aes(color = "Component Mu", linetype = "Fitted Model")) +
    geom_line(aes(x = offset_time_as_date(t, start_date), y = c1pred, color = "Component Mu"), data = ci_data) +
    geom_ribbon(aes(ymin = c1pred_lb, ymax = c1pred_ub, x = offset_time_as_date(t, start_date), fill = "Component Mu"), data = ci_data, alpha = 0.2) +
    geom_ribbon(aes(ymin = lwr, ymax = upr, x = t, fill = "Component Mu"), data = sim_pi_survreg_boot(df, fit = fitted_comp, alpha = 0.05, nSims = 10000) %>% offset_time_as_date_in_df(., start_date), alpha = 0.15) +
    scale_color_manual(breaks = c("Component Mu"), values = c(corresponding_color), labels = c(corresponding_label), name = "Component Mean") +
    scale_fill_manual(breaks = c("Component Mu"), values = c(corresponding_color), labels = c(corresponding_label), name = "Component Mean") +


    #scale_colour_gradientn(colours = c("purple", "orange")) +
    #ylim(plot_min - 0.5, plot_max + 0.5) +
    ggtitle(title) +
    xlab("Time") +
    # ylab(TeX(r'(MIC ($\mu$g/mL))')) +
    ylim(plot_min - 1, plot_max + 1) +
    scale_y_continuous(breaks = scales::breaks_extended((plot_max - plot_min)/1.5)) +
    #scale_x_continuous(breaks = scales::breaks_extended(6)) +
    theme_minimal()
  return(mean)
}

make_2C_mean_plot = function(output, df, start_date, title, plot_min, plot_max){
  ci_data = get_two_comp_ci(output)

  mean <- df %>%
    offset_time_as_date_in_df(., start_date) %>%
    ggplot(aes(x = t)) +
    scale_colour_gradient2(high = "#00BFC4", low = "#F8766D", mid = "green", midpoint = 0.5, name = "P(C=2|y,t)") +
    #geom_point(aes(x = t, y = mid, color = `P(C=c|y,t)`), data = df %>% filter(c == "2"), alpha = 0) +
    geom_segment(aes(x = t, xend = t, y = left_bound, yend = right_bound, color = `P(C=c|y,t)`), data = (df %>% filter(cens == "int" & c == "2") %>% offset_time_as_date_in_df(., start_date)), alpha = 0.3) +
    geom_segment(aes(x = t, xend = t, y = right_bound, yend = left_bound, color = `P(C=c|y,t)`), data = (df %>% filter(cens == "lc" & c == "2") %>% mutate(plot_min) %>% offset_time_as_date_in_df(., start_date)), arrow = arrow(length = unit(0.03, "npc")), alpha = 0.3) +
    geom_segment(aes(x = t, xend = t, y = left_bound, yend = right_bound, color = `P(C=c|y,t)`), data = (df %>% filter(cens == "rc" & c == "2") %>% mutate(plot_max) %>% offset_time_as_date_in_df(., start_date)), arrow = arrow(length = unit(0.03, "npc")), alpha = 0.3) +
    geom_point(aes(x = t, y = left_bound,  color = `P(C=c|y,t)`), data = df %>% filter(left_bound != -Inf & c == "2") %>% offset_time_as_date_in_df(., start_date), alpha = 0.3) +
    geom_point(aes(x = t, y = right_bound,  color = `P(C=c|y,t)`), data = df %>% filter(right_bound != Inf & c == "2") %>% offset_time_as_date_in_df(., start_date), alpha = 0.3) +
    ggnewscale::new_scale_color() +

    #geom_bar(aes(x = mid, fill = cens)) +
    #geom_function(fun = function(t){predict(output$mu_model[[1]], newdata = data.frame(t = as_offset_time(x = t, start_date)))}, aes(color = "Component 1 Mu", linetype = "Fitted Model")) +
    #geom_function(fun = function(t){predict(output$mu_model[[2]], newdata = data.frame(t = as_offset_time(x = t, start_date)))}, aes(color = "Component 2 Mu", linetype = "Fitted Model")) +
    geom_line(aes(x = offset_time_as_date(t, start_date), y = c1pred, color = "Component 1 Mu"), data = ci_data) +
    geom_line(aes(x = offset_time_as_date(t, start_date), y = c2pred, color = "Component 2 Mu"), data = ci_data) +
    geom_ribbon(aes(ymin = c1pred_lb, ymax = c1pred_ub, x = offset_time_as_date(t, start_date), fill = "Component 1 Mu"), data = ci_data, alpha = 0.25) +
    geom_ribbon(aes(ymin = c2pred_lb, ymax = c2pred_ub, x = offset_time_as_date(t, start_date), fill = "Component 2 Mu"), data = ci_data, alpha = 0.25)
  if(attr(df, "model") != "mgcv"){
    mean = mean +
      geom_ribbon(aes(ymin = lwr, ymax = upr, x = t, fill = "Component 1 Mu"), data = sim_pi_survreg_boot(df, fit = output$mu_model[[1]], alpha = 0.05, nSims = 10000) %>% offset_time_as_date_in_df(., start_date), alpha = 0.15) +
      geom_ribbon(aes(ymin = lwr, ymax = upr, x = t, fill = "Component 2 Mu"), data = sim_pi_survreg_boot(df, fit = output$mu_model[[2]], alpha = 0.05, nSims = 10000) %>% offset_time_as_date_in_df(., start_date), alpha = 0.15) +
      scale_fill_manual(breaks = c("Component 1 Mu", "Component 2 Mu"), values = c("#e4190b", "#00999d"), labels = c(TeX(r'(Component 1 Mean: $\hat{\mu}_{1,t}$)'), TeX(r'(Component 2 Mean: $\hat{\mu}_{2,t}$)')), name = "Component Means")
  }
  mean = mean + scale_color_manual(breaks = c("Component 1 Mu", "Component 2 Mu"), values = c("#e4190b", "#00999d"), labels = c(TeX(r'(Component 1 Mean: $\hat{\mu}_{1,t}$)'), TeX(r'(Component 2 Mean: $\hat{\mu}_{2,t}$)')),name = "Component Means") +
    ggtitle(title) +
    xlab("Time") +
    # ylab(TeX(r'(MIC ($\mu$g/mL))')) +
    ylim(plot_min - 1, plot_max + 1) +
    scale_y_continuous(breaks = scales::breaks_extended((plot_max - plot_min)/1.5)) +
    theme_minimal()

  return(mean)

}

