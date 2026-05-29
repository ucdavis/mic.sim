set_up_2C_mean_plot = function(){
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

}


add_splits_to_mean_plot = function(){
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
  labels = c(TeX(paste0("Susceptible Breakpoint: ", s_breakpoint,r'($\mu$)',"g/mL")), TeX(paste0("Resistant Breakpoint: ", r_breakpoint,r'($\mu$)',"g/mL")), TeX(paste0("ECOFF: ", ecoff,r'($\mu$)',"g/mL")), TeX(paste0("ECOFF: ", visual_split,r'($\mu$)',"g/mL"))),
  values = c("#7CAE00", "#C77CFF", "#ffd700", "#DF4601"),
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

}


adjust_scales_mean_plot = function(){

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
}
