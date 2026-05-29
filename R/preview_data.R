#' Plot Data
#'
#' Produces a basic plot of data generated from simulate_mics or import_mics_with_metadata
#'
#' @param data Tibble or data frame from simulate_mics or import_mics_with_metadata
#' @param title Title of plot
#' @param y_min Minimum value for plot, some extra space will be added below so the lowest tested concentration is a reasonable value
#' @param y_max Maximum value for plot, some extra space will be added below so the highest tested concentration is a reasonable value
#' @param ECOFF
#' @param ECOFF_scale defaults to "MIC", meaning the ECOFF provided will be in concentration directly, if you have already taken log2(ECOFF) then change this to "log"
#' @param covariate string, name of a column in data
#' @param covariate_title what to name the legend for the covariate
#'
#' @importFrom ggnewscale new_scale_color
#'
#' @return
#' @export
#'
#' @examples
preview_data = function(data, title = "", y_min = NULL, y_max = NULL, ECOFF = NULL, ECOFF_scale = "MIC", covariate = NULL, covariate_title = "Legend"){

  if(is.null(y_max)){
    y_max = max(data$high_con)
  }

  if(is.null(y_min)){
    y_min = min(data$low_con)
  }

  if(!is.null(ECOFF) && !is.na(ECOFF)){
    ECOFF = parse_number(as.character(ECOFF))
      if(ECOFF_scale %in% c("MIC", "mic", "concentration")){
    ECOFF = round(log2(ECOFF))
      }

  if(ECOFF > y_max){
    y_max = ECOFF
  }
  if(ECOFF < y_min){
    y_min = ECOFF
  }
}

  if(is.null(covariate)){
    plot = data %>%
      ggplot() +
      geom_segment(aes(x = t, xend = t, y = left_bound, yend = right_bound), color = "black", data = (. %>% filter(left_bound != -Inf & right_bound != Inf)), alpha = 0.2) +
      geom_segment(aes(x = t, xend = t, y = right_bound, yend = left_bound), color = "black", data = (. %>% filter(left_bound == -Inf) %>% mutate(left_bound = y_min - 2)), arrow = arrow(length = unit(0.03, "npc")), alpha = 0.2) +
      geom_segment(aes(x = t, xend = t, y = left_bound, yend = right_bound), color = "black", data = (. %>% filter(right_bound == Inf) %>% mutate(right_bound = y_max + 2)), arrow = arrow(length = unit(0.03, "npc")), alpha = 0.2) +
      geom_point(aes(x = t, y = left_bound), color = "black", data = . %>% filter(left_bound != -Inf), alpha = 0.2) +
      geom_point(aes(x = t, y = right_bound), color = "black", data = . %>% filter(right_bound != Inf), alpha = 0.2) +
      ggtitle(title) + ylab(TeX(r'(MIC ($\mu$g/mL))')) + xlab("Time")
  }else if(covariate %in% colnames(data)){
    data = data %>% rename(group = covariate)

    three_and_up = viridisLite::viridis(length(unique(data$group)) + 1 , option = "A")[-(length(unique(data$group)) + 1)]
    two = c("#0E3386", "#CC3433")

    hex_colors = case_when(length(unique(data$group)) > 2 ~ three_and_up, TRUE ~ two)

    plot = data %>%
      ggplot() +
      geom_segment(aes(x = t, xend = t, y = left_bound, yend = right_bound, color = group), data = (. %>% filter(left_bound != -Inf & right_bound != Inf)), alpha = 0.2) +
      geom_segment(aes(x = t, xend = t, y = right_bound, yend = left_bound, color = group), data = (. %>% filter(left_bound == -Inf) %>% mutate(left_bound = low_con - 2)), arrow = arrow(length = unit(0.03, "npc")), alpha = 0.2) +
      geom_segment(aes(x = t, xend = t, y = left_bound, yend = right_bound, color = group), data = (. %>% filter(right_bound == Inf) %>% mutate(right_bound = high_con + 2)), arrow = arrow(length = unit(0.03, "npc")), alpha = 0.2) +
      geom_point(aes(x = t, y = left_bound, color = group), data = . %>% filter(left_bound != -Inf), alpha = 0.2) +
      geom_point(aes(x = t, y = right_bound, color = group), data = . %>% filter(right_bound != Inf), alpha = 0.2) +
      ggtitle(title)  + xlab("Time") +
      scale_color_discrete(name = covariate_title, type = hex_colors)
  }else{
    errorCondition("covariate should be a string and the name of a column in data")
  }

  if(!is.null(ECOFF)){

    plot = plot +
      ggnewscale::new_scale_color() +
      geom_hline(aes(yintercept = ECOFF, color = "ECOFF")) +
      scale_color_manual(values = c("ECOFF" = "darkorange"), name = NULL) +
      ylim(y_min - 2, y_max + 2) %>% suppressWarnings()




  }else{
    plot = plot +
      ylim(y_min - 2, y_max + 2) %>% suppressWarnings()
  }

  plot = plot + #scale_y_continuous(labels = set_y_labels, breaks = function(limits) seq(floor(limits[1]), ceiling(limits[2]), by = 1), minor_breaks = NULL)
    scale_y_continuous( breaks = function(limits) seq(floor(limits[1]), ceiling(limits[2]), by = 1), minor_breaks = NULL, name = latex2exp::TeX(r'($\log_2$ MIC (Fold Scale))'),
                        sec.axis = sec_axis(
                          trans = ~ ., # Transformation: multiply by 5
                          labels = set_y_labels,
                          name = TeX(r'(MIC (Logarithmic Spacing) [$\mu$g/mL])'),
                          breaks = function(limits) seq(floor(limits[1]), ceiling(limits[2]), by = 1)
                        )) +
    theme_bw() +
    theme(legend.position = "bottom")

  plot %>% return()
}
