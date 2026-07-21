#' simulate_mics
#'
#' function that wraps together all the functions that determine t's distribution,
#' pi and its trends, trends in the mean, component draws, epsilon, covariates,
#' and censors the data
#'
#' @param n Number of observations
#' @param t_dist A function of n for drawing values of t
#' @param pi A function of time that returns a vector of weights that sum to 1.
#' @param `E[X|T,C]` A function of time and component that returns a value of mu (component mean) for any given time and component
#' @param sd_vector A vector with length equal to the number of components, with the elements named "1", "2",...
#' @param covariate_list List of covariates, each one has its own format, see examples of numeric and categorical covariates
#' @param covariate_effect_vector Vector of covariate effects corresponding to the covariates listed above
#' @param conc_limits_table If concentration limits vary by some covariate use this table to specify limits for each value of the covariate. Is right-joined to data by the covariate values
#' @param low_con If concentration limits are constant for all observations, used to set the lowest tested concentration on the log2(MIC) scale
#' @param high_con If concentration limits are constant for all observations, used to set the highest tested concentration on the log2(MIC) scale
#' @param scale What scale ("log" or "MIC") the data returned by simulate_mics is. Default is "log" which corresponds to log2(MIC)
#'
#' @return
#' @export
#'
#' @importFrom dplyr tibble mutate inner_join group_by case_when
#' @importFrom magrittr %>%
#' @importFrom purrr map as_vector map_chr map_dbl
#' @importFrom tibble tibble
#' @import dplyr
#'
#' @examples
#' #Covariate List
#' covariate_list = list(c("numeric", "normal", 0, 1), c("categorical", c(0.3, 0.4)), c("numeric", "uniform", 0, 5))
#' #Covariate Effect Vector
#' covariate_effect_vector = c(2, #intercept for all covariates combined
#'                                                                 10, #slope for covariate_1
#'                                                                 100, #effect of level b vs a of covariate 2
#'                                                                 3 #slope for covariate_3
#'                                                                 )
#' #Concentration Limits Table
#' conc_limits_table = tibble::as_tibble(rbind(c("a", -3, 3),
#'                                     c("b", -4, 4),
#'                                     c("c", -4, 4)),`.name_repair` = "unique"
#' ) %>% dplyr::rename("covariate_2" = 1, "low_cons" = 2, "high_cons" = 3)
#'
#' simulate_mics(
#' n = 300,
#' t_dist = function(n){runif(n, min = 0, max = 16)},
#' pi = function(t) {
#'   z <- 0.17 + 0.025 * t - 0.00045 * t ^ 2
#'   tibble::tibble("1" = 1 - z, "2" = z)
#' },
#' `E[X|T,C]` = function(t, c)
#' {
#'   dplyr::case_when(c == "1" ~ -4.0 + (0.24 * t) - (0.0055 * t ^ 2),
#'             c == "2" ~ 3 + 0.001 * t,
#'             TRUE ~ NaN)
#' },
#' sd_vector = c("1" = 1, "2" = 1.05),
#' covariate_list = covariate_list,
#' covariate_effect_vector = covariate_effect_vector,
#' conc_limits_table = conc_limits_table,
#' scale = "log")
#'
#'
simulate_mics <- function(n = 300,
                          t_dist = function(n){runif(n, min = 0, max = 16)},
                          pi = function(t) {
                            z <- 0.17 + 0.025 * t - 0.00045 * t ^ 2
                            tibble("1" = 1 - z, "2" = z)
                          },
                          `E[X|T,C]` = function(t, c)
                          {
                            case_when(c == "1" ~ -4.0 + (0.24 * t) - (0.0055 * t ^ 2),
                                      c == "2" ~ 3 + 0.001 * t,
                                      TRUE ~ NaN)
                          },
                          sd_vector = c("1" = 1, "2" = 1.05),
                          covariate_list = NULL,
                          covariate_effect_vector = c(0),
                          conc_limits_table = NULL,
                          low_con = -3,
                          high_con = 6,
                          scale = "log"){
                          # covariate_list = list(c("numeric", "normal", 0, 1), c("categorical", c(0.3, 0.4, 0.3))),
                          # covariate_effect_vector = c(0, #intercept for all covariates combined
                          #                             0.2, #slope for covariate_1
                          #                             -1, 0.2), #effect of level b vs a of covariate 2, and level c vs a of covariate 2
                          # conc_limits_table = as_tibble(rbind(c("a", -3, 3),
                          #                                     c("b", -4, 4),
                          #                                     c("c", -4, 4)),`.name_repair` = "unique"
                          # ) %>% rename("covariate_2" = 1, "low_cons" = 2, "high_cons" = 3),) {
                          if (is.null(covariate_list)) {
                            base_data <- draw_epsilon(n, t_dist, pi, `E[X|T,C]`, sd_vector)
                            simulated_obs <-
                              base_data %>% mutate(observed_value = epsilon + x)
                            simulated_obs <-
                              simulated_obs %>% mutate(low_cons = low_con, high_cons = high_con)
                            censored_obs <-
                              censor_values(simulated_obs, #simulated_obs$observed_value, #low_con, high_con,
                                            #tested_concentrations,
                                            scale)
                            df <-
                              inner_join(
                                simulated_obs,
                                censored_obs,
                                by = join_by(
                                  t,
                                  p,
                                  comp,
                                  x,
                                  sd,
                                  epsilon,
                                  observed_value,
                                  low_cons,
                                  high_cons
                                )
                              ) %>% mutate(low_con = as.numeric(low_cons),
                                           high_con = as.numeric(high_cons))
                            attr(df, "scale") <- scale
                            return(df)
                          } else{
                            base_data <- draw_epsilon(n, t_dist, pi, `E[X|T,C]`, sd_vector)
                            covariate_data <-
                              add_covariate(covariate_list = covariate_list, input = base_data$t)
                            merged_data <- tibble(base_data, covariate_data)
                            total_cov_effect <-
                              covariate_effect_total(merged_data, covariate_effect_vector)
                            simulated_obs <- tibble(merged_data, total_cov_effect) %>%
                              mutate(observed_value = epsilon + total_cov_effect + x)
                            if (!is.null(conc_limits_table)) {
                              simulated_obs <- left_join(simulated_obs, conc_limits_table)
                            } else{
                              simulated_obs <-
                                simulated_obs %>% mutate(low_cons = low_con, high_cons = high_con)
                            }

                            censored_obs <- censor_values(simulated_obs,
                                                          #simulated_obs$observed_value, #low_con, high_con,
                                                          #tested_concentrations,
                                                          scale)
                            df <-
                              inner_join(simulated_obs, censored_obs) %>% mutate(low_con = as.numeric(low_cons),
                                                                                 high_con = as.numeric(high_cons)) %>%
                              select(-c(low_cons, high_cons))
                            attr(df, "scale") <- scale
                            return(df)
                          }
                          }

###LATER NEED TO ADD AN OPTION TO CHANGE CONCENTRATIONS BY TIME PERIOD
