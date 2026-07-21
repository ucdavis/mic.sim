#' Import MICs and Covariates
#'
#' @param data Data frame containing the MICs and covariates
#' @param mic_column String, name of column in data corresponding to the MIC values
#' @param metadata_columns Vector of column names (as strings) for covariates to be included in the data frame produced
#' @param code_column String, name of column containing any symbols for MICs (if data is in a 2 column format)
#' @param combination_agent Numerical, if the MIC is not a combination agent, 0. If combination agent must be separated by '/', use 1 to select the value before the '/', or 2 for the value after. Recommend a value where log2(value) is an integer.
#' @param log_reg_value Logical, TRUE if a column for logistic regression model should be included in output (MICs with `>` will be doubled, MICs with `<`, `<=`, or `≤` are halved)
#' @param scale String, "log" if MIC value should be converted to log2 scale (preferred for subsequent fitting of the model using the EM algorithm)
#' @param round Set to true if log2(MIC values) are integers, but decimal MIC values are rounded (e.g. 0.12 in place of 0.125)
#' @param include_mic_bounds Logical, if TRUE includes left and right boundaries of interval on MIC scale (in addition to on log2 scale if scale is "log")
#' @param low_con Numerical, the lowest concentration tested, specify on the same scale as the data. If left null and concentration_by_covariate table is not supplied, will be set based on the data
#' @param high_con Numerical, the highest concentration tested, specify on the same scale as the data. If left null and concentration_by_covariate table is not supplied, will be set based on the data
#' @param concentration_by_covariate Data frame. Table, columns "low_con" and "high_con" are matched to data using covariates. Include one row for each combination of covariates.
#'
#' @return
#' @export
#'
#' @importFrom dplyr any_of
#'
#' @examples
#'import_mics_with_metadata(data = tibble::tibble(MIC_A = c("≤0.12", ">16", 4, 2), t = runif(4, 0, 10)),
#'                          mic_column = "MIC_A",
#'                          metadata_columns = "t",
#'                          log_reg_value = TRUE,
#'                          scale = "log",
#'                          round = TRUE)
#'
#'import_mics_with_metadata(data = tibble::tibble(MIC_A = c(0.125, 16, 4, 2), code_A = c("≤", ">", NA, NA), t = runif(4, 0, 10)),
#'                          mic_column = "MIC_A",
#'                          metadata_columns = "t",
#'                          code_column = "code_A",
#'                          log_reg_value = FALSE,
#'                          scale = "log",
#'                          round = FALSE,
#'                          include_mic_bounds = TRUE)
#'
#'import_mics_with_metadata(data = tibble::tibble(MIC_A = c("≤10/1", ">80/8", "40/4", "20/2"), t = runif(4, 0, 10)),
#'                          mic_column = "MIC_A",
#'                          metadata_columns = "t",
#'                          combination_agent = 2,
#'                          log_reg_value = FALSE,
#'                          scale = "log",
#'                          round = FALSE)
#'
import_mics_with_metadata = function(data, mic_column, metadata_columns = NULL, code_column = NULL, combination_agent = 0, log_reg_value = FALSE, scale = "log", round = FALSE, include_mic_bounds = FALSE,
                                     low_con = NULL,
                                     high_con = NULL,
                                     concentration_by_covariate = NULL){
  mic_col = data %>% select(all_of(mic_column)) %>% pull() %>% as.character()
  if(!is.null(metadata_columns)){
    metadata_col = data %>% select(all_of(metadata_columns))
  }else{
    metadata_col = NULL
  }

  if(!is.null(code_column)){
    code_col = data %>% select(all_of(code_column)) %>% pull()
  }else{
    code_col = NULL
  }

  df = import_mics(mic_column = mic_col, code_column = code_col, combination_agent = combination_agent, log_reg_value = log_reg_value, scale = scale, round = round, include_mic_bounds = include_mic_bounds) %>%
    mutate(metadata_col)
  if(!is.null(metadata_col)){
    attr(df, "metadata") = TRUE
    attr(df, "metadata_names") = metadata_columns
  }

  df = df %>% mutate(obs_id = row_number()) %>% relocate(obs_id, .before = everything())

  if(!is.null(low_con) & !is.null(high_con)){
    df = df %>% mutate(
      low_con = ifelse(scale == "log", log2(low_con), low_con),
      high_con = ifelse(scale == "log", log2(high_con), high_con)
    )
  }else if(!is.null(concentration_by_covariate)){
    names = colnames(concentration_by_covariate)
    cov = subset(names, !names %in% c("low_con", "high_con"))
    df = left_join(df, concentration_by_covariate, .by = cov)
  }

  present = df %>% select(any_of(c("low_con", "high_con"))) %>% colnames()

  if(!"low_con" %in% present){
    df = df %>% mutate(low_con = case_when(
      min(left_bound, na.rm = TRUE) == -Inf ~ min(right_bound, na.rm = TRUE),
      TRUE ~ min(left_bound, na.rm = TRUE)
    ))
  }

  if(!"high_con" %in% present){
    df = df %>% mutate(high_con = case_when(
      max(right_bound, na.rm = TRUE) == Inf ~ max(left_bound, na.rm = TRUE),
      TRUE ~ max(right_bound, na.rm = TRUE)
    ))
  }


  return(df)
}
