#' import_mics
#'
#' @param mic_column
#' @param code_column
#' @param combination_agent
#'
#' @return
#' @keywords internal
#'
#' @importFrom dplyr mutate case_when
#' @importFrom magrittr %>%
#' @importFrom stringr str_remove_all
#' @importFrom readr parse_number
#'
#' @examples
import_mics = function(mic_column, code_column = NULL, combination_agent = NULL, log_reg_value = FALSE, scale = "log", round = FALSE, include_mic_bounds = FALSE){

  if(!is.null(combination_agent) && combination_agent == 2 & is.null(code_column)){
    code_column = tibble(mic_column) %>% mutate(
      code_column = dplyr::case_when(
        grepl(pattern = "(\u2264)|(<=)|(=<)", x = mic_column) ~ "<=",
        grepl(pattern = ">", x = mic_column) ~ ">",
        TRUE ~ NA
      )
    ) %>% pull(code_column)
  }


  if( !is.null(combination_agent) && combination_agent %in% c(1,2)){
    mic_column = stringr::str_split_i(mic_column, "/", combination_agent)
  }



  df_temp <- dplyr::tibble(mic_column, code_column)



  if(is.null(code_column)){
    df <- df_temp %>%
      mutate(left_bound =
               dplyr::case_when(
                 grepl(pattern = "(\u2264)|(<=)|(=<)", x = mic_column) ~ 0,
                 grepl(pattern = ">", x = mic_column) ~ readr::parse_number(mic_column),
                 TRUE ~ readr::parse_number(mic_column)/2
               ),
             right_bound =
               dplyr::case_when(
                 grepl(pattern = "(\u2264)|(<=)|(=<)", x = mic_column) ~ stringr::str_remove_all(mic_column, "[\u2264<=]"),
                 grepl(pattern = ">", x = mic_column) ~ "Inf",
                 TRUE ~ mic_column
               ) %>% as.numeric()
      )

  } else{
    df <- df_temp %>%
      mutate(left_bound =
               dplyr::case_when(
                 grepl(pattern = "(\u2264)|(<=)|(=<)", x = code_column) ~ 0,
                 grepl(pattern = ">", x = code_column) ~ readr::parse_number(mic_column),
                 TRUE ~ readr::parse_number(mic_column)/2
               ),
             right_bound =
               dplyr::case_when(
                 grepl(pattern = ">", x = code_column) ~ Inf,
                 TRUE ~ readr::parse_number(mic_column)
               )
      )

  }

  if(scale == "log"){
    df = df %>% mutate(
      left_bound_mic = left_bound,
      right_bound_mic = right_bound,
      left_bound = log2(left_bound),
      right_bound = log2(right_bound)
    ) %>% relocate(all_of(c("left_bound", "right_bound")), .before = everything())
  }

  if(scale == "log" & round){
    df = df %>%
      mutate(
        left_bound = round(left_bound),
        right_bound = round(right_bound)
      ) %>%
      relocate(all_of(c("left_bound", "right_bound")), .before = everything())
  }

  if(!include_mic_bounds){
    df = df %>%
      select(-c(left_bound_mic,
             right_bound_mic)
             )
  }


  if(log_reg_value){
    df = df %>% mutate(
      mic_column = paste0(code_column, mic_column),
      lr_column =
        case_when(
          grepl(pattern = "(\u2264)|(<=)|(=<)", x = mic_column) ~ parse_number(mic_column),
          grepl(pattern = ">", x = mic_column) ~ parse_number(mic_column) * 2,
          TRUE ~ parse_number(mic_column)
        )
    )
  }

  attr(df, "source") <- "imported"
  attr(df, "lr_col") <- log_reg_value
  attr(df, "mic_class") <- "imported_mic_column"
  attr(df, "metadata") <- FALSE
  attr(df, "scale") <- scale
  return(df)
}






