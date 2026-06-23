#' Convert Elapsed Time to Date
#'
#' @param t time in decimal years since start_date
#' @param start_date start date in decimal years
#'
#' @return
#' @keywords internal
#'
#' @importFrom lubridate date_decimal is.Date
#'
#' @examples \dontrun{ offset_time_as_date(2.109589, 2007) }
offset_time_as_date = function(t, start_date){
  if(!lubridate::is.Date(t)){
    x = lubridate::date_decimal(t + start_date)
  }else{
    x = t
  }
  return(x)
}
