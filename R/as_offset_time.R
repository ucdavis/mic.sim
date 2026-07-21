#' Convert Date to Elapsed Time
#'
#' @param x Date
#' @param start_date Initial timepoint, in decimal years
#'
#' @return
#' @keywords internal
#'
#' @importFrom lubridate decimal_date is.Date
#'
#' @examples \dontrun{ as_offset_time(lubridate::mdy("2/10/2010"), 2010) }
#'
as_offset_time = function(x, start_date){
  if(lubridate::is.Date(x) | lubridate::is.POSIXct(x)){
    t = lubridate::decimal_date(x) - start_date
  }else{
    t = x
  }
return(t)
}
