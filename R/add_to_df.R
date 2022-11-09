#' add_to_df
#'
#' Adds a column with number of elements different than the number of rows of the data frame we want to fill the column. cbind does not work in this scenario.
#'
#'
#'
#' @param df  a data frame where we want to include a columns probably with dimensions different than the ones the data.frame allows
#' @param v   a vector containing values to be added as column
#'
#' @return
#' @export
#'
#' @examples
add_to_df <- function(df, v){
  nRow <- nrow(df)
  lngth <- length(v)
  if(nRow > lngth){
    length(v) <- nRow
  }else if(nRow < lngth){
    df[(nRow+1):lngth, ] <- NA
  }
  cbind(df,v)
}
