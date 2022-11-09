#' calculate_trig_dis
#'
#' Calculates trigonometric distances based on know correlations
#'
#' @param known_cor Values of known correlation to be transformed into trigonometric distances
#'
#' @return calculated trigometric distances
#' @export
#'
#' @examples  calculate_trig_dis(c(0.80,-0.99))
calculate_trig_dis=function(known_cor){
  trig_dist=asin(known_cor)
  return(trig_dist)

}
