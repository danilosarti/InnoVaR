#'mirror_string
#'
#'Makes a mirror of strings of the form a_b
#'
#' @param string a string of the form a_b which needs to be mirrored in the for b_a. both a_b and b_a form a complete mirroed vector in the end.
#'
#' @return a mirrored string
#' @export
#' @import stringr
#'
#' @examples mirror_string(c("2_1","1_3","10_2","100_4"))
#'
mirror_string=function(string){

  mirrored=c()
  for( i in 1:length(string)){
    row=str_extract(string[i], "[^_]+")
    j_mirror=row
    j=str_extract(string[i], "[^_]*$")
    row_mirror=j
    mirrored[i]=paste0(row_mirror, "_", j_mirror)
  }
  final_string=c(string,(mirrored))
  return(final_string)
}


