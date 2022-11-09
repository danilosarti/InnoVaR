#' sym_mtx
#'
#' Create a symetric matrix with n_var rows and n_var columns
#'
#' @param vector a vector containing values to be used in the matrix creation
#' @param n_var number of variables
#' @param var_names name of the variables in the symmetric matrix
#'
#' @return a symmetric matrix
#' @export
#' @import Matrix
#'
#' @examples sym_mtx(vector=c(1:25), n_var=5,var_names=c("v1","v2","v3","v4","v5"))
sym_mtx=function(vector=1:25, n_var=5,var_names=c("v1","v2","v3","v4","v5")){
  m=matrix(vector,n_var,n_var)
  colnames(m)=var_names
  rownames(m)=var_names
  m_sym=forceSymmetric(m)
  return(m_sym)
}
