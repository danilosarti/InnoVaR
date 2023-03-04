#' assure_corr
#'
#' @param corr a matrix with values between -1 and 1 to be transforme in a proper correlation matrix
#' @param n_var number of variables inside a module
#' @return
#' @export
#'
#' @examples
assure_corr<-function(n_var=n_var,corr=corr){
  n<-n_var
  # generate a 5 x 5 matrix
  ## how to add here the part of creating correlations based on distances.
  A <- corr
  A[!upper.tri(A)] <- 0

  ## add to its transpose
  B <- A + t(A)

  ## look at smallest eigenvalue of the matrix
  ## if > - 1, add the identity and it is fine
  ## if < - 1, multiply by 1/(|smallest eigenvalue| + epsilon)
  min_eigen <- min(eigen(B)$values)
  min_eigen

  D <- B * 1/(abs(min_eigen) + .01)

  ## now add the identity and we have a valid correlation matrix
  assured_corr_mat <- D + diag(rep(1, n))
  assured_corr_mat=as.matrix(assured_corr_mat)
  return(assured_corr_mat)
}






