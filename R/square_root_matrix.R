#'square_root_matrix
#'
#'Calculates the squared root of a matrix
#'
#' @param x a matrix
#' @references Hands-On Matrix Algebra Using R chapter 11 describes root matrix
#' @details the eigen vectors need of x need to  be positive.
#' @return x_sqrt the square root of a matrix.
#' @export
#'
#' @examples
#'
square_root_matrix <- function(x){

  # When Q = 1, x will be a scalar
  if (nrow(x) == 1) {return(sqrt(x))}

  # When Q > 1, then x will be a matrix
  if (nrow(x) > 1) {
    # Jordan normal form
    X = eigen(x)
    P = X$vectors
    A = diag(X$values)

    A_sqrt = diag(sqrt(X$values))
    P_inv = solve(P)
    x_sqrt = P %*% A_sqrt %*%  P_inv
    return(x_sqrt)
  }
}
