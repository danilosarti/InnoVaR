#'generate_gamma_delta
#'
#'
#'Generates parameters lambda and delta of Bilinear terms. Imported from earlier versions of AMBARTI package.
#'
#' @param INDEX describes number of genotypes I and environments J
#' @param Q number of componentes captured by the bi linear term. (number of lambdas)
#'
#' @return a matrix n_genotypes by Q or n_environments by Q
#' @export
#'
#' @examples
#' n_genotypes<-4
#' n_environments<-3
#' generate_gamma_delta(I=n_genotypes,Q=2)
#' generate_gamma_delta(I=n_environments,Q=2)
generate_gamma_delta <- function(INDEX, Q) {

  first_row = TRUE
  while(first_row) {
    raw_par  = matrix(rnorm(INDEX*Q), ncol=Q)
    par_mean = matrix(rep(apply(raw_par,2,mean), each = nrow(raw_par)), ncol=Q)
    par_aux  = raw_par - par_mean

    # Constraints ----
    # apply(par_aux,2,sum)
    parTpar = solve(t(par_aux)%*%(par_aux))
    A       = square_root_matrix(parTpar)
    samples = par_aux%*%A

    # Force the first to be positive
    for (i in 1:nrow(samples)){
      row1 = samples[1, ]
      if (all(samples[i, ] > 0)) {
        aux = samples[i, ]
        samples[1,] = aux
        samples[i,] = row1
        return(samples)
      }
    }
  }
}
