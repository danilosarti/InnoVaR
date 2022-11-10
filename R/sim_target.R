# sim_target
#'
#' Simulate univariate target variables
#'
#' @param X_gene matrix of genomic variables
#' @param X_env  matrix of environmental variables
#' @param method a method used to generate the variable. It can be Lasso, Gilberg model or AMMI
#' @param marginal_mean Marginal mean of the target variable
#' @param marginal_sd Marginal Standard deviation of the target variable
#' @param pars  a list of parameters.
#'
#' @return
#'
#' @import glmnet rmutil checkmate jmuOutlier
#' @export
#'
#'
#'
#' @examples
#'
#'
#'
#'
sim_target <- function(X_gene,
                       X_env,
                       method = c("lasso", "ammi", "gillberg"),
                       marginal_mean = 0,
                       marginal_sd = 1,
                       pars = list()) {

  # Make sure the method argument matches only one value from the list above
  method <- match.arg(method)

  # Get the function
  sim_fun <- get(paste0("sim_", method))

  # Call the simulation function
  target_raw <- scale(sim_fun(X_gene = X_gene, X_env = X_env, pars = pars))[, 1]

  # Re-scale
  target <- target_raw * marginal_sd + marginal_mean

  return(list(target=target,
              method=method,
              marginal_mean=marginal_mean,
              marginal_sd=marginal_sd,
              X_gene=X_gene,
              X_env=X_env,
              pars=pars))

}


sim_lasso <- function(X_gene, X_env, pars = list(lambda = 1, sigma = 1,k=3)) {
    # Check its a data frame
    assert_data_frame(X_gene)
    assert_data_frame(X_env)
    ## create a random subset of gens and envs ...
    X <- cbind(X_gene, X_env)## next column would be a matrix x_interactions
    k=pars$k
    # Turn X into a design matrix
    mat1 <- model.matrix(~.^2, data = X)
    names_mat1=colnames(mat1)
    which_int=grep(":",names_mat1)
    mat2 <- mat1[, which_int]
    mat=cbind(mat1[,- which_int],
              mat2[,sample(ncol(mat2),k)])


    ## create the pcik interaction... pick coliumsn x an y...
    # Get the parameters from the list object
    lambda <- pars$lambda
    sigma <- pars$sigma

    # Simulate the coefficients
    beta <- rlaplace(n = ncol(mat), s = lambda)

    # Simulate the target and return it
    target <- rnorm(nrow(X), mean = mat %*% beta, sd = sigma)
    return(target)
    # Note: you might want to return the lambda, sigma, and beta values too if performing model identification
  }

# Version to simulate from AMMI
sim_ammi <- function(X_gene,
                     X_env,
                     pars = list(
                       Q = 1,
                       sigma_lambda = 1,
                       sigma_g = 1,
                       sigma_e = 1,
                       sigma_delta = 1,
                       sigma_gamma = 1,
                       sigma = 1
                     )) {

  # Make sure the data is just two categorical variables
  assert_factor(X_gene)
  assert_factor(X_env, len = length(X_gene))

  # Get the number of genotypes and environments
  n_g <- length(table(X_gene))
  n_e <- length(table(X_env))

  # Simulate the values
  Q <- pars$Q
  lambda <- sort(rnorm(Q, mean = 0, sd = pars$sigma_lambda), decreasing = TRUE)
  #delta <- matrix(rnorm(n_g * Q, mean = 0, sd = pars$sigma_delta), ncol = Q, nrow = n_g)
  #gamma <- matrix(rnorm(n_e * Q, mean = 0, sd = pars$sigma_gamma), ncol = Q, nrow = n_e)

  delta <- matrix(generate_gamma_delta_sigmas(INDEX=n_g, Q=Q,sigma=pars$sigma_delta), ncol = Q, nrow = n_g)
  gamma <- matrix(generate_gamma_delta_sigmas(INDEX=n_e, Q=Q,sigma=pars$sigma_gamma), ncol = Q, nrow = n_e)

  g <- rnorm(n_g, mean = 0, sd = pars$sigma_g)
  e <- rnorm(n_e, mean = 0, sd = pars$sigma_e)

  # Simulate the target and return
  n <- length(X_gene)
  blin <- rep(0, n)
  for (i in 1:n) {
    for (j in 1:Q) {
      blin[i] <- blin[i] + lambda[j] * delta[X_gene[i], j] * gamma[X_env[i], j]
    }
  }
  target <- rnorm(n, mean = g[X_gene] + e[X_env] + blin, sd = pars$sigma)
  #return(list(target=target,
  #        delta=delta,
  #       gamma=gamma,
  #   bilinear=blin))
  return(list(target=target,
              delta=delta,
              gamma=gamma, lambdas=lambda))
}

# Version to simulate from Gillberg
sim_gillberg <- function(X_gene, X_env, pars = list(
  R = 1,
  lambda_g0 = 1,
  lambda_e0 = 1,
  lambda_g = 1,
  lambda_e = 1,
  sigma_g0 = 1,
  sigma_e0 = 1,
  sigma_g = 1,
  sigma_e = 1,
  sigma = 1
)) {

  # This method follows (I think) the Gillberg et al paper in Bioinformatics in 2019

  # Create the kernel matrices - first make sure they're all numeric
  X_gene <- model.matrix(~ . - 1, data = X_gene)
  X_env <- model.matrix(~ . - 1, data = X_env)
  K_gene <- as.matrix(dist(X_gene))
  K_env <- as.matrix(dist(X_env))

  # Create the genomic effects
  n <- nrow(X_gene)
  a_gene <- rnorm(n, 0, pars$lambda_g0)
  g <- rnorm(n, K_gene %*% a_gene, pars$sigma_g0)

  # Create the environment effects
  a_env <- rnorm(n, 0, pars$lambda_e0)
  e <- rnorm(n, K_env %*% a_env, pars$sigma_e0)

  # Interaction effects
  R <- pars$R
  A_gene <- matrix(rnorm(n * R, 0, pars$lambda_g), ncol = R, nrow = n)
  A_env <- matrix(rnorm(n * R, 0, pars$lambda_e), ncol = R, nrow = n)

  # Get the matrix errors
  E_gene <- matrix(rnorm(n * R, 0, pars$sigma_g), ncol = R, nrow = n)
  E_env <- matrix(rnorm(n * R, 0, pars$sigma_e), ncol = R, nrow = n)

  # And finally the H matrices
  H_g <- K_gene %*% A_gene + E_gene
  H_e <- K_env %*% A_env + E_env

  # Simulate the target and return
  blin <- rep(0, n)
  for (i in 1:n) {
    blin[i] <- blin[i] + t(H_g[i, ]) %*% H_e[i, ]
  }
  target <- rnorm(n, mean = g + e + blin, sd = pars$sigma)
  return(target)
}

