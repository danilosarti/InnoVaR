% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_target.R
\name{sim_target}
\alias{sim_target}
\title{simulate univariate target variables}
\usage{
sim_target(
  X_gene,
  X_env,
  method = c("lasso", "ammi", "gillberg"),
  marginal_mean = 0,
  marginal_sd = 1,
  pars = list()
)
}
\arguments{
\item{X_gene}{matrix of genomic variables}

\item{X_env}{matrix of environmental variables}

\item{method}{a method used to generate the variable. It can be Lasso, Gilberg model or AMMI}

\item{marginal_mean}{Marginal mean of the target variable}

\item{marginal_sd}{Marginal Standard deviation of the target variable}

\item{pars}{a list of parameters.}
}
\description{
simulate univariate target variables
}
