#' set_sim_copula
#'
#' @param d dimension of copula to be set. The same as the number of variables inside a module
#' @param lower_tri_corr lower triangular part of a valida correlation matrix
#' @param n_cont_var   number of continuous variables inside the module and which will be described by specific marginal distributions
#' @param cont_var_par a list of lists containing names and values of the parameters of the
#' @param n_unique number of sites for which the simulations should be made
#' @param mar_cont_dists a vector of strings containing strings with the names of the marginal distributions
#' @param var_names name of the variables inside the moduls
#'
#'@import copula
#' @return a list containing the copula settled for the simulations its mvdc component and the simulated continuos values. some of them will be siliced to generate categorical components
#' @export
#'
#' @examples
set_sim_copula=function(
  d=9,
  lower_tri_corr=soil_corr[lower.tri(soil_corr,diag=FALSE)],
  n_cont_var=ncol(soil_num),
  cont_var_par=list(list(lambda=40),list(df=30,ncp=3),
                    list(shape=40,scale=10),list(rate=30)) ,
  n_unique=length(unique(na.omit(soil_att_ex$siteid))),
  mar_cont_dists=c("pois","chisq","gamma","exp"),
  var_names=colnames(soil_att_ex)
){

  # d=ncol(gen_att)
  #                           lower_tri_corr = gen_i_corr[lower.tri(gen_i_corr,diag=FALSE)]
  #                           n_cont_var=0
  #                           cont_var_par = 0
  #                          n_sites=length(unique(na.omit(gen_att$siteid)))
  #
  #                         mar_cont_dists=c()
  #
  #          var_names=colnames(gen_att)

 n_sites=n_unique
  ## setting a copula.
  mycop<-normalCopula(lower_tri_corr,dim=d,dispstr="un")
  ### setting a list to contain marginals
  z<-vector("list", length = d)
  ## the i -n_cont_var categorical variables we set normal 0 distributions
  if(n_cont_var==0 & length(cont_var_par )== 0){
    for(i in 1:(length(z))){
      z[[i]]$mean<-0
      z[[i]]$sd<-1
    }
    mymvd<-mvdc(copula=mycop,margins   =c(rep("norm",d),mar_cont_dists),paramMargins=(z))
    r<-rMvdc(n_sites,mymvd) ## this will be our copulas simulated dataset which we will use to slice the categorical
    colnames(r)<-var_names

  }

  else{
    for(i in 1:(length(z)-n_cont_var)){
      z[[i]]$mean<-0
      z[[i]]$sd<-1
    }
    rest_of_list=length(z)-n_cont_var+1
    df=data.frame(j=rest_of_list:length(z),w=1:length(cont_var_par))
    for(k in 1:nrow(df)){
      z[df$j[k]]=cont_var_par[df$w[k]]
    }
    mymvd<-mvdc(copula=mycop,margins =c(rep("norm",d-n_cont_var),mar_cont_dists),paramMargins=(z))

    ## now we simulate from the copula
    r<-rMvdc(n_sites,mymvd) ## this will be our copulas simulated dataset which we will use to slice the categorical
    colnames(r)<-var_names
  }

  ## now we need to create a list containing name of       parameters and
  #head(r)
  return(list(my_cop=mycop,my_mvd=mymvd,simulated=data.frame(r)))
}
