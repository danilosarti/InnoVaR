#' cut_levels
#'
#'
#' This function creates a vector of levels based on a continuous vector. It allows to force all the levels desired for some of the variables respecting the correlation of the variable with other variable in a module.
#'
#' @param categorical categorical variables in a module
#' @param cat_inc_all_levels categorical variables which all levels must be in final simulations
#' @param levels_inc_all levels of above
#' @param continuous_to_cut a matrix containing values simulated from copulas
#' @param corrs a coorelation matrix
#'
#' @return
#' @export
#'
#' @examples
cut_levels=function(categorical=soil_cat,
                    cat_inc_all_levels=c("siteid","yes_no_sample"),
                    levels_inc_all=list(siteid=c(1:16),yes_no_sample=c(1:2))
                    , continuous_to_cut=sim_data_cop$simulated,
                    corrs=soil_corr
){

   # categorical=gen_att
   #  cat_inc_all_levels=c("genotype_id")
   #  levels_inc_all=list(genotype=c(1:n_genotypes))
   #  continuous_to_cut=sim_data_cop_gen_i$simulated
   #  corrs=gen_i_corr

  ## copulas stuff returned as a matrix so we need to create a df
  continuous_to_cut=data.frame(continuous_to_cut)
  ## we need consider variables which does not need to have all levels of factor included when the slicing process is made
  cat_not_all=categorical %>% select(-all_of(cat_inc_all_levels))
  # we need to consider variables which all levels need to be included
  cat_all=categorical %>% select(all_of(cat_inc_all_levels))
  #   ## we need to separate the continuous according to columns which need to correspond to cat_all
  #    sim_data_all=continuous_to_cut %>%
  #      select(all_of(cat_inc_all_levels))
  # # same for not cat_all
  #    sim_data_n_all=continuous_to_cut %>%
  #      select(-all_of(cat_inc_all_levels))
  #   # creating a empty data frame to store the cutted levels.
  ## create and empty dataframe to insert cut stuff
  df_fin_cutdown_n_all<-data.frame()
  ##
  df_final<-data.frame(continuous_to_cut)

  ## loop over the non_all_to_include
  for(i in 1:ncol(cat_not_all)){
    #i=1 just to check loop
    # we need to informe the number of breaks to be conisidered for the column i of data set cat_not_all
    brks=nrow(na.omit((cat_not_all[colnames(cat_not_all)[i]])))
    ## we need now to create the labels to be used when cutting. labels are the same as the levels of col i in the data set cat_not_all
    labls=data.frame((unique(na.omit(cat_not_all[colnames(cat_not_all)[i]]))))[,1]
    ## we need now to cut the simulated thing according to the levels in labls
    final_cutdown_n_all=cut(df_final[colnames(cat_not_all)[i]][,1], breaks=brks,
                            labels=labls)
    df_fin_cutdown_n_all=add_to_df(df_fin_cutdown_n_all,final_cutdown_n_all)
  }
  names(df_fin_cutdown_n_all)=colnames(cat_not_all)

  ## we need now to cut the other variables but making sure that all the levels are included. We will need to use sample but correlation with other things must be respected.
  ## first we discretize the variables into unique levels
  # for ( j in cat_inc_all_levels ){
  #   df_aus_x=data.frame("j"=(levels_inc_all[j]))
  #   df_fin_cutdown_n_all[j]=sample(df_aus_x[[j]],nrow(df_aus_x),replace = FALSE)
  #
  # }

  ## create the all to include respecting correlation
  for( j in cat_inc_all_levels ){
    #j="genotype_id" used just to check functions do not uncomment
    var_to_sort=names(sort(corrs[,j][order(corrs[,j])],decreasing = T)[2])
    orderer_by_j=continuous_to_cut %>% arrange(var_to_sort)
    df_aus_x=data.frame("j"=(levels_inc_all[j]))
    orderer_by_j[j]=sprintf(paste0(j,"_%s"),seq(1:nrow(df_aus_x)))
    df_fin_cutdown_n_all[j]=orderer_by_j[j]
  }
  df_fin_cutdows_comp=df_fin_cutdown_n_all

  return(df_fin_cutdows_comp=df_fin_cutdows_comp)
}







