#' module_coords
#'
#' creates a data set named coords that describes the association to be considered between variables of a module to be simulated, and details to map such information into a matrix.
#'
#' @param n_var  number of the variables in the module  to be simulated
#' @param which_neg string with coordinates of the variables which are negativelly correlated.
#' @param dis_neg_know_ass Distances describing negative assocaitions known a priori
#' @param coo_known_neg_ass Coordinates that describe the positive association between two variables ex 1_2 indicates association between first and second variable
#' @param dis_pos_know_ass Distances describing positive assocaitions known a priori
#' @param coo_known_pos_ass Coordinates that describe the positive association between two variables ex 1_2 indicates association between first and second variable
#' @param by  parameter like in the sample function used to generate randomly unknown positive and negative associations
#'
#' @return a list with an coords a data set named coords that describes the association to be considered between variables of a module to be simulated, and details to map such information into a matrix.
#' @export
#'
#' @examples module_coords(n_var=9,which_neg=c("2_1","6_2","1_9","1_6","2_3"),
#' dis_neg_know_ass=c(-0.9300003, -0.9310004, -0.8684074),
#' coo_known_neg_ass=c("2_1","6_2","1_9"),
#' dis_pos_know_ass=c(0.5235988, 0.6435011, 1.1197695),
#' coo_known_pos_ass=c("3_9","4_5","6_8"),by=0.001)
module_coords<-function(
  # n_var=9,
  #                       which_neg=c("2_1","6_2","1_9","1_6","2_3"),
  #                       dis_neg_know_ass=c(-0.9300003, -0.9310004, -0.8684074),
  #                       coo_known_neg_ass=c("2_1","6_2","1_9"),
  #                       dis_pos_know_ass=c(0.5235988, 0.6435011, 1.1197695),
  #                       coo_known_pos_ass=c("3_9","4_5","6_8"),
  #                       by=0.001,
  n_var=ncol(gen_att),
  which_neg = c("1_2","4_5","2_3","20_4","30_3"),
  dis_neg_know_ass =known_neg_dis_gen,
  coo_known_neg_ass=c("4_5","2_3","20_4","30_3"),
  dis_pos_know_ass =known_pos_dis_gen,
  coo_known_pos_ass = c("3_9","4_5","6_8","7_8"),
  by=0.0001

) {
  ## the user needs to inform which pairs should be mapped
  ## mirroring the informed coordinates
  mirroed_negs=mirror_string(which_neg)
  mirroed_coo_known_neg_ass=mirror_string(coo_known_neg_ass)
  mirroed_coo_known_pos_ass=mirror_string(coo_known_pos_ass)
  ## mirroring the informed know dis asso
  mirroed_dis_neg_know_ass=c(dis_neg_know_ass,dis_neg_know_ass)
  mirroed_dis_pos_know_ass=c(dis_pos_know_ass,dis_pos_know_ass)
  ## auxiliary dataframes for known negatives
  df_aux_neg=data.frame(coord_ij=mirroed_coo_known_neg_ass,aux_neg_dis=mirroed_dis_neg_know_ass)
  ## auxiliary dataframes for known positives
  df_aux_pos=data.frame(coord_ij=mirroed_coo_known_pos_ass,aux_pos_dis=mirroed_dis_pos_know_ass)
  #### creating coords data frame to store information about associations and their mapping into a matrix
  x<-c(1:n_var)
  coords<-expand.grid(i=x,j=x,pos_neg=NA)
  coords$coord_ij<-factor(paste0(coords$i,"_",coords$j))
  ## declaring which cartesian products will describe the negative correlations
  ## we need to be sure to include at least the combinations which will be
  ## in the lower triangular part withouth the diagonal of the matrix. ## copulas package requires this structure.
  ## including a columns to be completed later with trigonometric distances
  coords$distance<-NA
  ## create a mirror for the informed negatives.
  which_negative<-mirroed_negs##
  # declaring the negative associations
  coords[coords$coord_ij %in% which_negative,]$pos_neg="neg"
  coords$pos_neg=ifelse(is.na(coords$pos_neg),"pos","neg")
  # informing the known neg ass distances via its coord_ij
  # we have the coordinates
  #coords[coords$coord_ij %in%mirroed_coo_known_neg_ass,]$distance=mirroed_dis_neg_know_ass
  # updating known neg distances
  for( i in 1: nrow(df_aux_neg)){
    coords[coords$coord_ij == df_aux_neg$coord_ij[i],]$distance=df_aux_neg$aux_neg_dis[i]
  }
  # updating known pos distances
  for( i in 1: nrow(df_aux_pos)){
    coords[coords$coord_ij == df_aux_pos$coord_ij[i],]$distance=df_aux_pos$aux_pos_dis[i]
  }
  ####coords[coords$coord_ij %in% df_aux_pos$coord_ij,]
  ### mapping position in matrices
  coords$lower=ifelse(coords$i>coords$j,"low_tri","non_low_tri")
  coords$diag=ifelse(coords$i==coords$j,"diag","non_diag")
  coords$upper=ifelse(coords$i<coords$j,"uper_tri","non_uper_tri")
  ## dealing with not known distances.
  coords$known_priori=NA
  coords$known_priori=ifelse(is.na(coords$distance),"not_know","known")
  coords[(coords$known_priori=="not_know"&coords$pos_neg=="neg"),]
  ## ordering the coords data frame by i column

  #### generating negative now known a priori distances

  coords_fil_not_k_neg=coords[(coords$known_priori=="not_know"&coords$pos_neg=="neg"&coords$upper=="uper_tri"),]
  cfnkn_or=coords_fil_not_k_neg[order(coords_fil_not_k_neg$i),]
  size_unknou_neg_dis=nrow(cfnkn_or)
  cfnkn_or$distance=sample(seq(to=(3*(pi)/2), from=pi, by=by),size=size_unknou_neg_dis)
  # we need to mirror the strings of the ordered data set
  mirr_cfnkn_or_coord_ij=mirror_string( as.vector(cfnkn_or$coord_ij))
  # we need to mirror the distances of the ordered data set
  mirr_unk_neg_dis=c(cfnkn_or$distance,cfnkn_or$distance)
  # we need to create a data frame with the distances random and the mirroed coordinates to loop through
  aux_df_un_neg=data.frame(coords_ij=mirr_cfnkn_or_coord_ij,distance=mirr_unk_neg_dis)
  ## we need to update the coords object with the unknown a priori random distances generated
  # updating unknown random_neg distances
  for( i in 1: nrow(aux_df_un_neg)){
    coords[coords$coord_ij == aux_df_un_neg$coords_ij[i],]$distance=aux_df_un_neg$distance[i]
  }


  ### now we deal with the positive unknown distances that need to be sampled.
  coords_fil_not_k_pos=coords[(coords$known_priori=="not_know"&coords$pos_neg=="pos"&coords$upper=="uper_tri"),]
  cfnkp_or=coords_fil_not_k_pos[order(coords_fil_not_k_pos$i),]
  size_unknou_pos_dis=nrow(cfnkp_or)
  cfnkp_or$distance=sample(seq(to=((pi)/2), from=0, by=by),size=size_unknou_pos_dis)
  # we need to mirror the strings of the ordered data set
  mirr_cfnkp_or_coord_ij=mirror_string( as.vector(cfnkp_or$coord_ij))
  # we need to mirror the distances of the ordered data set
  mirr_unk_pos_dis=c(cfnkp_or$distance,cfnkp_or$distance)
  # we need to create a data frame with the distances random and the mirroed coordinates to loop through
  aux_df_un_pos=data.frame(coords_ij=mirr_cfnkp_or_coord_ij,distance=mirr_unk_pos_dis)
  ## we need to update the coords object with the unknown a priori random distances generated
  # updating unknown random_pos distances
  for( i in 1: nrow(aux_df_un_pos)){
    coords[coords$coord_ij == aux_df_un_pos$coords_ij[i],]$distance=aux_df_un_pos$distance[i]
  }

  ## now we need to insert trig distances in a way to fill the distances between x1 and x1 ,,,.
  coords$distance=ifelse(is.na(coords$distance),1.570796,coords$distance)
  coords[order(coords$i),][,c(3,4,5,9)]
  ## inserting a columns with the sins to be used in matrices later
  coords$sin=sin(coords$distance)
  return(list(coords=coords, which_negative=which_negative))
}
