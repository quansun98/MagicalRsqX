#' Calculate PC metric quantifying the distance between two cohorts
#'
#' This function is used to calculate PC dissimilarity metric to quantify the dissimilarity between two cohorts in terms of the sample genetic profiles.
#' This metric is not symmetric. It measures how well the reference cohort can represent the target cohort.
#' Harmonized PCs should be calculated first. Please see https://github.com/covid19-hg/ for more details.
#' This metric will be standardized accounting for the number of variants used in calculation due to some missing variants.
#'
#' @param ref_pc Required. The plink2 .sscore file containing PC information for each sample in the reference. 
#' @param target_pc Required. The plink2 .sscore containing PC information for each sample in the target cohort.
#' @param ref_pc_col Optional. A vector of column numbers containing the PCs in the reference dataframe. We recommend use 10 or 20 PCs. 
#' Default will be 4:23 using 20 PCs starting from the 4th column, following the format in our provided .sscore file.
#' @param target_pc_col Optional. A vector of column numbers containing the PCs in the target dataframe. The numbers of columns should match with the reference.
#' Default will be 4:23 using 20 PCs starting from the 4th column, following the format in our provided .sscore file.
#' @param ref_nvar_col Optional. The column number of `ALLELE_CT` from plink2 .sscore file. Default 2.
#' @param target_nvar_col Optional. The column number of `ALLELE_CT` from plink2 .sscore file. Default 2.
#' @param ind_out Optional. Whether the individual distance (for each individual in target) to a reference will be written or not. Default FALSE.
#'
#' @return If ind_out == T, it will output a vector of each individual distance to the reference.
#' The order of individuals will follow the input order in `target_pc`.
#' If ind_out == F (default), it will only output a number of the PC distance between two cohorts.
#'
#' @examples
#'
#' dist = pc_dist(ref_pc = "../data/test_cohort1.txt", 
#' target_pc = "../data/test_cohort2.txt", ind_out = F)
#'
#' @importFrom matrixStats colMedians
#' @importFrom data.table fread
#'
#' @export

pc_dist = function(ref_pc, target_pc, ref_pc_col = 4:23, target_pc_col = 4:23, ref_nvar_col = 2, target_nvar_col = 2, ind_out = F){

	target_pc = as.data.frame(fread(target_pc, header = T))
	ref_pc = as.data.frame(fread(ref_pc, header = T))

	n1 = dim(target_pc)[1]
	n2 = dim(ref_pc)[1]

	target_pc_std = as.data.frame(apply(target_pc[,target_pc_col], 2, function(x) (x/sqrt(target_pc[,target_nvar_col]))))
	ref_pc_std = as.data.frame(apply(ref_pc[,ref_pc_col], 2, function(x) (x/sqrt(ref_pc[,ref_nvar_col]))))

	ref_center = colMedians(as.matrix(ref_pc_std))

	dist_ind_to_ref = function(x, ref){
        	sum((as.vector(x)-ref)^2)
	}

	dist = apply(target_pc_std, 1, dist_ind_to_ref, ref = ref_center)

	if(ind_out){
		return(dist)
	}else{
		return(median(dist))
	}
}




