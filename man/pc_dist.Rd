% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pc_dist.R
\name{pc_dist}
\alias{pc_dist}
\title{Calculate PC metric quantifying the distance between two cohorts}
\usage{
pc_dist(
  ref_pc,
  target_pc,
  ref_pc_col = 4:23,
  target_pc_col = 4:23,
  ref_nvar_col = 2,
  target_nvar_col = 2,
  ind_out = F
)
}
\arguments{
\item{ref_pc}{Required. The plink2 .sscore file containing PC information for each sample in the reference.}

\item{target_pc}{Required. The plink2 .sscore containing PC information for each sample in the target cohort.}

\item{ref_pc_col}{Optional. A vector of column numbers containing the PCs in the reference dataframe. We recommend use 10 or 20 PCs.
Default will be 4:23 using 20 PCs starting from the 4th column, following the format in our provided .sscore file.}

\item{target_pc_col}{Optional. A vector of column numbers containing the PCs in the target dataframe. The numbers of columns should match with the reference.
Default will be 4:23 using 20 PCs starting from the 4th column, following the format in our provided .sscore file.}

\item{ref_nvar_col}{Optional. The column number of \code{ALLELE_CT} from plink2 .sscore file. Default 2.}

\item{target_nvar_col}{Optional. The column number of \code{ALLELE_CT} from plink2 .sscore file. Default 2.}

\item{ind_out}{Optional. Whether the individual distance (for each individual in target) to a reference will be written or not. Default FALSE.}
}
\value{
If ind_out == T, it will output a vector of each individual distance to the reference.
The order of individuals will follow the input order in \code{target_pc}.
If ind_out == F (default), it will only output a number of the PC distance between two cohorts.
}
\description{
This function is used to calculate PC dissimilarity metric to quantify the dissimilarity between two cohorts in terms of the sample genetic profiles.
This metric is not symmetric. It measures how well the reference cohort can represent the target cohort.
Harmonized PCs should be calculated first. Please see https://github.com/covid19-hg/ for more details.
This metric will be standardized accounting for the number of variants used in calculation due to some missing variants.
}
\examples{

dist = pc_dist(ref_pc = "../data/test_cohort1.txt", 
target_pc = "../data/test_cohort2.txt", ind_out = F)

}
