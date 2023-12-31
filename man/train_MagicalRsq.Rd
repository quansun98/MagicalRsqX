% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/train_MagicalRsq.R
\name{train_MagicalRsq}
\alias{train_MagicalRsq}
\title{Train MagicalRsq model}
\usage{
train_MagicalRsq(
  file,
  header = T,
  outfile = NULL,
  nvar_use = NULL,
  MAF_cate = NULL,
  p_train = 0.7,
  trueR2Col = 5,
  MAFCol = NULL,
  FeatureCols = 7:93,
  seed = 123,
  nrounds = 3000,
  max.depth = 6,
  verbose = 1,
  early_stopping_rounds = 50,
  print_every_n = 20
)
}
\arguments{
\item{file}{Required. A vector of strings with each element being one file name. Could be a single file or multiple files (if separated by chromosomes).}

\item{header}{Logical. The parameter is passing to the \code{fread} function. TRUE by default. Strongly recommended containing a header in your data file.}

\item{outfile}{Name of the MagicalRsq model output. Strings. Strongly recommended to specify the output file name to save the model.}

\item{nvar_use}{Number of variants to use for the model training and validation purpose. If not specified, all variants in the data file will be used.
Note that the model training will take much longer time if using full set of variants. Recommended values ranging from 10,000 to 1,000,000.}

\item{MAF_cate}{Minor allee frequency category. Must be chosen from c("common","lowfreq","rare). If not specified, it will use all the variants to train models.}

\item{p_train}{Proportion of variants used in training part. The rest of the variants would be used for validation purpose. Default value is 70\%.
Note that the total number of variants for validation is \code{nvar_use}*(1 - \code{p_train}).}

\item{trueR2Col}{Column number of true R2. Default value 5.}

\item{MAFCol}{Column number of minor allele frequency. Optional. Must be specified if \code{MAF_cate} is used.}

\item{FeatureCols}{Vector of column numbers of all the features used for model training. By default, it will be 7:93 including 87 features covering
population genetics, allele frequencies, LD scores across populations, etc., provided online. More user-defined features are also feasible.}

\item{seed}{Random seed. Default value 123.}

\item{nrounds}{Parameter passing to the \code{xgb.train} function: max number of boosting iterations. Default value 3000.}

\item{max.depth}{Parameter passing to the \code{xgb.train} function: maximum depth of a tree. Default value 6.}

\item{verbose}{Parameter passing to the \code{xgb.train} function: If 0, xgboost will stay silent. If 1, it will print information about performance. If 2, some additional information will be printed out.}

\item{early_stopping_rounds}{Parameter passing to the \code{xgb.train} function: training with a validation set will stop if the performance doesn't improve for \code{k} rounds. Default value 50.}

\item{print_every_n}{Parameter passing to the \code{xgb.train} function: Print each n-th iteration evaluation messages when \code{verbose > 0}. Default is 20.}
}
\value{
The trained MagicalRsq XGBoost model.
}
\description{
This function is used to train MagicalRsq models with data files containing all the features used.
}
\examples{

toy_model = train_MagicalRsq(file = "toy_chr22_50k_integrated.txt.gz", 
outfile = "toy_model", nvar_use = 500000)

toy_common = train_MagicalRsq(file = "toy_chr22_50k_integrated.txt.gz", 
MAF_cate = "common", MAFCol = 6, outfile = "toy_common")

## If data are separated by chromosomes 
test_model = train_MagicalRsq(file = paste0("test_chr",1:22,".txt.gz"), 
outfile = "test", nvar_use = 100000)

}
