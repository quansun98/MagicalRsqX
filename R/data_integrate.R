#' Integrate data for MagicalRsq model training and testing
#' 
#' This function is used to integrate user-input dataset containing Rsq information with provided recombination rate, SHIC features, allele frequencies, and LD scores in different populations.
#' Recombination rate can be downloaded at ftp://yunlianon:anon@rc-ns-ftp.its.unc.edu/MagicalRsqX/RecombRate/.
#' SHIC features can be downloaded at ftp://yunlianon:anon@rc-ns-ftp.its.unc.edu/MagicalRsq/SHIC/. 
#' Population-specific allele frequencies and LD scores can be downloaded at ftp://yunlianon:anon@rc-ns-ftp.its.unc.edu/MagicalRsqX/AFLD/.
#'
#' @param file Required. Input file name. The file should be chromosom-specific, i.e., all the variants should be on the same chromosome. 
#' It should contain minimum of position (in hg38), reference allele and alternative allele, with column names as POS, REF and ALT respectively. 
#' Strongly recommended to include Rsq from minimac (or info from IMPUTE) given by any imputation software to improve model performance. 
#' The file could also contain true R2 information for model training or evaluation purpose.
#' @param header Logical. The parameter is passing to the `fread` function. TRUE by default. Strongly recommended containing a header in your data file.
#' @param chr Chromosome number. Required, integer. 
#' @param pos_col GRCH38/hg38 position column number. Required, integer.
#' @param SHIC_dir Directory containing downloaded SHIC features. Optional.
#' @param AFLD_dir Directory containing downloaded population-specifc allele frequencies and LD scores. Optional.
#' @param rr_dir Directory containing downloaded recombination rate from 1000G. Optional.
#' @param outfile Output file name. Optional but strongly recommended.
#' 
#' @return Dataframe of integrated data that could be used as MagicalRsq training and testing input.
#'
#' @examples
#' 
#' toy_chr22 = data_integrate("toy_chr22_50k.txt.gz", chr = 22, 
#' pos_col = 2, SHIC_dir = "../SHIC/", AFLD_dir = "../AFLD/",
#' rr_dir = "../RecombRate/", outfile = "toy_chr22_integrated.txt.gz")
#' 
#' @importFrom data.table fread fwrite
#' @importFrom dplyr left_join
#'
#' @export


data_integrate = function(file, header = T, chr, pos_col, rr_dir = NULL, SHIC_dir = NULL, AFLD_dir = NULL, outfile = NULL){

dat = fread(file, header = header)

# merge with SHIC features

if(!is.null(SHIC_dir)){

nrow = dim(dat)[1]
ncol = dim(dat)[2]

dat = as.data.frame(cbind(dat, as.data.frame(matrix(0, nrow = nrow, ncol = 66)))) # result data frame

pops = c("CEU","GWD", "JPT", "LWK", "PEL", "YRI")
for (j in 1:6) {
  pop = pops[j]
  begin_col = ncol + 1 + (j-1)*11
  end_col = begin_col + 10
  
  ss = as.data.frame(fread(paste0(SHIC_dir, "/chr", chr, "_", pop, "_b38.txt")))
  
  for (i in 1:nrow) {
    index = (dat[i,pos_col] > ss[,2]) & (dat[i,pos_col] < ss[,3])
    index = which(index==1)
    if(length(index) == 0) {dat[i,begin_col:end_col] = rep(NA,11) }
    else{ 
      index = index[1]
      dat[i,begin_col:end_col] = as.numeric(ss[index,-c(1:4,6)]) 
    }
  }
  colnames(dat)[begin_col:end_col] = paste0(colnames(ss)[-c(1:4,6)],".",pop) 
  
}

drop <- apply(dat,1,anyNA)
dat <- dat[!drop,]

message(paste(dim(dat)[1], "variants left after merging with SHIC features."))
}

# merge with allele freq and LD scores

if(!is.null(AFLD_dir)){

nrow = dim(dat)[1]
ncol = dim(dat)[2]

AF = fread(paste0(AFLD_dir,"/chr",chr,".gz"))

## split into 2 chunks to speed up

num_fold = floor(nrow/2)
dat1 = dat[1:num_fold,]
dat2 = dat[(num_fold+1):nrow, ]

## chunk 1
res1 = left_join(dat1, AF[,-1], by = c("POS"="POS", "REF"="REF", "ALT"="ALT"))
drop = apply(res1, 1, anyNA)
res1 = res1[!drop,]

## chunk 2
res2 = left_join(dat2, AF[,-1], by = c("POS"="POS", "REF"="REF", "ALT"="ALT"))
drop = apply(res2, 1, anyNA)
res2 = res2[!drop,]

dat = rbind(res1, res2)

message(paste(dim(dat)[1], "variants left after merging with population specific AF and LD scores."))

}

# merge with recombination rate

num <- dim(dat)[1]
dat$Recomb_rate = rep(NA, num)
rr = fread(paste0(rr_dir,"/genetic_map_chr",chr, "_combined_b38.txt.gz"), header = T)
rr = as.data.frame(rr)

for (i in 1:num) {
	diff = as.vector(rr[,2] - dat[i,2])
	index <- which(diff >= 0)[1]
	if(is.na(index)){
		dat[i,"Recomb_rate"] = rr[dim(rr)[1],3]
	}else if(index == 1){
		dat[i,"Recomb_rate"] = rr[1,3]
	}else if(abs(diff[index-1]) > abs(diff[index])){
		dat[i,"Recomb_rate"] = rr[index,3]
	}else{
		dat[i,"Recomb_rate"] = rr[index-1,3]
	}
}

drop <- apply(dat,1,anyNA)
dat <- dat[!drop,]

message(paste(dim(dat)[1], "variants left after merging with recombination rate."))


if(!is.null(outfile)){
message(paste("Writing integrated data to", outfile))
fwrite(dat, file = outfile, quote = F, sep = "\t", col.names = T, row.names = F)
}

return(dat)

}


