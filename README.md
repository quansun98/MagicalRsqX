# MagicalRsqX

## Introduction

MagicalRsq-X is a cross-cohort machine learning based method to obtain better calibrated imputation quality scores. 
It was built based on our prevoius [MagicalRsq](https://github.com/quansun98/MagicalRsq), where the training and testing cohorts are homogenous.
We extended the MagicalRsq models to MagicalRsq-X with modified features, leading to a better generalized method applicable to more real life scenarios.
Training, testing and evaluating MagicalRsq-X are exactly the same as for MagicalRsq.
We also repeated some of the previous documentions for MagicalRsq here for convenience. 

For training MagicalRsq-X, we need true imputation quality (true R2) in training samples at a reasonable number of markers (>10k is recommended). 
Please check out the section **Calculate true R2 in training data** below for details. 
During training, XGBoost method is used to build models to predict true R2 with a number of variant-level features, 
including the standard quality metric (standard Rsq) among training samples, 
as well as various population genetics statistics (such as nucleotide diversity, Tajima’s D, and haplotype homozygosity statistic) 
and allele frequencies, recombination rate, and LD scores of major continental populations from external sources (e.g., the 1000 Genomes Project or the TOPMed Project). 
Users can also train MagicalRsq models themselves with other variant-level features as they view fit. 
Details regarding the features we used are under section **Variant-level features** below and 
model training procedure via XGBoost is detailed under section **Train MagicalRsq-X models** below.
We release the models we trained using data from the TOPMed cohorts under folder `model/`. 
To apply the pretrained model to target datasets, please refer to section **Select pre-trained models** and **Obtain MagicalRsq-X in target data**. 
 
## Data

* Population genetics features summarized by the [S/HIC](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005928) paper can be downloaded at: ftp://yunlianon:anon@rc-ns-ftp.its.unc.edu/MagicalRsq/SHIC/
* [TOP-LD](<https://www.cell.com/ajhg/fulltext/S0002-9297(22)00154-9>) population specific allele frequencies and LD scores can be downloaded at: ftp://yunlianon:anon@rc-ns-ftp.its.unc.edu/MagicalRsqX/AFLD/
* Recombination rate from 1000G can be downloaded at: ftp://yunlianon:anon@rc-ns-ftp.its.unc.edu/MagicalRsqX/RecombRate/.
 
## Installation

1. First, install the `devtools` package if you haven't already. 

		install.packages("devtools")

2. Load the devtools package and install through Github.

		library(devtools)
		install_github("quansun98/MagicalRsqX")


## Instructions

Below are the steps to train MagicalRsq models in training data and apply to target data. Note that it is not necessary to have the same variants in training and target datasets.

### Integrate variant-level features

Follow the steps below if you want to integrate the same population genetics features we used. 
You can skip these steps if your data already contains all the variant-level features you want to leverage.
The demonstration uses one chromosome (chr22) as an example. Note that:

* All the positions are in hg38. You need to liftover if you are using reference panels in hg19, e.g. 1000G, HRC, etc..
* If you choose to integrate the population genetic features we used, you need to run `data_integrate` function separately for each chromosome.
* Training data and testing data must have the same features. No need to keep the features in the same order, since you can specify the feature column numbers when training and calculating MagicalRsq, in training and target data respectively. But we highly recommend you keep the same feature order in training and target data.


1. Start from minimac output `.info.gz` file:

		zcat data/chr22.info.gz | sed '1d' | cut -f 1,5,7 | sed 's/:/\t/g' | sed 's/^chr//g' | sed '1iCHR POS REF ALT MAF Rsq' | sed 's/ /\t/g' | bgzip > chr22.imp.sumstat.txt.gz

You should also integrate true R2 information (detailed in the section **Calculate true R2 in training data** below) in this file if you want to use it for model training or evaluation purposes.

2. Download S/HIC, allele frequency and LD scores, and recombination rate features.

		wget -r ftp://yunlianon:anon@rc-ns-ftp.its.unc.edu/MagicalRsq/SHIC/
		wget -r ftp://yunlianon:anon@rc-ns-ftp.its.unc.edu/MagicalRsqX/AFLD/
		wget -r ftp://yunlianon:anon@rc-ns-ftp.its.unc.edu/MagicalRsqX/RecombRate/

3. Run `data_integrate` function in R to integrate imputation summary statistics with other variant-level features.  
Here we use population genetics statistics summarized in the S/HIC paper, population specific allele frequencies and LD scores, as well as recombination rate as included features.
You can include them and other variant-level features at your own discretion.

		library(MagicalRsqX)
		chr22_for_MRsqX = data_integrate("chr22.imp.sumstat.txt.gz", chr = 22, pos_col = 2, SHIC_dir = "../SHIC/", 
			rr_dir = "../RecombRate/", AFLD_dir = "../AFLD/", outfile = "test_chr22.txt.gz")

You will get integrated data for chr22 written in the file `test_chr22.txt.gz`.

### Calculate true R2 in training data

True R2 is defined as the squared Pearson correlation between imputed dosages and true genotypes. 
It can only be calculated if you have additional genotypes other than the array-genotype data that used as imputation target.
You can download the [script](https://yunliweb.its.unc.edu/software/doseR2_VCF.tar.gz) we used to calculate true R2.
Detailed instructions of this script is included in the downloaded zip folder.


### Train MagicalRsq/MagicalRsq-X models

In order to train MagicalRsq-X models, you need true R2 (true imputation quality) for a reasonable number (again, >10k is recommended) of variants among training individuals. 
Refer to the section **Calculate true R2 in training data** above for obtaining true R2.
Note that training individuals do not have to be independent of target individuals, if the variants used for training are independent of those imputed among target individuals.
The procedure will be demonstrated with the example data `toy_chr22_50k_integrated.txt.gz` in the `data/` folder.


1. Similarly, you need to integrate imputation summary statistics with other variant-level features.
Note that the input file should contain true R2 information.

		integrated = data_integrate("data/toy_chr22_50k.txt.gz", chr = 22, pos_col = 2, SHIC_dir = "../SHIC/",
			AFLD_dir = "../AFLD/", rr_dir = "../RecombRate/", outfile = "toy_chr22_50k_integrated.txt.gz") 

2. Train MagicalRsq-X models with the integrated data. 
You can choose how many variants you want to include in the training models. We recommend using 10k - 1million variants. 
Too many variants will increase computational burden, while leading to little performance gain.

		toy_model = train_MagicalRsq(file = "toy_chr22_50k_integrated.txt.gz", outfile = "toy_model", nvar_use = 100000) 

Note that:

(1) the input can contain markers across all chromosomes in *one* file (like the example above).
Alternatively, the `train_MagicalRsq` function also accepts input in multiple files, one for each chromosome.
In the latter case, input files would be a vector of file names as in the sample command line below.

	toy2_model = train_MagicalRsq(file = paste0("toy2_chr",1:22,".txt.gz"), outfile = "toy2", nvar_use = 500000)

(2) You can also specify MAF category (must be either "common", "lowfreq", or "rare") if providing the column number of MAF, to train models specific for one MAF category variants.

	toy_model_common = train_MagicalRsq(file = "toy_chr22_50k_integrated.txt.gz", outfile = "toy_model_common", MAF_cate = 'common', MAFCol = 6, nvar_use = 100000)
 
If not specified, it will train one model with all (or randomly subsetted if specifying `nvar_use`) variants disregarding MAF information.


### Select pre-trained models to calculate MagicalRsq-X

We developed a distance metric based on the harmonized PCs to quantify how well a reference can represent a target cohort.
Please calculate harmonized PCs first following the COVID-19 HGI instructions [here](<https://github.com/covid19-hg/pca_projection/blob/master/docs/prerequisites.md>).

Our pre-trained models and corresponding harmonized PCs are available to download [here]()

1. Calculate dissimilarity metrics for your target cohort with respect to the provided references.
Here we use AA cohort as a demonstration.

	all_dist = c()
	refs = c("BioMe_AA","MESA_AA","JHS","WHI_AA")
	for(cohort in refs){
		all_dist = c(all_dist, 
			pc_dist(ref_pc = paste0("reference/",cohort,"_harm_PC.sscore"),
				target_pc = "target_harm_PC.sscore", ind_out = F))
	}
	names(all_dist) = refs
	all_dist
	
	BioMe_AA     MESA_AA         JHS      WHI_AA 
	0.018375499 0.005504560 0.003868932 0.004081391	

We can see that the dissimilarity metrics < 0.02 for all the four references, with the smallest reference being JHS.
Therefore, JHS should be a good reference to go.

2. Download the pre-trained MagicalRsq-X models based on the chosen model.

	wget ftp://yunlianon:anon@rc-ns-ftp.its.unc.edu/MagicalRsqX/models/JHS_v3_common_100k_rep_1

For each training cohort and each variant MAF category, we provided two models trained with either 100k or 500k variants.

A special note that our pre-trained models have the following order of features: Rsq, SHIC features, AF features, recombination rate, long-range LD scores, short-range LD scores.
If you would like to directly apply our pre-trained models and your target cohort follows the above feature integration step,
please **specify `FeatureCols = c(6:83,92,84:91)`**. See the next section for more details.


### Obtain MagicalRsq-X in target data

1. Similarly, you need to integrate imputation summary statistics with other variant-level features.
Detailed in the section **Integrate variant-level features** above.

2. After getting the integrated data, you can calculate MagicalRsq values for each variant in the integrated data file

        chr22_MagicalRsq  = calc_MagicalRsq(file = "test_chr22.txt.gz", model = "toy_model",
                FeatureCols = 6:92, keptCols = 1:7, outfile = "test_MagicalRsq.txt.gz")

You will get variant-level MagicalRsq values, together with the original Rsq and estimated MAF written in the file `test_MagicalRsq.txt.gz`

Note that the testing data feature orders should be the same as the training models.
Please specify `FeatureCols = c(6:83,92,84:91)` if you directly use our pre-trained models with the provided data integration steps.

### Evaluate Rsq/MagicalRsq performance (true R2 required)

After calculation of MagicalRsq, you can evaluate the performance of Standard Rsq and/or MagicalRsq compared to true R2.
You should specify the column numbers of trueR2, Rsq and MagicalRsq (otherwise it will assume the same data format as the example data).
You can also choose to disregard extremely rare variants from evaluation by setting minimum MAF or (minimum MAC if you also specify the number of samples N).
You will need to specify the column number for MAF if you use this option.

	chr22_eval = eval_MagicalRsq(file = "test_MagicalRsq.txt.gz", minMAF = 0.001, trueR2Col = 8, MAFCol = 5, RsqCol = 6, MagicalRsqCol = 7)

Equivalently you can have (suppose your sample size for the evaluation data is 1,000):

	chr22_eval = eval_MagicalRsq(file = "test_MagicalRsq.txt.gz", minMAC = 2, N = 1000, trueR2Col = 8, MAFCol = 5, RsqCol = 6, MagicalRsqCol = 7)

All 22 chromosomes can also be evaluated together:

	all_eval = eval_MagicalRsq(file = paste0("test2_chr",1:22,".txt.gz"))



