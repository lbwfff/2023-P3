#####################################
#这里有三个function，都是为了把mvmr的输入输出给转换为可方便使用的形态

# library(remotes)
# install_github("WSpiller/MVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)

library(MVMR)
library(tidyr)
library(tibble)
library(dplyr)

tidy_mvmr_output <- function(mvmr_res) {
  #  tidy up MVMR returned output
  mvmr_res %>%
    as.data.frame() %>% 
    rownames_to_column("exposure") %>% 
    dplyr::rename(b='Estimate',  #会存在一个库错误
                  se="Std. Error",
                  pval="Pr(>|t|)") %>% 
    select(-c(`t value`)) %>% 
    TwoSampleMR::generate_odds_ratios()
}

tidy_pvals<-function(df){
  # round up output values and keep p-vals in scientific notation
  df %>% 
    mutate(pval= as.character(pval)) %>% 
    mutate_if(is.numeric, round, digits=4) %>% 
    mutate(pval=as.numeric(pval),
           pval=scales::scientific(pval, digits = 4),
           pval=as.numeric(pval))
}  #这个小数点位数搞那么少干嘛，用来画图的话不是很好

make_mvmr_input <- function(exposure_dat, outcome.id.mrbase=NULL, outcome.data=NULL){
  # provide exposure_dat created in the same way as for TwoSampleMR 
  # also specify the outcome argument [only ONE!] (MR-base ID or full gwas data in .outcome format)
  
  # extract SNPs for both exposures from outcome dataset
  # (for the selected option mr.base or local outcome data)
  if (!is.null(outcome.id.mrbase)) {
    # if mrbase.id is provided
    outcome_dat <- extract_outcome_data(snps = unique(exposure_dat$SNP),
                                        outcomes = outcome.id.mrbase)
  } else if (!is.null(outcome.data)){
    # if outcome df is provided
    outcome_dat <- outcome.data %>% filter(SNP %in% exposure_dat$SNP)
  }
  
  # harmonize datasets
  exposure_dat <- exposure_dat %>% mutate(id.exposure = exposure)
  outcome_harmonised <- mv_harmonise_data(exposure_dat, outcome_dat)
  
  exposures_order <- colnames(outcome_harmonised$exposure_beta)
  
  # Create variables for the analysis 
  
  ### works for many exposures
  no_exp = dim(outcome_harmonised$exposure_beta)[2] # count exposures
  # add beta/se names
  colnames(outcome_harmonised$exposure_beta) <- paste0("betaX", 1:no_exp)
  colnames(outcome_harmonised$exposure_se) <- paste0("seX", 1:no_exp)
  
  XGs <-left_join(as.data.frame(outcome_harmonised$exposure_beta) %>% rownames_to_column('SNP'), 
                  as.data.frame(outcome_harmonised$exposure_se)   %>%rownames_to_column('SNP'), 
                  by = "SNP")
  
  YG <- data.frame(beta.outcome = outcome_harmonised$outcome_beta,
                   se.outcome = outcome_harmonised$outcome_se) %>% 
    mutate(SNP = XGs$SNP)
  
  
  return(list(YG = YG,
              XGs = XGs,
              exposures = exposures_order))
}



