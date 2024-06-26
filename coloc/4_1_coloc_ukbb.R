##############################################
########## Run coloc on regions for qtl and gwas traits ##########
########## These traits are ukbb traits from nealelab ##########
##############################################
# load packages -----
rm(list = ls())
library(data.table)
library(tidyverse)
library(coloc)


# loop over all blood traits -----
file_blood_trait_info <- "/project2/xuanyao/llw/coloc/ukbb_coloc_blood_traits/ukbb_blood_traits.csv"
trait_info <- fread(file_blood_trait_info, sep = ",", header = TRUE)
pheno_seq <- trait_info$`GWAS ID`
n_gwas <- length(pheno_seq)


for(gwasPhenocode in pheno_seq){
  # I/O & paras -----
  qtlN <- 913
  qtlType <- "quant"
  gwasType <- "quant"
  pop <- "EUR"
  pp4Thre <- 0.75
  csThre <- 0.95
  
  ## coloc directory for the trait
  dir_coloc_gwas <- list.files(
    "/scratch/midway2/liliw1/coloc_MSigDB",
    pattern = paste0("ukbb_continuous_", gwasPhenocode),
    full.names = TRUE
  )
  file_qtlColocReg_gwas <- paste(dir_coloc_gwas, "data/qtlColocReg_gwas.txt.gz", sep = '/')
  file_gwasColocReg <- paste(dir_coloc_gwas, "data/gwasColocReg.txt.gz", sep = '/')
  file_gwasRegTruncPthre <- paste(dir_coloc_gwas, "data/gwasRegTruncPthre.txt", sep = '/')
  file_gwasTraitInfo <- "/project2/xuanyao/llw/GWAS/UKB_nealelab/phenotype_manifest.tsv.bgz"
  
  ## Add p-value and trait info
  file_qtlColocReg <- "/scratch/midway2/liliw1/coloc_MSigDB/qtlColocReg.txt.gz"
  file_blood_trait_info <- "/project2/xuanyao/llw/coloc/ukbb_coloc_blood_traits/ukbb_blood_traits.csv"
  
  ## output -----
  file_resColoc <- paste(dir_coloc_gwas, "data/resColoc.txt.gz", sep = '/')
  file_resColocSNP <- paste(dir_coloc_gwas, "data/resColocSNP.txt.gz", sep = '/')
  
  
  # read files -----
  qtlColocReg_gwas <- fread(file_qtlColocReg_gwas)
  gwasColocReg <- fread(file_gwasColocReg)
  gwasRegTruncPthre <- fread(file_gwasRegTruncPthre, header = FALSE, col.names = "Region")
  
  gwasTraitInfoCol <- c("trait_type", "phenocode", "pheno_sex", "description",
                        paste(c("n_cases", "n_controls"), pop, sep = "_") )
  gwasTraitInfo <- fread(cmd = paste("gunzip -c", file_gwasTraitInfo), select = gwasTraitInfoCol)
  
  qtlColocReg <- fread(file_qtlColocReg, header = TRUE)
  trait_info <- fread(file_blood_trait_info, sep = ",", header = TRUE)
  
  
  # GWAS trait info, sample size -----
  gwasTraitInfo_usedGWAS <- gwasTraitInfo %>% filter(phenocode == gwasPhenocode)
  if(nrow(gwasTraitInfo_usedGWAS) != 1) stop("There are more than one trait with the same given phenocode. Check more. e.g. The trait is not quantitative and have more than two phenotype categories. Or, e.g. The phenocode corresponds to multiple traits with different detailed descriptions.")
  n_cases <- gwasTraitInfo_usedGWAS[[paste0("n_cases_", pop)]]
  n_controls <- gwasTraitInfo_usedGWAS[[paste0("n_controls_", pop)]]
  gwasN <- sum(n_cases, n_controls, na.rm = TRUE)
  
  
  # remove SNPs whose se(\beta) equals 0, as coloc needs to use 1/se(\beta) -----
  if(any(gwasColocReg[[paste0("se_", pop)]] == 0)){
    gwasColocReg <- gwasColocReg[gwasColocReg[[paste0("se_", pop)]] != 0, ]
    qtlColocReg_gwas <- qtlColocReg_gwas[qtlColocReg_gwas$SNP_ID %in% gwasColocReg$SNP_ID, ]
  }
  
  
  # prepare coloc files and run coloc -----
  resColoc <- NULL
  resColocSNP <- NULL
  nRegion <- length(gwasRegTruncPthre$Region)
  for(reg in gwasRegTruncPthre$Region){
    
    ## construct coloc dataset D1 & D2 -----
    tmpqtlColocReg_gwas <- qtlColocReg_gwas %>% filter(Region == reg)
    tmpgwasColocReg <- gwasColocReg %>% filter(Region == reg)
    
    D1 <- list("pvalues" = tmpqtlColocReg_gwas$Pval,
               "N" = qtlN,
               "MAF" = tmpqtlColocReg_gwas$MAF,
               "type" = qtlType,
               "snp" = tmpqtlColocReg_gwas$SNP_ID)
    D2 <- list("pvalues" = tmpgwasColocReg[[paste0("pval_", pop)]],
               "beta" = tmpgwasColocReg[[paste0("beta_", pop)]],
               "varbeta" = (tmpgwasColocReg[[paste0("se_", pop)]])^2,
               "type" = gwasType,
               "s" = if(gwasType=="cc") n_cases/gwasN else NULL,
               "snp" = tmpgwasColocReg$SNP_ID,
               "MAF" = tmpgwasColocReg[[paste0("af_", pop)]],
               "N" = gwasN)
    
    if(gwasType=="cc"){
      D2 <- within(D2, rm(MAF))
    }else{
      D2 <- within(D2, rm(s))
    }
    
    
    ## do coloc -----
    coloc_res <- coloc.abf(D1, D2)
    
    ## credible set
    if(coloc_res$summary["PP.H4.abf"] > pp4Thre){
      o <- order(coloc_res$results$SNP.PP.H4,decreasing=TRUE)
      cs <- cumsum(coloc_res$results$SNP.PP.H4[o])
      w <- which(cs > csThre)[1]
      resColocSNP <- rbind(resColocSNP, data.table("Region" = reg,
                                                   "SNP_ID" = as.character(coloc_res$results[o,][1:w,]$snp),
                                                   "SNP.PP.H4" = coloc_res$results[o,][1:w,]$SNP.PP.H4))
    }
    
    ## aggregate results across regions -----
    resColoc <- rbind(resColoc, data.table("Region" = reg, t(coloc_res$summary)) )
    
    
    cat(grep(reg, gwasRegTruncPthre$Region, fixed = TRUE),
        "-th region:", reg, "(out of", nRegion, "regions)", "is done!", "\n")
  }
  
  
  # Add p-value and trait info -----
  resColoc <- qtlColocReg %>% select(c(Signal, Pval)) %>%
    right_join(y = resColoc, by = c("Signal" = "Region")) %>%
    rename("Region" = "Signal") %>%
    mutate("Phenocode" = gwasPhenocode)
  resColoc <- trait_info %>% select(c("GWAS ID", "Trait Abbreviation")) %>%
    right_join(y = resColoc, by = c("GWAS ID" = "Phenocode")) %>%
    rename("Phenocode" = "GWAS ID", "trait" = "Trait Abbreviation") %>%
    arrange(desc(PP.H4.abf))
  
  
  # print out key message or write out -----
  fwrite(resColoc, file_resColoc, quote = FALSE, sep = "\t")
  if(!is.null(resColocSNP)) fwrite(resColocSNP, file_resColocSNP, quote = FALSE, sep = "\t")
  
  cat(which(gwasPhenocode == pheno_seq), '-th out of', n_gwas, "GWASs is done. \n\n")
}

