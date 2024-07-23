rm(list = ls())
library(MASS)
source("transCCA.R")
source("transCCA_sv.R")
source("transCCA_mv.R")
source("transCCA_utils.R")
source("W_calculation.R")
source("PMD_dependence.R")
source("transCCA_plot.R")
#setwd("~/Library/CloudStorage/GoogleDrive-xwcaochn@gmail.com/My Drive/transQTL/transCCA/Simulation/releastic")
setwd("D:/TransPCO")
files <- list.files("data/", full.names = TRUE)
X <- Y_cis <- true_variants <- list()
for (i.file in 1:length(files)){
    data <- readRDS(files[i.file])
    X[[i.file]] <- data$X
    Y_cis[[i.file]] <- data$Y
    true_variants[[i.file]] <- data$variant
}
X_comb <- do.call(cbind, X)

coexpression <- list()
for (i in 1:21){
    coexpression[[i]] <- readRDS(paste0("Sigma_module", i, ".rds"))
}

library(MASS)
effect_trans = 0.5
Y_trans <- list()
for (i in 1:length(Y_cis)){
    # sigma <- generate_realistic_correlation_matrix(n)
    sigma <- coexpression[[i]]
    mu <- Y_cis[[i]] * effect_trans
    Y_trans[[i]] <- mu + mvrnorm(n = length(Y_cis[[i]]), mu = rep(0,nrow(sigma)), Sigma = sigma)
}
Y_trans_comb <- do.call(cbind, Y_trans)

Z <- lapply(1:ncol(Y_trans_comb), function(j){
    res <- susieR::univariate_regression(X_comb, Y_trans_comb[,j])
    res$betahat / res$sebetahat
})
Z_comb <- do.call(cbind, Z)
MAF <- colMeans(X_comb)/2

############## RUN ARCHIE and TransCCA ####################
# For one simualtion with 20 TAD regions: K=50
# 1. run ARCHIE - X_comb and Y_trans_comb ---- save ARCHIE results: original (SigmaGG = cor(X_comb); SigmaEE = block sigma you have) & diag (both SigmaGG and SigmaEE are identity)
# 2. run transCCCA - X_comb and Y_trans_comb ---- save transCCA results
# 3. generate null data based on effect_trans=0 (above is 0.5), then calculate Z_null and run ONLY transCCA for 100 replicates.
# 4. follow ARCHIE paper, generate null data 1000 times and run ONLY transCCA for 1000 null data. data.list = list("Sigma_GE" = Sigma_GE_null)


# - load transCCA code
#setwd("~/Library/CloudStorage/GoogleDrive-xwcaochn@gmail.com/My Drive/transQTL/transCCA/code/v2")
source("transCCA.R")
source("transCCA_sv.R")
source("transCCA_mv.R")
source("transCCA_utils.R")
source("W_calculation.R")
source("PMD_dependence.R")
source("transCCA_plot.R")
data.list <- list(list("Z" = Z_comb, "N" = nrow(X_comb)))
res_transcca <- transCCA(data.list = data.list, K=50)


#gc()
#sigmaGG=cor(X_comb)



effect_trans = 0.5
Y_trans <- list()
for (i in 1:length(Y_cis)){
    # sigma <- generate_realistic_correlation_matrix(n)
    sigma <- coexpression[[i]]
    mu <- Y_cis[[i]] * effect_trans
    Y_trans[[i]] <- mu + mvrnorm(n = length(Y_cis[[i]]), mu = rep(0,nrow(sigma)), Sigma = sigma)
}
Y_trans_comb <- do.call(cbind, Y_trans)

Z <- lapply(1:ncol(Y_trans_comb), function(j){
    res <- susieR::univariate_regression(X_comb, Y_trans_comb[,j])
    res$betahat / res$sebetahat
})
Z_comb <- do.call(cbind, Z)
MAF <- colMeans(X_comb)/2

############## RUN ARCHIE and TransCCA ####################
# For one simualtion with 20 TAD regions: K=50
# 1. run ARCHIE - X_comb and Y_trans_comb ---- save ARCHIE results: original (SigmaGG = cor(X_comb); SigmaEE = block sigma you have) & diag (both SigmaGG and SigmaEE are identity)
# 2. run transCCCA - X_comb and Y_trans_comb ---- save transCCA results
# 3. generate null data based on effect_trans=0 (above is 0.5), then calculate Z_null and run ONLY transCCA for 100 replicates.
# 4. follow ARCHIE paper, generate null data 1000 times and run ONLY transCCA for 1000 null data. data.list = list("Sigma_GE" = Sigma_GE_null)


# - load transCCA code
#setwd("~/Library/CloudStorage/GoogleDrive-xwcaochn@gmail.com/My Drive/transQTL/transCCA/code/v2")
source("transCCA.R")
source("transCCA_sv.R")
source("transCCA_mv.R")
source("transCCA_utils.R")
source("W_calculation.R")
source("PMD_dependence.R")
source("transCCA_plot.R")
data.list <- list(list("Z" = Z_comb, "N" = nrow(X_comb)))
res_transcca <- transCCA(data.list = data.list, K=50)



#-------------------------------------------------------------------------------------
res_transcca_null=list()
for (i in 1:10) {
  cat(paste("Simulation:",i,"start"))
  effect_trans = 0
  Y_trans <- list()
  for (i in 1:length(Y_cis)){
    # sigma <- generate_realistic_correlation_matrix(n)
    sigma <- coexpression[[i]]
    mu <- Y_cis[[i]] * effect_trans
    Y_trans[[i]] <- mu + mvrnorm(n = length(Y_cis[[i]]), mu = rep(0,nrow(sigma)), Sigma = sigma)
  }
  Y_trans_comb <- do.call(cbind, Y_trans)
  
  Z <- lapply(1:ncol(Y_trans_comb), function(j){
    res <- susieR::univariate_regression(X_comb, Y_trans_comb[,j])
    res$betahat / res$sebetahat
  })
  Z_comb <- do.call(cbind, Z)
  MAF <- colMeans(X_comb)/2

  data.list <- list(list("Z" = Z_comb, "N" = nrow(X_comb)))
  res_transcca_null[[i]] <- transCCA(data.list = data.list, K=50)
  saveRDS(res_transcca_null,file = "res_transcca_null.rds")
  cat(paste("Simulation:",i,"end"))
}
saveRDS(res_transcca_null,file = "res_transcca_null.rds")

