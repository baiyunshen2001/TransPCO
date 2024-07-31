rm(list = ls())
library(MASS)
#setwd("~/Library/CloudStorage/GoogleDrive-xwcaochn@gmail.com/My Drive/transQTL/transCCA/Simulation/releastic")
setwd("D:/TransPCO")
source("transCCA.R")
source("transCCA_sv.R")
source("transCCA_mv.R")
source("transCCA_utils.R")
source("W_calculation.R")
source("PMD_dependence.R")
source("transCCA_plot.R")
set.seed(123456)

files <- list.files("data/", full.names = TRUE)
X <- Y_cis <- true_variants <- list()
for (i.file in 1:length(files)){
    data <- readRDS(files[i.file])
    X[[i.file]] <- data$X
    Y_cis[[i.file]] <- data$Y
    true_variants [[i.file]] <- data$variant
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
rm(list = ls())
library(MASS)
#setwd("~/Library/CloudStorage/GoogleDrive-xwcaochn@gmail.com/My Drive/transQTL/transCCA/Simulation/releastic")
setwd("D:/TransPCO")
source("transCCA.R")
source("transCCA_sv.R")
source("transCCA_mv.R")
source("transCCA_utils.R")
source("W_calculation.R")
source("PMD_dependence.R")
source("transCCA_plot.R")
set.seed(123456)

files <- list.files("data/", full.names = TRUE)
X <- Y_cis <- true_variants <- list()
for (i.file in 1:length(files)){
  data <- readRDS(files[i.file])
  X[[i.file]] <- data$X
  Y_cis[[i.file]] <- data$Y
  true_variants [[i.file]] <- data$variant
}
X_comb <- do.call(cbind, X)

coexpression <- list()
for (i in 1:21){
  coexpression[[i]] <- readRDS(paste0("Sigma_module", i, ".rds"))
}

res_transcca_null=list()
Y_trans_comb_null=list()
MAF_null=list()
for (n.sim in 1:10) {
  cat(paste("Simulation:",n.sim,"start \n"))
  effect_trans = 0
  Y_trans <- list()
  for (i in 1:length(Y_cis)){
    # sigma <- generate_realistic_correlation_matrix(n)
    sigma <- coexpression[[i]]
    mu <- Y_cis[[i]] * effect_trans
    Y_trans[[i]] <- mu + mvrnorm(n = length(Y_cis[[i]]), mu = rep(0,nrow(sigma)), Sigma = sigma)
  }
  Y_trans_comb <- do.call(cbind, Y_trans)
  Y_trans_comb_null[[n.sim]]=Y_trans_comb
  file_Y=paste0("Y_trans_comb_null.rds")
  saveRDS(Y_trans_comb_null,file = file_Y)
  
  Z <- lapply(1:ncol(Y_trans_comb), function(j){
    res <- susieR::univariate_regression(X_comb, Y_trans_comb[,j])
    res$betahat / res$sebetahat
  })
  Z_comb <- do.call(cbind, Z)
  MAF <- colMeans(X_comb)/2
  MAF_null[[n.sim]]=MAF
  
  file_Z=paste0("Z_comb_null",n.sim,".rds")
  saveRDS(Z_comb,file= file_Z)
  saveRDS(MAF_null,file="MAF_null.rds")

  data.list <- list(list("Z" = Z_comb, "N" = nrow(X_comb)))
  res_transcca_null[[n.sim]] <- transCCA(data.list = data.list, K=50)
  saveRDS(res_transcca_null,file = "res_transcca_null.rds")
  cat(paste("Simulation:",n.sim,"end \n"))
}
saveRDS(res_transcca_null,file = "res_transcca_null.rds")

res_transcca_alt=readRDS("res_transcca_alt.rds")
res_transcca_null=readRDS("res_transcca_null.rds")
for (i in 1:10) {
  file_n=paste0("./plots/cov_v_null_",i,".png")
  png(filename =file_n )
  pheatmap::pheatmap(cor(res_transcca_null[[i]]$res$v))
  dev.off()
}



for (i in 1:10) {
  file_n=paste0("./plots/cov_v_alt_",i,".png")
  png(filename =file_n )
  pheatmap::pheatmap(cor(res_transcca_alt[[i]]$res$v))
  dev.off()
}

Level=c()
for (i in 1:10) {
  res_eig=eigen(cor(res_transcca_null[[i]]$res$v))
  prop=cumsum(res_eig$values)/sum(res_eig$values)
  Level[i]=which(prop>=0.8)[1]
}
Level

Level=c()
for (i in 1:9) {
  res_eig=eigen(cor(res_transcca_alt[[i]]$res$v))
  prop=cumsum(res_eig$values)/sum(res_eig$values)
  Level[i]=which(prop>=0.8)[1]
}
Level

q_matrix01_null=list()
for (i in 1:10) {
  File_Z=paste0("Z_comb_null",i,".rds")
  Z_matrix=readRDS(file = File_Z)
  p_matrix=t(apply(Z_matrix,1,pnorm))
  q_matrix=t(apply(p_matrix,1,function(x){
    qvalue::qvalue(x,pi0 = 1)$qvalue
  }))
  q_matrix01_null[[i]]=(q_matrix<0.05)*1
}
saveRDS(q_matrix01_null,file = "q_matrix01_null.rds")

q_matrix=list()
q_matrix01_alt=list()
for (i in 1:10) {
  File_Z=paste0("Z_comb_alt",i,".rds")
  Z_matrix=readRDS(file = File_Z)
  p_matrix=t(apply(Z_matrix,1,pnorm))
  q_matrix[[i]]=t(apply(p_matrix,1,function(x){
    qvalue::qvalue(x,pi0 = 1)$qvalue
  }))
  q_matrix01_alt[[i]]=(q_matrix[[i]]<0.05)*1
}
saveRDS(q_matrix01_alt,file = "q_matrix01_alt.rds")

q_matrix01_alt=readRDS(file = "q_matrix01_alt.rds")
q_matrix01_null=readRDS(file= "q_matrix01_null.rds")

sum(rowSums(q_matrix01_null[[1]]==1)!=0)
sum(rowSums(q_matrix01_alt[[1]]==1)!=0)

Y_comb_null=readRDS("Y_trans_comb_null.rds")
Y_comb_alt=readRDS("Y_trans_comb_alt.rds")

Y_comb_ano=Y_comb_alt[[1]]%*%res_transcca_alt[[1]]$res$v#[,1:K(1-50)]
Z <- lapply(1:ncol(Y_comb_ano), function(j){
  res <- susieR::univariate_regression(X_comb, Y_comb_ano[,j])
  res$betahat / res$sebetahat
})
Z_comb <- do.call(cbind, Z)



  Z_matrix=Z_comb
  p_matrix=t(apply(Z_matrix,1,pnorm))
  q_matrix=t(apply(p_matrix,1,function(x){
    qvalue::qvalue(x,pi0 = 1)$qvalue
  }))
  EV01=(q_matrix<0.05)*1
  sum(rowSums(EV01==1)!=0)

  
  
#---------------------------------------------------------------------------------------
true_vari=c()  
num_blo=as.vector(do.call(cbind,lapply(X, ncol)))
num_blo=c(0,num_blo[1:(length(num_blo)-1)])
index=cumsum(num_blo)
for (i.files in 1:length(files)) {
  vari=index[i.files]+true_variants[[i.files]]
  true_vari=c(true_vari,vari)
}
  
Y_comb_null=readRDS("Y_trans_comb_null.rds")
Y_comb_alt=readRDS("Y_trans_comb_alt.rds")
res_transcca_alt=readRDS("res_transcca_alt.rds")
res_transcca_null=readRDS("res_transcca_null.rds")

res_pred=list() 
power=c()
rowsum_Ev=list()
for (K in 1:50) {
  Y_comb_ano=Y_comb_alt[[1]]%*%res_transcca_alt[[1]]$res$v[,1:K]#[,1:K(1-50)]
  Z <- lapply(1:ncol(Y_comb_ano), function(j){
    res <- susieR::univariate_regression(X_comb, Y_comb_ano[,j])
    res$betahat / res$sebetahat
  })
  Z_comb <- do.call(cbind, Z)
  p_matrix=t(apply(Z_comb,1,pnorm))
  q_matrix=t(apply(p_matrix,1,function(x){
    qvalue::qvalue(x,pi0 = 1)$qvalue
  }))
  EV01=(q_matrix<0.05)*1
  res_pred[[K]]=which(rowSums(EV01==1)!=0)
  rowsum_Ev[[K]]=rowSums(EV01==1)
  power[K]=sum(res_pred[[K]]%in%true_vari)/length(true_vari)
  cat(paste0("K:",K," done \n"))
} 

saveRDS(rowsum_Ev,"rowsum_Ev.rds")
saveRDS(res_pred,"res_pred_transcca.rds")
saveRDS(power,"res_power_transcca.rds")


res_pred=readRDS(file ="res_pred_transcca.rds")
power=readRDS(file ="res_power_transcca.rds")
q_matrix01_alt=readRDS(file = "q_matrix01_alt.rds")
res_alt_pred=which(rowSums(q_matrix01_alt[[1]]==1)!=0)
power0=sum(res_alt_pred%in%true_vari)/length(true_vari)
res_power=c(power0,power)
power_table=data.frame(K=0:50,Power=res_power)


png(filename = "res_power.png")
ggplot(data = power_table,aes(x=K,y=Power))+
  geom_line()
dev.off()


#---------------------------------------------------------------------------------------
res_pred=readRDS("res_pred_transcca.rds")
FDR=c()
for (K in 1:50) {
  FDR[K]=(length(res_pred[[K]])-sum(res_pred[[K]]%in%true_vari))/length(res_pred[[K]])
}


q_matrix01_alt=readRDS(file = "q_matrix01_alt.rds")
res_alt_pred=which(rowSums(q_matrix01_alt[[1]]==1)!=0)
FDR0=(length(res_alt_pred)-sum(res_alt_pred%in%true_vari))/length(res_alt_pred)
res_FDR=c(FDR0,FDR)
FDR_table=data.frame(K=0:50,FDR=res_FDR)

png(filename = "res_FDR.png")
ggplot(data = FDR_table,aes(x=K,y=FDR))+
  geom_line()
dev.off()

#------------------------------------------------------------------------------------


