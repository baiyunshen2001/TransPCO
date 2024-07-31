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


res_transcca_null=list()
Y_trans_comb_null=list()
MAF_null_true=list()
for (n.sim in 1:1) {
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
  file_Y=paste0("Y_trans_comb_null_true.rds")
  saveRDS(Y_trans_comb_null,file = file_Y)
  
  Z <- lapply(1:ncol(Y_trans_comb), function(j){
    res <- susieR::univariate_regression(X_comb, Y_trans_comb[,j])
    res$betahat / res$sebetahat
  })
  Z_comb_true <- do.call(cbind, Z)
  MAF_true <- colMeans(X_comb)/2
  MAF_null_true[[n.sim]]=MAF_true
  
  file_Z=paste0("Z_comb_null",n.sim,"_true.rds")
  saveRDS(Z_comb_true,file= file_Z)
  saveRDS(MAF_null_true,file="MAF_null_true.rds")
  
  data.list <- list(list("Z" = Z_comb_true, "N" = nrow(X_comb)))
  res_transcca_null[[n.sim]] <- transCCA(data.list = data.list, K=8)
  saveRDS(res_transcca_null,file = "res_transcca_null_true.rds")
  cat(paste("Simulation:",n.sim,"end \n"))
}
saveRDS(res_transcca_null,file = "res_transcca_null_true.rds")


res_transcca_alt=list()
Y_trans_comb_alt=list()
MAF_alt_true=list()
for (n.sim in 1:1) {
  cat(paste("Simulation:",n.sim,"start \n"))
  effect_trans = 0.5
  Y_trans <- list()
  for (i in 1:length(Y_cis)){
    # sigma <- generate_realistic_correlation_matrix(n)
    sigma <- coexpression[[i]]
    mu <- Y_cis[[i]] * effect_trans
    Y_trans[[i]] <- mu + mvrnorm(n = length(Y_cis[[i]]), mu = rep(0,nrow(sigma)), Sigma = sigma)
  }
  Y_trans_comb <- do.call(cbind, Y_trans)
  Y_trans_comb_alt[[n.sim]]=Y_trans_comb
  file_Y=paste0("Y_trans_comb_alt_true.rds")
  saveRDS(Y_trans_comb_alt,file = file_Y)
  
  Z <- lapply(1:ncol(Y_trans_comb), function(j){
    res <- susieR::univariate_regression(X_comb, Y_trans_comb[,j])
    res$betahat / res$sebetahat
  })
  Z_comb_true <- do.call(cbind, Z)
  MAF_true <- colMeans(X_comb)/2
  MAF_alt_true[[n.sim]]=MAF_true
  
  file_Z=paste0("Z_comb_alt",n.sim,"_true.rds")
  saveRDS(Z_comb_true,file= file_Z)
  saveRDS(MAF_alt_true,file="MAF_alt_true.rds")
  
  data.list <- list(list("Z" = Z_comb_true, "N" = nrow(X_comb)))
  res_transcca_alt[[n.sim]] <- transCCA(data.list = data.list, K=8)
  saveRDS(res_transcca_alt,file = "res_transcca_alt_true.rds")
  cat(paste("Simulation:",n.sim,"end \n"))
}
saveRDS(res_transcca_alt,file = "res_transcca_alt_true.rds")


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
res_transcca_alt_true=readRDS("res_transcca_alt_true.rds")
res_transcca_null_true=readRDS("res_transcca_null_true.rds")


res_pred=list() 
power=c()
rowsum_Ev=list()
for (K in 1:8) {
  Y_comb_ano=Y_comb_alt[[1]]%*%res_transcca_alt_true[[1]]$res$v[,1:K]#[,1:K(1-50)]
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

saveRDS(rowsum_Ev,"rowsum_Ev_true.rds")
saveRDS(res_pred,"res_pred_transcca_true.rds")
saveRDS(power,"res_power_transcca_true.rds")


res_pred=readRDS(file ="res_pred_transcca_true.rds")
power=readRDS(file ="res_power_transcca_true.rds")

q_matrix=list()
q_matrix01_alt_true=list()
for (i in 1:1) {
  File_Z=paste0("Z_comb_alt",i,"_true.rds")
  Z_matrix=readRDS(file = File_Z)
  p_matrix=t(apply(Z_matrix,1,pnorm))
  q_matrix[[i]]=t(apply(p_matrix,1,function(x){
    qvalue::qvalue(x,pi0 = 1)$qvalue
  }))
  q_matrix01_alt_true[[i]]=(q_matrix[[i]]<0.05)*1
}

res_alt_pred=which(rowSums(q_matrix01_alt_true[[1]]==1)!=0)
power0=sum(res_alt_pred%in%true_vari)/length(true_vari)
res_power=c(power0,power)
power_table=data.frame(K=0:8,Power=res_power)


png(filename = "res_power_true.png")
ggplot(data = power_table,aes(x=K,y=Power))+
  geom_line()
dev.off()