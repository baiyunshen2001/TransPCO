rm(list = ls())
library(MASS)
setwd("D:/TransPCO")
source("ARCHIE/codes/helper.R")
source("ARCHIE/codes/main.R")
set.seed(123456)


gene_data=function(n,p,A_snp_cis,A_cis_trans,A_trans,MAF){
  res.list=list()
  
  G=lapply(1:p,function(x) rbinom(n,2,sample(MAF[MAF>=0.1&MAF<=0.4],1)))
  G=do.call(cbind,G)
  res.list$maf=colMeans(G)
  
  res.list$G_unscaled=G
  G=apply(G,2,scale)
  res.list$G=G
  
  
  E_cis=matrix(nrow = n,ncol=ncol(A_snp_cis))
  for (j in 1:ncol(A_snp_cis)) {
    asso=which(A_snp_cis[,j]==1)
    beta=rnorm(length(asso),mean=0.1,sd=sqrt(0.04))
    mean_cis=G[,asso]%*%as.matrix(beta)
    Ecis=lapply(mean_cis,function(x) rnorm(1,mean = x,sd=1))
    E_cis[,j]=as.vector(do.call(cbind,Ecis))
  }
  res.list$E_cis=E_cis
  
  E_trans=matrix(nrow = n,ncol=ncol(A_cis_trans))
  for (j in 1:ncol(A_cis_trans)) {
    asso_cis=which(A_cis_trans[,j]==1)
    asso_trans=which(A_trans[,j]==1)
    if(length(asso_cis)!=0){
      beta=rnorm(length(asso_cis),mean=0.1,sd=sqrt(0.1/length(asso_cis)))#heritablity
      mean_trans=E_cis[,asso_cis]%*%as.matrix(beta)
      Etrans=lapply(mean_trans,function(x) rnorm(1,mean=x,sd=1))
      E_trans[,j]=as.vector(do.call(cbind,Etrans))
    }else{
      beta=rnorm(length(asso_trans),mean=0.1,sd=sqrt(0.1/length(asso_trans)))
      mean_trans=E_cis[,asso_trans]%*%as.matrix(beta)
      Etrans=lapply(mean_trans,function(x) rnorm(1,mean=x,sd=1))
      E_trans[,j]=as.vector(do.call(cbind,Etrans))
    }
  }
  res.list$E_trans=E_trans
  
  return(res.list)
}

gene_data_null=function(n,p,A_snp_cis,A_cis_trans,A_trans,MAF){
  res.list=list()
  
  G=lapply(1:p,function(x) rbinom(n,2,sample(MAF[MAF>=0.1&MAF<=0.4],1)))
  G=do.call(cbind,G)
  res.list$maf=colMeans(G)
  
  res.list$G_unscaled=G
  G=apply(G,2,scale)
  res.list$G=G
  
  
  E_cis=matrix(nrow = n,ncol=ncol(A_snp_cis))
  for (j in 1:ncol(A_snp_cis)) {
    asso=which(A_snp_cis[,j]==1)
    beta=rnorm(length(asso),mean=0.1,sd=sqrt(0.04))
    mean_cis=G[,asso]%*%as.matrix(beta)
    Ecis=lapply(mean_cis,function(x) rnorm(1,mean = x,sd=1))
    E_cis[,j]=as.vector(do.call(cbind,Ecis))
  }
  res.list$E_cis=E_cis
  
  E_trans=matrix(nrow = n,ncol=ncol(A_cis_trans))
  for (j in 1:ncol(A_cis_trans)) {
    asso_cis=which(A_cis_trans[,j]==1)
    asso_trans=which(A_trans[,j]==1)
    if(length(asso_cis)!=0){
      E_trans[,j]=rnorm(n)
    }else{
      beta=rnorm(length(asso_trans),mean=0.1,sd=sqrt(0.1/length(asso_trans)))
      mean_trans=E_cis[,asso_trans]%*%as.matrix(beta)
      Etrans=lapply(mean_trans,function(x) rnorm(1,mean=x,sd=1))
      E_trans[,j]=as.vector(do.call(cbind,Etrans))
    }
  }
  res.list$E_trans=E_trans
  
  return(res.list)
}





n=1000 #sample size
p=40 #SNPs
MAF_sample=readRDS("MAF.rds")

A_snp_cis=lapply(1:8,function(x) rep(1,5))
A_snp_cis=as.matrix(do.call(Matrix::bdiag,A_snp_cis))
A_cis_trans=cbind(c(1,1,0,0,0,0,0,0),
                  c(0,0,1,1,0,0,0,0),
                  c(0,0,0,0,1,1,0,0),
                  c(0,0,0,0,0,0,1,1),0,0,0,0,0)
A_trans=cbind(0,0,0,0,c(1,rep(0,8)),
              c(0,1,1,rep(0,6)),
              c(0,0,0,1,rep(0,5)),
              c(0,0,0,0,1,1,0,0,0),
              c(rep(0,5),1,1,0,0)
)
#true_associations <- matrix(0, 40, 9)
#true_associations[1:10, c(1,5,8)] <- 1
#true_associations[11:20, c(2,6,8,9)] <- 1
#true_associations[21:30, c(3,6,8,9)] <- 1
#true_associations[31:40, c(4,7,9)] <- 1

#res_type_I_error=list()
#res_Power=list()
#res_FDR=list()


n.sim=10^4
K=9
standardmap=0
q2_matrix=matrix(NA,nrow = n.sim, ncol=K)
for(N.sim in 1:n.sim){
  data=gene_data_null(n,p,A_snp_cis,A_cis_trans,A_trans,MAF_sample)
  data_ind=gene_data_null(700,p,A_snp_cis,A_cis_trans,A_trans,MAF_sample)
  Z_matrix=matrix(nrow = p,ncol = ncol(A_trans))
  for (i in 1:p) {
    for (j in 1:ncol(A_trans)) {
      res_z=susieR::univariate_regression(data$G[,i],data$E_trans[,j])
      Z_matrix[i,j]=res_z$betahat/res_z$sebetahat
    }
  }
  pmat=2*(1-pnorm(abs(Z_matrix)))
  if(sum(pmat<0.001)>0)
    standardmap=standardmap+1
  
  maf=data$maf
  
  SigmaGG=diag(1,nrow = p)
  SigmaEE=cor(data_ind$E_trans)
  SigmaGE=Z_matrix/sqrt(2*n*maf*(1-maf))
  
  res_ARCHIE=archie_work(Sigma_GE = SigmaGE,Sigma_GG = SigmaGG,Sigma_EE = SigmaEE,K=K)
  q2=res_ARCHIE$qs
  q2_matrix[N.sim,]=q2
}
TIE_standardmap=standardmap/n.sim

q2_null_matrix=matrix(NA,nrow = n.sim, ncol=K)
for(N.sim in 1:n.sim){
  data=gene_data_null(n,p,A_snp_cis,A_cis_trans,A_trans,MAF_sample)
  data_ind=gene_data_null(700,p,A_snp_cis,A_cis_trans,A_trans,MAF_sample)
  
  SigmaGG_null=diag(1,nrow = p)
  SigmaEE_null=cor(data_ind$E_trans)
  Kronecker_GE=kronecker(X=SigmaEE_null,Y=SigmaGG_null)
  vec_GE=mvrnorm(n=1,mu=rep(0,nrow(SigmaEE_null)*nrow(SigmaGG_null)),Sigma = Kronecker_GE)
  SigmaGE_null=matrix(vec_GE,byrow = TRUE,nrow = nrow(SigmaGG_null),ncol = nrow(SigmaEE_null))
  
  res_ARCHIE_null=archie_work(Sigma_GE = SigmaGE_null,Sigma_GG = SigmaGG_null,Sigma_EE = SigmaEE_null,K=K)
  q2_null=res_ARCHIE_null$qs
  q2_null_matrix[N.sim,]=q2_null
}


pvalue_TIE=matrix(NA,nrow = n.sim, ncol=K)
for (i in 1:n.sim) {
  for (j in 1:K) {
    pvalue_TIE[i,j]=sum(q2_matrix[i,j]>q2_null_matrix[,j])/n.sim
  }
}

TIE=mean(apply(pvalue_TIE, 1, function(x) sum(x<(10^(-3)))))
TIE


q2_matrix=matrix(NA,nrow = n.sim/10, ncol=K)
for(N.sim in 1:(n.sim/10)){
  data=gene_data(n,p,A_snp_cis,A_cis_trans,A_trans,MAF_sample)
  data_ind=gene_data(700,p,A_snp_cis,A_cis_trans,A_trans,MAF_sample)
  Z_matrix=matrix(nrow = p,ncol = ncol(A_trans))
  for (i in 1:p) {
    for (j in 1:ncol(A_trans)) {
      res_z=susieR::univariate_regression(data$G[,i],data$E_trans[,j])
      Z_matrix[i,j]=res_z$betahat/res_z$sebetahat
    }
  }
  maf=data$maf
  SigmaGE=Z_matrix/sqrt(2*n*maf*(1-maf))
  SigmaGG=diag(1,nrow = p)
  SigmaEE=cor(data_ind$E_trans)
  res_ARCHIE=archie_work(Sigma_GE = SigmaGE,Sigma_GG = SigmaGG,Sigma_EE = SigmaEE,K=K)
  q2=res_ARCHIE$qs
  q2_matrix[N.sim,]=q2
}


q2_null_matrix=matrix(NA,nrow = n.sim, ncol=K)
for(N.sim in 1:n.sim){
  data=gene_data_null(n,p,A_snp_cis,A_cis_trans,A_trans,MAF_sample)
  data_ind=gene_data_null(700,p,A_snp_cis,A_cis_trans,A_trans,MAF_sample)
  
  SigmaGG_null=diag(1,nrow = p)
  SigmaEE_null=cor(data_ind$E_trans)
  Kronecker_GE=kronecker(X=SigmaEE_null,Y=SigmaGG_null)
  vec_GE=mvrnorm(n=1,mu=rep(0,nrow(SigmaEE_null)*nrow(SigmaGG_null)),Sigma = Kronecker_GE)
  SigmaGE_null=matrix(vec_GE,byrow = TRUE,nrow = nrow(SigmaGG_null),ncol = nrow(SigmaEE_null))
  
  res_ARCHIE_null=archie_work(Sigma_GE = SigmaGE_null,Sigma_GG = SigmaGG_null,Sigma_EE = SigmaEE_null,K=K)
  q2_null=res_ARCHIE_null$qs
  q2_null_matrix[N.sim,]=q2_null
}


pvalue_power=matrix(NA,nrow = n.sim/10, ncol=K)
for (i in 1:(n.sim/10)) {
  for (j in 1:K) {
    pvalue_power[i,j]=sum(q2_matrix[i,j]>q2_null_matrix[,j])/n.sim
  }
}

res_power=mean(apply(pvalue_power, 1, function(x) sum(x<(10^(-3)))))
res_power
TIE
TIE_standardmap

