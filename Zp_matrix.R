library(WGCNA)
library(flashClust)
library(mvtnorm)
library(data.table)
library(tidyverse)

setwd("D:/Jupyter_notebook/TransPCO")
dir_pco="./simulation/script_lambda0.1/"
source(paste0(dir_pco, "ModifiedPCOMerged.R"))
source(paste0(dir_pco, "liu.R"))
source(paste0(dir_pco, "liumod.R"))
source(paste0(dir_pco, "davies.R"))
dyn.load(paste0(dir_pco, "qfc.so"))
source(paste0(dir_pco, "ModifiedSigmaOEstimate.R"))
dir_ARCHIE="./ARCHIE/codes/"
source(paste0(dir_ARCHIE,"main.R"))
source(paste0(dir_ARCHIE,"helper.R"))

set.seed(123456)
options(stringsAsFactors = FALSE);
enableWGCNAThreads()
load("oed.RData")
gene.names=rownames(oed)
trans.oed=t(oed)
n=500;
datExpr=trans.oed[,1:n]
SubGeneNames=gene.names[1:n]
powers = c(c(1:10), seq(from = 12, to=20, by=2));
sft=pickSoftThreshold(datExpr,dataIsExpr = TRUE,powerVector = powers,corFnc = cor,corOptions = list(use = 'p'),networkType = "unsigned")

softPower = 7;

#calclute the adjacency matrix
adj= adjacency(datExpr,type = "unsigned", power = softPower);

#turn adjacency matrix into topological overlap to minimize the effects of noise and spurious associations
TOM=TOMsimilarityFromExpr(datExpr,networkType = "unsigned", TOMType = "unsigned", power = softPower);
colnames(TOM) =rownames(TOM) =SubGeneNames
dissTOM=1-TOM

geneTree = flashClust(as.dist(dissTOM),method="average");
minModuleSize = 20;

# Module identification using dynamic tree cut

dynamicMods = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize);
#dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method="hybrid", deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);

#the following command gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
table(dynamicMods)
n_sim <- 10^3
var.b <- 0.001
caus <- 0.5
#Z_matrix=matrix(NA,nrow = 9000,ncol = 500)
#p_matrix=matrix(NA,nrow = 9000,ncol = 9)

N.seq <- c(200, 400, 600, 800)
models <- paste0("N=", N.seq)
N.sim=10
col=1



Z_matrix=NULL
p_matrix=NULL

for(model in models){
  for (i in 1:N.sim) {
    Z_matrix[[model]][[i]]=matrix(NA,nrow = 9000,ncol = 500)
    p_matrix[[model]][[i]]=matrix(NA,nrow = 9000,ncol = 9)
  }
}

col=1
  for (sel_module in 0:8) {
    module=datExpr[,dynamicMods==sel_module]
    Sigma=as.matrix(cor(module))
    SigmaO=ModifiedSigmaOEstimate(Sigma,method = "liu")
    K=dim(Sigma)[1]
    B <- matrix(rep(NA, K*length(models)), ncol = length(models),
                dimnames = list(NULL, models))
    #beta.n <- as.matrix(rnorm(K, sd = sqrt(var.b)))
    beta.n <- rbind(as.matrix(rnorm(floor(caus*K), sd = sqrt(var.b))), as.matrix(rep(0, K-floor(caus*K))))
    B[, models] <- beta.n %*% sqrt(N.seq)
    
    for (i in 1:N.sim) {
      
    for (model in models) {
      for (j in 1:9) {
        if(j==sel_module+1){
          Z.mat= rmvnorm(1000, B[,model], Sigma)
          p.mat= as.matrix(ModifiedPCOMerged(Z.mat=Z.mat, Sigma=Sigma, SigmaO=SigmaO,method = "liu"))
        }
        else{
          Z.mat= rmvnorm(1000, rep(0, K), Sigma)
          p.mat=as.matrix(ModifiedPCOMerged(
            Z.mat = Z.mat, Sigma = Sigma, SigmaO = SigmaO, method = "liu"
          ))
        }
        Z_matrix[[model]][[i]][(1000*(j-1)+1):(1000*j),col:(col+K-1)]=Z.mat
        p_matrix[[model]][[i]][(1000*(j-1)+1):(1000*j),sel_module+1]=p.mat
      }
    }
    
    cat("module",sel_module," Simulation: ", i, "\n")
    }
    col=col+K
}



saveRDS(Z_matrix,file ="Z_matrix.rds")
saveRDS(p_matrix,file ="p_matrix.rds")

Z_matrix=readRDS("Z_matrix.rds")
p_matrix=readRDS("p_matrix.rds")

library(PMA)
a = PMD(Z_matrix$`N=800`[[1]], K=9, sumabsu = 30, sumabsv = 5, type = "standard", center = FALSE)
View(a$u)
View(a$v)

res_error=NULL
for (sel_module in 0:8) {
  for(model in models){
    for (i in 1:N.sim) {
      res_error[[paste0("module=",sel_module)]][[model]][i]=c(rep(NA,N.sim))
    }
  }
}

for (sel_module in 0:8) {
  for(model in models){
    for (i in 1:N.sim) {
    res_error[[paste0("module=",sel_module)]][[model]][i]=sum(qvalue::qvalue(
      as.numeric(p_matrix[[model]][[i]][((sel_module)*1000+1):((sel_module+1)*1000),-(sel_module+1)]),pi0 = 1
    )$qvalue<0.05)/8000
    }
  }
}

res_power=NULL
for(model in models){
  for (i in 1:10) {
    res.pow=c()
    for (sel_module in 0:8) {
      res.pow=c(res.pow,as.numeric(p_matrix[[model]][[i]][((sel_module)*1000+1):((sel_module+1)*1000),(sel_module+1)]))
    }
    res_power[[model]][i]=sum(qvalue::qvalue(res.pow,pi0 = 1)$qvalue<0.05)/length(res.pow)
  }
}



sum(is.na(Z_matrix[["N=800"]][[10]]))

library(profvis)

p=profvis({
  
  Sigma_GG <- diag(1,nrow = 9000,ncol=9000) #9000*9000
  Sigma_EE <- diag(1,nrow=500,ncol=500) ##500*500 #block diag #identity
  Sigma_GE <- Z_matrix[["N=800"]][[1]] / sqrt(800) #9000*500
  res.alt2 <- archie_work(Sigma_GE, Sigma_GG, Sigma_EE,K=9, verbose = FALSE)
  saveRDS(res.alt2,file = "res_archie_test2.rds")
})
htmlwidgets::saveWidget(p, "profile_archie.html")

#Sigma_GG <- diag(1,nrow = 9000,ncol=9000) #9000*9000
#Sigma_EE <- diag(1,nrow=500,ncol=500) ##500*500 #block diag #identity
#Sigma_GE <- Z_matrix[["N=800"]][[1]] / sqrt(800) #9000*500
#res.alt <- archie_work(Sigma_GE, Sigma_GG, Sigma_EE,K=9, verbose = FALSE)


saveRDS(res.alt,file = "res_archie.rds")
res.alt=readRDS("res_archie.rds")

for(i in 4:5){
  Sigma_GG <- diag(1,nrow = 9000,ncol=9000) #9000*9000
  Sigma_EE <- diag(1,nrow=500,ncol=500) ##500*500 #block diag #identity
  Sigma_GE <- Z_matrix[["N=800"]][[i]] / sqrt(800) #9000*500
  res.alt <- archie_work(Sigma_GE, Sigma_GG, Sigma_EE,K=9, verbose = FALSE)
  file.name=paste0("res_archie",i,".rds")
  saveRDS(res.alt,file = file.name)
}

pheatmap::pheatmap(cor(Z_matrix[["N=800"]][[2]]),cluster_rows = FALSE,cluster_cols = FALSE)

pheatmap::pheatmap(Z_matrix[["N=800"]][[1]],cluster_rows = FALSE,cluster_cols = FALSE)
module_exp=datExpr[,dynamicMods==0]
cor(module_exp)

res.archie=readRDS("res_archie.rds")

View(res.archie$vs)


for (i in 0:8) {
  module=datExpr[,dynamicMods==i]
  Sigma=as.matrix(cor(module))
  saveRDS(Sigma,file = paste0("Sigma_module",i,".rds"))
}

num_nz=c()
for (i in 1:9) {
  num_nz[i]=sum(res.archie$us[,i]!=0)
}

ind=apply(res.archie$us, 2, function(x){
  index=which(x!=0)
  block=floor(index/1000)
  names(which.max(table(block)))
})

acc=c()
fdr=c()
for (i in 1:8) {
  index=as.numeric(ind[[i]])
  acc[i]=sum(res.archie$us[,i][(index*1000+1):((index+1)*1000)]!=0)/1000
  fdr[i]=sum(res.archie$us[,i][-((index*1000+1):((index+1)*1000))]!=0)/num_nz[i]
}

res_archie=readRDS("res_archie.rds")
res_archie_test=readRDS("res_archie_test2.rds")
res_us=res_archie_test$us
res_us[res_us!=0]=1
sum(res_us!=res_archie$us)


View(res_archie_test$us)
