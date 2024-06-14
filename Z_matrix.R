library(WGCNA)
library(flashClust)
library(mvtnorm)
library(data.table)
library(tidyverse)
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
Z_matrix=matrix(0,nrow = 9000,ncol = 500)
col=1
for (i in 0:8) {
  module=datExpr[,dynamicMods==i]
  Sigma=as.matrix(cor(module))
  K=dim(Sigma)[1]
  for (j in 1:9) {
    
    if(j==i+1){
      beta.n <- rbind(as.matrix(rnorm(floor(caus*K), sd = sqrt(var.b))), as.matrix(rep(0, K-floor(caus*K))))
      B=beta.n *sqrt(800)
      Z.mat= rmvnorm(1000, B, Sigma)
    }
    else{
      Z.mat= rmvnorm(1000, rep(0, K), Sigma)
    }
    Z_matrix[(1000*(j-1)+1):(1000*j),col:(col+K-1)]=Z.mat
  }
  col=col+K
}

saveRDS(Z_matrix,fil="Z_matrix.rds")
write.table(Z_matrix,file = "Z_matrix.csv")
