library(WGCNA)
library(flashClust)
library(mvtnorm)
library(data.table)
library(tidyverse)

dir_pco="./simulation/script_lambda0.1/"
source(paste0(dir_pco, "ModifiedPCOMerged.R"))
source(paste0(dir_pco, "liu.R"))
source(paste0(dir_pco, "liumod.R"))
source(paste0(dir_pco, "davies.R"))
dyn.load(paste0(dir_pco, "qfc.so"))
source(paste0(dir_pco, "ModifiedSigmaOEstimate.R"))

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

for (sel_module in 8) {
  


module=datExpr[,dynamicMods==sel_module]



Sigma=as.matrix(cor(module))
K=dim(Sigma)[1]

n_sim <- 10^3
p_null_all  <- list()
z_null <- rmvnorm(n_sim, rep(0, K), Sigma)

SigmaO=ModifiedSigmaOEstimate(Sigma,method = "liu")
eigen_res <- eigen(Sigma)
eigen_lamb <- eigen_res$values
eigen_vec <- eigen_res$vectors

p_null_all$'p.null.PCO' <- ModifiedPCOMerged(
  Z.mat = z_null, Sigma = Sigma, SigmaO = SigmaO, method = "liu"
) |> as.numeric()

PC1 <- z_null %*% eigen_vec[, 1]
p_null_all$'p.null.PC1' <- 2*pnorm(-abs(PC1/sqrt(eigen_lamb[1])))|> as.numeric()

p_null_all$'p.null.minp' <- apply(z_null, 1, function(x) min(1-pchisq(x^2, 1))*length(x) )



oracle_thre <- 0.1

var.b <- 0.001
caus <- 0.5
N <- 500

N.sample <- 10^4
N.sim <- 10



N.seq <- c(200, 400, 600, 800)
models <- paste0("N=", N.seq)
#c("u1", "uPCO", "uK", paste0("N=", N.seq), "100%", "70%", "30%")


res.alt <- NULL
z.alt.mat=NULL
p_alt_all <- array(dim = c(N.sample, 4, length(models), N.sim),
                   dimnames = list(NULL, c("Oracle", "PC1", "minp", "PCO"), models, NULL))

B <- matrix(rep(NA, K*length(models)), ncol = length(models),
            dimnames = list(NULL, models))
#beta.n <- as.matrix(rnorm(K, sd = sqrt(var.b)))
beta.n <- rbind(as.matrix(rnorm(floor(caus*K), sd = sqrt(var.b))), as.matrix(rep(0, K-floor(caus*K))))
B[, models] <- beta.n %*% sqrt(N.seq)
for(i in 1:N.sim){
  
  
  for(model in models){
    Z.alt <- rmvnorm(N.sample, B[, model], Sigma)
    z.alt.mat[[model]][[i]]=Z.alt
    # PC1
    method.tmp <- "PC1"
    PC1 <- Z.alt %*% eigen_vec[, 1]
    p_alt_all[, method.tmp, model, i] <- as.numeric(2*pnorm(-abs(PC1/sqrt(eigen_lamb[1]))))
    
    # univariate minp
    method.tmp <- "minp"
    p_alt_all[, method.tmp, model, i] <- apply(Z.alt, 1, function(x) min(1-pchisq(x^2, 1)) )
    
    # PCO
    method.tmp <- "PCO"
    p_alt_all[, method.tmp, model, i] <- as.numeric(ModifiedPCOMerged(Z.mat=Z.alt, Sigma=Sigma, SigmaO=SigmaO,method = "liu"))
  }
  cat("Simulation: ", i, "\n")
}
  
  

fdr_level =0.05
N.sample <- dim(p_alt_all)[1]
C <- length(p_null_all[[1]])/N.sample
res <- NULL
for (method.tmp in c("PCO", "PC1", "minp")) {
  for (n_seq in paste0("N=",N.seq)) {
  
    for (i in 1:N.sim) {
      q=qvalue::qvalue(p_alt_all[,method.tmp,n_seq,i],pi0 = 1)$qvalue
      power=sum(q < fdr_level)/N.sample
      res[[method.tmp]][[n_seq]][i]=power
    }
  }
}


res_power=c()
res_power_sd=c()
for (method in c("PCO","PC1","minp")) {
  for (sample_size in c(200,400,600,800)) {
    res_power=c(res_power,mean(res[[method]][[paste0("N=",sample_size)]]))
    res_power_sd=c(res_power_sd,sd(res[[method]][[paste0("N=",sample_size)]]))
  }
}
  
res.alt=data.frame(power=res_power,
                   power_sd=res_power_sd,
                   Method=rep(c("PCO","PC1","MinP"),each=4),
                   sample_size=rep(c(200,400,600,800),3))

file_name=paste0("power_simulation_samplesize_module",sel_module,"test",".png")
png(filename = file_name,width = 500,height = 300)

ggplot(data=res.alt,aes(x=sample_size,y=power,color=Method,group=Method))+
  geom_line(size=1)+
  geom_point()+
  geom_errorbar(aes(ymin=power-power_sd,ymax=power+power_sd),width = 50)+
  labs(x="Sample Size",y="Power",title = paste0("Power of three methods across different sample sizes module",sel_module))

dev.off()

#image_name=paste0("full",sel_module,".Rdata")
#save.image(file=image_name)
}

pheatmap::pheatmap(Z_matrix[1:1000,],cluster_rows = FALSE,cluster_cols = FALSE)


p_val=runif(50)
qvalue::qvalue(p_val,fdr.level = fdr_level)$qvalues


p_alt_all[,"PCO","N=800",8]
z.alt.mat[["N=800"]][[8]]
