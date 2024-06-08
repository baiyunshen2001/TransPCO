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


module=datExpr[,dynamicMods==0]

dir_pco="./simulation/script_lambda0.1/"
source(paste0(dir_pco, "ModifiedPCOMerged.R"))
source(paste0(dir_pco, "liu.R"))
source(paste0(dir_pco, "liumod.R"))
source(paste0(dir_pco, "davies.R"))
dyn.load(paste0(dir_pco, "qfc.so"))
source(paste0(dir_pco, "ModifiedSigmaOEstimate.R"))

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
N.sim <- 100



N.seq <- c(200, 400, 600, 800)
models <- paste0("N=", N.seq)
#c("u1", "uPCO", "uK", paste0("N=", N.seq), "100%", "70%", "30%")


res.alt <- NULL
p_alt_all <- array(dim = c(N.sample, 4, length(models), N.sim),
                   dimnames = list(NULL, c("Oracle", "PC1", "minp", "PCO"), models, NULL))
for(i in 1:N.sim){
  B <- matrix(rep(NA, K*length(models)), ncol = length(models),
              dimnames = list(NULL, models))
  #beta.n <- as.matrix(rnorm(K, sd = sqrt(var.b)))
  beta.n <- rbind(as.matrix(rnorm(floor(caus*K), sd = sqrt(var.b))), as.matrix(rep(0, K-floor(caus*K))))
  B[, models] <- beta.n %*% sqrt(N.seq)
  
  for(model in models){
    Z.alt <- rmvnorm(N.sample, B[, model], Sigma)
    
    # PC1
    method.tmp <- "PC1"
    PC1 <- Z.alt %*% eigen_vec[, 1]
    p_alt_all[, method.tmp, model, i] <- as.numeric(2*pnorm(-abs(PC1/sqrt(eigen_vec[1]))))
    
    # univariate minp
    method.tmp <- "minp"
    p_alt_all[, method.tmp, model, i] <- apply(Z.alt, 1, function(x) min(1-pchisq(x^2, 1))*length(x) )
    
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
  p.null <- data.table(paste0("null", 1:length(p_null_all[[paste0("p.null.", method.tmp)]])), p_null_all[[paste0("p.null.", method.tmp)]])
  
  MyArray <- p_alt_all[ , method.tmp, , , drop = FALSE]
  res[[method.tmp]] <- apply(MyArray, 3, function(x) {
    tmp_p = apply(x, 3, function(y) list(as.numeric(y)))
    
    pbmcapply::pbmclapply(tmp_p, function(z) {
      p.obs = data.table(paste0("obs", 1:N.sample), unlist(z) )
      p.obs.rank = frank(p.obs, V2)
      names(p.obs.rank) = p.obs$V1
      
      all.rank = frank(rbindlist(list(p.obs, p.null)), V2)
      names(all.rank) = c(p.obs$V1, p.null$V1)
      q = pmin((all.rank[names(p.obs.rank)]/p.obs.rank-1)/C, 1)
      
      sum(q < fdr_level)/N.sample
    }
    , mc.cores = 20) %>%
      list_c() %>%
      list()
  }) %>%
    setNames(dimnames(MyArray)[[3]])
  
  

}

res_power=c()
for (method in c("PCO","PC1","minp")) {
  for (sample_size in c(200,400,600,800)) {
    res_power=c(res_power,mean(res[[method]][[paste0("N=",sample_size)]][[1]]))
  }
}
  
res.alt=data.frame(power=res_power,
                   Method=rep(c("PCO","PC1","MinP"),each=4),
                   sample_size=rep(c(200,400,600,800),3))

png(filename = "power_simulation_samplesize.png",width = 500,height = 300)
ggplot(data=res.alt,aes(x=sample_size,y=power,color=Method,group=Method))+
  geom_line(size=1)+
  geom_point()+
  labs(x="Sample Size",y="Power",title = " Power of three methods across different sample sizes")
dev.off()
