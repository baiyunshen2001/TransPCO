file.ex = snakemake@input[["file_ex"]]
file.gene.annotation = snakemake@input[['file_gene_annotation']]
file.mappability = snakemake@input[['file_mappability']]
file.cross.mappability = snakemake@input[['file_cross_mappability']]
file.covariates = snakemake@input[["file_covariates"]]

file.gene.meta = snakemake@output[["file_gene_meta"]]
file.ex.var.regressed = snakemake@output[["file_ex_var_regressed"]]


require(data.table)

## Part 1: filter genes

### read expression, meta and mappability files
ex = readRDS(file.ex)
gene.annotation = fread(file.gene.annotation, header = TRUE)
low.mapp = fread(file.mappability)
cross.map = fread(file.cross.mappability)

### decide whether the gene names are gene symbols or EMSEMBLE ids.
gnames = sapply(colnames(ex), function(x) strsplit(x, "\\|")[[1]][1])
if(length(grep("^ENSG", gnames)) ==  length(gnames) ){
  cat("Gene names are ENS ids!", "\n")
  
  gnames.meta = gnames[gnames %in% gene.annotation$Geneid]
  rem.pseu = !(gene.annotation[match(gnames.meta, gene.annotation$Geneid), Class] %in% c("protein_coding", "lincRNA"))
  rem.low.mapp = gnames.meta %in% low.mapp[low.mapp$V2 < 0.9, V1]
  rem.cross.map = gnames.meta %in% unique(c(cross.map$V1, cross.map$V2))
  rem.auto = !(gene.annotation[match(gnames.meta, gene.annotation$Geneid), Chromosome] %in% paste0("chr", 1:22))
  
  GeneNameConv = gene.annotation[match(gnames.meta, gene.annotation$Geneid), GeneSymbol]
  GeneName = gnames.meta
  gnames.ensb = GeneName

}else{
  cat("Gene names are gene symbols!", "\n")
  
  ### remove genes not in annotation file
  gnames.meta = gnames[gnames %in% gene.annotation$GeneSymbol]
  
  ### remove genes that are not protein coding or lincRNA, including all types of pseudogenes.
  rem.pseu = !(gene.annotation[match(gnames.meta, gene.annotation$GeneSymbol), Class] %in% c("protein_coding", "lincRNA"))
  
  ### remove poorly mapped genes, i.e. with mappability < 0.9.
  rem.low.mapp = gnames.meta %in% gene.annotation[match(low.mapp[low.mapp$V2 < 0.9, V1], gene.annotation$Geneid), GeneSymbol]
  
  ### remove cross-mappable genes, i.e. gene pairs with cross mappability > 0
  rem.cross.map = gnames.meta %in% gene.annotation[match(unique(c(cross.map$V1, cross.map$V2)), gene.annotation$Geneid), GeneSymbol]
  
  ### remove genes that are not autosomal.
  rem.auto = !(gene.annotation[match(gnames.meta, gene.annotation$GeneSymbol), Chromosome] %in% paste0("chr", 1:22))
  
  GeneName = colnames(ex)[match(gnames.meta, gnames)]
  GeneNameConv = gene.annotation[match(gnames.meta, gene.annotation$GeneSymbol), Geneid]
  gnames.ensb = GeneNameConv

}

### write out genes meta info and removal info
res = data.frame("gene" = GeneName,
                 "chr" = gene.annotation[match(gnames.ensb, gene.annotation$Geneid), Chromosome],
                 "start" = gene.annotation[match(gnames.ensb, gene.annotation$Geneid), Start],
                 "end" = gene.annotation[match(gnames.ensb, gene.annotation$Geneid), End],
                 "GeneNameConv" = GeneNameConv,
                 "Remove" = (rem.pseu | rem.low.mapp | rem.cross.map | rem.auto),
                 "RemovePseu" = rem.pseu,
                 "RemoveLowmapp" = rem.low.mapp,
                 "RemoveCrossmapp" = rem.cross.map,
                 "RemoveAuto" = rem.auto)
write.table(as.matrix(res),
            file.gene.meta,
            sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)


## Part 2: regress out covariates

### Extract genes
gnames.use = GeneName[!(rem.pseu | rem.low.mapp | rem.cross.map | rem.auto)]
datExpr = ex[, gnames.use]

### read covariates
cov_all = t(as.matrix(read.table(file.covariates,
                                 sep = "\t", header = TRUE, row.names = 1,
                                 stringsAsFactors = FALSE, check.names = FALSE)))

### regress out covariates
extract_residual <- function(y, x){
  return(lm(y ~ x)$residuals)
}
ex_cov_regressed = apply(datExpr, 2, function(y) extract_residual(y, cov_all))

### save expression matrix with covariates regressed out
saveRDS(ex_cov_regressed, file.ex.var.regressed)


### print out message
rem.N = sum(rem.pseu | rem.low.mapp | rem.cross.map | rem.auto)
cat("All genes:", length(gnames), '\n',
    "With meta:", length(gnames.meta), '\n',
    "Final use:", length(gnames.meta) - rem.N, '\n',
    "All removed:", rem.N+length(gnames)-length(gnames.meta), '\n',
    "including:", c(length(gnames)-length(gnames.meta), sum(rem.auto), sum(rem.pseu), sum(rem.low.mapp), sum(rem.cross.map)), '\n')
