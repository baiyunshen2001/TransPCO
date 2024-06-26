# prepare gwas data
{
  rm(list = ls())

  dir_data = "./p/"
  chr_seq = 1:22


  pvalue_pc = NULL;
  chr_pc = NULL; bp = NULL; snp = NULL;
  for(chr in chr_seq){
    tmp_pc = readRDS(file = paste0(dir_data, "p.module1.chr", chr, ".rds"))
    #names(tmp_pc) = paste0("C", gene_cluster_id, ":", names(tmp_pc))
    pvalue_pc = c(pvalue_pc, tmp_pc)

    chr_pc = c(chr_pc, rep(chr, length(tmp_pc)))

    bp = c(bp, sapply(strsplit(names(tmp_pc), ":"), function(x) as.numeric(x[[2]])))
    snp = c(snp, names(tmp_pc))

    cat("Chr:", chr, "\n")
  }

  gwasResults = data.frame("CHR"=chr_pc, "BP"=bp, "SNP"=snp, "P"=pvalue_pc)
  saveRDS(gwasResults, file = "manhattan.rds")

}




# prepare signals and pvalues
{
    PCO_trans = read.table("signals.chr.module.perm10.txt", header = FALSE, col.names = c("module.snp", "pvalue", "qvalue"), row.names = 1, stringsAsFactors = FALSE)
    PCO_trans$SNP = sapply(strsplit(rownames(PCO_trans), ":"), function(x) paste(x[[2]], x[[3]], sep = ":"))

    PCO_trans_uniq = NULL
    for(i in unique(PCO_trans$SNP)){
      tmp = PCO_trans[PCO_trans$SNP == i, ]
      PCO_trans_uniq = rbind(PCO_trans_uniq, tmp[which.min(tmp$pvalue), ])
    }

    write.table(PCO_trans_uniq, file = "PCO.trans.uniq.txt", quote = FALSE)

}



# prepare annotation files
{
  require(data.table)

  PCO_trans_indep = read.table("indep.signals.chr.module.perm10.txt", stringsAsFactors = FALSE, col.names = "SNP")
  PCO_trans_indep$PCO = PCO_trans_indep$SNP; PCO_trans_indep$elife = "no"; PCO_trans_indep$both = "no"
  PCO_trans_indep[PCO_trans_indep$SNP=="3:56849749", c("elife", "both")] = "3:56849749"
  #PCO_trans_indep[PCO_trans_indep$SNP=="9:126985334", c("elife", "both")] = "9:126985334"

  gene_pos = fread(paste0("/project2/xuanyao/data/mappability/gencode.v19.annotation.table.txt"),
                        sep = '\t', header = T)
  gene_pos = as.data.frame(gene_pos[gene_pos$Class %in% c('protein_coding', 'lincRNA'), ])

  PCO_trans_indep$neargene = "NA"
  for(i in PCO_trans_indep$SNP){
    chr = as.numeric(strsplit(i, ":")[[1]][1])
    snpos = as.numeric(strsplit(i, ":")[[1]][2])

    {
      ind2 = gene_pos$Chromosome == paste0("chr", chr) & abs(gene_pos$Start - snpos)<5*10^5
      dis_snp_gene2 = gene_pos[ind2, "Start"]-snpos
      tmp2 = cbind(gene_pos[ind2, ], 'dis_snp_gene' = dis_snp_gene2)

      if(nrow(tmp2)>0){
        PCO_trans_indep[PCO_trans_indep$SNP==i, "neargene"] = tmp2[order(abs(tmp2$dis_snp_gene)), ][1, "GeneSymbol"]
      }
    }
  }

  #gene_pos[match(PCO_trans_indep$neargene, gene_pos$GeneSymbol), ]

  write.table(PCO_trans_indep, "PCO.trans.indep.txt", quote = FALSE, row.names = FALSE)


}




### plot manhattan plot
cd ./scratch/manhattan_tmp/
module load R/3.6.1
R

rm(list = ls())
library(qqman)
library(ggrepel)
library(dplyr)

gwasResults = readRDS("manhattan.rds")
PCO_trans_uniq = read.table(file = "PCO.trans.uniq.txt", stringsAsFactors = FALSE, row.names = 1, header = T)
PCO_trans_indep = read.table("PCO.trans.indep.txt", header = TRUE, stringsAsFactors = FALSE)
PCO_trans_uniq[PCO_trans_uniq$pvalue==0, "pvalue"] = 1e-18 #pvalues<1e-17 is displayed as 0
gwasResults$P = 1; gwasResults[match(PCO_trans_uniq$SNP, gwasResults$SNP), "P"] = PCO_trans_uniq$pvalue
PCO_trans_indep[25, 'neargene']='IKZF1'
{
  don <- gwasResults %>%

    # Compute chromosome size
    group_by(CHR) %>%
    summarise(chr_len=max(BP)) %>%

    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%

    # Add this info to the initial dataset
    left_join(gwasResults, ., by=c("CHR"="CHR")) %>%

    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( chromosome=BP+tot) %>%

    # Add highlight and annotation information
    mutate( is_highlight=ifelse(SNP %in% PCO_trans_indep$SNP, "yes", "no")) %>%
    mutate( is_annotate=ifelse(SNP %in% PCO_trans_indep$SNP, "yes", "no")) %>%
    mutate( is_elife=ifelse(SNP %in% PCO_trans_indep$elife, "yes", "no") )


  axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(chromosome) + min(chromosome) ) / 2 )
  label = PCO_trans_indep[match(don[don$is_annotate=="yes", "SNP"], PCO_trans_indep$SNP), "neargene"]

  manhattan_plot <-   ggplot(don, aes(x=chromosome, y=-log10(P))) +

    # Show all points
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("grey", "grey"), 22 )) +

    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = c(5, 20) ) +     # remove space between plot area and x axis
    # Add highlighted points
    geom_point(data=subset(don, is_highlight=="yes"), color="red", size=1.5) +
    geom_point(data=subset(don, is_elife=="yes"), color="red", size=3, shape=2) +
    # Add label using ggrepel to avoid overlapping
    geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=label), size=3, fill = NA) +
    # Custom the theme:
    theme_bw() +
    theme(
      text = element_text(size=15),
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  ggsave("manhattan.pdf", manhattan_plot)
}

