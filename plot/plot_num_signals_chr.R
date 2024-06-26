######################################################
#################### Plot #unique SNPs or #indep SNPs ####################
#################### v.s. chr or module ####################
######################################################
rm(list = ls())
library(data.table)
library(tidyverse)


### read the indep_SNP file and construct df for plotting
file_signal <- '~/xuanyao_llw/DGN_no_filter_on_mappability/FDR/signals.chr.module.perm10.fdr10.txt'
file_indep_SNP <- '/scratch/midway2/liliw1/figures/fig0/fig0_data/num_indep_SNP_chr.txt'

signal <- fread(file_signal, col.names = c("module_SNP", "p", "q") )
df_indep_SNP <- fread(file_indep_SNP, header = FALSE, col.names = c('pair', 'num_indep_SNP') )

# indep loci
df_indep_SNP <- df_indep_SNP %>%
  separate(pair, c(NA, "chr", NA), sep = '[.]', convert = TRUE)

# uniq SNPs
signal <- signal %>%
  separate(module_SNP, c("module", "chr", "pos"), ":", remove = FALSE, convert = TRUE) %>%
  unite("SNP_ID", chr, pos, sep = ":", remove = FALSE) %>%
  select(module, chr, SNP_ID)
df_uniq_SNP <- signal %>% group_by(chr) %>% distinct(SNP_ID) %>% summarise(num_uniq_SNP = n())

# combine two types of SNPs for plotting
df_plot_signal <- left_join(df_uniq_SNP, df_indep_SNP, by = "chr") %>%
  pivot_longer(c(num_uniq_SNP, num_indep_SNP), "signal_type", values_to = "num_SNP")

df_plot_signal$chr <- factor(df_plot_signal$chr, levels = 22:1)
df_plot_signal$signal_type <- factor(df_plot_signal$signal_type,
                                     levels = c("num_uniq_SNP", "num_indep_SNP"),
                                     labels = c("SNP", "Independent Loci") )

# only label modules with many indep loci, as too many modules
fig <- ggplot(df_plot_signal, aes(x = num_SNP, y = chr)) +
  geom_line(aes(x = num_SNP-1, group = chr),
            size = 0.4, color = "grey") +
  geom_point(aes(color = signal_type),
             size = 4, alpha = 0.8) +
  geom_text(data = subset(df_plot_signal, signal_type == "Independent Loci"),
            aes(label = num_SNP),
            color = "#2171B5", size = 4) +
  labs(x = "Number of Loci/SNPs", y = "Chromosome", color = "Type") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        
        legend.position = "top",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10, face = "bold"),
        legend.background = element_rect(color = "black", linetype = "dashed"),
        legend.key.size= unit(0.5, "cm"),
        
        axis.line = element_line(colour="black", size = 0.7),
        
        plot.margin=unit(c(10,5,5,5),"mm"),
        axis.text.x = element_text(size = 10, color = "black",
                                   vjust = 0.5, hjust = 1, angle = 90),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.title.x = element_text(vjust = -0.2, size=14),
        axis.title.y = element_text(angle=90,vjust =2, size=14)
  ) +
  scale_color_manual(values = c("#df9aa1", "#9ECAE1"))

fig

### save figure and plotting object for further editing
saveRDS(fig, "signal_chr.rds")

ggsave('signal_chr.pdf', fig, width = 6, height = 5)

