library(tidyverse)
library(dplyr)
library(magrittr)

#Read in files
coding_genes <- as.list(read.table("Protein_coding_genes", header = T))
gene_exp <- read.csv("nuclear_deseq_results_notrna.csv")
gene_geno <- read.delim("../Psass/psass_gene_metrics.tsv") %>% subset(ID %in% c(coding_genes$ID))
gene_info <- inner_join(gene_geno, gene_exp, by = "ID")
write.csv(gene_info, "gene_geno_exp.csv")
write.csv(gene_exp, "nuclear_deseq_results_notrna.csv")
DEGs <- subset(gene_exp, padj<0.05)

#Filter genomic regions without 3x depth
reseq <- read.delim("../Psass/psass_window_1kb.tsv")
filter_reseq <- subset(reseq, c(Abs_depth_females>=33,Abs_depth_males>=33))
filter_reseq_nm <- filter_reseq %>% drop_na(Fst)
rm(reseq)
rm(filter_reseq)

range(filter_reseq_nm$Fst)
high_fst <- subset(filter_reseq_nm, Fst>=0.1)
rm(filter_reseq_nm)

contig_list <- c(high_fst$Contig) %>% unique()

gene_info_match <- subset(gene_info, Contig %in% c(contig_list))

library(fuzzyjoin)
combined_reseq_gene <- fuzzy_left_join(high_fst, gene_info_match,
           by = c(
             "Contig" = "Contig",
             "Position" = "Start",
             "Position" = "End"), 
           match_fun = list(`==`, `>=`, `<=`))

write.csv(combined_reseq_gene, "psass_combined_high_FST_wgene_info_1kb.csv")


highfst_express <- read.csv("psass_combined_high_FST_wgene_info_1kb.csv", header = T)
unique_genes <- c(highfst_express$ID) %>% unique()
sum(!is.na(highfst_express$ID))
write(unique_genes, "high_FST_genes.txt")

highfst_exp_peaks <- highfst_express %>% subset(!is.na(highfst_express$log2FoldChange), .keep.all=T)

distinct(highfst_exp_peaks, ID, .keep.all = TRUE)

#Check for twin peaks
ggplot()+
  geom_density(aes(highfst_exp_peaks$log2FoldChange), size = 0.7)+
  xlim(-10, 10)+
  ylab("Density")+
  xlab("Log2 Fold-Change")+
  theme_classic()

#Check DGE and high fst
DEG_High_FST <- inner_join(DEGs, highfst_exp_peaks)

#Get list of biased genes
kb_highfst <- read.csv("psass_combined_high_FST_wgene_info.csv", header = T)
kb_highfst <- kb_highfst %>% drop_na(gene_id)
range(kb_highfst$log2.fold_change)

female_genes <- subset(kb_highfst, log2.fold_change. >= 1) 
unique_female_genes <- c(female_genes$gene_id) %>% unique()
write(unique_female_genes, "female_biased_genes.txt")

male_genes <- subset(kb_highfst, log2.fold_change. <= -1) 
unique_male_genes <- c(male_genes$gene_id) %>% unique()
write(unique_male_genes, "male_biased_genes.txt")

library(ggplot2)

