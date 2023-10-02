library(magrittr)
library(tidyverse)
###SCRIPT FOR SDR-like SNPS
vcf<- read.table("Pstr_RNA_snps_freebayes_Q30_BA_NS_no_head.vcf", header=TRUE, sep="\t", stringsAsFactors = FALSE)

#Columns corresponding to each species/sex
Males <-c("PohiBra081_RNA", "PohiBra084_RNA", "PohiBra034","PohiBra037","PohiBra042","PohiBra043","PohiBra064","PohiBra065","PstrBra078","PstrBra079","PohiBra081","PohiBra084","PohiBra085", "PohiBra085_RNA")
Females <-c("PohiBra083_RNA", "PohiBra086_RNA", "PstrBra093_RNA", "PstrBra094_RNA", "PstrBra098_RNA", "PohiBra044","PohiBra045","PohiBra066","PstrBra076","PstrBra077","PohiBra080","PohiBra080_RNA", "PohiBra082","PohiBra083","PohiBra086","PohiBra088","PohiBra089")

# Convert from vcf format to 0,1,2 genotype format (where "1" is heterozygous)
procfunc = function(x) {
  t1 = sub(":.*","",x)
  if(t1 == "0/0") return(0)
  else if(t1 == "0/1" | t1 == "1/0") return(1)
  else if(t1 == "1/1") return(2)
  else return(-1)
}

pos = vcf$POS 
geno = data.frame(apply(vcf[,10:ncol(vcf)],1:2,procfunc),stringsAsFactors = FALSE)
names(geno) = names(vcf[,10:ncol(vcf)])
geno$POS<-vcf$POS
geno$CHR<-vcf$CHR
write.table(geno,file="Potamilus_RNA.geno",append=FALSE,quote=FALSE,sep="\t",row.names=FALSE) #Save genotype format file

geno<- read.table("Potamilus.geno",header=TRUE,sep="\t",stringsAsFactors = FALSE)
geno_RNA <- read.table("Potamilus_RNA.geno",header=TRUE,sep="\t",stringsAsFactors = FALSE)

geno_combined <- full_join(geno, geno_RNA)
geno_combined <- geno_combined %>% mutate_all(funs(replace_na(.,-1)))
head(geno_combined)
rm(geno)
rm(geno_RNA)
#Remove SNPs with all missing data for each sex
geno.nodels<-geno_combined[apply(geno_combined[,c(Males)],1,sum)!=-14 & apply(geno_combined[,c(Females)],1,sum)!=-17,]
geno.pw<-geno.nodels[,c("CHR","POS",Males,Females)]

#Need at least two males & two females to identify sex-linked SNPs 
geno.pw<-geno.pw[apply(geno.pw[,c(Males)],1,function(x) length(which(x=="-1")))<13,]
geno.pw<-geno.pw[apply(geno.pw[,c(Females)],1,function(x) length(which(x=="-1")))<16,]

#For each SNP, calculate frequencies of possible genotypes for each species/sex
geno.pw$males.hetfreq<-apply(geno.pw[,c(Males)],1,function(x) length(which(x=="1")))/apply(geno.pw[,c(Males)],1,function(x) length(which(x!="-1")))
geno.pw$males.00freq<-apply(geno.pw[,c(Males)],1,function(x) length(which(x=="0")))/apply(geno.pw[,c(Males)],1,function(x) length(which(x!="-1")))
geno.pw$males.11freq<-apply(geno.pw[,c(Males)],1,function(x) length(which(x=="2")))/apply(geno.pw[,c(Males)],1,function(x) length(which(x!="-1")))
geno.pw$females.hetfreq<-apply(geno.pw[,c(Females)],1,function(x) length(which(x=="1")))/apply(geno.pw[,c(Females)],1,function(x) length(which(x!="-1")))
geno.pw$females.00freq<-apply(geno.pw[,c(Females)],1,function(x) length(which(x=="0")))/apply(geno.pw[,c(Females)],1,function(x) length(which(x!="-1")))
geno.pw$females.11freq<-apply(geno.pw[,c(Females)],1,function(x) length(which(x=="2")))/apply(geno.pw[,c(Females)],1,function(x) length(which(x!="-1")))

###XY
#All females must be homozygous for same allele
geno.filter_XY <-geno.pw[geno.pw$females.00freq==1 | geno.pw$females.11freq==1, ]

#All males must be heterozygous or homozygous for same allele as females
geno.filter_XY <-geno.filter_XY[geno.filter_XY$males.hetfreq==1 & geno.filter_XY$females.hetfreq==0,]
#geno.filter_XY_candidates <- geno.filter_XY %>% group_by(CHR) %>% filter(n() >= 25)
gene_info <- read.csv("gene_geno_exp.csv")
contig_list <- c(geno.filter_XY$CHR) %>% unique()
gene_info_match <- subset(gene_info, Contig %in% c(contig_list))

library(fuzzyjoin)
genes_XY_SNPs <- fuzzy_left_join(geno.filter_XY, gene_info_match,
                                       by = c(
                                         "CHR" = "Contig",
                                         "POS" = "Start",
                                         "POS" = "End"), 
                                       match_fun = list(`==`, `>=`, `<=`))

write.table(geno.filter_XY,file="Potamilus_XY_sites.txt",append=FALSE,quote=FALSE,sep="\t",row.names=FALSE) #Save genotype format file
write.table(genes_XY_SNPs,file="Potamilus_XY_gene_sites.txt",append=FALSE,quote=FALSE,sep="\t",row.names=FALSE) #Save genotype format file

###ZW
#All males must be homozygous for same allele
geno.filter_ZW <-geno.pw[geno.pw$males.00freq==1 | geno.pw$males.11freq==1, ]

#All females must be heterozygous or homozygous for same allele as females
geno.filter_ZW <-geno.filter_ZW[geno.filter_ZW$females.hetfreq==1 & geno.filter_ZW$males.hetfreq==0, ]
contig_list_f <- c(geno.filter_ZW$CHR) %>% unique()
gene_info_match <- subset(gene_info, Contig %in% c(contig_list_f))
genes_ZW_SNPs <- fuzzy_left_join(geno.filter_ZW, gene_info_match,
                                 by = c(
                                   "CHR" = "Contig",
                                   "POS" = "Start",
                                   "POS" = "End"), 
                                 match_fun = list(`==`, `>=`, `<=`))
geno.filter_ZW_candidates <- geno.filter_ZW %>% group_by(CHR) %>% filter(n() >= 25)
geno.filter_ZW_candidates %>% group_by(CHR) %>% count()

write.table(geno.filter_ZW,file="Potamilus_ZW_sites.txt",append=FALSE,quote=FALSE,sep="\t",row.names=FALSE) #Save genotype format file
write.table(geno.filter_ZW,file="Potamilus_gene_ZW_sites.txt",append=FALSE,quote=FALSE,sep="\t",row.names=FALSE) #Save genotype format file