#!/usr/bin/env Rscript

#============================#
#============================#
# ------ Load libraries ----
#============================#
#============================#
library(IRanges)
library(vcfR)
library(ggplot2)
library(adegenet)
library(hierfstat)
library(ggbreak)
library(IRanges)
library(stringr)
library(basetheme)
library(devtools)
library(GenomicRanges)

system2("git clone https://github.com/EliseGAY/Package_VCF2PopStructure.git")
load_all("/Package_VCF2PopStructure/")

#===============================#
#===============================#
# ------ PREPARE YOUR DATA ----
#===============================#
#===============================#

#--------------#
# Metadata pop
#--------------#

# the pop table has to be ordered in the same way as all the VCF header
#'''
#	pop	samples
#	pop1 sample_1
#	pop1 sample_2
#	pop1 sample_3
#	pop2 sample_4
#	pop2 sample_5
#'''

# read metada
metadata=read.table("metadata/metadata.txt", header = TRUE)
metadata=as.data.frame(metadata)
pop=unique(metadata$Population)

# read chr length
table_chr=read.table("Scaff_length_sorted.list", header = T)
#-----------------------------------------------------------#
# Generate Genotype tables needed in different R packages
#-----------------------------------------------------------#

# Read the VCF with vcfR :
VCFR_data=read.vcfR("data/VCF_example.vcf.gz")

# create a pop sorted by VCF colnames
metadata_sorted <- metadata[match(colnames(VCFR_data@gt)[-1], metadata$GT_sample),]
pop_list = split(metadata_sorted$GT_sample, metadata_sorted$Population)

# current_chr
args <- commandArgs(trailingOnly = TRUE)
chr <- as.character(args[1])
#===================================#
#===================================#
# ------ ACP Sliding windows-------
#===================================#
#===================================#

# INPUT
#------#
# VCF input : VCF filtered for Na and MAF

# Methods
#----------------#
# Create a table geno_table  with pos (col), ind (rows) and genotype as 0,1,2
# window : the size of the window (int)
# slide : the size of the slide (int)
# min_n_snp: minimum Nb of SNPs to compute PCA (100 should be enough, it is just to avoid error at the bounds of the scaffold)

# OUTPUT
#---------#
# The Get_sw_abs return a table of sliding windows coordinates made from 1 to the last SNPs position :
# Fine-tune the wind and slide 
#'  chr        start   end mid_point nb_snp
#'<chr>      <int> <int>     <dbl>  <int>
#'1 ptg000007l     1 10000     5000.     21
#'2 ptg000007l  2001 12000     7000.     60
#'3 ptg000007l  4001 14000     9000.     73
#'4 ptg000007l  6001 16000    11000.     87

# The final table is 
# PCs_table<-cbind(low_bound,upper_bound,pca_asse1, pca_asse2)	

# in "PCs" you have a matrix of four columns: 
# col1 = lower value of the window
# col2 = upper value of the window
# col3 = mid value of the window
# col3 = % of variance of the first axe
# col4 = % of variance of the second axe

# Create directory
if (!dir.exists("sliding_PCA")) {
  dir.create("sliding_PCA")
}

# get genotype table
loci_table = extract.gt(VCFR_data, element = "GT")
loci_table = as.data.frame(loci_table)
geno_table = t(Convert_GT(GT_table = loci_table))

# geno_table[c(1:10),c(1:10)]
#'                ptg000007l_8628 ptg000007l_8631 ptg000007l_8633 ptg000007l_8639 ptg000007l_9618 ptg000007l_9622 ptg000007l_9625 ptg000007l_9626
#'1622W1_S254               0               0               0               0              NA              NA              NA              NA
#'1623W1_S276               0               0               0               0              NA              NA              NA              NA
#'1624W1_S392               0               0               0               0              NA               0               0               0
#'1625Q_S431                0               0               0               0              NA              NA              NA              NA
#'1625W1_S384               0               0               0               0               1               1               1               1

# get SW table :
sw_data = Get_SW_abs(VCFR_data, slide = 5000, window = 10000)

# Run PCA 
min_snp=100
for(i in 1:nrow(sw_data)){
  if(sw_data[i,]$nb_snp >= min_snp){
    pos_start = sw_data[i, "start"]
    pos_end = sw_data[i, "end"]
    indexes_col = which(getPOS(VCFR_data) >= pos_start & getPOS(VCFR_data) < pos_end)
                  
    sub_geno = as.data.frame(geno_table[,indexes_col])
    pca = indpca(sub_geno, scale = T)
    PCs = pca$ipca$eig / sum(pca$ipca$eig) * 100
    sw_data[i,"PC1"] = PCs[1]
    sw_data[i,"PC2"] = PCs[2]
  }
  else{print("no enough snp")}
}

# wite table :
sw_chr_name = paste("PCs_Sliding" , chr, "table", sep = "_")
write.table(sw_data, sw_chr_name)

# Plot PC1 variance along chr
#-----------------------------#
# read the table if needed
# sw_data=read.table(paste("Sliding_PCA/PCs_Sliding", , chr, "table" sep = "_"), header=TRUE)

# get current chr length
length_chr=table_chr[which(table_chr$scaffold==chr),"length"]
  
# plot 
p=ggplot() +
  geom_point(aes(x=sw_data$mid_point, y=sw_data$PC1),
             color = "darkgreen",
             alpha = 0.5,
             size=1)+
  
  labs(x="", y = "PC 1") +
  
  ggtitle(chr) +
  theme(text = element_text(size = 10), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text = element_text(size = 10),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5)) +
  
  scale_x_continuous(breaks = seq(0, length_chr, 
                                  by = round(length_chr/15))) +
  
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(1,100))

p

ggsave(p, path="sliding_PCA/", width = 12, height = 10, device = "pdf", filename = paste(chr, "PCA_scans.pdf",sep = "_"))

#===================================#
#===================================#
# ------ ACP on SNPs range -------
#===================================#
#===================================#

# INPUT
#---------#
# VCF input : vcfR vcf

# Methods :
#---------#
# create a sliding windows table based on number of SNP (and not chr length)
# Segment genotype table according to windows

# OUTPUT
#---------#
# Return the final table with windows chosen and their corresponding PCA values

# In 'table_PCA' you have a matrix of four col: 
# col1 = start position of the window
# col2 = end position of the window
# col3 = Mid size of the windows
# col3 = % of variance of the first axe
# col4 = % of variance of the second axe

#----------------------------#
# Run PCA on range SNP
#----------------------------#

# Get SNP range sliding windows 
# Set the value of nb snp windows and nb snp slide 
sw_snp_data = Get_SW_SNPrange(VCF = VCFR_data)

# compute PCA
for(i in 1:nrow(sw_snp_data)){
    pos_start = sw_snp_data[i, "start"]
    pos_end = sw_snp_data[i, "end"]

    sub_geno = as.data.frame(geno_table[,c(pos_start : pos_end)])
    pca = indpca(sub_geno, scale = T)
    PCs = pca$ipca$eig / sum(pca$ipca$eig) * 100
    sw_snp_data[i,"PC1"] = PCs[1]
    sw_snp_data[i,"PC2"] = PCs[2]
}

sw_snp_data
#      chr         start   end mid_point  PC1       PC2
#1   ptg000007l     1  1000     500.5 18.21882  9.548866
#2   ptg000007l   501  1500    1000.5 14.38647  7.604225
#3   ptg000007l  1001  2000    1500.5 12.09975  7.025071


# wite table :
sw_chr_name = paste("PCs_SNP_Sliding" , chr, "table", sep = "_")
write.table(sw_snp_data, sw_chr_name)
# Plot PC1 variance along chr
#-----------------------------#
# read the table if needed
# sw_data=read.table(paste("Sliding_PCA/PCs_Sliding", , chr, "table" sep = "_"), header=TRUE)

# get current chr length
length_chr=table_chr[which(table_chr$scaffold==chr),"length"]

# plot 
p=ggplot() +
  geom_point(aes(x=sw_snp_data$mid_point, y=sw_snp_data$PC1),
             color = "darkred",
             alpha = 0.5,
             size=1)+
  
  labs(x="", y = "PC 1") +
  
  ggtitle(chr) +
  theme(text = element_text(size = 10), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text = element_text(size = 10),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5)) +
  
  scale_x_continuous(breaks = seq(0, length(getPOS(VCFR_data)), 
                                  by = round(length(getPOS(VCFR_data))/15))) +
  
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(1,100))

p

ggsave(p, path="sliding_PCA/", width = 12, height = 10, device = "pdf", filename = paste(chr, "PCA_scans.pdf",sep = "_"))
