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

system("git clone https://github.com/EliseGAY/Package_VCF2PopStructure.git")
load_all("../Package_VCF2PopStructure/")

#===============================#
#===============================#
# ------ PREPARE YOUR DATA ----
#===============================#
#===============================#
# read arg 
chr = arg[2]
metadata = arg[3]
length_table = arg[4]
vcf = arg[1]

#===============================#
#===============================#
# ------ PREPARE YOUR DATA ----
#===============================#
#===============================#
getwd()
dir()
if (!dir.exists("PCA_scans")) {
  dir.create("PCA_scans")
}


#----------------------------#
# to test in the toy example
#----------------------------#

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
table_chr=read.table("metadata/table_chr_length.txt", header = T)

# set chr name
chr = "CHR1"

# Read the VCF with vcfR :
VCFR_data=read.vcfR(paste("data/",chr, "_Example.vcf.gz", sep =""))

# create a pop sorted by VCF colnames
metadata_sorted <- metadata[match(colnames(VCFR_data@gt)[-1], metadata$sample),]
pop_list = split(metadata_sorted$sample, metadata_sorted$Social_Morph)

# read chr langth table
chr_len = table_chr[which(table_chr$Chr == chr),]$length

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
# Create a table geno_table  with pos (col), pos (rows) and genotype as 0,1,2
# window : the size of the window (int)
# slide : the size of the slide (int)
# min_n_snp: minimum Nb of SNPs to compute PCA (100 should be enough, it is just to avoid error at the bounds of the scaffold)

# OUTPUT
#---------#
# The Get_sw_abs() return a table of sliding windows coordinates made from 1 to the last SNP position :
# Fine-tune the wind and slide 
#'  chr        start   end mid_point nb_snp
#'<chr>      <int> <int>     <dbl>  <int>
#'1 ptg000007l     1 10000     5000.     21
#'2 ptg000007l  2001 12000     7000.     60
#'3 ptg000007l  4001 14000     9000.     73
#'4 ptg000007l  6001 16000    11000.     87

# in "PCs" you have a matrix of four columns: 
# col1 = lower value of the window
# col2 = upper value of the window
# col3 = mid value of the window
# col3 = % of variance of the first axe
# col4 = % of variance of the second axe


# get genotype table
loci_table = extract.gt(VCFR_data, element = "GT")
loci_table = as.data.frame(loci_table)
geno_table = t(Convert_GT(GT_table = loci_table))

geno_table[c(1:10),c(1:10)]
#                   CHR1_159039 CHR1_219063 CHR1_254450 CHR1_359471 CHR1_376366 CHR1_416191 CHR1_464593 CHR1_486240
#10_Q_S_A_S9            0           0           0           0           0           1           0           0
#12_Q_S_A_S11           1           0           0           1           0           2           0           0
#13_Q_A_S12             0           0           1           0           0           1           0           1
#14_Q_A_S13             1           0           0           1           0           1           0           0
#15_Q_A_S14             0           0           0           0           0           1           0           0
#16_Q_A_S15             0           0           0           1           0           2           0           0
#17_Q_A_S16             0           0           0           0           0           2           0           0
#18_Q_A_S17             0           0           0           1           0           2           0           0
#19_Q_A_S18             0           0           0           1           0           1           0           0
#1_Q_S_A_S1             0           0           0           2           0           2           0           0

# get SW table :
sw_data = Get_SW_abs(VCFR_data, slide = 100000, window = 200000)

# Run PCA 
min_snp=10
# loop over the windows :
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
sw_chr_name = paste("PCA_scans/PCs_Sliding" , chr, "table", sep = "_")
write.table(sw_data, sw_chr_name)

# Plot PC1 variance along chr
#-----------------------------#
# read the table if needed
# sw_data=read.table(paste("Sliding_PCA/PCs_Sliding", , chr, "table" sep = "_"), header=TRUE)

# plot 
p=ggplot() +
  geom_point(aes(x=sw_data$mid_point, y=sw_data$PC1),
             alpha = 0.5,
             fill = "lightblue3",
             color = "lightcyan4",
             size = 2,
             pch=21) +
  
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
  
  scale_x_continuous(breaks = seq(0, chr_len, 
                                  by = round(chr_len/15))) +
  
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(1,100))

p

ggsave(p, path="PCA_scans/", width = 10, height = 5, device = "pdf", filename = paste(chr, "PCA_scans.pdf",sep = "_"))

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
# create a sliding windows table based on number of SNP using the Get_SW_SNPrange() function (and not chr length)
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
sw_snp_data = Get_SW_SNPrange(VCF = VCFR_data,nb_snp_wind = 20, nb_snp_slide = 10)

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
sw_chr_name = paste("PCA_scans/PCs_SNP_Sliding" , chr, "table", sep = "_")
write.table(sw_snp_data, sw_chr_name)

# Plot PC1 variance along chr
#-----------------------------#
# read the table if needed
# sw_data=read.table(paste("Sliding_PCA/PCs_Sliding", , chr, "table" sep = "_"), header=TRUE)

# plot 
p=ggplot() +
  geom_point(aes(x=sw_snp_data$mid_point, y=sw_snp_data$PC1),
             alpha = 0.5,
             fill = "lightblue3",
             color = "lightcyan4",
             size = 2,
             pch=21) +
  
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

ggsave(p, path="PCA_scans/", width = 10, height = 5, device = "pdf", filename = paste(chr, "PCA_BySNP_scans.pdf",sep = "_"))
