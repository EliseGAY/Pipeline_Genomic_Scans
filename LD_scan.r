#!/usr/bin/env Rscript

#============================#
#============================#
# ------ Load libraries ----
#============================#
#============================#
library(IRanges)
library(vcfR)
library(pegas)
library(ggplot2)
library(adegenet)
library(hierfstat)
library(reshape2)
library(pcadapt)
library(ggrepel)
library(gridExtra)
library(stringr)
library(devtools)
library(GenomicRanges)
library(LDheatmap)
library(RColorBrewer)
library(snpStats)
library(tibble)
system("git clone https://github.com/EliseGAY/Package_VCF2PopStructure.git")
load_all("Package_VCF2PopStructure/")

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
metadata=read.table("data/metadata.txt", header = TRUE)
metadata=as.data.frame(metadata)
pop=unique(metadata$Population)
pop_list = split(metadata$sample, metadata$Population)

# read chr length
table_chr=read.table("Scaff_length_sorted.list", header = T)

#-----------------------------------------------------------#
# Generate Genotype tables needed in different R packages
#-----------------------------------------------------------#

# Read the VCF with vcfR :
VCFR_data=read.vcfR("data/CHR1_Flowqual_Noindels_Norepeat_DP10_50_Na0_SNP.vcf.gz", verbose = TRUE)

# current_chr
chr = arg[1]
# to test in local example 
chr = "chr1"
chr_len = table_chr[which(table_chr$scaffold == chr),]$length

#=========================#
#=========================#
# ------ LD heatmap ----
#=========================#
#=========================#
# bin data
data_bin <- bin_snps(vcf = VCFR_pop1, bin_size = 200)
data_bin

# run LDheatmap
snp_matrix <- vcfR2SnpMatrix(data_bin, phased = NULL, subjects = NULL)
# explore snp_matrix object 
snp_matrix$data # 00 are NA and the genotypes are 01, 02, 03
snp_matrix$genetic.distances # these are the snps positions

# get LD and plot
colors <- brewer.pal(9,"Reds")
colors <- colors[c(9,8,7,6,5,4,3,2,1)]

nompng=paste0("LDheatmap_2",chr,".png")
png(nompng, width = 11.69*1000, height = 8.27*1000, res = 1000)
LD <- LDheatmap(snp_matrix$data, snp_matrix$genetic.distances, flip = T, color = colors)
dev.off()

# explore LD object
dim(LD$LDmatrix)
LD$LDmatrix[c(1:10), c(1:10)]
LD$genetic.distances # these are the snps positions

# extract SNP with correlation superior to "threshold"
indices_table = as.data.frame(which(LD$LDmatrix >= 0.9, arr.ind=TRUE))
index_snp1 = indices_table$col
index_snp2 = indices_table$row
snp1 = LD$genetic.distances[index_snp1]
snp2 = LD$genetic.distances[index_snp2]

length(snp1) = max(length(snp1), length(snp2))
length(snp2) = max(length(snp1), length(snp2))
snp_table = data.frame("snp1" = as.numeric(snp1), "snp2" = as.numeric(snp2))

# write result
write.table(snp_table, paste0("Coord_",chr,".txt"), sep = " ", row.names = TRUE, col.names = TRUE)
