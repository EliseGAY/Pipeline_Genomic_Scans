#!/usr/bin/env Rscript

#============================#
#============================#
# ------ Load libraries ----
#============================#
#============================#
library(vcfR)
library(ggplot2)
library(reshape2)
library(stringr)
library(devtools)
# system("git clone https://github.com/EliseGAY/Package_VCF2PopStructure.git")
load_all("../Package_VCF2PopStructure/")

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
table_chr=read.table("metadata/table_chr_length.txt", header = T)

#-----------------------------------------------------------#
# Generate Genotype tables needed in different R packages
#-----------------------------------------------------------#

# Read the VCF with vcfR :
VCFR_data=read.vcfR("../data/VCF_example.NoNa.SNP.vcf.gz")

# create a pop sorted by VCF colnames
metadata_sorted <- metadata[match(colnames(VCFR_data@gt)[-1], metadata$GT_sample),]
pop_list = split(metadata_sorted$GT_sample, metadata_sorted$Population)

# current_chr
chr = arg[1]
# to test in local 
chr = "ptg000007l"
chr_len = table_chr[which(table_chr$scaffold == chr),]$length

#===============================================#
#===============================================#
# ------ Nb of freq of chosen genotype ----
#===============================================#
#===============================================#
# Create directory
getwd()
dir()
if (!dir.exists("Geno_Freq")) {
  dir.create("Geno_Freq")
}

table = GetGenoFreqByPop(geno_table = loci_table, chr_name = chr, ploidy = "diploid", pop_list = pop_list)

# write final table
write.table(table, "Geno_Freq/table_geno_freq.txt", quote=F, row.names = F)

# read file from local folder
geno_file=read.table("Geno_Freq/table_geno_freq.txt")

# get size of the current chromosome
length_chr=table_chr[which(table_chr$scaffold==chr),"length"]

# plot
ggplot() +
  geom_point(data=table,
             aes(x=table$position,y=table$Cesseras_het_freq), color="red", pch=19,alpha=0.1) +
  geom_point(data=table,
             aes(x=table$position,y=table$Termes_het_freq), color="blue",pch=19, alpha=0.1) +
  
  ylab("Heterozygous rate by position in two pop") +
  
  xlab("pb") +
  
  ggtitle(paste(chr, "Freq het in both pop", sep = " ")) +
  
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
  scale_color_manual(values=c("red","blue"),
                     name = "pop", 
                     label=c(names(pop_list[1]), names(pop_list[2]))) +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1,
                                   size = 12),
        text = element_text(size=12),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5)) +
  scale_x_continuous(breaks = seq(0, length_chr, by = length_chr/10), 
                     limits = c(0, length_chr))
# save
ggsave(p, path="Geno_Freq/", width = 12, height = 10, device = "pdf", filename = paste(chr, "Freq_scans.pdf",sep = "_"))       
