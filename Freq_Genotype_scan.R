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
# When runing loop over all chr vcf : 
# chr = arg[2]
# metadata = arg[3]
# length_table = arg[4]
# vcf = arg[1]

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

# get genotype table
loci_table = extract.gt(VCFR_data, element = "GT")
loci_table = as.data.frame(loci_table)
geno_table = Convert_GT(GT_table = loci_table)

# geno_table[c(1:10),c(1:10)]
#                10_Q_S_A_S9 12_Q_S_A_S11 13_Q_A_S12 14_Q_A_S13 15_Q_A_S14 16_Q_A_S15 17_Q_A_S16 18_Q_A_S17 19_Q_A_S18 1_Q_S_A_S1
#CHR1_159039           0            1          0          1          0          0          0          0          0          0
#CHR1_219063           0            0          0          0          0          0          0          0          0          0
#CHR1_254450           0            0          1          0          0          0          0          0          0          0
#CHR1_359471           0            1          0          1          0          1          0          1          1          2
#CHR1_376366           0            0          0          0          0          0          0          0          0          0
#CHR1_416191           1            2          1          1          1          2          2          2          1          2
#CHR1_464593           0            0          0          0          0          0          0          0          0          0
#CHR1_486240           0            0          1          0          0          0          0          0          0          0
#CHR1_493425           0            0          0          0          0          0          1          0          0          0
#CHR1_514070           0            0          0          0          1          0          0          0          0          0
 
table = GetGenoFreqByPop(geno_table = loci_table, chr_name = chr, ploidy = "diploid", pop_list = pop_list)

head(table)

# write final table
write.table(table, "Geno_Freq/table_geno_freq.txt", quote=F, row.names = F)

# read file from local folder
geno_file=read.table("Geno_Freq/table_geno_freq.txt")

# get size of the current chromosome
length_chr=table_chr[which(table_chr$Chr==chr),"length"]

# plot
p = ggplot() +
  geom_point(data=table,
             aes(x=table$position,y=table$monogyne_het_freq), color="red", pch=19,alpha=0.5) +
  geom_point(data=table,
             aes(x=table$position,y=table$polygyne_het_freq), color="blue",pch=19, alpha=0.5) +
  
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


#=====================================================#
#=====================================================#
# ------ Make average frequencies by individuals ----
# 
# IF you want to get get the average rate of each genotype 
# over one specific region (typically supergene) within each samples
#
#=====================================================#
#=====================================================#

# define your region
start=10000000
end=25000000

# format table :
position_vector=as.numeric(str_remove(rownames(geno_table), paste(chr, "_", sep ="")))
geno_table$position = position_vector
geno_table_region = geno_table[c(geno_table$position > start & geno_table$position < end),]

# get average genotype by individual
table_by_sample = apply(geno_table_region[,-ncol(geno_table_region)], 2, table) / nrow(geno_table_region)
table_by_sample_t = as.data.frame(t(table_by_sample))
table_by_sample_t$sample <- rownames(table_by_sample_t)

# Add pop name from pop_list 
table_by_sample_t$sociality <- ifelse(table_by_sample_t$sample %in% pop_list$monogyne,
                       "monogyne", "polygyne")

# format for plot
table_by_sample_t_melt = melt(table_by_sample_t)

freq = ggplot(data = table_by_sample_t_soc_melt) +
  geom_violin(aes(x = vec_mono_poly, y = value, colour = variable, fill = variable),
              alpha = 0.3) +
  ylab("Genotype frequencies") +
  xlab("groups") +
  
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
  
  scale_colour_manual(values=RColorBrewer::brewer.pal(3, "Set1"),
                     name = "genotype",
                     label = c("Homozygous to ref", "heterozygous", "homozygous alt")) +
  scale_fill_manual(values=RColorBrewer::brewer.pal(3, "Set1"), guide = "none") +
  theme(text = element_text(size=12),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5))
freq
ggsave(freq, path="Geno_Freq/", width = 10, height = 10, device = "pdf", filename = paste(chr, "Group_Genotype.pdf",sep = "_"))       
