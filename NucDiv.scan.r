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
VCFR_data=read.vcfR("data/Chr1_Example.vcf.gz")

# create a fake vector of called and valid position (more or less same filter than the SNP.VCF used to compute the theta)
Pos_Tot = sort(round(runif(n = length(getPOS(VCFR_data)) * 100 , min = 1, max = max(getPOS(VCFR_data)))))

# create a pop sorted by VCF colnames
metadata_sorted <- metadata[match(colnames(VCFR_data@gt)[-1], metadata$sample),]
pop_list = split(metadata_sorted$sample, metadata_sorted$Social_Morph)

# current_chr
chr = arg[1]
# to test in local example 
chr = "CHR1"
chr_len = table_chr[which(table_chr$Chr == chr),]$length

#===============================================#
#===============================================#
# ------ FST on sliding windows on whole chr ----
#===============================================#
#===============================================#
# Create directory
getwd()
dir()
if (!dir.exists("Sliding_Div")) {
  dir.create("Sliding_Div")
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

# get SW table :
sw_data = Get_SW_abs(VCFR_data, slide = 200000, window = 500000)

sw_data[, "Diversity"] <- NA_real_
colcount = which((colnames(sw_data) == "Diversity"))
head(sw_data)
# chr   start     end mid_point nb_snp              Diversity
#1 CHR1       1  500000  250000.5      9                NA
#2 CHR1  200001  700000  450000.5     22                NA
#3 CHR1  400001  900000  650000.5     34                NA
#4 CHR1  600001 1100000  850000.5     45                NA
#5 CHR1  800001 1300000 1050000.5     52                NA
#6 CHR1 1000001 1500000 1250000.5     56                NA

min_snp = 20

# compute Theta for each windows
for(i in 1:nrow(sw_data)){
  if(sw_data[i,]$nb_snp >= min_snp){
    pos_start = sw_data[i, "start"]
    pos_end = sw_data[i, "end"]
    indexes_raw = which(getPOS(VCFR_data) >= pos_start & getPOS(VCFR_data) < pos_end)
    print(pos_start)
    # Allele ALT count with col = position and rows = samples
    sub_geno = t(geno_table[indexes_raw,])
    # tot pos in the windows :
    nb_pos_tot = length(which(Pos_Tot >= pos_start & Pos_Tot < pos_end))
    print("nb tot position")
    print(nb_pos_tot)
    # computing thetaPi with HierFstat. 
    Theta = theta.Watt.dosage(sub_geno,L=nb_pos_tot)
  
    # create pairs column
    sw_data[i, "Diversity"] = Theta
  }
}

# save
write.table(sw_data, paste("Sliding_Div/", chr, "_Theta_Pi_Scans.txt", sep = ""))

#----------#
# Plot Div
#----------#
# read the Div indices
sw_data=read.table(paste("Sliding_Div/", chr, "_Theta_Pi_Scans.txt"), header = TRUE)

# Format table : discard negative value
sw_data[sw_data < 0] = 0

# plot
plot_name=paste(chr, "_Div", sep="")
head(sw_data)

p = ggplot(data = sw_data) + 
  geom_point(aes_string(x="mid_point",y="Diversity"),
             alpha = 0.5,
             fill = "lightblue3",
             color = "lightcyan4",
             size = 2,
             pch=21) +

  labs(y="Div", x = chr)+
  ggtitle(plot_name) +
  
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
  
  scale_x_continuous(breaks = seq(0, chr_len, 
                                  by = round(chr_len/15))) +
  scale_y_continuous(breaks = seq(0, 0.1, 
                                  by = 0.01), limits = c(0,0.1))
p
ggsave(p, path="Sliding_Fst/", width = 10, height = 5, device = "pdf", filename = paste(chr, "Fst_scans.pdf",sep = "_"))
