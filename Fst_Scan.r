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
if (!dir.exists("Sliding_Fst")) {
  dir.create("Sliding_Fst")
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

# Prepare table to dynamically adapt to the number of pairwise comparisons
pop_pairs_table = getPairwisePop(pop_vec = sort(unique(metadata_sorted$Social_Morph)))
pop_pairs_table
#   V1       V2
# monogyne polygyne

# initialize your variable
pop_pairs_names = paste(pop_pairs_table$V1,pop_pairs_table$V2, sep = "_" )
# [1] "monogyne_polygyne"

sw_data[, pop_pairs_names] <- NA_real_
colcount = which((colnames(sw_data) %in% pop_pairs_names))
head(sw_data)
# chr   start     end mid_point nb_snp monogyne_polygyne
#1 CHR1       1  500000  250000.5      9                NA
#2 CHR1  200001  700000  450000.5     22                NA
#3 CHR1  400001  900000  650000.5     34                NA
#4 CHR1  600001 1100000  850000.5     45                NA
#5 CHR1  800001 1300000 1050000.5     52                NA
#6 CHR1 1000001 1500000 1250000.5     56                NA

# compute Fst for each windows
for(i in 1:nrow(sw_data)){
  if(sw_data[i,]$nb_snp >= min_snp){
    pos_start = sw_data[i, "start"]
    pos_end = sw_data[i, "end"]
    indexes_raw = which(getPOS(VCFR_data) >= pos_start & getPOS(VCFR_data) < pos_end)
    print(pos_start)
    sub_geno = geno_table[indexes_raw,]
    Fst_all = getGlobalFst_Count(loci_table_T = sub_geno, pop_table = metadata_sorted)
    # create pairs column
    sw_data[i, colcount] = Fst_all
  }
}
# fix colnames
colnames(sw_data) = c(colnames(sw_data)[1:5], names(Fst_all))

# save
write.table(sw_data, paste("Sliding_Fst/", chr, "_Hudson_Fst_scans.txt", sep = ""))

#----------#
# Plot FST
#----------#
# pair to plot
pairs_FST=colnames(sw_data)[6]

# read the FST indices
sw_data=read.table(paste("Sliding_Fst/", chr, "_Hudson_Fst_scans.txt"), header = TRUE)

# Format table : discard negative value
sw_data[sw_data < 0] = 0

# plot
plot_name=paste(chr, "_FST_",pairs_FST, sep="")
head(sw_data)

p = ggplot(data = sw_data) + 
  geom_point(aes_string(x="mid_point",y=pairs_FST),
             alpha = 0.5,
             fill = "lightblue3",
             color = "lightcyan4",
             size = 2,
             pch=21) +

  labs(y="FST", x = chr)+
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
  
  scale_x_continuous(breaks = seq(0, length_chr, 
                                  by = round(length_chr/15))) +
  scale_y_continuous(breaks = seq(0, 1, 
                                  by = 0.1), limits = c(0,1))
p
ggsave(p, path="Sliding_Fst/", width = 10, height = 5, device = "pdf", filename = paste(chr, "Fst_scans.pdf",sep = "_"))

#========================#
#========================#
# ------ FST by SNP ----
#========================#
#========================#

# Create directory
getwd()
dir()
if (!dir.exists("Sliding_Fst")) {
  dir.create("Sliding_Fst")
}

# 1. Run getFstBySNP_Count function 
Fst_by_SNP = getFstBySNP_Count(loci_table_T = geno_table, pop_table = metadata_sorted, Na_rate = 0.20, MAF_threshold = 0.05)

# 2. Create a global table with all position and the Fst pairs results :

# get all intersect of all position in the pairwise comparisons
position=sort(unique(unlist(Fst_by_SNP[which(grepl(pattern = "position", names(Fst_by_SNP)))])))

# initiate dataframe
all_fst = data.frame("position" = position)

# 3. loop over Fst pair:
for(name in names(Fst_by_SNP)){
  if(grepl(pattern = "position", name)){
    next}
  # get intersect for each pop
  pairs_pos = paste(name, "_position", sep = "")
  # get indexes in position present in current pair
  vec_pos_pair = as.numeric(unlist(Fst_by_SNP[pairs_pos]))
  # Add a column corresponding to the current Fst pair
  all_fst[which(position %in% vec_pos_pair),name] = Fst_by_SNP[name]
}

# check for restuls :
#position polygyne_monogyne
#1   159039       0.071784891
#2   359471      -0.014948461
#3   416191      -0.041384288
#4   493425      -0.073682407

head(Fst_by_SNP$polygyne_monogyne_position)
head(Fst_by_SNP$polygyne_monogyne)
# [1] 0.50375940 0.50375940 0.02957393 0.71929825 0.17293233 0.71929825
# [1]  9678 10382 10484 10982 11771 13963

# save the table 
write.table(all_fst, paste("Sliding_Fst/", chr, "_Hudson_Fst_bySNPs.txt", sep =""))

#-----------------------------#
# Plot FST
#-----------------------------#
# read chr length
length_chr=table_chr[which(table_chr$Chr==chr),"length"]

# Format table : discard negative value
all_fst[all_fst < 0] = 0
summary(all_fst)

# plot
all_fst_melt = melt(all_fst, id.vars = c("position"))
plot_name=paste(chr, "_FST_",pairs_FST, sep="")
head(all_fst_melt)
p=ggplot() + 
  geom_point(aes(x=all_fst_melt$position,y=all_fst_melt$value, 
                 group = all_fst_melt$variable, 
                 colour = all_fst_melt$variable,
                 fill = all_fst_melt$variable),
             alpha = 0.4,
             size = 1,
             pch=21) +
  
  labs(y = "FST", x = chr, colour = "Pop pairs", fill = "Pop pairs") +
  
  scale_colour_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  
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
  scale_y_continuous(breaks = seq(0, 1, 
                                  by = 0.1), limits = c(0,1))


ggsave(p, path="Sliding_Fst/", width = 10, height = 5, device = "pdf", filename = paste(chr, "Fst_BySNP_scans.pdf",sep = "_"))
