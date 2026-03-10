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
metadata=read.table("metadata/metadata.txt", header = TRUE)
metadata=as.data.frame(metadata)
pop=unique(metadata$Population)

# read chr length
table_chr=read.table("Scaff_length_sorted.list", header = T)

#-----------------------------------------------------------#
# Generate Genotype tables needed in different R packages
#-----------------------------------------------------------#

# Read the VCF with vcfR :
VCFR_data=read.vcfR("/data/VCF_example.NoNa.SNP.vcf.gz")

# create a pop sorted by VCF colnames
metadata_sorted <- metadata[match(colnames(VCFR_data@gt)[-1], metadata$GT_sample),]
pop_list = split(metadata_sorted$GT_sample, metadata_sorted$Population)

# current_chr
chr = arg[1]
# to test in local example 
chr = "ptg000007l"
chr_len = table_chr[which(table_chr$scaffold == chr),]$length

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
#'                ptg000007l_8628 ptg000007l_8631 ptg000007l_8633 ptg000007l_8639 ptg000007l_9618 ptg000007l_9622 ptg000007l_9625 ptg000007l_9626
#'1622W1_S254               0               0               0               0              NA              NA              NA              NA
#'1623W1_S276               0               0               0               0              NA              NA              NA              NA
#'1624W1_S392               0               0               0               0              NA               0               0               0
#'1625Q_S431                0               0               0               0              NA              NA              NA              NA
#'1625W1_S384               0               0               0               0               1               1               1               1

# get SW table :
sw_data = Get_SW_abs(VCFR_data, slide = 5000, window = 10000)
min_snp = 50

# Prepare table to dynamically adapt to the number of pairwise comparisons
pop_pairs_table= getPairwisePop(pop_vec = sort(unique(metadata_sorted$Population)))
pop_pairs_names = paste(pop_pairs_table$V1,pop_pairs_table$V2, sep = "_" )
sw_data[, pop_pairs_names] <- NA_real_
colcount = which((colnames(sw_data) %in% pop_pairs_names))

# compute Fst
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

# save
write.table(sw_data, "Sliding_Fst/Hudson_Fst_scans.txt")

#----------#
# Plot FST
#----------#
# pair to plot
pairs_FST="Cesseras_Cucugnan"


# read the FST indices
sw_data=read.table("SlidingFst/Hudson_Fst_scans.txt", header = TRUE)
# read chr length
length_chr=table_chr[which(table_chr$Chr==chr),"length"]
  
# Format table : discard negative value
sw_data[sw_data < 0] = 0

# plot
plot_name=paste(chr, "_FST_",pairs_FST, sep="")

p = ggplot(data = sw_data) + 
  geom_point(aes_string(x="mid_point",y=pairs_FST),
             alpha = 0.6,
             color = "darkred",
             pch=1) +

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

ggsave(p, path="Sliding_Fst/", width = 12, height = 10, device = "pdf", filename = paste(chr, "Fst_scans.pdf",sep = "_"))

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

# get all intersect of all position
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
#position Termes_Cesseras Termes_Cucugnan Cesseras_Cucugnan
#1       9678     0.503759398              NA       0.590062112
#2       9728              NA    -0.043478261       0.092627599
#3      10382     0.503759398              NA       0.590062112
#4      10434              NA    -0.043478261       0.092627599
#5      10484     0.029573935     0.306666667                NA
#6      10982     0.719298246     0.733333333                NA

# head(Fst_by_SNP$Termes_Cesseras)
# head(Fst_by_SNP$Termes_Cesseras_position)
# [1] 0.50375940 0.50375940 0.02957393 0.71929825 0.17293233 0.71929825
# [1]  9678 10382 10484 10982 11771 13963

# head(Fst_by_SNP$Termes_Cucugnan)
# [1] -0.04347826 -0.04347826  0.30666667  0.73333333 -0.09565217 -0.04347826
# head(Fst_by_SNP$Termes_Cucugnan_position)
# [1]  9728 10434 10484 10982 11761 11931

# head(Fst_by_SNP$Cesseras_Cucugnan)
# [1] 0.59006211 0.09262760 0.59006211 0.09262760 0.04725898 0.31677019
# head(Fst_by_SNP$Cesseras_Cucugnan_position)
# [1]  9678  9728 10382 10434 11761 11771

# save the table 
write.table(all_fst, "Sliding_Fst/Hudson_Fst_bySNPs.txt")

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

ggplot() + 
  geom_point(aes(x=all_fst_melt$position,y=all_fst_melt$value, 
                 group = all_fst_melt$variable, 
                 colour = all_fst_melt$variable),
             alpha = 0.9,
             pch=1) +
  
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
  
  scale_x_continuous(breaks = seq(0, chr_len, 
                                  by = round(chr_len/15))) +
  scale_y_continuous(breaks = seq(0, 1, 
                                  by = 0.1), limits = c(0,1))

ggsave(p, path="Sliding_Fst/", width = 12, height = 10, device = "pdf", filename = paste(chr, "Fst_BySNP_scans.pdf",sep = "_"))
