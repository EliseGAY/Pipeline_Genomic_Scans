### Genomic Scans in Sliding Windows


---

### Aim

Provide example of genomic scan to run :

`PCA_scans.r`
  - PCA by SW
  - PCA by SNPs

`Fst_scans.r`

  - Fst by SW
  - Fst by SNPs

`Freq_Genotype_scan.R`

- Get overall freq of het and hom genotype
- Get distribution of genotypes (0,1,2) on specific region over individual or group
  
`LD_scan.r` :  (Contributors: Hugo Deshayes, PhD, EPHE)

- Use LDHeatmap to compute LD on binned VCF


`NucDiv_scan.r`

  - theta by SW

### TO DO  (deleted for now because too old)

- SFS in sliding windows: Euclidean distance of SFS in sliding windows vs. global SFS on chromosome
 
---

### Input

1. **metadata summary table**
```
pop pop
pop1  sample1
pop1  sample2
pop1  sample3
pop2  sample4
pop2  sample5
```

3. **VCF files**

4. **Position files** (two columns: scaffold and position — adapt to your file names):

5. **CHR name in bash list**

chr_name=("CHR1" "CHR2" "CHR3")

7. **Scaffold size table**:
```
Chr length
CHR1 274107945
CHR2 239769335
CHR3 233268383
```
---

## Methods

  - Genomic scans can be run independantly in order to be adapted easly within HPC (one SLURM job by chr)

  - Fill the arguments value according to the R script

```
vcf = arg[1] #VCF path of the corresponding chromosome
chr = arg[2] #chr name
metadata = arg[3] #path to metadata (col1 = population , col2 = samples ID)
length_table = arg[4] #table of chr length (see metadata)
```


`**example**` 

  - SLURM SCRIPT "run_Scan.sh"
    
```
#!/bin/sh
#SBATCH -p std
#SBATCH -n 10
#SBATCH --time=12:00:00
#SBATCH --job-name=Fst_SW_%a
#SBATCH --time=02:00:00
#SBATCH --array=1-20
#SBATCH -o scan_Fst_W_AnoNa_%a.o
#SBATCH -e scan_Fst_W_AnoNa_%a.e

# SET THE VARIABLE TO FEED TO RScript
##########################################

# iterate over chr name
idx=$((SLURM_ARRAY_TASK_ID-1))
curren_sc=chr_name[$idx]

# get VCf path
current_vcf=$("data/"${curren_sc}"_Example.vcf.gz")
# get chr length table
chr_lentgh="table_chr_length.txt"
# metadata
metadata="metadata/metadata.txt"
# Load r module in your HPC
module load R/4.3.3

# run the Rscript
################

Rscript Fst_scans.r $current_vcf $chr $curren_sc $metadata $length_table 

```

  - Run slurm script :
```sh
sbatch run_Scan.sh
```
