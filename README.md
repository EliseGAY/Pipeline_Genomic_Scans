### Genomic Scans in Sliding Windows

**Contributors:** 

---

### Aim

Provide example of genomic scan to run :

`PCA_scans.r`
  - PCA by SW
  - PCA by SNPs

### TO DO  (deleted for now because too old)
`Fst_scans.r`

  - Fst by SW
  - Fst by SNPs
    
`theta_Scan.r`
  - Fst by SW

`Fixed_Genotype_scans.r`

  - Get fixed genotype in x pops vs y pops

`Freq_genotype_scans.r`

- Get overall freq of het and hom genotype

- Global SFS by chromosome  
- SFS in sliding windows: Euclidean distance of SFS in sliding windows vs. global SFS on chromosome
- LD (to do)  
---

### Input

1. **Sample summary table**
```
samples pop
sample_1 pop1
sample_2 pop1
sample_3 pop1
sample_4 pop2
sample_5 pop2
```

3. **VCF files**

4. **Position files** (two columns: scaffold and position — adapt to your file names):


6. **Scaffold size table**:

Chr length
SUPER_1 274107945
SUPER_2 239769335
SUPER_3 233268383

---

## Methods

Genomic scans can be run independantly in order to be adapted easly within HPC (one job by chr)

`**example**` 

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

chr_name=/MYPATH/chr_name.txt
sc_list < <(cat $chr_name)
idx=$((SLURM_ARRAY_TASK_ID-1))
curren_sc="${sc_list[$idx]//[[:space:]]/}"

module load R/4.3.3
Rscript Fst_scans.r $curren_sc
```
- And then in `Fst_scans.r`:

```r
args <- commandArgs(trailingOnly = TRUE)
chr <- as.character(args[1])
```
