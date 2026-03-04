### Genomic Scans in Sliding Windows

**Contributors:** 

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


### TO DO  (deleted for now because too old)

`theta_Scan.r`
  - Fst by SW

- Global SFS by chromosome  
- SFS in sliding windows: Euclidean distance of SFS in sliding windows vs. global SFS on chromosome
- LD (to do)  
---

### Input

1. **Sample summary table**
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


6. **Scaffold size table**:
```
Chr length
SUPER_1 274107945
SUPER_2 239769335
SUPER_3 233268383
```
---

## Methods

Genomic scans can be run independantly in order to be adapted easly within HPC (one SLURM job by chr)

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

# iterate over chr name
chr_name=/MYPATH/chr_name.txt
sc_list < <(cat $chr_name)
idx=$((SLURM_ARRAY_TASK_ID-1))
curren_sc="${sc_list[$idx]//[[:space:]]/}"

# Load r module in your HPC
module load R/4.3.3
Rscript Fst_scans.r $curren_sc
```

  - Prepare the iterative variable in your R script  `Fst_scans.r`:

```r
args <- commandArgs(trailingOnly = TRUE)
chr <- as.character(args[1])
```

  - Run slurm script :
```sh
sbatch run_Scan.sh
```
