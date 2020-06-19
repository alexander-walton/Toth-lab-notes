
#featurecounts2
```
#!/bin/bash

#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 24:00:00
#SBATCH -J Hisat2
#SBATCH -o Hisat2.o%j
#SBATCH -e Hisat2.e%j
#SBATCH --mail-user=awalton@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
set -o xtrace


cd /work/LAS/amytoth-lab/awalton/03_hitsat/
ulimit -s unlimited
ANNOT="/work/LAS/amytoth-lab/awalton/03_hitsat/fuscatus.gff"
mkdir -p counts
ODIR="counts"


module purge

module load subread
module load parallel

parallel -j 4 "featureCounts -B 100 --multiMapping -T 4 -s 2 -p -t gene -g ID -a $ANNOT -o $ODIR/{/.}.gene.txt {}" :::  hisatOut/*.bam
scontrol show job $SLURM_JOB_ID
```
convert into a single file, to be manipulated in R
```
paste <(awk 'BEGIN {OFS="\t"} {print $1,$7}' PFUS_1.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_2.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_3.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_4.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_5.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_6.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_7.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_8.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_9.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_10.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_11.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_12.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_13.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_14.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_15.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_16.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_17.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_18.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_19.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_20.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_21.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_22.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_23.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_24.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_25.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_26.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_27.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_28.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_29.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_30.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_31.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_32.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_33.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_34.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_35.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_36.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_37.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_38.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_39.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_40.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_41.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_42.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_43.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_44.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_45.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_46.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_47.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_48.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_49.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_50.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_51.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_52.gene.txt) | grep -v '^\#' > Pf_count.txt
```

#Then this file needs to be cleaned up, removing excess text from column names.

```
's/hisatOuty\///g' Pf_count.txt
perl -i -pe 's/\.bam//g' Pf_count.txt
```
#Also need to go into vim and remove the "geneid" from the first column, top cell.



```
paste <(awk 'BEGIN {OFS="\t"} {print $1,$7}' AW1724_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1754_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1761_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1841_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1887_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1933_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1726_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1756_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1816_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1853_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1906_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1949_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1728_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1758_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1820_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1863_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1909_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1742_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1759_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1830_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1872_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1919_1.fq.gz.gene.txt) | grep -v '^\#' > Foundress_count.txt

```
```
perl -i -pe 's/hisatOut\///g' Foundress_count.txt
perl -i -pe 's/\.bam//g' Foundress_count.txt
perl -i -pe 's/\_1//g' Foundress_count.txt
```

#6-3-20 working directory /Users/alexanderwalton/Box Sync/bionformatics training with Severin/combined data from me and Legan
#combining data from my foundresses and Andrew Legan's RNAseq data, below is just the files from Andrew that are heads (he also had abdomens, antennae, and thoraces).
```
paste <(awk 'BEGIN {OFS="\t"} {print $1,$7}' PFUS_1.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_3.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_4.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_5.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_6.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_7.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_8.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_10.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_12.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_31.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_33.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' PFUS_35.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1724_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1754_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1761_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1841_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1887_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1933_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1726_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1756_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1816_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1853_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1906_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1949_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1728_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1758_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1820_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1863_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1909_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1742_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1759_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1830_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1872_1.fq.gz.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' AW1919_1.fq.gz.gene.txt) | grep -v '^\#' > waltonlegan_heads_count.txt


perl -i -pe 's/hisatOuty\///g' waltonlegan_heads_count.txt
perl -i -pe 's/\.bam//g' waltonlegan_heads_count.txt
perl -i -pe 's/\_1.fq.gz//g' waltonlegan_heads_count.txt
