### 3-5-20 trying HiSat2 for mapping
#working directory


#formatting .gff files
```
#grabbed header
head -n 1 fuscatusHIC.gff> fuscatus.gff

#double checked that the number of scaffolds in the new .gff file matched the old one.
 grep Scaffold fuscatusHIC.gff | awk '{print $1}' | sort -k 1,1n | uniq| grep -v ">Scaffold" | wc
    168     168    1908
 grep Scaffold fuscatus.gff | awk '{print $1}' | sort -k 1,1n | uniq| grep -v ">Scaffold" | wc
    168     168    1908

### creating a new genome file
more GCA_010416935.1_CU_Pfus_HIC_genomic.fna | perl -pe 's/\>.*Scaffold/\>Scaffold/g' | perl -pe 's/, whole genome shotgun sequence//g' > fuscatus.fasta
```


hisat2 database build
file name = hisatmap
```
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -t 24:00:00
#SBATCH -J HI_build
#SBATCH -o HI_build.o%j
#SBATCH -e HI_build.e%j
#SBATCH --mail-user=awalton@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
set -o xtrace
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited

module load hisat2

GENOME="/work/LAS/amytoth-lab/awalton/03_hitsat/fuscatus.fasta"
hisat2-build $GENOME ${GENOME%.*}
scontrol show job $SLURM_JOB_ID
```



HiSat2 mapping script
hisatmap

```
#!/bin/bash
set -o xtrace
# set the reference index:
GENOME="/work/LAS/amytoth-lab/awalton/03_hitsat/fuscatus"
# make an output directory to store the output aligned files
mkdir -p hisatOuty
# set that as the output directory
ODIR="hisatOuty"


p=8 # use 8 threads
R1_FQ="$1" # first argument
R2_FQ="$2" # second argument

# purge and load relevant modules.
module purge

module load hisat2
module load samtools

OUTPUT=$(basename ${R1_FQ} |cut -f 1,2 -d "_");

hisat2 \
  -p ${p} \
  -x ${GENOME} \
  -1 ${R1_FQ} \
  -2 ${R2_FQ} \
  -S $ODIR\/${OUTPUT}.sam &> ${OUTPUT}.log
samtools view --threads 8 -bS -o $ODIR\/${OUTPUT}.bam $ODIR\/${OUTPUT}.sam

rm $ODIR\/${OUTPUT}.sam
```


loop_hisat2.sh
script that will run through all fasq files and map for every paired read of every sample

```
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -t 24:00:00
#SBATCH -J Hisat2
#SBATCH -o Hisat2.o%j
#SBATCH -e Hisat2.e%j
#SBATCH --mail-user=awalton@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
set -o xtrace
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
scontrol show job $SLURM_JOB_ID

for fq1 in *1.*gz;
do
fq2=$(echo $fq1 | sed 's/1/2/g');
/work/LAS/amytoth-lab/awalton/03_hitsat/hisatmap ${fq1} ${fq2};
done >& hisat2_1.log
```
to submit the job: sbatch loop_hisat2

####I also needed to make both the hisatmap and the loop_hisat2 files executeable with the code
```chmod 755 hisatmap```
#and
```chmod 755 loop_hisat2```
