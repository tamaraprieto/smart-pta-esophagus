#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu 100G

threads=4

# Load the input and create variables 
CHR=$1
patient=$2
rna_dir=/gpfs/commons/groups/landau_lab/ResolveOME/StartDir/all_cell_rna/
bams_rna=$(ls ${rna_dir}/${patient}*bam | \
  awk '{print "-I "$0}' | tr -s "\n" " ")

module load GATK/4.5.0.0-GCCcore-12.3.0-Java-17
module load BCFtools/1.19-GCC-13.2.0
vcf=/gpfs/commons/groups/landau_lab/ResolveOME/StartDir/${patient}/manual_collate/sequoia/myselectedvars.snvswithsingletons.LR.vcf.gz
vcf_region=/gpfs/commons/groups/landau_lab/ResolveOME/StartDir/${patient}/manual_collate/sequoia/myselectedvars.snvswithsingletons.LR.${CHR}.vcf.gz
bcftools view -r $CHR $vcf | awk 'BEGIN{prev1=""; prev2=""}{if ($1 != prev1 || $2 != prev2) print; prev1=$1; prev2=$2}' | bgzip  > $vcf_region
tabix -p vcf $vcf_region

gatk ASEReadCounter \
  -R /gpfs/commons/groups/landau_lab/ResolveOME/Resources/references/hg38/Homo_sapiens_assembly38.fasta \
  ${bams_rna} \
  -L $CHR \
  --disable-read-filter HasReadGroupReadFilter \
  -V $vcf_region \
  -O ${rna_dir}/allelic_counts_atSNVs/${patient}.${CHR}.alleliccounts_RNA.tsv

rm $vcf_region
rm ${vcf_region}.tbi
