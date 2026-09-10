#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 4
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu 1G

threads=4

# Load the input and create variables 
INPUT=$1
WORKDIR=$(dirname $INPUT)
PATIENT=$(basename $INPUT | sed 's/_.*//')
CHR=$2

OUTDIR=${WORKDIR}/phasing/
mkdir -p $OUTDIR


LARGE_CHUNKS=/gpfs/commons/groups/landau_lab/ResolveOME/Resources/shapeit5/resources/chunks/b38/UKB_WGS_200k/large_chunks_25cM/
MAP=/gpfs/commons/groups/landau_lab/ResolveOME/Resources/shapeit5/resources/maps/b38/chr${CHR}.b38.gmap.gz
CHUNKS=${LARGE_CHUNKS}/chunks_chr${CHR}.txt
# reference in the globus for 100G project: 
# Collection= 'EMBL-EBI Public Data'
# Path= '/1000g/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/'
REFERENCE="/gpfs/commons/groups/landau_lab/ResolveOME/Resources/shapeit5/resources/1kGP_high_coverage_Illumina/b38/1kGP_high_coverage_Illumina.chr${CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

module load shapeit/5.1.0
module load bcftools/1.23.1

while read LINE; do
    echo $LINE
    REG=$(echo $LINE | awk '{ print $3; }')
    CHUNK_NBR=$(echo $LINE | awk '{ print $1; }')
    OUT=${PATIENT}_chr${CHR}.chunk_${CHUNK_NBR}.shapeit5_common.bcf
    LOG=${PATIENT}_chr${CHR}.chunk_${CHUNK_NBR}.shapeit5_common.log
phase_common_static --input $INPUT \
  --map $MAP \
  --reference $REFERENCE \
  --output ${OUTDIR}${OUT} \
  --thread $threads \
  --log ${OUTDIR}${LOG} \
  --filter-maf 0.001 \
  --region $REG \
  && bcftools index -f ${OUTDIR}${OUT} --threads $threads
done < ${CHUNKS}
