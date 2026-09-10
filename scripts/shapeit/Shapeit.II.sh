#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu 1G

threads=1

# Load the input and create variables 
INPUT=$1
WORKDIR=$(dirname $INPUT)
PATIENT=$(basename $INPUT | sed 's/_.*//')
CHR=21
CHR=$2
LARGE_CHUNKS=/gpfs/commons/groups/landau_lab/ResolveOME/Resources/shapeit5/resources/chunks/b38/UKB_WGS_200k/large_chunks_25cM/
OUTDIR=${WORKDIR}/phasing/

IN=${OUTDIR}${PATIENT}_chr${CHR}.chunk_*.shapeit5_common.bcf
OUT=${OUTDIR}${PATIENT}_chr${CHR}.shapeit5_common_ligate.bcf
module purge
module load shapeit/5.1.0
ls -1v $IN > ${OUTDIR}list_ligate.chr${CHR}.txt && ligate_static --input ${OUTDIR}list_ligate.chr${CHR}.txt --output ${OUT} --thread $threads --index


