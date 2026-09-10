#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH -t 150:00:00
#SBATCH --mem 5G


PATIENT=$1
WORKDIR=$2
SAMPLELIST=$3
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${SAMPLELIST})
echo $SAMPLE

# Create directory if it does not exist
GINKGO_INST_DIR=/gpfs/commons/groups/landau_lab/tprieto/apps/ginkgo/uploads/
mkdir -p ${GINKGO_INST_DIR}${PATIENT}

# Prepare input
module purge
module load SAMtools/1.21
module load BEDTools/2.31.0-GCC-12.3.0

samtools view -bq 40 ${WORKDIR}${SAMPLE}.bam | \
	bedtools bamtobed -i stdin | gzip > \
	${GINKGO_INST_DIR}${PATIENT}/${SAMPLE}.bed.gz
