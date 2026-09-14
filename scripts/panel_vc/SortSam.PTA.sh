#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH -t 100:00:00
#SBATCH --mem 6G

source ReadConfig.sh $1
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${ORIDIR}/${IDLIST})
echo $SAMPLE
module purge
module load picard/3.0.0-Java-17

java -jar $EBROOTPICARD/picard.jar \
	SortSam \
	I=${WORKDIR}/${SAMPLE}.bam \
	TMP_DIR=${WORKDIR} \
	O=${WORKDIR}/${SAMPLE}.sorted.bam \
	CREATE_INDEX=true \
	SORT_ORDER=coordinate
