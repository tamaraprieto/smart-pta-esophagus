#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 4
#SBATCH -t 20:00:00
#SBATCH --mem 7G

source ReadConfig.sh $1
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${ORIDIR}/${IDLIST})
echo $SAMPLE
module purge
module load bwa/0.7.18-GCCcore-12.3.0
module load samtools/1.21

ID=${SAMPLE}
SM=$(echo ${SAMPLE} | sed 's/-.*//')
PL=${PLATFORM}    #Allowable options are ILLUMINA,SLX,SOLEXA,SOLID,454,LS454,COMPLETE,PACBIO,IONTORRENT,CAPILLARY,HELICOS,UNKNOWN
LB=${LIBRARY}
PU=`zcat ${WORKDIR}${SAMPLE}.trimmed.R1.fastq.gz | head -1 | sed 's/[:].*//' | sed 's/@//' | awk '{print $1}'`

echo "SAMPLE: "${SAMPLE}" ID: "${ID}" SM: "${SM}
RG="@RG\\tID:${ID}\\tSM:${SM}\\tPL:${PL}\\tLB:${LB}\\tPU:${PU}"
echo $RG

bwa mem -t 4 \
	-R ${RG} \
	${RESDIR}/${REF}.fasta \
	${WORKDIR}/${SAMPLE}.trimmed.R1.fastq.gz \
	${WORKDIR}/${SAMPLE}.trimmed.R2.fastq.gz | samtools view -bS - > ${WORKDIR}/${SAMPLE}.bam
