#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 2G
#SBATCH -t 20:00:00

source ReadConfig.sh $1
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${ORIDIR}/${SAMPLELIST})
echo $SAMPLE
module purge
module load gatk/4.5.0.0-GCCcore-12.3.0-Java-17

suffix=".dedup"

gatk BaseRecalibrator \
	-R ${RESDIR}/${REF}.fasta \
	-I ${WORKDIR}/${SAMPLE}${suffix}.bam \
	--known-sites ${RESDIR}/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
	--known-sites ${RESDIR}Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
	-O ${WORKDIR}/RecalibrationReportI_${SAMPLE}.grp
