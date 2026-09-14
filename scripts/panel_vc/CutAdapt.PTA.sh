#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH -t 10:00:00
#SBATCH --mem 1G


source ReadConfig.sh $1
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${ORIDIR}/${IDLIST})
echo "SAMPLE: "$SAMPLE
module purge
module load cutadapt/4.4-GCCcore-12.2.0

number_rows=$(awk -v lib=${LIBRARY} '$1 == lib {print $0}' OFS='\t' ${SCRIPTDIR}/${ADAPTERS} | wc -l | awk '{print $1}')

if [ $number_rows -gt 1 ]; then
	echo "ERROR: Two lines with adapters for ${LIBRARY} were found in ${SCRIPTDIR}/${ADAPTERS}. Please, remove one"
	exit 1
elif [ $number_rows -eq 0 ];then
	echo "None adapter specified for library ${LIBRARY}"
	exit
else
	number_adapters=$(awk -v lib=${LIBRARY} '$1 == lib {print $0}' OFS='\t' ${SCRIPTDIR}/${ADAPTERS} | head -1 | awk -F'\t' '{print NF - 1; exit}')
	if [ $number_adapters -gt 2 ]; then
		echo "ERROR: Only possible to remove two adapters."
        	exit 1
	else 
		Adapter1=$(awk -v lib=${LIBRARY} '$1 == lib {print $2}' ${SCRIPTDIR}/${ADAPTERS})
		Adapter2=$(awk -v lib=${LIBRARY} '$1 == lib {print $3}' ${SCRIPTDIR}/${ADAPTERS})
		cutadapt --minimum-length=20 \
			-a AdapterA=$Adapter1 -A AdapterB=$Adapter2 \
			-o ${WORKDIR}/${SAMPLE}.trimmed.R1.fastq.gz \
			-p ${WORKDIR}/${SAMPLE}.trimmed.R2.fastq.gz \
			${ORIDIR}/${SAMPLE}_R1.fastq.gz ${ORIDIR}/${SAMPLE}_R2.fastq.gz > \
			${WORKDIR}/Cutadapt.${SAMPLE}.txt
	fi
fi
