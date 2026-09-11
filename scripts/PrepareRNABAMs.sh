#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 4
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu 1G

threads=4

patient=$1
rna_dir=/gpfs/commons/groups/landau_lab/ResolveOME/StartDir/all_cell_rna/

sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${rna_dir}/${patient}.rna.cells.txt)


module load all/Java/17.0.6
module load picard/3.0.0-Java-17
module load samtools
echo $sample
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
      I=${rna_dir}darkshore/output/star/${sample}.Aligned.sortedByCoord.out.bam \
      O=${rna_dir}/${sample}.bam \
      RGID=${sample} \
      RGLB=${sample} \
      RGPL=ILLUMINA \
      RGPU=unit1 \
      RGSM=${sample}
samtools index ${rna_dir}/${sample}.bam
