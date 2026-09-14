#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH -t 10:00:00
#SBATCH --mem-per-cpu 5G

module purge
module load samtools 
module load bedtools
module load bedops

# read arguments
patient=$1
bed=$1
newdir=$2
samplelist=$3
binSize=$4

# select sample
bam=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $samplelist)
cell="$(basename -- $bam | sed 's/.dedup//' | sed 's/.recal//' | sed 's/.bam//')"
echo "Analyzing alignment in $bam"


echo "Creating bam file containing just reads in genes..."
samtools view -b -L $bed $bam > ${newdir}${cell}.mygenes.bam
samtools view -b -F 1024 ${newdir}${cell}.mygenes.bam | samtools sort > ${newdir}${cell}.mygenes.nondup.bam
samtools index ${newdir}${cell}.mygenes.nondup.bam
bedtools coverage -a $bed -b ${newdir}${cell}.mygenes.nondup.bam -d > ${newdir}${cell}.mygenes.perbasecoverage.txt

# NON DUPL
bedtools genomecov -ibam ${newdir}${cell}.mygenes.nondup.bam -bg > ${newdir}${cell}.mygenes.nondup.bedgraph

echo "Creating coverage by windows..."
# transform bedgraph into file with intervals (same distance)
sort-bed ${newdir}${cell}.mygenes.nondup.bedgraph | awk -vOFS="\t" '{ print $1, $2, $3, ".", $4 }' - > ${newdir}${cell}.mygenes.${binSize}.nondup.bed
awk '{print $1"\t"0"\t"$2}' /gpfs/commons/groups/landau_lab/ResolveOME/Resources/references/hg38/Homo_sapiens_assembly38.fasta.fai | \
       head -n24 | sort-bed - > ${newdir}${cell}.hg38.bed
bedops --chop ${binSize} ${newdir}${cell}.hg38.bed | bedmap --echo --mean --skip-unmapped --indicator - ${newdir}${cell}.mygenes.${binSize}.nondup.bed | \
       sed 's/|/\t/g' | awk '{if ($5=="1"){print $0}}' | \
       bedtools intersect  -wao -a - -b $bed | \
       awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$9}' | \
       awk -v cell=$cell '{if ($5!="."){print $0"\t"cell}}' > \
	${newdir}${cell}.mygenes.${binSize}.hg38.nondup.bed
