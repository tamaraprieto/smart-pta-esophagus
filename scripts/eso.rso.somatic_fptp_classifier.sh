#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH -t 100:00:00
#SBATCH --mem 10G

module purge
module load bedtools/2.31.0-GCC-12.3.0
module load samtools/1.9 
module load bcftools/1.9

patient=$1
tag=$2
indir=/gpfs/commons/groups/landau_lab/ResolveOME/StartDir/${patient}/manual_collate/sequoia/
outdir=${indir}/fpr/
mkdir $outdir
bamdir=/gpfs/commons/groups/landau_lab/ResolveOME/StartDir/${patient}/dna/merged/
prefix=""
reference=/gpfs/commons/groups/landau_lab/ResolveOME/Resources/references/hg38/Homo_sapiens_assembly38.fasta
vcfname=${patient}_${tag}.sequoia.vaf.eso
somaticvars=${indir}/myselectedvars.snvswithsingletons.vcf.gz 
chrom_sizes="/gpfs/commons/groups/landau_lab/ResolveOME/Resources/references/hg38/Homo_sapiens_assembly38.fasta.fai"
germlinevars=${indir}/${patient}_${tag}.sequoia.vaf.eso.germline.newfilters.het.vcf.gz
samplelist=${bamdir}../Samples.${patient}.txt
cell=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $samplelist)
echo $cell
echo "bamfile:"${bamdir}${cell}${prefix}.bam

  # create a bed file with somatic  per cell
  bcftools view --samples ${cell} $somaticvars |\
	bcftools view -i 'VAF>0' | \
	grep -v "^#" | awk -v cell=$cell '{print $1"\t"$2-1"\t"$2"\t"$1":"$2"\t"cell"\t"$4"\t"$5}' \
  	> ${outdir}${cell}.bed

  # create a germline file with het germline sites present in the cell
  bcftools view --samples ${cell} $germlinevars | \
	bcftools view -i 'GT=="het"' -O z > ${outdir}${cell}.germline.het.vcf.gz  

  # pad the bed file intervals  
  bedtools slop \
   -i ${outdir}${cell}.bed \
   -g $chrom_sizes \
   -header \
   -b 100 \
   > ${outdir}${cell}.expanded.bed

  # intersect somatic expanded bed with hetSNPs 
  bedtools intersect \
            -wa \
            -wb \
            -a ${outdir}${cell}.expanded.bed \
            -b ${outdir}${cell}.germline.het.vcf.gz \
            -header \
            > ${outdir}${cell}.expanded_with_hetSNPs.bed

  # pad the bed file intervals even more for downsampling bam files
  bedtools slop \
   -i ${outdir}${cell}.expanded_with_hetSNPs.bed \
   -g $chrom_sizes \
   -header \
   -b 2000 \
   > ${outdir}${cell}.extra-expanded_with_hetSNPs.bed


   myregions=$(awk '{print $1":"$2"-"$3}' ${outdir}${cell}.extra-expanded_with_hetSNPs.bed | sort | uniq -c | awk '{print $2}' | tr -s "\n" " ")  
   
   # reduce the size of the bam file   
   samtools view -h ${bamdir}${cell}${prefix}.bam $myregions | samtools sort -  > ${outdir}${cell}.reduced.bam   
   samtools index ${outdir}${cell}.reduced.bam 

 
# count linked reads
echo "Analyzing: "${outdir}${cell}.reduced.bam
module purge
  # first arg: the bam file
  # second arg: the bed file
module load python/3.9.6
python fpr/CountReads.Intersect.py \
    --bamfile ${outdir}${cell}.reduced.bam \
    --bedfile ${outdir}${cell}.expanded_with_hetSNPs.bed \
    > ${outdir}${cell}.linkedreadcounts.txt


# classify somatic into FP/TPs based on linked read counts
python fpr/ClassifyVariants.py \
    --countsfile ${outdir}${cell}.linkedreadcounts.txt \
    --bedfile ${outdir}${cell}.expanded_with_hetSNPs.bed \
    > ${outdir}${cell}.linkedreadcounts.classified.txt
