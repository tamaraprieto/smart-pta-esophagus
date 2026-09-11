#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tamara.prieto.fernandez@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 100:00:00
#SBATCH --mem 1G

module purge
module load samtools/1.9 
module load bcftools/1.9

patient=$1
tag=$2

indir=/gpfs/commons/groups/landau_lab/ResolveOME/StartDir/${patient}/manual_collate/sequoia/
outdir=${indir}fnr/
mkdir $outdir
vcfname=${patient}_${tag}.sequoia.vaf.eso
germlinevars=${indir}/${vcfname}.germline.newfilters.het.vcf.gz
samplelist=/gpfs/commons/groups/landau_lab/ResolveOME/StartDir/${patient}/dna/Samples.${patient}.txt


cell=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $samplelist)
echo $cell

# create a file with read counts for the reference and alternative alleles to load in R after merging
bcftools view --samples $cell $germlinevars | \
	 bcftools query  \
        -f  '%CHROM:%POS:%REF>%ALT\t%avsnp150\t[%AD{0}\t%AD{1}]\n' |\
	awk -v cell=$cell 'BEGIN{print cell"-ref\t"cell"-alt\t"cell"-sum\t"cell"-vaf"}{printf $3"\t"$4"\t"$4+$3"\t";if ($4+$3>0){print $4/($4+$3)}else{print "NA"}}'> \
	${outdir}${cell}.germline.het.counts_recall.txt 

awk -v name=$cell 'BEGIN{sum=0}{if ($4>0){sum+=1}}END{print sum"\t"NR"\t"sum/NR"\t"name}' ${outdir}${cell}.germline.het.counts_recall.txt > \
        ${outdir}${cell}.recall.txt 
