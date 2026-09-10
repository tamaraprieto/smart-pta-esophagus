#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 5
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu 20G


# variable to declare
patient=$1
tag=$2
# path and name of a bgzipped vcf
workdir=/gpfs/commons/groups/landau_lab/ResolveOME/StartDir/${patient}/manual_collate/sequoia/
vcfname=${patient}_${tag}.sequoia.vaf.eso
VCF=${workdir}/${vcfname}
resources=/gpfs/commons/groups/landau_lab/ResolveOME/Resources/
reference=references/hg38/Homo_sapiens_assembly38.fasta

module purge
module load bcftools/1.9

#####################
# Select a in avsnp #
#####################

echo "Generating germline set..."
bcftools filter -i 'avsnp150!="." && Germline_qval > 0.00005 &&  Mean_Depth > 3' \
	${VCF}.overlap.vcf.gz -O z > ${workdir}/${vcfname}.germline.newfilters.vcf.gz

##############################################
# Create a set of heterozygous germline SNPs #
##############################################

echo "Select heterozygous..."
# include only sites with one or more heterozygous sites (het) and VAF betwen 0.25 and 0.75
bcftools view -i 'GT=="het" && AF > 0.25 && AF < 0.75' ${workdir}/${vcfname}.germline.newfilters.vcf.gz -O z > \
       ${workdir}/${vcfname}.germline.newfilters.het.vcf.gz

