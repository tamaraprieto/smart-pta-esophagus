#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type END
#SBATCH --cpus-per-task 4
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu 80G


##################################
# Define arguments and set paths #
##################################

patient=$1
tag=$2
resources=/gpfs/commons/groups/landau_lab/ResolveOME/Resources/
reference=references/hg38/Homo_sapiens_assembly38.fasta
vcf_input_dir=/gpfs/commons/groups/landau_lab/ResolveOME/StartDir/${patient}/manual_collate/sequoia/
vcf_input_name=${patient}_${tag}.sequoia.vaf.eso
myvcf=${vcf_input_dir}${vcf_input_name}.vcf.gz
samples_to_keep=/gpfs/commons/groups/landau_lab/ResolveOME/StartDir/cellswithgoodqual.txt
# define output
out_dir=$vcf_input_dir

###############################################################
# Normalize the vcf and remove cells & sites with low quality #
###############################################################

module purge
module load bcftools
bcftools view --samples-file $samples_to_keep --force-samples --threads $SLURM_CPUS_PER_TASK $myvcf | \
	bcftools view  --threads $SLURM_CPUS_PER_TASK -e 'QUAL<20' |\
	bcftools norm  --threads $SLURM_CPUS_PER_TASK -m -any --fasta-ref ${resources}/${reference} \
	-O z --output ${out_dir}/${vcf_input_name}.normalized.vcf.gz

echo "Indexing file"
bcftools index ${out_dir}/${vcf_input_name}.normalized.vcf.gz
