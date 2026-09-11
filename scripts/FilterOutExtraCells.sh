#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 4
#SBATCH -t 01:00:00
#SBATCH --mem-per-cpu 4G


patient=$1
classifier_dir=/gpfs/commons/groups/landau_lab/ResolveOME/StartDir/${patient}/manual_collate/sequoia/classifier/

module purge
module load bcftools/1.23.1

grep -v -F -x -f ${classifier_dir}contaminated_cells_manualinspection.txt \
	${classifier_dir}noncontaminated_cells.txt > ${classifier_dir}noncontaminated_cells2.txt 

# Extract only selected cells and remove variants with no alternative allele in a single step
bcftools view -S ${classifier_dir}noncontaminated_cells2.txt -e 'COUNT(VAF>0)<1' ${classifier_dir}/myselectedvars.snvswithsingletons.classifierfiltered.nocontaminatedcells.vcf.gz -Oz -o ${classifier_dir}/myselectedvars.snvswithsingletons.classifierfiltered.nocontaminatedcells2.vcf.gz
mv ${classifier_dir}/myselectedvars.snvswithsingletons.classifierfiltered.nocontaminatedcells2.vcf.gz ${classifier_dir}/myselectedvars.snvswithsingletons.classifierfiltered.nocontaminatedcells.vcf.gz
bcftools index -f ${classifier_dir}/myselectedvars.snvswithsingletons.classifierfiltered.nocontaminatedcells.vcf.gz
