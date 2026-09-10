#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 4
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu 20G


module purge
module load bcftools/1.15.1

main_dir=/gpfs/commons/groups/landau_lab/ResolveOME/StartDir/
out_dir=${main_dir}callset_intersection
file_a=${main_dir}eso29/manual_collate/sequoia/eso29_full.sequoia.vaf.eso.normalized.vcf.gz
file_b=${main_dir}cu01/manual_collate/sequoia/cu01_full.sequoia.vaf.eso.normalized.vcf.gz
file_c=${main_dir}cu02/manual_collate/sequoia/cu02_993cells.updated.sequoia.vaf.eso.normalized.vcf.gz

bcftools index $file_a
bcftools index $file_b
bcftools index $file_c


# generate a file with the variant coordinates which overlap
bcftools isec -p $out_dir -n +2 \
	$file_a $file_b $file_c 

overlap_info=${out_dir}/sites

echo "> Create the tsv with annotations and index it..."
awk 'BEGIN{print "#CHROM\tPOS\tEND\tREF\tALT\tOverlap"}{print $1"\t"$2"\t"$2+length($4)-1"\t"$3"\t"$4"\t"$5}' ${overlap_info}.txt > \
         ${overlap_info}.tab
bgzip -c ${overlap_info}.tab > \
        ${overlap_info}.tab.gz
tabix -s 1 -b 2 -e 3 -f \
        --comment '#' \
        ${overlap_info}.tab.gz
