#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH -t 10:00:00
#SBATCH --mem-per-cpu 1G


patient=$1
classifier_all_dir=/gpfs/commons/groups/landau_lab/ResolveOME/StartDir/classifier/
classifier_dir=/gpfs/commons/groups/landau_lab/ResolveOME/StartDir/${patient}/manual_collate/sequoia/classifier/

module purge
module load bcftools/1.23.1

CELL=$(bcftools query -l ${classifier_dir}/myselectedvars.snvswithsingletons.classifierfiltered.vcf.gz | sed "${SLURM_ARRAY_TASK_ID}q;d")
echo $CELL
contamination_file=${classifier_dir}"contamination/Contamination_${CELL}.txt"
mkdir ${classifier_dir}contamination/

    proportion=$(bcftools view -s $CELL -c 1 ${classifier_dir}/myselectedvars.snvswithsingletons.classifierfiltered.vcf.gz | \
    awk -F'\t' -v cell=$CELL 'BEGIN {total=0; rs_count=0}
        /^#/ {next}
        {
            total++;
            if ($8 ~ /avsnp150=rs[0-9]/) rs_count++;
        }
        END {
            if (total > 0) print rs_count/total;
            else print "NA";
        }')

echo -e "$CELL\t$proportion" > $contamination_file
