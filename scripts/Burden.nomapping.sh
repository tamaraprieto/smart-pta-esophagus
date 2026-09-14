#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 4
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu 50G

module purge
module load GATK/4.5.0.0-GCCcore-12.3.0-Java-17
module load BCFtools/1.19-GCC-13.2.0 


VCF=$1
echo $VCF
gatk IndexFeatureFile \
     -I $VCF


OUT=$(echo $VCF | sed 's/.vcf/.4round3burden/' | sed 's/.gz//')

info_fields=$(bcftools view --header-only   $VCF | grep "^##INFO=<" | sed 's/##INFO=<ID=//' | sed 's/,.*//' | grep -v "^$" | awk '{print "-F "$0}' | tr -t "\n" " " | sed 's/ -F $//' )
genotype_fields=$(bcftools view --header-only  $VCF | grep "^##FORMAT=<" | sed 's/##FORMAT=<ID=//' | sed 's/,.*//' | \
        tr -t "\n" " " | sed 's/ / -GF /g' | awk 'BEGIN{printf "-GF "}{print $0}' | sed 's/ -GF $//')

echo ""
echo ${info_fields}
echo ${genotype_fields}

genotype_fields="-GF VAF"
gatk VariantsToTable \
        --add-output-vcf-command-line \
        --show-filtered \
        -V ${VCF} \
        -F CHROM -F POS -F REF -F ALT -F FILTER -F QUAL -F INFO \
        ${genotype_fields} \
        -O ${OUT}.table \
        -LE

