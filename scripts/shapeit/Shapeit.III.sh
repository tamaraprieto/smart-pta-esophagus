#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu 80G

threads=1

# Load the input and create variables 
INPUT=$1
WORKDIR=$(dirname $INPUT)
PATIENT=$(basename $INPUT | sed 's/_.*//')
OUTDIR=${WORKDIR}/phasing/

my_chromosome=chr${2}
IN=${OUTDIR}${PATIENT}_${my_chromosome}.chunk_*.shapeit5_common.bcf
OUT=${OUTDIR}${PATIENT}_${my_chromosome}.shapeit5_common_ligate.bcf


module purge
module load bcftools/1.9

echo '##FORMAT=<ID=PhasedGT,Number=1,Type=String,Description="Phased genotypes from ShapeIT">' > ${OUTDIR}annots.${my_chromosome}.hdr
# Create the files with genotypes to annotated the original VCF
myvcf_samples=$(zgrep -m 1 "^#CHR" $INPUT | cut -f10- | tr -s "\t" ",")
num_cells=$(echo $myvcf_samples | awk -F ',' '{print NF}')
bcftools view ${OUTDIR}/${PATIENT}_${my_chromosome}.shapeit5_common_ligate.bcf | \
  grep -v "^#" | \
  awk -v numcells=$num_cells 'BEGIN{printf "#CHROM\tSTART\tEND\tREF\tALT\t"; for(i=1; i <= numcells; i++) {printf "PhasedGT\t"}; print ""}{printf $1"\t"$2-1"\t"$2"\t"$4"\t"$5"\t"; for(i=10; i<=NF; i++){printf $i"\t"}; print ""}' | \
  awk '{sub(/\t$/, ""); print}'  > ${OUTDIR}phasing_annotation.${my_chromosome}.tsv
bgzip -c ${OUTDIR}phasing_annotation.${my_chromosome}.tsv > ${OUTDIR}phasing_annotation.${my_chromosome}.tsv.gz
tabix -c 1 -b 2 -e 3 -0 -f --comment '#' ${OUTDIR}phasing_annotation.${my_chromosome}.tsv.gz  
bcftools view -r $my_chromosome $INPUT | \
bcftools annotate \
 -a ${OUTDIR}phasing_annotation.${my_chromosome}.tsv.gz \
 -h ${OUTDIR}annots.${my_chromosome}.hdr \
 -c CHROM,FROM,TO,-,-,FMT/PhasedGT \
 -s $myvcf_samples \
 --output-type z \
 --output ${OUTDIR}/${PATIENT}_${my_chromosome}.shapeit5.${my_chromosome}.vcf.gz

echo "Done!"


# remove intermediate files
rm ${OUTDIR}annots.${my_chromosome}.hdr ${OUTDIR}phasing_annotation.${my_chromosome}.tsv ${OUTDIR}phasing_annotation.${my_chromosome}.tsv.gz


##############################################
# CREATE A TABLE FROM THE VCF TO LOAD INTO R #
##############################################

module purge
module load gatk/4.6.1.0-GCCcore-13.3.0-Java-17
module load bcftools/1.9
VCF=${OUTDIR}/${PATIENT}_${my_chromosome}.shapeit5.${my_chromosome}.vcf.gz
gatk IndexFeatureFile \
     -I $VCF
OUT=$(echo $VCF | sed 's/.vcf//' | sed 's/.gz//')
info_fields=$(bcftools view --header-only $VCF | grep "^##INFO=<" | sed 's/##INFO=<ID=//' | sed 's/,.*//' | grep -v "^$" | awk '{print "-F "$0}' | tr -t "\n" " " | sed 's/ -F $//' )
genotype_fields=$(bcftools view --header-only $VCF | grep "^##FORMAT=<" | sed 's/##FORMAT=<ID=//' | sed 's/,.*//' | \
        tr -t "\n" " " | sed 's/ / -GF /g' | awk 'BEGIN{printf "-GF "}{print $0}' | sed 's/ -GF $//')

echo ""
echo ${info_fields}
echo ${genotype_fields}

genotype_fields="-GF AD -GF VAF -GF PhasedGT -GF DP"
gatk VariantsToTable \
        --add-output-vcf-command-line \
        --show-filtered \
        -V ${OUT}.vcf.gz \
        -F CHROM -F POS -F REF -F ALT -F FILTER -F QUAL -F INFO \
        -F EVENTLENGTH \
        -F MULTI-ALLELIC \
        -F TRANSITION \
        -F HOM-REF \
        -F HET \
        -F HOM-VAR \
        -F VAR \
        -F NO-CALL \
        -F NCALLED \
        -F NSAMPLES \
        -F TYPE \
        ${info_fields} ${genotype_fields} \
        -O ${OUT}.table \
        -LE

