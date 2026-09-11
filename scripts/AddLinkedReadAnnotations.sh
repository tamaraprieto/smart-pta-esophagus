#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 4
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu 50G


VCF_PATH=$1
VCF=$(echo $VCF_PATH | sed 's/.vcf//' | sed 's/.gz//')
workdir=$(dirname $VCF)
outdir=${workdir}/fpr

echo "VCF: "$VCF
echo "Outputing to: "$outdir

#####################################################################################
# Annotate the vcf file with the linked read annotations obtained within an R chunk #
#####################################################################################

module purge
module load BCFtools/1.19-GCC-13.2.0
module load vcftools
module load bedtools


echo "> Create the tsv tab with annotations"
awk -F'\t' 'BEGIN{OFS = FS;printf "a"}{print $0}'  ${outdir}/linked-reads-annotations.bed  | sort --version-sort  | \
	sed 's/^achrom/#chrom/' > ${outdir}/linked-reads-annotations.tab
bgzip -c ${outdir}/linked-reads-annotations.tab > ${outdir}/linked-reads-annotations.tab.gz
rm ${outdir}/linked-reads-annotations.tab.gz.tbi
tabix -c 1 -f -b 2 -0 -e 3 \
        --comment '#' \
        ${outdir}/linked-reads-annotations.tab.gz

grep "^chrom" ${outdir}/linked-reads-annotations.bed | sed 's/#//' | sed 's/\t/\n/g' > ${outdir}/newfilters.txt
# 1-:chrom   
# 2-:start
# 3-:end
# 4-a:LR_closesthetSNP
# 5-b:LR_class        
# 6-c:LR_class_stringent      
# 7-d:LR_shanon_div_median    
# 8-e:LR_matchingreads_sum
# 9-f:LR_nonmatchingreads_sum
# 10-g:LR_numcellstested  
# 11-h:LR_cellstested   
# 12-i:LR_ratiotp
# 13-j:LR_ratiofp
a=$(head -n 4 ${outdir}/newfilters.txt | tail -n 1)
b=$(head -n 5 ${outdir}/newfilters.txt | tail -n 1)
c=$(head -n 6 ${outdir}/newfilters.txt | tail -n 1)
d=$(head -n 7 ${outdir}/newfilters.txt | tail -n 1)
e=$(head -n 8 ${outdir}/newfilters.txt | tail -n 1)
f=$(head -n 9 ${outdir}/newfilters.txt | tail -n 1)
g=$(head -n 10 ${outdir}/newfilters.txt | tail -n 1)
h=$(head -n 11 ${outdir}/newfilters.txt | tail -n 1)
i=$(head -n 12 ${outdir}/newfilters.txt | tail -n 1)
j=$(head -n 13 ${outdir}/newfilters.txt | tail -n 1)
cat ${outdir}/newfilters.txt
rm ${outdir}/newfilters.txt

bcftools annotate \
       --annotations ${outdir}/linked-reads-annotations.tab.gz \
       -h <(echo '##INFO=<ID='${a}',Number=1,Type=Integer,Description="Position of the heterozygous SNP used for the linked-read analysis">\n') \
       -c CHROM,FROM,TO,${a},-,-,-,-,-,-,-,-,- \
       --output  ${VCF}.v1.vcf \
       ${VCF}.vcf.gz
bcftools annotate \
       --annotations ${outdir}/linked-reads-annotations.tab.gz \
       -h <(echo '##INFO=<ID='${b}',Number=1,Type=String,Description="Classification of the site after the linked-read approach">\n') \
       -c CHROM,FROM,TO,-,${b},-,-,-,-,-,-,-,- \
       --output  ${VCF}.v2.vcf \
       ${VCF}.v1.vcf
bcftools annotate \
       --annotations ${outdir}/linked-reads-annotations.tab.gz \
       -h <(echo '##INFO=<ID='${c}',Number=1,Type=String,Description="Stringent classification of the site after the linked-read approach">\n') \
       -c CHROM,FROM,TO,-,-,${c},-,-,-,-,-,-,- \
       --output  ${VCF}.v3.vcf \
       ${VCF}.v2.vcf
bcftools annotate \
       --annotations ${outdir}/linked-reads-annotations.tab.gz \
       -h <(echo '##INFO=<ID='${d}',Number=1,Type=Float,Description="Shanon entropy of the linked-read counts for the four potential single-cell haplotypes (median across cells when detected in more than 1 cell). Higher value means more haplotypes have been detected which is unexpected after the infinite sites assumption">\n') \
       -c CHROM,FROM,TO,-,-,-,${d},-,-,-,-,-,- \
       --output  ${VCF}.v4.vcf \
       ${VCF}.v3.vcf
bcftools annotate \
       --annotations ${outdir}/linked-reads-annotations.tab.gz \
       -h <(echo '##INFO=<ID='${e}',Number=1,Type=Integer,Description="Number of linked-read counts which support any of the 4 single-cell haplotypes added across cells">\n') \
       -c CHROM,FROM,TO,-,-,-,-,${e},-,-,-,-,- \
       --output  ${VCF}.v5.vcf \
       ${VCF}.v4.vcf 
bcftools annotate \
       --annotations ${outdir}/linked-reads-annotations.tab.gz \
       -h <(echo '##INFO=<ID='${f}',Number=1,Type=Integer,Description="Number of linked-read counts which support an haplotype different from the four determined by the reference and alternative alleles at the somatic and heterozygous germline site interrogated for the linked-read analysis">\n') \
       -c CHROM,FROM,TO,-,-,-,-,-,${f},-,-,-,- \
       --output  ${VCF}.v6.vcf \
       ${VCF}.v5.vcf  
bcftools annotate \
       --annotations ${outdir}/linked-reads-annotations.tab.gz \
       -h <(echo '##INFO=<ID='${g}',Number=1,Type=Integer,Description="Number of cells in which the linked-read approach has been tested (with enough depth)">\n') \
       -c CHROM,FROM,TO,-,-,-,-,-,-,${g},-,-,- \
       --output  ${VCF}.v7.vcf \
       ${VCF}.v6.vcf
bcftools annotate \
       --annotations ${outdir}/linked-reads-annotations.tab.gz \
       -h <(echo '##INFO=<ID='${h}',Number=1,Type=String,Description="Name of the cells in which the linked-read approach has been tested">\n') \
       -c CHROM,FROM,TO,-,-,-,-,-,-,-,${h},-,- \
       --output ${VCF}.v8.vcf \
       ${VCF}.v7.vcf
bcftools annotate \
       --annotations ${outdir}/linked-reads-annotations.tab.gz \
       -h <(echo '##INFO=<ID='${i}',Number=1,Type=Float,Description="Proportion of cells in which the site has been called as true positive based on the linked-read approach">\n') \
       -c CHROM,FROM,TO,-,-,-,-,-,-,-,-,${i},- \
       --output ${VCF}.v9.vcf \
       ${VCF}.v8.vcf
bcftools annotate \
       --annotations ${outdir}/linked-reads-annotations.tab.gz \
       -h <(echo '##INFO=<ID='${j}',Number=1,Type=Float,Description="Proportion of cells in which the site has been called as false positive based on the linked-read approach">\n') \
       -c CHROM,FROM,TO,-,-,-,-,-,-,-,-,-,${j} \
       --output ${VCF}.LR.vcf \
       ${VCF}.v9.vcf


echo "Remove intermediate files"
rm ${VCF}.v9.vcf ${VCF}.v8.vcf ${VCF}.v7.vcf ${VCF}.v6.vcf ${VCF}.v5.vcf ${VCF}.v4.vcf ${VCF}.v3.vcf ${VCF}.v2.vcf ${VCF}.v1.vcf
bgzip -f ${VCF}.LR.vcf
tabix -p vcf ${VCF}.LR.vcf.gz


##############################################################
# Select our training dataset and create a tabular structure #
##############################################################

bcftools filter --include 'LR_class=="tp" || LR_class=="fp" || LR_class=="uncertain"' ${VCF}.LR.vcf.gz -O z > \
        ${VCF}.linkedreadstested.vcf.gz
        

module purge
module load GATK/4.5.0.0-GCCcore-12.3.0-Java-17
module load BCFtools/1.19-GCC-13.2.0

gatk IndexFeatureFile \
     -I ${VCF}.LR.vcf.gz
OUT=$(echo ${VCF}.LR.vcf.gz | sed 's/.vcf//' | sed 's/.gz//')
info_fields=$(bcftools view --header-only ${VCF}.LR.vcf.gz | grep "^##INFO=<" | sed 's/##INFO=<ID=//' | sed 's/,.*//' | grep -v "^$" | awk '{print "-F "$0}' | tr -t "\n" " " | sed 's/ -F $//' )
genotype_fields=$(bcftools view --header-only ${VCF}.LR.vcf.gz | grep "^##FORMAT=<" | sed 's/##FORMAT=<ID=//' | sed 's/,.*//' | \
        tr -t "\n" " " | sed 's/ / -GF /g' | awk 'BEGIN{printf "-GF "}{print $0}' | sed 's/ -GF $//')


genotype_fields="-GF AD"
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



VCF=${VCF}.linkedreadstested.vcf.gz

gatk IndexFeatureFile \
     -I $VCF
OUT=$(echo $VCF | sed 's/.vcf//' | sed 's/.gz//')
info_fields=$(bcftools view --header-only $VCF | grep "^##INFO=<" | sed 's/##INFO=<ID=//' | sed 's/,.*//' | grep -v "^$" | awk '{print "-F "$0}' | tr -t "\n" " " | sed 's/ -F $//' )
genotype_fields=$(bcftools view --header-only $VCF | grep "^##FORMAT=<" | sed 's/##FORMAT=<ID=//' | sed 's/,.*//' | \
        tr -t "\n" " " | sed 's/ / -GF /g' | awk 'BEGIN{printf "-GF "}{print $0}' | sed 's/ -GF $//')

echo ""
echo ${info_fields}
echo ${genotype_fields}

genotype_fields="-GF AD"
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



