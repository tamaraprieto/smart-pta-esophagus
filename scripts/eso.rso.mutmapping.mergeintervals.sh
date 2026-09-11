#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 4
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu 10G

module purge
module load bcftools/1.9

VCF_PATH=$1
VCF=$(echo $VCF_PATH | sed 's/.vcf//' | sed 's/.gz//')
workdir=$(dirname $VCF)

echo "VCF: "$VCF
#echo "Phylogeny: "$tree
echo "Outputing to: "$workdir

SCRIPTDIR=/gpfs/commons/groups/landau_lab/ResolveOME/Scripts/
  
intermediate_dir="${workdir}/mutmapping/"
maxIt_pval=100
maxlines=$(wc -l ${intermediate_dir}/genes.bed | awk '{print $1}')

###########################
# Combine the rds objects #
###########################

module purge
module load Miniconda3/23.10.0-1
source activate /gpfs/commons/groups/landau_lab/tprieto/conda/myR4
Rscript ${SCRIPTDIR}/rscripts/mutmapping.mergeintervals.R $intermediate_dir $maxIt_pval $maxlines 
source deactivate


#######################
# Merge the bed files #
#######################

outdir=${intermediate_dir}maxIt.${maxIt_pval}
bed_files=$(ls ${outdir}/intervals/MutationMappingAnnotations.*.bed | sort --version-sort)
first_file=$(ls $bed_files | head -n 1)
head -n 1 $first_file > ${outdir}/MutationMappingAnnotations.bed
for file in `ls ${bed_files[@]}`
do
grep -v "^chrom" $file 
done | sort -k1 --version-sort >> ${outdir}/MutationMappingAnnotations.bed


###########################################
# Annotate the vcf file with the bed file #
###########################################


module purge
module load BCFtools/1.19-GCC-13.2.0
module load VCFtools/0.1.16-GCC-12.3.0
module load BEDTools/2.31.0-GCC-12.3.0

echo "> Create the tsv tab with annotations"
awk -F'\t' 'BEGIN{OFS = FS;printf "#"}{print $0}' ${outdir}/MutationMappingAnnotations.bed  > \
         ${outdir}/mutmapping.tab
bgzip -c ${outdir}/mutmapping.tab > ${outdir}/mutmapping.tab.gz
tabix -c 1 -f -b 2 -0 -e 3 \
        --comment '#' \
        ${outdir}/mutmapping.tab.gz

grep "^chrom" ${outdir}/MutationMappingAnnotations.bed | sed 's/#//' | sed 's/\t/\n/g' > ${outdir}/newfilters.txt
#1:chrom,2:start,3:end,4:branch_filter,5:mapping_filter,6:branch,7:pvalue_branch,8:pvalue_otherbranches,9:mapped
a=$(head -n 4 ${outdir}/newfilters.txt | tail -n 1)
b=$(head -n 5 ${outdir}/newfilters.txt | tail -n 1)
c=$(head -n 6 ${outdir}/newfilters.txt | tail -n 1)
d=$(head -n 7 ${outdir}/newfilters.txt | tail -n 1)
e=$(head -n 8 ${outdir}/newfilters.txt | tail -n 1)
f=$(head -n 9 ${outdir}/newfilters.txt | tail -n 1)
rm ${outdir}/newfilters.txt

bcftools annotate \
       --annotations ${outdir}/mutmapping.tab.gz \
       -h <(echo '##INFO=<ID='${a}',Number=1,Type=String,Description="Germline or somatic annotation based on mapping">\n') \
       -c CHROM,FROM,TO,${a},-,-,-,-,- \
       --output  ${VCF}.treemut1.vcf \
       ${VCF}.vcf.gz
bcftools annotate \
       --annotations ${outdir}/mutmapping.tab.gz \
       -h <(echo '##INFO=<ID='${b}',Number=1,Type=String,Description="Congruence of the variant based on mapping">\n') \
       -c CHROM,FROM,TO,-,${b},-,-,-,- \
       --output  ${VCF}.treemut2.vcf \
       ${VCF}.treemut1.vcf
tip_order=$(cat ${outdir}/PhylogenyMappingTipOrder.txt | tr -t "\n" " ")
bcftools annotate \
       --annotations ${outdir}/mutmapping.tab.gz \
       -h <(echo '##INFO=<ID='${c}',Number=1,Type=String,Description="Branch to which the mutation is mapped '${tip_order}'">\n') \
       -c CHROM,FROM,TO,-,-,${c},-,-,- \
       --output  ${VCF}.treemut3.vcf \
       ${VCF}.treemut2.vcf
bcftools annotate \
       --annotations ${outdir}/mutmapping.tab.gz \
       -h <(echo '##INFO=<ID='${d}',Number=1,Type=Float,Description="p-value of the mutation being mapped to the branch specified">\n') \
       -c CHROM,FROM,TO,-,-,-,${d},-,- \
       --output  ${VCF}.treemut4.vcf \
       ${VCF}.treemut3.vcf 
bcftools annotate \
       --annotations ${outdir}/mutmapping.tab.gz \
       -h <(echo '##INFO=<ID='${e}',Number=1,Type=Float,Description="p-value of the mutation being mapped to any other branch than the specified">\n') \
       -c CHROM,FROM,TO,-,-,-,-,${e},- \
       --output  ${VCF}.treemut5.vcf \
       ${VCF}.treemut4.vcf  
bcftools annotate \
       --annotations ${outdir}/mutmapping.tab.gz \
       -h <(echo '##INFO=<ID='${f}',Number=1,Type=Integer,Description="Mutations being significantly mapped (1) or not (0)">\n') \
       -c CHROM,FROM,TO,-,-,-,-,-,${f} \
       --output  ${VCF}.treemut6.vcf.gz \
       ${VCF}.treemut5.vcf
rm ${VCF}.treemut5.vcf ${VCF}.treemut4.vcf ${VCF}.treemut3.vcf ${VCF}.treemut2.vcf ${VCF}.treemut1.vcf

bcftools filter --include 'outbranch_filter=="germline"' ${VCF}.treemut6.vcf.gz -O z > \
        ${VCF}.treemut.germline.vcf.gz

bcftools filter --include 'outbranch_filter=="somatic" && mapping_filter=="congruent"' ${VCF}.treemut6.vcf.gz -O z > \
        ${VCF}.treemut.vcf.gz

bcftools filter --include 'outbranch_filter=="somatic" && mapping_filter=="inconsistent"' ${VCF}.treemut6.vcf.gz -O z > \
        ${VCF}.treemut.inconsistent.vcf.gz
        
        
bcftools query \
        --print-header -f '%CHROM\t%POS-%REF-%ALT\t%Gene.refGene\t%Func.refGene\t%ExonicFunc.refGene\t%OGAINI\t%QUAL[\t%VAF:%DP]\n' \
        ${VCF}.treemut.vcf.gz | \
        sed 's/^# //' | \
        awk '{if ($1 ~ /^\[/){gsub(/\[[0-9]+\]/,"",$0); print $0}else{print $0}}' > \
        ${VCF}.treemut.tsv
    
higher_sample_index=$(grep -m 1 "^CHR" ${VCF}.treemut.tsv | awk '{print NF}')
myfile=${VCF}.treemut.tsv
    
    for index in `seq 8 $higher_sample_index`
      do
        samplename=$(grep -m 1 "^CHR" $myfile | cut -f $index | sed 's/:.*//')
        if [[ $index = 8 ]]
          then
          grep -m 1 "^CHROM" $myfile | cut -f 1-7  | awk '{print $0"\tVAF\tdepth\tsample"}'
        fi
        grep -v "^CHROM" $myfile | awk \
          -v myindex=$index -v mysample=$samplename '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$myindex"\t"mysample}'  |\
        grep -v "0:" | sed 's/:/\t/'
    done | sort --version-sort > ${VCF}.treemut.longformat.tsv

rm ${VCF}.treemut6.vcf.gz              
