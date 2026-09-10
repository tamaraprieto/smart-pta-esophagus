#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type FAIL
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu 10G
#SBATCH --cpus-per-task 4

##################################
# Define arguments and set paths #
##################################

patient=$1
tag=$2
intervals=$3
sex=$4
workdir=/gpfs/commons/groups/landau_lab/ResolveOME/StartDir/${patient}/manual_collate/
vcfname=${patient}_${tag}.glnexus
input_vcf=${workdir}/annovar/${tag}/${vcfname}.annovar.vcf.gz

scriptdir=/gpfs/commons/groups/landau_lab/ResolveOME/Scripts/
resources=/gpfs/commons/groups/landau_lab/ResolveOME/Resources/
reference=references/hg38/Homo_sapiens_assembly38.fasta

module purge
module load bedops/2.4.41
module load BCFtools/1.9
echo "nthreads: ${SLURM_CPUS_PER_TASK}; mem-per-thread: "${SLURM_MEM_PER_CPU}
binSize=$(awk '{print $1"\t"0"\t"$2}' ${resources}${reference}.fai | head -n24 | awk -v desired_intervals=$intervals '{sum+=$3}END{print int(sum/desired_intervals)}')
echo "bin size: "$binSize
interval=$(awk '{print $1"\t"0"\t"$2}' ${resources}${reference}.fai | head -n24 | sort-bed - | bedops --chop ${binSize} - | sort --version-sort | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{print $1":"$2"-"$3}')
echo "interval: "$interval

############################
# Create input for sequoia #
############################

echo "Create downsampled vcf..."
echo ""
# Define the intermediate directory
intermediate_dir="${workdir}sequoia"
# Ensure the intermediate directory exists
mkdir -p "$intermediate_dir"
# Define the VCF file path
VCF_PATH=${intermediate_dir}/${patient}_${tag}.${interval}.vcf.gz
bcftools view --threads $SLURM_CPUS_PER_TASK $input_vcf \
	$interval -O z > $VCF_PATH

# Check if the file doesn't contain lines starting with "chr"
if ! zgrep -q "^chr" "$VCF_PATH"; then
    # Rename the file
    mv "$VCF_PATH" "${intermediate_dir}/intervals/${patient}_${tag}.sequoia.vaf.eso.${SLURM_ARRAY_TASK_ID}.vcf.gz"
    tabix -p vcf "${intermediate_dir}/intervals/${patient}_${tag}.sequoia.vaf.eso.${SLURM_ARRAY_TASK_ID}.vcf.gz"
    echo "Empty VCF created."
    exit
else
    echo "File contains lines starting with 'chr'."
fi


echo "Create alternative and total read count tables..."
sample_names=$(bcftools query -l "$VCF_PATH")

# Create the header with sample names as column names for total depth
total_depth_header="CHROM_POS_REF_ALT"
for sample_name in $sample_names; do
    # Replace any dashes with underscores
    sample_name=${sample_name//-/_}
    total_depth_header="$total_depth_header\t${sample_name}"
done

# Print the total depth header
echo -e "$total_depth_header" | sed 's/___//g' | sed 's/\t$//' > "${intermediate_dir}/total_depth_${SLURM_ARRAY_TASK_ID}.txt"

# Process each variant and calculate the total depth
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' "$VCF_PATH" | \
awk -F"\t" -v intermediate_dir="$intermediate_dir" -v sample_names="$sample_names" '
BEGIN {
    OFS = "\t"
}
{
    chrom_pos_ref_alt = $1 "_" $2 "_" $3 "_" $4
    printf "%s", chrom_pos_ref_alt

    # Calculate total depth for each sample
    split($0, fields, "\t")
    for (i = 5; i <= NF; i++) {
        n = split(fields[i], ad, ",")
        total_depth = 0
        for (j = 1; j <= n; j++) {
            total_depth += ad[j]
        }
        printf "\t%s", total_depth
    }
    print ""
}' | uniq  >> "${intermediate_dir}/total_depth_${SLURM_ARRAY_TASK_ID}.txt"

# Create the header with sample names as column names for alt depth
alt_depth_header="CHROM_POS_REF_ALT"
for sample_name in $sample_names; do
    # Replace any dashes with underscores
    sample_name=${sample_name//-/_}
    alt_depth_header="$alt_depth_header\t${sample_name}"
done

# Print the alt depth header
echo -e "$alt_depth_header" | sed 's/___//g' | sed 's/\t$//' > "${intermediate_dir}/alt_depth_${SLURM_ARRAY_TASK_ID}.txt"

# Process each variant and print the results for alt depth
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' "$VCF_PATH" | \
awk -F"\t" -v intermediate_dir="$intermediate_dir" -v sample_names="$sample_names" '
BEGIN {
    OFS = "\t"
}
{
    chrom_pos_ref_alt = $1 "_" $2 "_" $3 "_" $4
    printf "%s", chrom_pos_ref_alt

    # Print the allele depths for each sample for alt depth
    split($0, fields, "\t")
    for (i = 5; i <= NF; i++) {
        n = split(fields[i], ad, ",")
        second_allele = (n >= 2 ? ad[2] : 0)
        printf "\t%s", second_allele
    }
    print ""
}' | uniq >> "${intermediate_dir}/alt_depth_${SLURM_ARRAY_TASK_ID}.txt" # there were some duplicated lines for some intervals, not sure why



echo "Matrix total: "${intermediate_dir}/total_depth_${SLURM_ARRAY_TASK_ID}.txt
echo "Matrix alt: "${intermediate_dir}/alt_depth_${SLURM_ARRAY_TASK_ID}.txt


############################
# Run sequoia per interval #
############################

echo "Run sequoia Rscript..."
echo ""
# matrices of read depth and variant read depth 
# matrix format: variants × samples
module purge
module load Miniconda3/23.10.0-1
source activate /gpfs/commons/groups/landau_lab/tprieto/conda/myR4
which Rscript

Rscript ${scriptdir}/sequoia/sequoia_filters_germline.R \
	-r ${intermediate_dir}/total_depth_${SLURM_ARRAY_TASK_ID}.txt \
	-v ${intermediate_dir}/alt_depth_${SLURM_ARRAY_TASK_ID}.txt \
	-o ${intermediate_dir}/intervals/ \
	--sex=${sex} \
	--maxIt_pval=25 \
	--ncores=${SLURM_CPUS_PER_TASK} 
	
conda deactivate


###################################################
# Add annotations to VCF so I can filter later on #
###################################################

module purge
module load bcftools
# Load sequoia output
sequoia_filters=${intermediate_dir}/intervals/Filter_annotations${SLURM_ARRAY_TASK_ID}

echo "> Create the tsv with annotations and index it..."
awk 'BEGIN{printf "#"}{print $0}' ${sequoia_filters}.txt > \
         ${sequoia_filters}.tab
bgzip -c ${sequoia_filters}.tab > \
        ${sequoia_filters}.tab.gz
tabix -c 1 -f -b 2 -0 -e 3 \
        --comment '#' \
        ${sequoia_filters}.tab.gz

grep "^chrom"  ${sequoia_filters}.txt | sed 's/#//' | sed 's/\t/\n/g' > ${sequoia_filters}.names.txt
# chrom start end Ref Alt Mean_Depth Depth_filter Germline_qval Germline Rho Beta_binomial
sequoia_meandepth=$(head -n 6 ${sequoia_filters}.names.txt | tail -n 1)
sequoia_germlinepval=$(head -n 8 ${sequoia_filters}.names.txt | tail -n 1)
sequoia_rhobetabinomial=$(head -n 10 ${sequoia_filters}.names.txt | tail -n 1)
randtopo_branch_filter=$(head -n 12 ${sequoia_filters}.names.txt | tail -n 1)
randtopo_pvalextbranch=$(head -n 14 ${sequoia_filters}.names.txt | tail -n 1)
randtopo_pvalotherbranches=$(head -n 15 ${sequoia_filters}.names.txt | tail -n 1)

rm ${sequoia_filters}.names.txt
echo "> Annotate the vcf files..."
bcftools annotate --threads $SLURM_CPUS_PER_TASK \
	--annotations ${sequoia_filters}.tab.gz \
       -h <(echo '##INFO=<ID='${sequoia_meandepth}',Number=1,Type=Float,Description="Mean depth across cells calculated from AD in sequoia">\n') \
       -c CHROM,FROM,TO,-,-,${sequoia_meandepth},-,-,-,-,-,-,-,-,- \
	-O z \
       --output "${intermediate_dir}/intervals/${patient}_${tag}.intermediate1.${SLURM_ARRAY_TASK_ID}.vcf.gz" \
       $VCF_PATH
bcftools annotate --threads $SLURM_CPUS_PER_TASK \
        --annotations ${sequoia_filters}.tab.gz \
       -h <(echo '##INFO=<ID='${sequoia_germlinepval}',Number=1,Type=Float,Description="P-value of a site being germline from sequoia">\n') \
       -c CHROM,FROM,TO,-,-,-,-,${sequoia_germlinepval},-,-,-,-,-,-,- \
	-O z \
       --output "${intermediate_dir}/intervals/${patient}_${tag}.intermediate2.${SLURM_ARRAY_TASK_ID}.vcf.gz" \
       "${intermediate_dir}/intervals/${patient}_${tag}.intermediate1.${SLURM_ARRAY_TASK_ID}.vcf.gz"
bcftools annotate --threads $SLURM_CPUS_PER_TASK \
        --annotations ${sequoia_filters}.tab.gz \
       -h <(echo '##INFO=<ID='${sequoia_rhobetabinomial}',Number=1,Type=Float,Description="Rho calculated by maximum likelihood for the overdispersion in sequoia">\n') \
       -c CHROM,FROM,TO,-,-,-,-,-,-,${sequoia_rhobetabinomial},-,-,-,-,- \
       -O z \
       --output "${intermediate_dir}/intervals/${patient}_${tag}.intermediate3.${SLURM_ARRAY_TASK_ID}.vcf.gz" \
       "${intermediate_dir}/intervals/${patient}_${tag}.intermediate2.${SLURM_ARRAY_TASK_ID}.vcf.gz"
bcftools annotate --threads $SLURM_CPUS_PER_TASK \
        --annotations ${sequoia_filters}.tab.gz \
       -h <(echo '##INFO=<ID='${randtopo_branch_filter}',Number=1,Type=String,Description="Germline or somatic based on mapping to a random topology">\n') \
       -c CHROM,FROM,TO,-,-,-,-,-,-,-,-,${randtopo_branch_filter},-,-,- \
       -O z \
       --output "${intermediate_dir}/intervals/${patient}_${tag}.intermediate4.${SLURM_ARRAY_TASK_ID}.vcf.gz" \
       "${intermediate_dir}/intervals/${patient}_${tag}.intermediate3.${SLURM_ARRAY_TASK_ID}.vcf.gz"
bcftools annotate --threads $SLURM_CPUS_PER_TASK \
        --annotations ${sequoia_filters}.tab.gz \
       -h <(echo '##INFO=<ID='${randtopo_pvalextbranch}',Number=1,Type=Float,Description="Pvalue external branch from treemut on a random topology">\n') \
       -c CHROM,FROM,TO,-,-,-,-,-,-,-,-,-,-,${randtopo_pvalextbranch},- \
       -O z \
       --output "${intermediate_dir}/intervals/${patient}_${tag}.intermediate5.${SLURM_ARRAY_TASK_ID}.vcf.gz" \
       "${intermediate_dir}/intervals/${patient}_${tag}.intermediate4.${SLURM_ARRAY_TASK_ID}.vcf.gz"
bcftools annotate --threads $SLURM_CPUS_PER_TASK \
        --annotations ${sequoia_filters}.tab.gz \
       -h <(echo '##INFO=<ID='${randtopo_pvalotherbranches}',Number=1,Type=Float,Description="Pvalue any other non-external branch from treemut on a random topology">\n') \
       -c CHROM,FROM,TO,-,-,-,-,-,-,-,-,-,-,-,${randtopo_pvalotherbranches} \
       -O z \
       --output "${intermediate_dir}/intervals/${patient}_${tag}.sequoia.${SLURM_ARRAY_TASK_ID}.vcf.gz" \
       "${intermediate_dir}/intervals/${patient}_${tag}.intermediate5.${SLURM_ARRAY_TASK_ID}.vcf.gz"

################
# Annotate VAF #
################

echo "> Annotating the VCF with VAF per sample"
# The FORMAT/VAF will be 0 instead of NA for those samples with no reads
bcftools +fill-tags "${intermediate_dir}/intervals/${patient}_${tag}.sequoia.${SLURM_ARRAY_TASK_ID}.vcf.gz" \
        --threads ${SLURM_CPUS_PER_TASK} \
        -O z -o "${intermediate_dir}/intervals/${patient}_${tag}.sequoia.vaf.${SLURM_ARRAY_TASK_ID}.vcf.gz" \
        -- \
        --tags FORMAT/VAF

#######################################
# Add esophageal-specific annotations #
#######################################

echo "> Create the tsv with annotations for esophageal mutations"
awk 'BEGIN{printf "#"}{print $0}' ${resources}/dna_mutations/Ogawa-Inigo-obs_mutations.liftover.hg38.bed > \
         ${intermediate_dir}/intervals/${SLURM_ARRAY_TASK_ID}.sc.annotations.tab
bgzip -c   ${intermediate_dir}/intervals/${SLURM_ARRAY_TASK_ID}.sc.annotations.tab > \
         ${intermediate_dir}/intervals/${SLURM_ARRAY_TASK_ID}.sc.annotations.tab.gz
bcftools annotate --threads $SLURM_CPUS_PER_TASK \
       --annotations  ${intermediate_dir}/intervals/${SLURM_ARRAY_TASK_ID}.sc.annotations.tab.gz \
       -h <(echo '##INFO=<ID='OGAINI',Number=1,Type=String,Description="Gene annotation from the Ogawa and Martincorena papers">\n') \
       -c CHROM,FROM,TO,OGAINI \
       --output "${intermediate_dir}/intervals/${patient}_${tag}.sequoia.vaf.eso.${SLURM_ARRAY_TASK_ID}.vcf.gz" \
       "${intermediate_dir}/intervals/${patient}_${tag}.sequoia.vaf.${SLURM_ARRAY_TASK_ID}.vcf.gz"
tabix -p vcf "${intermediate_dir}/intervals/${patient}_${tag}.sequoia.vaf.eso.${SLURM_ARRAY_TASK_ID}.vcf.gz"

echo "> Remove intermediate data..."
rm ${VCF_PATH}
rm ${intermediate_dir}/total_depth_${SLURM_ARRAY_TASK_ID}.txt
rm ${intermediate_dir}/alt_depth_${SLURM_ARRAY_TASK_ID}.txt
rm ${intermediate_dir}/intervals/${SLURM_ARRAY_TASK_ID}.sc.annotations*
rm "${intermediate_dir}/intervals/${patient}_${tag}.sequoia.vaf.${SLURM_ARRAY_TASK_ID}.vcf.gz"
rm "${intermediate_dir}/intervals/${patient}_${tag}.sequoia.${SLURM_ARRAY_TASK_ID}.vcf.gz"
rm "${intermediate_dir}/intervals/${patient}_${tag}.intermediate1.${SLURM_ARRAY_TASK_ID}.vcf.gz"
rm "${intermediate_dir}/intervals/${patient}_${tag}.intermediate2.${SLURM_ARRAY_TASK_ID}.vcf.gz"
rm "${intermediate_dir}/intervals/${patient}_${tag}.intermediate3.${SLURM_ARRAY_TASK_ID}.vcf.gz"
rm "${intermediate_dir}/intervals/${patient}_${tag}.intermediate4.${SLURM_ARRAY_TASK_ID}.vcf.gz"
rm "${intermediate_dir}/intervals/${patient}_${tag}.intermediate5.${SLURM_ARRAY_TASK_ID}.vcf.gz"
rm  ${sequoia_filters}.tab
rm  ${sequoia_filters}.tab.gz
rm  ${sequoia_filters}.tab.gz.tbi
rm ${sequoia_filters}.txt
