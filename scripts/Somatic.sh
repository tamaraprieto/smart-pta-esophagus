#!/bin/bash
#SBATCH -N 1
#SBATCH --partition io
#SBATCH -n 1
#SBATCH --mail-type FAIL,COMPLETED
#SBATCH --cpus-per-task 4
#SBATCH -t 50:00:00
#SBATCH --mem-per-cpu 20G


##################################
# Define arguments and set paths #
##################################

patient=$1
tag=$2
echo ${patient} ${tag}
resources=/gpfs/commons/groups/landau_lab/ResolveOME/Resources/
reference=references/hg38/Homo_sapiens_assembly38.fasta
# define input
vcf_input_dir=/gpfs/commons/groups/landau_lab/ResolveOME/StartDir/${patient}/manual_collate/sequoia/
vcf_input_name=${patient}_${tag}.sequoia.vaf.eso
myvcf=${vcf_input_dir}${vcf_input_name}.vcf.gz

# define output
out_dir=${vcf_input_dir}/
mkdir $out_dir
vcf_output_name=myselectedvars

module purge
module load bcftools

#############################
# Add overlap between files #
#############################

overlap_info=/gpfs/commons/groups/landau_lab/ResolveOME/StartDir/callset_intersection/sites


grep "^#CHROM"  ${overlap_info}.tab | sed 's/#//' | sed 's/\t/\n/g' > ${overlap_info}.names.txt
# chrom start end Ref Alt Mean_Depth Depth_filter Germline_qval Germline Rho Beta_binomial
overlap=$(head -n 6 ${overlap_info}.names.txt | tail -n 1)
echo $overlap
rm ${overlap_info}.names.txt
echo "> Annotate the vcf files..."
bcftools annotate --threads $SLURM_CPUS_PER_TASK \
        --annotations ${overlap_info}.tab.gz \
       -h <(echo '##INFO=<ID='${overlap}',Number=1,Type=String,Description="Presence absence of mutations string across patients">\n') \
       -c CHROM,FROM,TO,-,-,${overlap} \
        -O z \
       --output ${vcf_input_dir}/${vcf_input_name}.overlap.vcf.gz \
       ${vcf_input_dir}/${vcf_input_name}.normalized.vcf.gz

################
# Annotate VAF #
################

filename=${VCF}.${chrom}
echo "> Annotating the VCF with VAF per sample"
# The FORMAT/VAF will be 0 instead of NA for those samples with no reads
bcftools +fill-tags ${vcf_input_dir}/${vcf_input_name}.overlap.vcf.gz \
        --threads $SLURM_CPUS_PER_TASK  \
        -O z -o ${vcf_input_dir}/${vcf_input_name}.overlap.vaf.vcf.gz \
        -- \
        --tags FORMAT/VAF


#######################################
# Select our set of somatic mutations #
#######################################

echo "Select the somatic mutations..."
bcftools filter  --threads $SLURM_CPUS_PER_TASK \
	-e 'INFO/avsnp150!="." && branch_filter == "germline" && Germline_qval > 0.00005' \
	${vcf_input_dir}/${vcf_input_name}.overlap.vaf.vcf.gz | \
	bcftools filter  --threads $SLURM_CPUS_PER_TASK -i \
	'(INFO/avsnp150=="." && Rho>0.5 &&  Mean_Depth<30 && Overlap==".") || (INFO/avsnp150!="." && Rho>0.65 && Mean_Depth<30 && Overlap==".")' \
        -O z -o ${out_dir}${vcf_output_name}.temp.vcf.gz 

# Explore discarded germline
#bcftools filter  --threads $SLURM_CPUS_PER_TASK -i 'INFO/avsnp150!="." && branch_filter == "germline" && Germline_qval > 0.00005' \
#        ${vcf_input_dir}/${vcf_input_name}.overlap.vaf.vcf.gz \
#        -O z -o ${vcf_input_dir}${vcf_output_name}.germline_pluspotentialresidualsomatic.vcf.gz


#Explore somatic and noise 
bcftools filter  --threads $SLURM_CPUS_PER_TASK -e 'INFO/avsnp150!="." && branch_filter == "germline" && Germline_qval > 0.00005' \
        ${vcf_input_dir}/${vcf_input_name}.overlap.vaf.vcf.gz \
	-O z -o ${vcf_input_dir}${vcf_output_name}.somaticplusnoise.vcf.gz

##############################################
# Select just snvs as input for cloudcellphy #
##############################################

echo "Select somatic mutations for the phylogenetic reconstruction..."
# this is for mapping
bcftools view  --threads $SLURM_CPUS_PER_TASK -i 'TYPE="snp" && Mean_Depth > 5' ${out_dir}${vcf_output_name}.temp.vcf.gz -O z \
	-o ${out_dir}${vcf_output_name}.snvswithsingletons.vcf.gz
# this is for cloudcellphy
bcftools view  --threads $SLURM_CPUS_PER_TASK -e 'COUNT(VAF>0)<2 || F_MISSING > 0.5 || Mean_Depth <= 9' ${out_dir}${vcf_output_name}.snvswithsingletons.vcf.gz -O z | \
	zgrep -v "^chr9" | bcftools view  --threads $SLURM_CPUS_PER_TASK -O z > \
        ${out_dir}${vcf_output_name}.snvs.vcf.gz

##################################
# Select not noisy indels shared #
##################################

echo "Select indels..."
bcftools view  --threads $SLURM_CPUS_PER_TASK -i 'TYPE="indel" && Mean_Depth > 5' ${out_dir}${vcf_output_name}.temp.vcf.gz \
	-O z -o	${out_dir}${vcf_output_name}.indels.vcf.gz

rm ${out_dir}${vcf_output_name}.temp.vcf.gz

############################
# Count number of variants #
############################

zgrep -v "^#" ${out_dir}${vcf_output_name}.snvs.vcf.gz | wc -l | awk '{print $0"\tsnvs_shared"}' > ${out_dir}${vcf_output_name}.count.txt
zgrep -v "^#"  ${out_dir}${vcf_output_name}.snvswithsingletons.vcf.gz | wc -l | awk '{print $0"\tall_snvs"}' >> ${out_dir}${vcf_output_name}.count.txt
zgrep -v "^#" ${out_dir}${vcf_output_name}.indels.vcf.gz | wc -l | awk '{print $0"\tindels"}' >> ${out_dir}${vcf_output_name}.count.txt
