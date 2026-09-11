#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 4
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu 50G

module purge
module load bcftools

VCF_PATH=$1
VCF=$(echo $VCF_PATH | sed 's/.vcf//' | sed 's/.gz//')
workdir=$(dirname $VCF)

echo "VCF: "$VCF
echo "Outputing to: "$workdir

# depth I want to assign to the fake outgroup (called zeros)
depth=15 # we sequence at ~15x we changed it
SCRIPTDIR=/gpfs/commons/groups/landau_lab/ResolveOME/Scripts/

# Define the intermediate directory
intermediate_dir="${workdir}/mutmapping/"
# Ensure the intermediate directory exists
mkdir -p "$intermediate_dir"


#########################################
# Create the alt and total count tables #
#########################################


echo "Create alternative and total read count tables..."
sample_names=$(bcftools query -l "$VCF_PATH")

# Create the header with sample names as column names for total depth
total_depth_header="CHROM_POS_REF_ALT"
for sample_name in $sample_names; do
    total_depth_header="$total_depth_header\t${sample_name}"
done

# Print the total depth header
echo -e "$total_depth_header" | sed 's/___//g' | sed 's/\t$//' > "${intermediate_dir}/total_depth.txt"

# Process each variant and calculate the total depth
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' "$VCF_PATH" | \
awk -F"\t" -v intermediate_dir="$intermediate_dir" -v sample_names="$sample_names" '
BEGIN {
    OFS = "\t"
}
{
    chrom_pos_ref_alt = $1 ":" $2 "_" $3 ">" $4
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
}' | uniq  >> "${intermediate_dir}/total_depth.txt"

# Create the header with sample names as column names for alt depth
alt_depth_header="CHROM_POS_REF_ALT"
for sample_name in $sample_names; do
    alt_depth_header="$alt_depth_header\t${sample_name}"
done

# Print the alt depth header
echo -e "$alt_depth_header" | sed 's/___//g' | sed 's/\t$//' > "${intermediate_dir}/alt_depth.txt"

# Process each variant and print the results for alt depth
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' "$VCF_PATH" | \
awk -F"\t" -v intermediate_dir="$intermediate_dir" -v sample_names="$sample_names" '
BEGIN {
    OFS = "\t"
}
{
    chrom_pos_ref_alt = $1 ":" $2 "_" $3 ">" $4
    printf "%s", chrom_pos_ref_alt

    # Print the allele depths for each sample for alt depth
    split($0, fields, "\t")
    for (i = 5; i <= NF; i++) {
        n = split(fields[i], ad, ",")
        second_allele = (n >= 2 ? ad[2] : 0)
        printf "\t%s", second_allele
    }
    print ""
}' | uniq >> "${intermediate_dir}/alt_depth.txt"



echo "Matrix total: "${intermediate_dir}/total_depth.txt
echo "Matrix alt: "${intermediate_dir}/alt_depth.txt


bcftools query -f '%CHROM %POS %POS %Func.refGene %Gene.refGene:%AAChange.refGene %OGAINI %ExonicFunc.refGene\n' \
        ${VCF}.vcf.gz | \
        awk '{if ($4 !~ /(^exonic$)/ ){$5="NA"}; print $0}' | \
        awk '{if ($6=="."){$6="NA"}; print $0}' | \
	sed -e "s/:NM.*:p\./\./g" | \
        awk '{$2=$2-1;print $0}' > \
        ${workdir}/mutmapping/genes.bed

