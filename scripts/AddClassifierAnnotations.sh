#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu 120G

module purge
module load VCFtools/0.1.16-GCC-12.3.0
module load bedtools/2.31.0-GCC-12.3.0
module load BCFtools/1.9

VCF_PATH=$1
PATIENT=$2
VCF_TABLE=$(echo $VCF_PATH | sed 's/.vcf.gz/.table/')
VCF=$(echo $VCF_PATH | sed 's/.vcf//' | sed 's/.gz//')
workdir=$(dirname $VCF)
outdir=${workdir}/classifier/
mkdir -p $outdir

echo "VCF: "$VCF
echo "Outputing to: "$outdir

###############################
# OBTAIN BED WITH ANNOTATIONS #
###############################

# There are some mutations occuring in the same position (multiple alleles)
# The annotations will be the same because annovar was run before normalization
bcftools view ${VCF}.vcf.gz | \
        grep  -v "#" | sort -k1,1 -k2,2n | awk \
        'BEGIN{printf "#chrom\tstart\tend\n"}
        {OFS="\t"; print $1, $2-1, $2}' | uniq  > ${outdir}/coordinates.bed

#################################
# Obtain the nucleotide context #
#################################

CONTEXT_LEN=10
REFERENCE=/gpfs/commons/groups/landau_lab/ResolveOME/Resources/references/hg38/Homo_sapiens_assembly38.fasta
awk -v var=$CONTEXT_LEN \
        '{OFS="\t"; print $1, $2-var, $3+var}' ${outdir}/coordinates.bed | \
        grep -v "^#" | \
        bedtools getfasta -fi $REFERENCE -bed - | \
        grep -v "^>" | awk 'BEGIN{print "PTA_20bpCONTEXTwithREF"}{print toupper($0)}' > ${outdir}/context.bed
        
#########################################################################
# Obtain PER SITE GENEBODY, SIMPLEREPEAT, STRAND and REPLICATION TIMING #
#########################################################################

# closest 
resources=/gpfs/commons/groups/landau_lab/ResolveOME/Resources/
ptato_resources=${resources}/ptato/
genebodystrand_bed=${ptato_resources}genebody.hg38.bed
simplerepeats_bed=${ptato_resources}simpleRepeats.hg38.bed
replicationtime_bed=${ptato_resources}all_RepliSeq_median.hg38.bed

# Create resource files only if they do not already exist
if [ ! -f "$genebodystrand_bed" ]; then
    awk '{if ($1~/([0-9]$|X$|Y$)/){print $0}}' \
        ${ptato_resources}/genebody.sorted.bed | \
        awk '{print "chr"$0}' | \
        sort -k1,1 -k2,2n > "$genebodystrand_bed"
fi
if [ ! -f "$genebodyonly_bed" ]; then
    awk '{print $1"\t"$2"\t"$3"\t"$4}' \
        "$genebodystrand_bed" > "$genebodyonly_bed"
fi
if [ ! -f "$strandonly_bed" ]; then
    awk '{print $1"\t"$2"\t"$3"\t"$6}' \
        "$genebodystrand_bed" > "$strandonly_bed"
fi
if [ ! -f "$simplerepeats_bed" ]; then
    awk '{if ($1~/([0-9]$|X$|Y$)/){print $0}}' \
        ${ptato_resources}/simpleRepeats.feature.sorted.bed | \
        awk '{print "chr"$0}' | \
        sort -k1,1 -k2,2n > "$simplerepeats_bed"
fi
if [ ! -f "$replicationtime_bed" ]; then
    awk '{if ($1~/([0-9]$|X$|Y$)/){print $0}}' \
        ${ptato_resources}/all_RepliSeq_median.sorted.bed | \
        awk '{print "chr"$0}' | \
        sort -k1,1 -k2,2n > "$replicationtime_bed"
fi

bedtools \
    closest \
    -a ${outdir}/coordinates.bed \
    -b ${ptato_resources}/genebodyonly.hg38.bed \
    -D a | bedtools merge -i - -d -1 -c 8 -o distinct | \
    awk 'BEGIN{print "PTA_GENEBODY"}{print $4}' \
    > ${outdir}/genebody.bed

bedtools \
    closest \
    -a ${outdir}/coordinates.bed \
    -b ${ptato_resources}/strandonly.hg38.bed | \
	bedtools merge -i - -d -1 -c 7 -o distinct | \
    awk 'BEGIN{print "PTA_STRAND"}{print $4}' \
    > ${outdir}/strand.bed

bedtools \
    closest \
    -a ${outdir}/coordinates.bed \
    -b ${simplerepeats_bed} \
    -D a | bedtools merge -i - -d -1 -c 7 -o distinct | \
    awk 'BEGIN{print "PTA_SIMPLEREPEAT"}{print $4}' \
    > ${outdir}/simplerepeat.bed
    
bedtools \
    closest \
    -a ${outdir}/coordinates.bed \
    -b ${replicationtime_bed} | \
    bedtools merge -i - -d -1 -c 7 -o median | \
    awk 'BEGIN{print "PTA_REPLISEQ"}{print $4}' \
    > ${outdir}/replicationtime.bed    

echo "make sure the length of all the files is the same: "
wc -l ${outdir}/*bed

##########################################
# CREATE A BED FILE WITH ALL ANNOTATIONS #
##########################################

paste ${outdir}/coordinates.bed ${outdir}/context.bed \
        ${outdir}/genebody.bed \
        ${outdir}/strand.bed \
        ${outdir}/simplerepeat.bed \
        ${outdir}/replicationtime.bed \
        -d'\t' > ${outdir}/classifier-annotations-unsorted.bed
  
  
# sort all except first row  
head -n 1 ${outdir}/classifier-annotations-unsorted.bed > ${outdir}/classifier-annotations.bed
tail -n +2 ${outdir}/classifier-annotations-unsorted.bed | sort --version-sort >> ${outdir}/classifier-annotations.bed
rm ${outdir}/classifier-annotations-unsorted.bed


#########################################################
# Annotate the vcf file with the classifier annotations #
#########################################################


echo "> Create the tsv tab with annotations"
cp ${outdir}/classifier-annotations.bed ${outdir}/classifier-annotations.tab
bgzip -c ${outdir}/classifier-annotations.tab > ${outdir}/classifier-annotations.tab.gz
rm ${outdir}/classifier-annotations.tab.gz.tbi
tabix -c 1 -f -b 2 -0 -e 3 \
        --comment '#' \
        ${outdir}/classifier-annotations.tab.gz

grep "^#chrom" ${outdir}/classifier-annotations.bed | sed 's/#//' | sed 's/\t/\n/g' > ${outdir}/newfilters.txt
# 1-:chrom   
# 2-:start
# 3-:end
# 4-a:PTA_20bpCONTEXTwithREF
# 5-b:PTA_GENEBODY        
# 6-c:PTA_STRAND      
# 7-d:PTA_SIMPLEREPEAT    
# 8-e:PTA_REPLISEQ
a=$(head -n 4 ${outdir}/newfilters.txt | tail -n 1)
b=$(head -n 5 ${outdir}/newfilters.txt | tail -n 1)
c=$(head -n 6 ${outdir}/newfilters.txt | tail -n 1)
d=$(head -n 7 ${outdir}/newfilters.txt | tail -n 1)
e=$(head -n 8 ${outdir}/newfilters.txt | tail -n 1)
cat ${outdir}/newfilters.txt
rm ${outdir}/newfilters.txt

bcftools annotate \
       --annotations ${outdir}/classifier-annotations.tab.gz \
       -h <(echo '##INFO=<ID='${a}',Number=1,Type=String,Description="10bp sequence context upstream and downstream of the SNV. The reference base is included in the middle. The context has not been reverse complemented for any mutation type.">\n') \
       -c CHROM,FROM,TO,${a},-,-,-,-,- \
       --output  ${VCF}.v1.vcf \
       ${VCF}.vcf.gz
bcftools annotate \
       --annotations ${outdir}/classifier-annotations.tab.gz \
       -h <(echo '##INFO=<ID='${b}',Number=1,Type=Integer,Description="Distance to the closest gene. When the mutation is present in the gene then the distance is zero">\n') \
       -c CHROM,FROM,TO,-,${b},-,-,-,- \
       --output  ${VCF}.v2.vcf \
       ${VCF}.v1.vcf
bcftools annotate \
       --annotations ${outdir}/classifier-annotations.tab.gz \
       -h <(echo '##INFO=<ID='${c}',Number=1,Type=String,Description="Strand of the DNA which is transcribed in the closest gene to the SNV">\n') \
       -c CHROM,FROM,TO,-,-,${c},-,-,- \
       --output  ${VCF}.v3.vcf \
       ${VCF}.v2.vcf
bcftools annotate \
       --annotations ${outdir}/classifier-annotations.tab.gz \
       -h <(echo '##INFO=<ID='${d}',Number=1,Type=Integer,Description="Distance to the closest simple repeat">\n') \
       -c CHROM,FROM,TO,-,-,-,${d},-,- \
       --output  ${VCF}.v4.vcf \
       ${VCF}.v3.vcf
bcftools annotate \
       --annotations ${outdir}/classifier-annotations.tab.gz \
       -h <(echo '##INFO=<ID='${e}',Number=1,Type=Float,Description="Replication time of the region in which the SNV occurs">\n') \
       -c CHROM,FROM,TO,-,-,-,-,${e},- \
       --output  ${VCF}.v5.vcf \
       ${VCF}.v4.vcf 

########################
# ADD REDI ANNOTATIONS #
########################

# In this case the reference and alternative alleles have to match to add the annotation
redi_database=${resources}rna_editing/REDIPortal_RNAEditSites.hg38.txt
if [ ! -f "${resources}rna_editing/REDIPortal_RNAEditSites.hg38.tsv" ]; then
    awk '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$6}' "$redi_database" | \
        grep -v "^Region" | \
        sort -k1,1 -k2,2n | \
        awk 'BEGIN{print "#CHROM\tSTART\tEND\tREF\tALT\tPTA_REDI"}{OFS="\t"; print $0}' \
        > "${resources}rna_editing/REDIPortal_RNAEditSites.hg38.tsv"
fi

cp ${resources}rna_editing/REDIPortal_RNAEditSites.hg38.tsv ${outdir}/classifier-annotations.tab
bgzip -c ${outdir}/classifier-annotations.tab > ${outdir}/classifier-annotations.tab.gz
rm ${outdir}/classifier-annotations.tab.gz.tbi
tabix -c 1 -f -b 2 -e 3 -0  \
        --comment '#' \
        ${outdir}/classifier-annotations.tab.gz

bcftools annotate -a ${outdir}/classifier-annotations.tab.gz \
  -c CHROM,FROM,TO,REF,ALT,INFO/PTA_REDI \
  -h <(echo "##INFO=<ID=PTA_REDI,Number=1,Type=String,Description='RNA editing database in which the change has been found'>") \
  --output ${VCF}.redi.vcf \
  ${VCF}.v5.vcf


#############################
# REMOVE INTERMEDIATE FILES #
#############################

echo "Remove intermediate files"
rm ${VCF}.v5.vcf ${VCF}.v4.vcf ${VCF}.v3.vcf ${VCF}.v2.vcf ${VCF}.v1.vcf


###########################################################
# OBTAIN ANNOTATIONS FOR THE CELLS CARRYING THE MUTATIONS #
###########################################################


# MAX VAF IN CELLS WITH MUTATION #
awk 'BEGIN {FS="\t"} 
NR==1 { 
    # Identify columns with .VAF suffix
    for (i=1; i<=NF; i++) {
        if ($i ~ /\.VAF$/) vaf_cols[i] = 1
    }
} 
NR > 1 { 
    max_vaf = 0
    # Iterate over the identified VAF columns
    for (i=1; i<=NF; i++) {
        if (vaf_cols[i] && $i != 0) {
            # Update max_vaf if the current value is greater
            if ($i > max_vaf) {
                max_vaf = $i
            }
        }
    }
    # Print the maximum VAF value
    print max_vaf
}' $VCF_TABLE | awk 'BEGIN{print "PTA_MAXVAF"}{print $0}' > ${outdir}/MaxVAFVariable.bed


# MAX ALT
awk 'BEGIN {FS="\t"} 
NR==1 { 
    # Identify columns with .VAF suffix
    for (i=1; i<=NF; i++) {
        if ($i ~ /\.AD$/) vaf_cols[i] = 1
    }
    # Print a new header
    print "PTA_MAXALT"
} 
NR > 1 { 
    max_alt = 0
    # Iterate over the identified VAF columns
    for (i=1; i<=NF; i++) {
        if (vaf_cols[i] && $i != 0) {
            # Split the VAF value into REF and ALT counts
            split($i, vaf_values, ",")
            alt_count = vaf_values[2] # Extract ALT count
            # Update max_alt if the current ALT count is greater
            if (alt_count > max_alt) {
                max_alt = alt_count
            }
        }
    }
    # Print the maximum ALT count
    print max_alt
}' $VCF_TABLE > ${outdir}/MaxAltVariable.bed

# MEAN MAD FOR CELLS WITH MUTATION #
MAD_file=/gpfs/commons/groups/landau_lab/ResolveOME/StartDir/MAD_scores.txt
zgrep -m 1 "^#CHR" $VCF_PATH | cut -f10- | \
  tr -s "\t" "\n" > ${outdir}/cellnames.txt  
# Create a file with MAD scores in the same order than the cells
while read cell
do
  grep -m 1 -e $cell"\s" $MAD_file
done < ${outdir}/cellnames.txt  | awk '{print $2}' | \
  paste -sd'\t' > ${outdir}/MADs.txt
# Create a file with VAFs
awk 'NR==1 { 
    for (i=1; i<=NF; i++) 
        if ($i ~ /\.VAF$/) vaf_cols[i]=1; 
    next 
} 
{ 
    for (i=1; i<=NF; i++) 
        if (i in vaf_cols) printf "%s ", $i; 
    printf "\n"
}' $VCF_TABLE > ${outdir}/zeros_ones.txt  
# Substitute each VAF of a cell higher than zero by the MAD score of the cell
# Select the MAX cell MAD per row
awk 'NR==FNR {for (i=1; i<=NF; i++) vec[i]=$i; next} {for (i=1; i<=NF; i++) $i = ($i > 0 ? vec[i] : $i)} 1' ${outdir}/MADs.txt ${outdir}/zeros_ones.txt | \
  awk '{ 
    min = 9999999; 
    for (i=1; i<=NF; i++) { 
        if ($i != 0) { 
            if ($i < min) min = $i; 
        } 
    } 
    print min; 
}' | awk 'BEGIN{print "PTA_MINMAD"}{print $0}' > ${outdir}/MinMADVariable.txt


# GENE EXPRESSION
# Get coverage of RNA at the somatic positions
rna_coverage=/gpfs/commons/groups/landau_lab/ResolveOME/StartDir/all_cell_rna/${PATIENT}.alleliccounts_RNA.tsv
awk '{print $1"\t"$2-1"\t"$2"\t"$8"\t"$7}' $rna_coverage > ${outdir}/RNACoverage.txt

# CONCATENATE BED FILES
awk '{print $1"\t"$2"\t"$3"\t"$4}' $VCF_TABLE | \
  paste - ${outdir}/MinMADVariable.txt ${outdir}/MaxVAFVariable.bed ${outdir}/MaxAltVariable.bed | \
  awk '{if ($1 ~ /^CHROM/){gsub($1, "#chrom\tSTART", $1);gsub($2, "END", $2); print $0} else {print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}}' OFS='\t' | \
  sed 's/NaN/\./g' | sed 's/NA/\./' > ${outdir}/classifier-annotations2.bed

ls ${outdir}/classifier-annotations2.bed

head -n 1 ${outdir}/classifier-annotations2.bed | awk '{printf $0}' | \
awk '{print $0"\tPTA_RNACOV\tPTA_RNAALT"}' > ${outdir}/classifier-annotations-extra.bed
bedtools intersect -a ${outdir}/classifier-annotations2.bed \
  -b ${outdir}/RNACoverage.txt -loj | cut -f 1-8,12-13 >> ${outdir}/classifier-annotations-extra.bed


###############################################################
# Annotate the vcf file with the EXTRA classifier annotations #
###############################################################

echo "> Create the tsv tab with annotations"
cp ${outdir}/classifier-annotations-extra.bed ${outdir}/classifier-annotations.tab
bgzip -c ${outdir}/classifier-annotations.tab > ${outdir}/classifier-annotations.tab.gz
rm ${outdir}/classifier-annotations.tab.gz.tbi
tabix -c 1 -f -b 2 -0 -e 3 \
        --comment '#' \
        ${outdir}/classifier-annotations.tab.gz

grep "^#chrom" ${outdir}/classifier-annotations-extra.bed | sed 's/#//' | sed 's/\t/\n/g' > ${outdir}/newfilters.txt
# 1-:chrom   
# 2-:START
# 3-:END
# 4-:REF
# 5-:ALT        
# 6-a:PTA_MINMAD      
# 7-b:PTA_MAXVAF    
# 8-c:PTA_MAXALT
# 9-d:PTA_RNACOV
# 10-e:PTA_RNAALT
# 11-f
a=$(head -n 6 ${outdir}/newfilters.txt | tail -n 1)
b=$(head -n 7 ${outdir}/newfilters.txt | tail -n 1)
c=$(head -n 8 ${outdir}/newfilters.txt | tail -n 1)
d=$(head -n 9 ${outdir}/newfilters.txt | tail -n 1)
e=$(head -n 10 ${outdir}/newfilters.txt | tail -n 1)

cat ${outdir}/newfilters.txt
rm ${outdir}/newfilters.txt

bcftools annotate \
       --annotations ${outdir}/classifier-annotations.tab.gz \
       -h <(echo '##INFO=<ID='${a}',Number=1,Type=Float,Description="Best (lowest) MAD score (coverage dispersion) of a cell carrying an alternative allele.">\n') \
       -c CHROM,FROM,TO,REF,ALT,${a},-,-,-,-,- \
       --output  ${VCF}.v1.vcf \
       ${VCF}.redi.vcf
bcftools annotate \
       --annotations ${outdir}/classifier-annotations.tab.gz \
       -h <(echo '##INFO=<ID='${b}',Number=1,Type=Float,Description="Highest cell variant-allele frequency at this SNV.">\n') \
       -c CHROM,FROM,TO,REF,ALT,-,${b},-,-,-,- \
       --output  ${VCF}.v2.vcf \
       ${VCF}.v1.vcf
bcftools annotate \
       --annotations ${outdir}/classifier-annotations.tab.gz \
       -h <(echo '##INFO=<ID='${c}',Number=1,Type=Integer,Description="Highest cell alternative allele count at this SNV">\n') \
       -c CHROM,FROM,TO,REF,ALT,-,-,${c},-,-,- \
       --output  ${VCF}.v3.vcf \
       ${VCF}.v2.vcf
bcftools annotate \
       --annotations ${outdir}/classifier-annotations.tab.gz \
       -h <(echo '##INFO=<ID='${d}',Number=1,Type=Integer,Description="Total coverage at the RNA bam files for this SNV">\n') \
       -c CHROM,FROM,TO,REF,ALT,-,-,-,${d},-,- \
       --output  ${VCF}.v4.vcf \
       ${VCF}.v3.vcf
bcftools annotate \
       --annotations ${outdir}/classifier-annotations.tab.gz \
       -h <(echo '##INFO=<ID='${e}',Number=1,Type=Float,Description="Alternative counts at the RNA bam files for this SNV">\n') \
       -c CHROM,FROM,TO,-,-,-,-,${e},- \
       --output ${VCF}.v5.vcf \
       ${VCF}.v4.vcf 


#############################
# REMOVE INTERMEDIATE FILES #
#############################

echo "Remove intermediate files"
rm ${VCF}.v4.vcf ${VCF}.v3.vcf ${VCF}.v2.vcf ${VCF}.v1.vcf

#######################
# ADD ADO ANNOTATIONS #
#######################

# ADO
ado_information=/gpfs/commons/groups/landau_lab/ResolveOME/StartDir/${PATIENT}/manual_collate/sequoia/phasing/PTA_ADO.txt
awk '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' $ado_information | sed 's/^CHROM\t-1/#CHROM\tSTART/' > ${outdir}/ADO.txt

cp ${outdir}/ADO.txt ${outdir}/classifier-annotations.tab
bgzip -c ${outdir}/classifier-annotations.tab > ${outdir}/classifier-annotations.tab.gz
rm ${outdir}/classifier-annotations.tab.gz.tbi
tabix -c 1 -f -b 2 -0 -e 3 \
        --comment '#' \
        ${outdir}/classifier-annotations.tab.gz

grep "^#CHROM" ${outdir}/classifier-annotations.tab | sed 's/#//' | sed 's/\t/\n/g' 

bcftools annotate -a ${outdir}/classifier-annotations.tab.gz \
  -c CHROM,FROM,TO,ALT,PTA_ADO200kbinterval,-,- \
  -h <(echo "##INFO=<ID=PTA_ADO200kbinterval,Number=1,Type=Integer,Description='200kb window identifier for PTA ADO analysis'>") \
  --output ${VCF}.v6.vcf \
  ${VCF}.v5.vcf
  
bcftools annotate -a ${outdir}/classifier-annotations.tab.gz \
  -c CHROM,FROM,TO,ALT,-,PTA_ADOResidualLoessPseudobulk,- \
  -h <(echo "##INFO=<ID=PTA_ADOResidualLoessPseudobulk,Number=1,Type=Float,Description='Residual of the Loess regression for hetSNPs in the window'>") \
  --output ${VCF}.v7.vcf \
  ${VCF}.v6.vcf  

bcftools annotate -a ${outdir}/classifier-annotations.tab.gz \
  -c CHROM,FROM,TO,ALT,-,-,PTA_ADOmaxbinompvalue \
  -h <(echo "##INFO=<ID=PTA_ADOmaxbinompvalue,Number=1,Type=Float,Description='Pvalue of a binomial test for the alternative counts of the somatic SNV given the predicted probability of alternative alleles by modelling the germline hetSNPs for PTA in 200kb windows'>") \
  --output ${VCF}.4classifier.vcf \
  ${VCF}.v7.vcf 

rm ${VCF}.v5.vcf ${VCF}.v6.vcf ${VCF}.v7.vcf 
bgzip -f ${VCF}.4classifier.vcf
tabix -p vcf ${VCF}.4classifier.vcf.gz

##############################################################
# Select our training dataset and create a tabular structure #
##############################################################

module purge
module load gatk/4.6.1.0-GCCcore-13.3.0-Java-17
module load BCFtools/1.9
VCF=${VCF}.4classifier.vcf.gz

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
