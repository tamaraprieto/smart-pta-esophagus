#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 4
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu 10G


patient=$1
classifier_all_dir=/gpfs/commons/groups/landau_lab/ResolveOME/StartDir/classifier/
classifier_dir=/gpfs/commons/groups/landau_lab/ResolveOME/StartDir/${patient}/manual_collate/sequoia/classifier/

module load BCFtools/1.23.1
rm  ${classifier_dir}${patient}_CLASSIFIER.bed
zcat ${classifier_all_dir}${patient}_v3_predictions.csv.gz | \
  awk -F ',' '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$(NF-2)"\t"$(NF-1)"\t"$NF}' | \
  awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$7"\t"$8}' | \
  sed 's/^CHROM\t-1/#CHROM\tSTART/' > ${classifier_dir}${patient}_CLASSIFIER.bed

# ANNOTATE ORIGINAL VCF FILES WITH THE CLASSIFIER PROBABILITY
bgzip -c ${classifier_dir}${patient}_CLASSIFIER.bed > ${classifier_dir}/classifier-annotations.tab.gz
rm ${classifier_dir}/classifier-annotations.tab.gz.tbi
tabix -f -b 2 -e 3 -0  \
        --comment '#' \
        ${classifier_dir}/classifier-annotations.tab.gz
bcftools annotate -a ${classifier_dir}/classifier-annotations.tab.gz \
  -c CHROM,-,POS,REF,ALT,INFO/CLASSIFIER_PREDICTION \
  -h <(echo '##INFO=<ID=CLASSIFIER_PREDICTION,Number=1,Type=Float,Description="Classifier probability of being a SMART-PTA error">') \
  --output-type z \
  --output ${classifier_dir}/myselectedvars.snvswithsingletons.classifier.vcf.gz \
  /gpfs/commons/groups/landau_lab/ResolveOME/StartDir/${patient}/manual_collate/sequoia/myselectedvars.snvswithsingletons.LR.4classifier.vcf.gz
  
classifier_dir=/gpfs/commons/groups/landau_lab/ResolveOME/StartDir/${patient}/manual_collate/sequoia/classifier/

if [ "$patient" == "cu08" ]; then
# bioskryb v2
bcftools filter -e '(COUNT(VAF>0)<2 && CLASSIFIER_PREDICTION > 0.5) || (COUNT(VAF>0)<2 && (PTA_MAXVAF>0.9 || PTA_MAXVAF<0.25) && CHROM !~ "chr[XY]") || (COUNT(VAF>0)>=2 && CLASSIFIER_PREDICTION > 0.5)' ${classifier_dir}/myselectedvars.snvswithsingletons.classifier.vcf.gz \
  -O z -o ${classifier_dir}/myselectedvars.snvswithsingletons.classifierfiltered.vcf.gz
else
# bioskryb v1
bcftools filter -e '(COUNT(VAF>0)<2 && CLASSIFIER_PREDICTION > 0.5) || (COUNT(VAF>0)<2 && (PTA_MAXVAF>0.9 || PTA_MAXVAF<0.25) && CHROM !~ "chr[XY]") || (COUNT(VAF>0)>=2 && CLASSIFIER_PREDICTION > 0.9)' ${classifier_dir}/myselectedvars.snvswithsingletons.classifier.vcf.gz \
  -O z -o ${classifier_dir}/myselectedvars.snvswithsingletons.classifierfiltered.vcf.gz 
fi

bcftools index -f ${classifier_dir}/myselectedvars.snvswithsingletons.classifierfiltered.vcf.gz

