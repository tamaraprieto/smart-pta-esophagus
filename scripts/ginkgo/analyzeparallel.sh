#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH -t 150:00:00
#SBATCH --mem 50M

# ------------------------------------------------------------------------------
# -- Variables
# ------------------------------------------------------------------------------

home=/gpfs/commons/groups/landau_lab/tprieto/apps/ginkgo
dir=${home}/uploads/${1}
source ${dir}/config
distMet=$distMeth
#touch $dir/index.html

inFile=list
#statFile=status.xml

if [ ${chosen_genome} == "hg19" ];
then
  genome=${home}/genomes/${chosen_genome}/original
else
  genome=${home}/genomes/${chosen_genome} # use hg38
fi


# if [ "$rmpseudoautosomal" == "1" ];
# then
#   genome=${genome}/pseudoautosomal
# else
#   genome=${genome}/original
# fi

# ------------------------------------------------------------------------------
# -- Map Reads & Prepare Samples For Processing
# ------------------------------------------------------------------------------


file=$(sed -n "${SLURM_ARRAY_TASK_ID}p"  ${dir}/${inFile})
firstLineChr=$(zcat "${dir}/${file}" | head -n 1 | cut -f1)

if [[ "${firstLineChr}" != chr* ]]; then
  echo "Adding chr prefix to $file"

  zcat "${dir}/${file}" \
  | awk '{print "chr"$0}' \
  | ${home}/scripts/binUnsorted "${genome}/${binMeth}" \
      "$(wc -l < ${genome}/${binMeth})" \
      /dev/stdin \
      "$(basename "${file}" .bed.gz)" \
      "${dir}/${file}_mapped"

else
  ${home}/scripts/binUnsorted "${genome}/${binMeth}" \
    "$(wc -l < ${genome}/${binMeth})" \
    <(zcat -cd "${dir}/${file}") \
    "$(basename "${file}" .bed.gz)" \
    "${dir}/${file}_mapped"
fi
