#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH -t 50:00:00
#SBATCH --mem 10G



# -------------------------------
# set vars
# -------------------------------
VCF=$1
BINSIZE=$2

WORKDIR=$(dirname $VCF)
WORKDIR=${WORKDIR}/
GINKGO_INST_DIR=/gpfs/commons/groups/landau_lab/tprieto/apps/ginkgo/

CELL=$(zcat ${VCF} | grep -v "^#" -B 1 | head -1 | cut -f 10- | tr "\t" "\n" | sed "${SLURM_ARRAY_TASK_ID}q;d")

module purge
#module load bcftools/1.15
module load BCFtools/1.9

# clear directory if needed
mkdir -p ${WORKDIR}/cnv
mkdir -p ${WORKDIR}/cnv/${BINSIZE}
rm -f ${WORKDIR}/cnv/${BINSIZE}/${CELL}_${BINSIZE}_avg.vaf.csv



echo "VCF: $VCF" 
echo "working directory: $WORKDIR"
echo "binsize: $BINSIZE"
echo "cell: $CELL"

# --------------------------------------------------------------
# for each region (bin), calculate average VAF per cell 
# pipe region to query VAF to awk and calculate average, write as a single col
# --------------------------------------------------------------
while read region; do 
        bcftools view \
                -Oz ${VCF} \
                -r ${region} \
                -s ${CELL} | \
        bcftools query -f '[ %VAF]\n' | 
                awk '{ VAFsum += sqrt(($1 - 0.5)^2) } END  { if (NR == 0) {print NA} else { print VAFsum/NR }}' \
                >> ${WORKDIR}cnv/${BINSIZE}/${CELL}_${BINSIZE}_avg.vaf.csv
done < ${GINKGO_INST_DIR}uploads/${BINSIZE}_regions.csv

