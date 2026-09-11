#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH -t 50:00:00
#SBATCH --mem 2G


# -------------------------------
# set vars
# -------------------------------
VCF=$1
BINSIZE=$2

WORKDIR=$(dirname $VCF)
WORKDIR=${WORKDIR}/
GINKGO_INST_DIR=/gpfs/commons/groups/landau_lab/tprieto/apps/ginkgo/


module purge

echo "VCF: $VCF" 
echo "working directory: $WORKDIR"
echo "binsize: $BINSIZE"


# -------------------------------
# concat the columns as one table: 
# -------------------------------
# remove existing file if needed
rm -f ${WORKDIR}cnv/${BINSIZE}_avg.vaf.concat.csv

# paste the names of the cells
zcat ${VCF} | grep -v "^#" -B 1 | head -1 | cut -f 10- | \
        tr "\t" "," | \
        sed 's/^/CHR,START,END,/' \
        >> ${WORKDIR}cnv/${BINSIZE}_avg.vaf.concat.csv

# paste in the chr cols + vaf info
paste <(cat ${GINKGO_INST_DIR}uploads/${BINSIZE}_regions.csv | tr ":" "," | tr "-" ",") \
        $(zcat ${VCF} | grep -v "^#" -B 1 | head -1 | cut -f 10- | \
                tr "\t" "\n" | sed 's@^@'"$WORKDIR"'/cnv\/'"$BINSIZE"'\/@' | sed 's/$/_'"$BINSIZE"'_avg.vaf.csv/'| tr "\n" " ") -d "," \
        >> ${WORKDIR}cnv/${BINSIZE}_avg.vaf.concat.csv
