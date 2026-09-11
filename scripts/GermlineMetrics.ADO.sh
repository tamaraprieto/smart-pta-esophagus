#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu 80G

threads=1

# Load the input and create variables 
INPUT=$1
WORKDIR=$(dirname $INPUT)
PATIENT=$(basename $INPUT | sed 's/_.*//')
OUTDIR=${WORKDIR}/phasing/

my_chromosome=chr${2}

###################################
# CALCULATE DIFFERENT ADO METRICS #
###################################

VCF=${OUTDIR}/${PATIENT}_${my_chromosome}.shapeit5.${my_chromosome}.vcf.gz
OUT=$(echo $VCF | sed 's/.vcf//' | sed 's/.gz//')
somatic_table=/gpfs/commons/groups/landau_lab/ResolveOME/StartDir/${PATIENT}/manual_collate/sequoia/myselectedvars.snvswithsingletons.LR.table
germline_table=${OUT}.table
module load Miniconda3/23.10.0-1
conda activate /gpfs/commons/groups/landau_lab/tprieto/conda/myR4
which Rscript
SCRIPTDIR=/gpfs/commons/groups/landau_lab/ResolveOME/Scripts/
/gpfs/commons/groups/landau_lab/tprieto/conda/myR4/bin/Rscript ${SCRIPTDIR}/rscripts/ADO.simplified.R $somatic_table $germline_table $my_chromosome
conda deactivate

