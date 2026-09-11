#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 4
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu 4G


mutmappedtree=$1
patient=$2
workdir=$(dirname $mutmappedtree)
treeRDS=$(basename $mutmappedtree)

SCRIPTDIR=/gpfs/commons/groups/landau_lab/ResolveOME/Scripts/
age=$(awk -v donor=$patient '{if ($1==donor){print $2}}' /gpfs/commons/groups/landau_lab/ResolveOME/StartDir/donor.metadata.txt)

module purge
module load Miniconda3/23.10.0-1
conda activate /gpfs/commons/groups/landau_lab/tprieto/conda/myR4
Rscript ${SCRIPTDIR}/rscripts/megacc.input.R ${workdir}"/" $treeRDS $age $patient
