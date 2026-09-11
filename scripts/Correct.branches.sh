#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 4
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu 4G


mutmappedtree=$1
patient=$2
# the recall object contains at least three columns called cell, recall and donor. 
# The column names have to be exactly the same. Example of values for these three columns:  eso29_p1A4 0.656 eso29
recall_object=$3
workdir=$(dirname $mutmappedtree)
treeRDS=$(basename $mutmappedtree)

SCRIPTDIR=/gpfs/commons/groups/landau_lab/ResolveOME/Scripts/

module purge
source activate /gpfs/commons/groups/landau_lab/tprieto/conda/myR4
Rscript ${SCRIPTDIR}/rscripts/megacc.input.recallcorrection.all.R ${workdir}"/" $treeRDS $patient $recall_object
