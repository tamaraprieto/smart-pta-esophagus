#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu 20G

module purge
module load bcftools/1.9

tree=$1
patient=$(echo $1 | sed 's/.*CellPhy.//' | sed 's/.newick//')
workdir=/gpfs/commons/groups/landau_lab/ResolveOME/StartDir/${patient}/panel/preprocessing/coverage/

echo "Phylogeny: "$tree
echo "Outputing to: "$workdir

depth=15
SCRIPTDIR=/gpfs/commons/groups/landau_lab/ResolveOME/Scripts/
# Define the intermediate directory
intermediate_dir="${workdir}/mutmapping/"
maxIt_pval=100
interval=$(echo $SLURM_ARRAY_TASK_ID)
maxlines=$(wc -l ${workdir}/mutmapping/genes_filtered.bed | awk '{print $1}')

module purge
source deactivate
conda deactivate
module load Miniconda3/23.10.0-1
source activate /gpfs/commons/groups/landau_lab/tprieto/conda/myR4

echo "OK"

Rscript ${SCRIPTDIR}/panel_vc/mutmapping.interval.R $intermediate_dir $tree $depth $maxIt_pval $maxlines $interval 
source deactivate



