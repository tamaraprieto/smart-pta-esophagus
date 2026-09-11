#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu 20G
#SBATCH --exclude=ne1dc3-[001-026]

module purge
module load bcftools/1.9

VCF_PATH=$1
VCF=$(echo $VCF_PATH | sed 's/.vcf//' | sed 's/.gz//')
tree=$2
workdir=$(dirname $VCF)

echo "VCF: "$VCF
echo "Phylogeny: "$tree
echo "Outputing to: "$workdir

# depth I want to assign to the fake outgroup (called zeros)
depth=15 # we sequence at ~15x
SCRIPTDIR=/gpfs/commons/groups/landau_lab/ResolveOME/Scripts/

# Define the intermediate directory
intermediate_dir="${workdir}/mutmapping/"
maxIt_pval=100

interval=$(echo $SLURM_ARRAY_TASK_ID)

maxlines=$(wc -l ${workdir}/mutmapping/genes.bed | awk '{print $1}')

module load Miniconda3/23.10.0-1
conda activate /gpfs/commons/groups/landau_lab/tprieto/conda/myR4

/gpfs/commons/groups/landau_lab/tprieto/conda/myR4/bin/Rscript ${SCRIPTDIR}/rscripts/mutmapping.interval.R $intermediate_dir $tree $depth $maxIt_pval $maxlines $interval 



