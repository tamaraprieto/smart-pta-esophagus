#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu 1G

module purge
module load bcftools

PATIENT=$1
WORKDIR=/gpfs/commons/groups/landau_lab/ResolveOME/StartDir/${PATIENT}/panel/preprocessing/coverage/
SAMPLELIST=${WORKDIR}Bams.txt
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $SAMPLELIST)
CELL="$(basename -- ${SAMPLE} | sed 's/.dedup//' | sed 's/.recal//' | sed 's/.bam//')"
BAM=${WORKDIR}${CELL}.mygenes.nondup.bam
BED=/gpfs/commons/groups/landau_lab/ResolveOME/Resources/gene_panel/implemented/Targets-XGEN.87EF358171F54BE7A9133DE926D12A03.chrprefix.merged.allpositionsallmutations.annovar.bed
REF=/gpfs/commons/groups/landau_lab/ResolveOME/Resources/references/hg38/Homo_sapiens_assembly38.fasta
SNV=/gpfs/commons/groups/landau_lab/ResolveOME/Resources/gene_panel/implemented/Targets-XGEN.87EF358171F54BE7A9133DE926D12A03.chrprefix.merged.allpositionsallmutations.avinput

echo "BAM: "${BAM}
OUTDIR="${WORKDIR}/mutmapping/"
mkdir -p "$OUTDIR"

#########################################
# Create the alt and total count tables #
#########################################

module load Miniconda3/23.10.0-1
source activate /gpfs/commons/groups/landau_lab/tprieto/conda/myR4
echo "Run the python script"
# RUN PYTHON SCRIPT TO CREATE DEPTH TABLES
/nfs/sw/easybuild/software/Miniconda3/23.10.0-1/bin/python /gpfs/commons/groups/landau_lab/ResolveOME/Scripts/panel_vc/BaseCountsBAMBedFile.py \
	--bam $BAM \
	--snv $SNV \
	--ref $REF \
	--outdir $OUTDIR
source deactivate

echo "Matrix total: "${OUTDIR}${WORKDIR}${CELL}.mygenes.nondup_total_depth.tsv
echo "Matrix alt: "${OUTDIR}${WORKDIR}${CELL}.mygenes.nondup_total_depth.tsv
cp /gpfs/commons/groups/landau_lab/ResolveOME/Resources/gene_panel/implemented/Targets-XGEN.87EF358171F54BE7A9133DE926D12A03.chrprefix.merged.allpositionsallmutations.annovar.bed ${OUTDIR}genes.bed
echo "Bed file: "${OUTDIR}genes.bed
