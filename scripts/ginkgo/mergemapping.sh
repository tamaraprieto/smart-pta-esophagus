#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH -t 150:00:00
#SBATCH --mem 50G


# ==============================================================================
# == Launch analysis
# ==============================================================================

module purge
module load php/8.1.0
module load miniconda3/23.10.0-1
source activate /gpfs/commons/groups/landau_lab/tprieto/conda/myR4

# ------------------------------------------------------------------------------
# -- Variables
# ------------------------------------------------------------------------------

home=/gpfs/commons/groups/landau_lab/tprieto/apps/ginkgo
dir=${home}/uploads/${1}
source ${dir}/config
distMet=$distMeth
touch $dir/index.html

inFile=list
statFile=status.xml

if [ ${chosen_genome} == "hg19" ];
then
  genome=${home}/genomes/${chosen_genome}/original
else
  genome=${home}/genomes/${chosen_genome} # use hg38
fi


# ------------------------------------------------------------------------------
# -- Error Check and Reformat User Files
# ------------------------------------------------------------------------------

if [ "$f" == "0" ]; then
  touch ${dir}/ploidyDummy.txt
  facs=ploidyDummy.txt
else 
  # In case upload file with \r instead of \n (Mac, Windows)
  tr '\r' '\n' < ${dir}/${facs} > ${dir}/quickTemp
  mv ${dir}/quickTemp ${dir}/${facs}
  # 
  sed "s/.bed//g" ${dir}/${facs} | sort -k1,1 | awk '{print $1"\t"$2}' > ${dir}/quickTemp 
  mv ${dir}/quickTemp ${dir}/${facs}
fi

if [ "$init" == "1" ];
then

  # Concatenate binned reads to central file  
  paste ${dir}/*_mapped > ${dir}/data
  rm -f ${dir}/*_mapped ${dir}/*_binned

fi

# ------------------------------------------------------------------------------
# -- Map User Provided Reference/Segmentation Sample
# ------------------------------------------------------------------------------

if [ "$segMeth" == "2" ]; then
    ${home}/scripts/binUnsorted ${genome}/${binMeth} `wc -l < ${genome}/${binMeth}` ${dir}/${ref} Reference ${dir}/${ref}_mapped
else
    ref=refDummy.bed
    touch ${dir}/${ref}_mapped
fi

# ------------------------------------------------------------------------------
# -- Run Mapped Data Through Primary Pipeline
# ------------------------------------------------------------------------------

if [ "$process" == "1" ]; then
  echo "Launching process.R $genome $dir $statFile data $segMeth $binMeth $clustMeth $distMet $color ${ref}_mapped $f $facs $sex $rmbadbins"
  ${home}/scripts/process.R $genome $dir $statFile data $segMeth $binMeth $clustMeth $distMet $color ${ref}_mapped $f $facs $sex $rmbadbins
fi

# ------------------------------------------------------------------------------
# -- Recreate Clusters/Heat Maps (With New Parameters)
# ------------------------------------------------------------------------------

if [ "$fix" == "1" ]; then
  echo "Launching reclust.R $genome $dir $statFile $binMeth $clustMeth $distMet $f $facs $sex"
  ${home}/scripts/reclust.R $genome $dir $statFile $binMeth $clustMeth $distMet $f $facs $sex
fi

# ------------------------------------------------------------------------------
# -- Create CNV profiles
# ------------------------------------------------------------------------------

nbCols=$(awk '{ print NF; exit; }' $dir/SegCopy)
for (( i=1; i<=$nbCols; i++ ));
do
  currCell=$(cut -f$i $dir/SegCopy | head -n 1 | tr -d '"')
  if [ "$currCell" == "" ]; then
    continue;
  fi
  cut -f$i $dir/SegCopy | tail -n+2 | awk '{if(NR==1) print "1,"$1; else print NR","prev"\n"NR","$1;prev=$1; }' > $dir/$currCell.cnv
done

# ------------------------------------------------------------------------------
# -- Call CNVs
# ------------------------------------------------------------------------------

echo "Launching ${home}/scripts/CNVcaller ${dir}/SegCopy ${dir}/CNV1 ${dir}/CNV2"
${home}/scripts/CNVcaller ${dir}/SegCopy ${dir}/CNV1 ${dir}/CNV2

# ------------------------------------------------------------------------------
# -- Email notification of completion
# ------------------------------------------------------------------------------

if [ "$email" != "" ]; then
	echo -e "Your analysis on Ginkgo is complete! Check out your results at $permalink" | mail -s "Your Analysis Results" $email -- -F "Ginkgo"
fi
