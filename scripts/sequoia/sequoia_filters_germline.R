
#----------------------------------
# Load packages 
#----------------------------------

library("optparse", character.only=T,quietly = T, warn.conflicts = F)

options(stringsAsFactors = F)
cran_packages=c("ggplot2","ape","seqinr","stringr","data.table","tidyr","dplyr","VGAM","MASS","devtools")
bioconductor_packages=c("Rsamtools","GenomicRanges")

for(package in cran_packages){
    library(package, character.only=T,quietly = T, warn.conflicts = F)
}
for(package in bioconductor_packages){
    library(package, character.only=T,quietly = T, warn.conflicts = F)
}
library("treemut",character.only=T,quietly = T, warn.conflicts = F)

library(ape)  
library(castor)  
library(TreeTools)

#----------------------------------
# Input options
#----------------------------------
option_list = list(
  make_option(c("-v", "--input_nv"), action="store", default=NULL, type='character', help="Input NV matrix (rows are variants, columns are samples)"),
  make_option(c("-r", "--input_nr"), action="store", default=NULL, type='character', help="Input NR matrix (rows are variants, columns are samples)"),
  make_option(c("-o", "--output_dir"), action="store", default="", type='character', help="Output directory for files"),
  make_option(c("-n", "--ncores"), action="store", default=1, type='numeric', help="Number of cores to use for the beta-binomial step"),
  make_option(c("--sex"), action="store", default='female', type='character', help="Patient sex: male or female as indicated by sex chromosomes"),
  make_option(c("--snv_rho"), action="store", default=0.1, type='numeric', help="Rho value threshold for SNVs"),
  make_option(c("--indel_rho"), action="store", default=0.15, type='numeric', help="Rho value threshold for indels"),
  make_option(c("--min_cov"), action="store", default=10, type='numeric', help="Lower threshold for mean coverage across variant site"),
  make_option(c("--max_cov"), action="store", default=500, type='numeric', help="Upper threshold for mean coverage across variant site"),
  make_option(c("-x","--exclude_samples"), action="store", default=NULL, type='character', help="Option to manually exclude certain samples from the analysis, separate with a comma"),
  make_option(c("--cnv_samples"), action="store", default=NULL, type='character', help="Samples with CNVs, exclude from germline/depth-based filtering, separate with a comma"),
  make_option(c("--pval_thres"), action="store", default=0.05, type='numeric', help="P-value to reject or accept the mapping of a mutations to a branch of a random topology"),
  make_option(c("--maxIt_pval"), action="store", default=100, type='numeric', help="Number of iterations for the treemut::assign_to_tree function"),
  make_option(c("--germline_cutoff"), action="store", default=-5, type='numeric', help="Log10 of germline qval cutoff")
)
opt = parse_args(OptionParser(option_list=option_list, add_help_option=T))

print(opt)

nv_path=opt$input_nv
nr_path=opt$input_nr
output_dir=opt$output_dir
ncores=opt$ncores
gender=opt$sex
snv_rho=opt$snv_rho
indel_rho=opt$indel_rho
min_cov=opt$min_cov
max_cov=opt$max_cov
pval_thres=opt$pval_thres
maxIt_pval=opt$maxIt_pval
if(is.null(opt$exclude_samples)) {samples_exclude=NULL} else {samples_exclude=unlist(strsplit(x=opt$exclude_samples,split = ","))}
if(is.null(opt$cnv_samples)) {samples_with_CNVs=NULL} else {samples_with_CNVs=unlist(strsplit(x=opt$cnv_samples,split = ","))}
germline_cutoff=opt$germline_cutoff




#----------------------------------
# Functions
#----------------------------------

exact.binomial=function(gender,NV,NR,adjust=F){
  # Function to filter out germline variants based on unmatched
  # variant calls of multiple samples from same individual (aggregate coverage
  # ideally >150 or so, but will work with less). NV is matrix of reads supporting 
  # variants and NR the matrix with total depth (samples as columns, mutations rows, 
  # with rownames as chr_pos_ref_alt or equivalent). Returns a logical vector, 
  # TRUE if mutation is likely to be germline.
  
  XY_chromosomal = grepl("X|Y",rownames(NR))
  autosomal = !XY_chromosomal
  NV <- apply(NV, 2, as.numeric) 
  NR <- apply(NR, 2, as.numeric) 
  
  if(gender=="female"){
    NV_vec = rowSums(NV, na.rm = T) 
    NR_vec = rowSums(NR, na.rm = T) 
    pval = rep(1,length(NV_vec))
    for (n in 1:length(NV_vec)){
      if(NR_vec[n]>0){
        pval[n] = binom.test(x=NV_vec[n],
                             n=NR_vec[n],
                             p=0.5,alt='less')$p.value
      }
    }
  }
  # For male, split test in autosomal and XY chromosomal part
  if(gender=="male"){
    pval=rep(1,nrow(NV))
    # for autosomatic positions
    NV_vec = rowSums(NV, na.rm = T)[autosomal] 
    NR_vec = rowSums(NR, na.rm = T)[autosomal] 
    pval_auto = rep(1,sum(autosomal))
    pval_XY = rep(1,sum(XY_chromosomal))
    if (length(NV_vec)>0){
    for (n in 1:sum(autosomal)){
      if(NR_vec[n]>0){
        pval_auto[n] = binom.test(x=NV_vec[n],
                                  n=NR_vec[n],
                                  p=0.5,alt='less')$p.value
      }
    }}else{pval_auto=NULL}
    
    # for X and Y positions
    NV_vec = rowSums(NV, na.rm = T)[XY_chromosomal] 
    NR_vec = rowSums(NR, na.rm = T)[XY_chromosomal] 
    if (length(NV_vec)>0){
      for (n in 1:sum(XY_chromosomal)){
        if(NR_vec[n]>0){
          pval_XY[n] = binom.test(x=NV_vec[n],
                                  n=NR_vec[n],
                                  p=0.95,alt='less')$p.value
        }      
    }} else {pval_XY=NULL}
    
    pval[autosomal]=pval_auto
    pval[XY_chromosomal]=pval_XY
  }
  # multiple correction (dependent on number of variants)
  if(adjust){
    pval = p.adjust(pval,method="BH")
  }
  return(pval)
}

estimateRho_gridml = function(NV_vec,NR_vec) {
  # Function to estimate maximum likelihood value of rho for beta-binomial
  rhovec = 10^seq(-6,-0.05,by=0.05) # rho will be bounded within 1e-6 and 0.89
  mu=sum(NV_vec)/sum(NR_vec)
  ll = sapply(rhovec, function(rhoj) sum(dbetabinom(x=NV_vec, size=NR_vec, rho=rhoj, prob=mu, log=T), na.rm = T))
  
  return(rhovec[ll==max(ll)][1])
}

beta.binom.filter = function(NR,NV){
  # Function to apply beta-binomial filter for artefacts. Works best on sets of
  # clonal samples (ideally >10 or so). As before, takes NV and NR as input. 
  # Optionally calculates pvalue of likelihood beta-binomial with estimated rho
  # fits better than binomial. This was supposed to protect against low-depth variants,
  # but use with caution. Returns logical vector with good variants = TRUE
  
  rho_est = pval = rep(NA,nrow(NR))
  # calculate rho per site
  for (k in 1:nrow(NR)){
    rho_est[k]=estimateRho_gridml(NV_vec = as.numeric(NV[k,]),
                                  NR_vec=as.numeric(NR[k,]))
  }
  return(rho_est)
}



#----------------------------------
# Read in data
#----------------------------------

print("Reading in data...")

  if(!is.null(nr_path)&!is.null(nv_path)){
    NR = fread(nr_path,data.table=F)
    rownames(NR)=NR[,1]
    NR=NR[,-1]
    NR=NR[,!colnames(NR)%in%samples_exclude]
    NV = fread(nv_path,data.table=F)
    rownames(NV)=NV[,1]
    NV=NV[,-1]
    NV=NV[,!colnames(NV)%in%samples_exclude]
    samples=colnames(NV)
    Muts=rownames(NV)
  }else{
    print("Please provide either NV and NR files or a path to CGPVaf output")
    break
  }

# Determine if analyzing indels, snvs or both
Muts_coord=matrix(ncol=4,unlist(strsplit(Muts,split="_")),byrow = T)
if(all(nchar(Muts_coord[,3])==1&nchar(Muts_coord[,4]))==1){
  mut_id="snv"
} else{
  if(all(nchar(Muts_coord[,3])>1|nchar(Muts_coord[,4])>1)){
    mut_id="indel"
  } else{
    mut_id="both"
  }
}
print(paste0("Mutations in data (snv/indel):", mut_id))

# Count how many mutations are in autosomes and/or sex chromosomes
XY_chromosomal = grepl("X|Y",Muts)
print(paste0("Number of variants in sex chromosome: ",as.character(sum(XY_chromosomal))))
autosomal = !XY_chromosomal
print(paste0("Number of variants in autosomes: ",as.character(sum(autosomal))))
print(paste0("Gender: ", gender))

noCNVs=!samples%in%samples_with_CNVs

#----------------------------------
# Create the file with filter annotations
#----------------------------------
if(output_dir!="") system(paste0("mkdir -p ",output_dir))
print("Start creating a file with filter annotations...")

filter_df=as.data.frame(matrix(ncol=4,unlist(strsplit(rownames(NV),split="_")),byrow = T)) %>%
  dplyr::mutate(start=as.character(as.integer(as.numeric(V2)-1))) %>%
  dplyr::select(V1,start,V2,V3,V4)
rownames(filter_df)=rownames(NV)
colnames(filter_df)=c("chrom","start","end","Ref","Alt")

#----------------------------------
# Calculate mean depth per cell and add it to the table
#----------------------------------
filter_df$Mean_Depth=rowMeans(NR[,noCNVs], na.rm = T)
# Add tag to sites with high and low depth across samples taking into account chromosomal sex 
if(gender=='male'){
  filter_df$Depth_filter = (rowMeans(NR[,noCNVs])>min_cov&rowMeans(NR[,noCNVs])<max_cov&autosomal)|
    (rowMeans(NR[,noCNVs])>(min_cov/2)&rowMeans(NR[,noCNVs])<(max_cov/2)&XY_chromosomal)
}else{
  filter_df$Depth_filter = rowMeans(NR)>min_cov&rowMeans(NR)<max_cov
}

#----------------------------------
# Add pvalue of germline to the table as well as well as another column with germline vs non-germline based on threshold
#----------------------------------
germline_qval=exact.binomial(gender=gender,NV=NV[,noCNVs],NR=NR[,noCNVs],adjust=F) 
filter_df$Germline_qval=germline_qval
filter_df$Germline=as.numeric(log10(germline_qval)<germline_cutoff)


#----------------------------------
# Calculate beta-binomial dispersion
#----------------------------------


  print("Calculating beta-binomial...")
  
  if(ncores>1){
    print("...in parallel")
    rho_est=unlist(mclapply(1:nrow(NR),function(x){
      estimateRho_gridml(NR_vec=as.numeric(NR[x,]),NV_vec=as.numeric(NV[x,]))
    },mc.cores=ncores))
  }else{
    rho_est=beta.binom.filter(NR=NR, NV=NV)
  }
  
# add rho value to table with filter annotations
# add a 0 or 1 depending on threshold
  filter_df$Rho=rho_est
  if(mut_id=="snv")filter_df$Beta_binomial=as.numeric(rho_est>snv_rho&!is.na(rho_est))
  if(mut_id=="indel")filter_df$Beta_binomial=as.numeric(rho_est>indel_rho&!is.na(rho_est))
  if(mut_id=="both"){
    is.indel=nchar(filter_df$Ref)>1|nchar(filter_df$Alt)>1
    filter_df$Beta_binomial=as.numeric(((rho_est>indel_rho&is.indel)|(rho_est>snv_rho&!is.indel))&!is.na(rho_est))
  }

#----------------------------------------  
#  Map mutations to a random topology 
#----------------------------------------  
   
mytree <- ape::rcoal(n=length(colnames(NR)), tip.label = c(colnames(NR)))
myroot <- castor::find_root(mytree)
# add a tip but also create an internal node where we can map germline mutations
# the line that they suggest in the code to add a zeros tip does not add an internal branch
#mytree=ape::bind.tree(mytree,ape::read.tree(text="(zeros:0);"))
mytree <- TreeTools::AddTip(tree = mytree, where = myroot, label = "zeros")
myroot <- castor::find_root(mytree)

# Create node names (very important for downstream processing)
tree <- ape::makeNodeLabel(mytree, method = "number", prefix = "N")  
 
# FUNCTION TO MAP THE MUTATIONS TO THE TREE 
# MUTATIONS ARE ONLY ASSIGNED ONCE (hardassigned after pval calculations), get NEW LENGTHS
# The code calls another function called "assign_to_df" which I have to download from the github if I need to debug
# It is working unless there are some NAs in the depth.
# Try different mtr and depth variants if it doesn't work!!!!
res=treemut::assign_to_tree(tree,as.matrix(NV),as.matrix(NR),error_rate=0, maxits = maxIt_pval)

# Code to map the mutations to the tree and visualize it
table_tree_withmappinginfo <- cbind(res$tree$edge, original.branch.length=tree$edge.length,
                    res$df$df) %>%
  dplyr::full_join(cbind(res$summary,mutid=rownames(as.matrix(NV))) %>%
                     dplyr::mutate(all=ifelse(pval>pval_thres,1,1)) %>%
                     # pval: a heuristic pvalue assessing the hypothesis that the mutation is consistent with the provided tree topology
                     # null hypothesis: pval>pval_threshold -> consistent
                     # reject hypothesis: pval<pval_threshold -> inconsistent
                     dplyr::mutate(signif=ifelse(as.numeric(pval)>=pval_thres | as.numeric(p_else_where)==0,1,0)) %>%
                     dplyr::rename(id=edge_ml), by="id") %>%
  # if pval close to 1 then create another list of mutations
  dplyr::mutate(mutid2=ifelse(signif==1,mutid,NA)) %>%
  dplyr::rename(label2=label) %>%
  dplyr::rename(parent=`1`) %>%
  dplyr::rename(node=`2`)

# germline snps will be mapped to the outgroup branch
profile_order <- res$df$samples 
germline_branch <- paste0(c(rep(1, length(profile_order)-1),0), collapse = "")

# Load mutation matrices 
# Create a bed file to annotate the original file
mutation_mapping_classification <- table_tree_withmappinginfo %>%
  dplyr::select(mutid, profile.x, label2, edge_length , expected_edge_length, original.branch.length, mut_count, 
                pval, p_else_where, signif) %>%
  dplyr::filter(!is.na(mutid)) %>%
  dplyr::group_by(label2) %>%
  dplyr::mutate(mut_mapped_2branch=n()) %>%
  as.data.frame() %>%
  dplyr::mutate(mapping_filter=ifelse(signif==1, "congruent","inconsistent"))  %>%
  dplyr::mutate(branch_filter=ifelse(profile.x==germline_branch,"germline","somatic")) %>%
  dplyr::mutate(branch=stringr::str_sub(label2,1,-2)) %>%
  dplyr::mutate(chrom=gsub(":.*","",gsub(";.*","",mutid))) %>%
  dplyr::mutate(end=as.numeric(gsub(".*:","",gsub(".*;","",mutid)))) %>%
  dplyr::mutate(start=end-1) %>%
  dplyr::select(chrom, start, end, branch_filter, mapping_filter, branch, pval, p_else_where, signif) %>%
  dplyr::arrange(chrom,end)
colnames(mutation_mapping_classification) <- c("chrom","start","end","branch_filter","mappingrandomtopo_filter","branch","pvaluerandomtopo_branch","pvaluerandomtopo_otherbranches","mapped")

print(paste0("filter_df dims: ",dim(filter_df)))
print(paste0("mutmapping dims:",dim(mutation_mapping_classification)))

filter_df$branch_filter <- mutation_mapping_classification$branch_filter
filter_df$mapping_filter <- mutation_mapping_classification$mappingrandomtopo_filter
filter_df$pvalue_branch <- mutation_mapping_classification$pvaluerandomtopo_branch
filter_df$pvalue_otherbranches <- mutation_mapping_classification$pvaluerandomtopo_otherbranches

  
#-----------------------------------------
# Write down file for annotating the vcf  
#-----------------------------------------  
  
interval <- gsub(".txt","",gsub(".*_","",nv_path))
write.table(filter_df,paste0(output_dir,"Filter_annotations",interval,".txt"), quote = F, row.names = F, sep = "\t")



