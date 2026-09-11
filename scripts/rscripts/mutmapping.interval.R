#!/usr/bin/env Rscript
library(treemut)
library(ape)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggtree)
library(data.table)
library(vroom)
library(castor)
library(TreeTools)
library(treeio)

# treemut was designed to enable users to map clonal somatic mutations to a specified tree.  The
# method requires:
#   - phylogeny
#   - the mutant read counts and depth for each sample/loci in the form of matrices
#   - error rate per single cell

####################################################
# Obtain the patient name and phylogeny directory #
####################################################

args = commandArgs(trailingOnly=TRUE)
mydir <- args[1]
treefile <- args[2]
seqdepth <- args[3]
maxIt_pval <- as.integer(args[4]) 
max_line <- as.integer(args[5])
print(max_line)
interval <- as.integer(args[6])
interval_size <- 10000
pval_thres <- 0.05  # I am also selecting all sites for which the hypothesis of being consistent is rejected (pval<0.05) but with 0 probability anywhere else in the phylogeny (p_else_where=0)
print(treefile)

#######################
# Load the dataframes #
#######################

print("Get start and end for the interval...")

# create an interval start and end coordinates
# Function to get start and end coordinates given an interval number
get_coordinates <- function(interval, maxlines=max_line) {
  start_coordinate <- (interval - 1) * interval_size + 1
  end_coordinate <- interval * interval_size
  if (start_coordinate>maxlines){
    stop(paste0("there are not ",interval," intervals of size ", interval_size," if the total length is ",max_line))
  }
  if (end_coordinate>maxlines){
    end_coordinate=maxlines
  }
  return(c(start_coordinate, end_coordinate))
}
coordinates <- get_coordinates(interval)
start_coordinate <- coordinates[1]
end_coordinate <- coordinates[2]

print("Load the files...")

nr_path=paste0(mydir,"total_depth.txt")
nv_path=paste0(mydir,"alt_depth.txt")
command=paste0('awk -v start=',start_coordinate,' -v end=',end_coordinate,' \'{if (NR >= start && NR <= end){print $0}}\' ')
depth <- data.table::fread(paste0(command, nr_path),data.table=F, check.names = F)
colnames(depth) <- data.table::fread(nr_path,data.table=F ,nrows = 1, header = F) %>% unlist() %>% as.vector()
mtr <- data.table::fread(paste0(command, nv_path),data.table=F, check.names = F)
colnames(mtr) <- data.table::fread(nv_path,data.table=F ,nrows = 1, header = F) %>% unlist() %>% as.vector()
# Bed of the genomic mutations
bed <- data.table::fread(paste0(command,mydir,'genes.bed'), col.names = c("CHROM","START","END","Type","Gene","Ogaini","ExonFun")) %>%
  # add a counter to each gene (same gene with several mutations)
  dplyr::group_by(Gene,Type) %>% mutate(counter = dplyr::row_number(Gene)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Gene=ifelse(Type=="exonic",paste0(Gene,"-exon",ExonFun),Gene)) %>%
  dplyr::mutate(Gene=ifelse(!is.na(Ogaini),paste0(Gene,"-ogacor",Ogaini),Gene)) %>%
  dplyr::mutate(pos=paste0(CHROM,":",END)) %>%
  dplyr::mutate(Gene=gsub(".*:","",Gene)) %>%
  dplyr::select(pos,Type,Gene) %>% unique()


# change mutation names
# add gene name at the beginning of the gene
depth <- depth %>%
  dplyr::mutate(pos=gsub("_.*","",CHROM_POS_REF_ALT)) %>%
  dplyr::left_join(bed) %>%
  dplyr::mutate(CHROM_POS_REF_ALT=ifelse(!is.na(Gene),paste0(Gene,";",CHROM_POS_REF_ALT),CHROM_POS_REF_ALT)) %>%
  dplyr::select(-pos,-Type,-Gene)
mtr <- mtr %>%
  dplyr::mutate(pos=gsub("_.*","",CHROM_POS_REF_ALT)) %>%
  dplyr::left_join(bed) %>%
  dplyr::mutate(CHROM_POS_REF_ALT=ifelse(!is.na(Gene),paste0(Gene,";",CHROM_POS_REF_ALT),CHROM_POS_REF_ALT)) %>%
  dplyr::select(-pos,-Type,-Gene)


rownames(depth)=depth[,1]
depth=depth[,-1]
rownames(mtr)=mtr[,1]
mtr=mtr[,-1]
mtr2<-mtr
depth2<-depth

# Transform into matrices (avoid issues later on)
mtr <- as.matrix(sapply(mtr, as.numeric))
mtr[is.na(mtr)] <- 0
rownames(mtr) <- rownames(mtr2)
depth <- as.matrix(sapply(depth, as.numeric))
rownames(depth) <- rownames(depth2)

#######################################
# load the tree and run the functions #
#######################################

print("Load tree...")
mytree <- ape::read.tree(treefile)

# Drop tips not in VCF
badqual_cells <- setdiff(mytree$tip.label, colnames(mtr))
mytree<- treeio::drop.tip(mytree, tip = badqual_cells)

myroot <- castor::find_root(mytree)
# add a tip but also create an internal node where we can map germline mutations
mytree <- TreeTools::AddTip(tree = mytree, where = myroot, label = "zeros")
myroot <- castor::find_root(mytree)
# Create node names (very important for downstream processing)
tree <- ape::makeNodeLabel(mytree, method = "number", prefix = "N")

#################################
# take only samples in the tree #
#################################

mtr <- mtr[,tree$tip.label[!grepl("zeros$",tree$tip.label, fixed = F)]]
depth <- depth[,tree$tip.label[!grepl("zeros$",tree$tip.label, fixed = F)]]


print("> Mapping the mutations...")

# FUNCTION TO MAP THE MUTATIONS TO THE TREE
# MUTATIONS ARE ONLY ASSIGNED ONCE (hardassigned after pval calculations), get NEW LENGTHS
res=treemut::assign_to_tree(tree,mtr,depth,error_rate=0, maxits = maxIt_pval)

# Code to map the mutations to the tree and visualize it
mydataRaw1 <- cbind(res$tree$edge, original.branch.length=tree$edge.length,
                    res$df$df) %>%
  dplyr::full_join(cbind(res$summary,mutid=rownames(mtr)) %>%
                     dplyr::mutate(all=ifelse(pval>pval_thres,1,1)) %>%
                     # pval: a heuristic pvalue assessing the hypothesis that the mutation is consistent with the provided tree topology
                     # null hypothesis: pval>0.05 -> consistent
                     # reject hypothesis: pval<0.05 -> inconsistent
                     dplyr::mutate(signif=ifelse(as.numeric(pval)>=pval_thres | as.numeric(p_else_where)==0,1,0)) %>%
                     dplyr::rename(id=edge_ml), by="id") %>%
  # if pval close to 1 then create another list of mutations
  dplyr::mutate(mutid2=ifelse(signif==1,mutid,NA)) %>%
  dplyr::rename(label2=label) %>%
  dplyr::rename(parent=`1`) %>%
  dplyr::rename(node=`2`)


############################################
# Create a bed file for annotating the vcf #
############################################

# germline snps will be mapped to the outgroup branch
profile_order <- res$df$samples
germline_branch <- paste0(c(rep(1, length(profile_order)-1),0), collapse = "")

# Load mutation matrices
# Create a bed file to annotate the original file
mutation_mapping_classification <- mydataRaw1 %>%
  dplyr::select(mutid, profile.x, label2, edge_length , expected_edge_length, original.branch.length, mut_count,
                pval, p_else_where, signif) %>%
  dplyr::filter(!is.na(mutid)) %>%
  dplyr::group_by(label2) %>%
  dplyr::mutate(mut_mapped_2branch=dplyr::n()) %>%
  as.data.frame() %>%
  dplyr::mutate(mapping_filter=ifelse(signif==1, "congruent","inconsistent"))  %>%
  dplyr::mutate(branch_filter=ifelse(profile.x==germline_branch,"germline","somatic")) %>%
  dplyr::mutate(branch=ifelse(nchar(label2)==nchar(germline_branch),stringr::str_sub(label2,1,-2),label2)) %>%
  dplyr::mutate(chrom=gsub(":.*","",gsub(".*;","",mutid))) %>%
  dplyr::mutate(end=as.numeric(gsub("_.*","",gsub(".*:","",gsub(".*;","",mutid))))) %>%
  dplyr::mutate(start=as.character(as.integer(end-1))) %>%
  dplyr::select(chrom, start, end, branch_filter, mapping_filter, branch, pval, p_else_where, signif) %>%
  dplyr::arrange(chrom,end)

colnames(mutation_mapping_classification) <- c("chrom","start","end","outbranch_filter","mapping_filter","branch","pvalue_branch","pvalue_otherbranches","mapped")

mydataRaw <- mydataRaw1 %>%
  dplyr::left_join(res$tree %>% as_tibble(), by=c("node","parent"))

vaf=mtr/depth

  sites_to_plot <- unique(mydataRaw$mutid2[!is.na(mydataRaw$mutid2) &
                                             !grepl("^chr",mydataRaw$mutid2) #&
                                           ])



mydata <- mydataRaw %>%
  # Create a list of mutations to map to each node
  dplyr::group_by(parent,node,label) %>%
  dplyr::mutate(mutList = paste0(coalesce(as.character(mutid2),""), collapse = ",")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(mutList=gsub(",$","",gsub("^,","",gsub("(,)\\1+", "\\1", mutList)))) %>%
  # Create a third list with a few mutations
  dplyr::mutate(mutid3=ifelse(mutid2 %in% sites_to_plot,mutid2,NA)) %>%
  dplyr::group_by(parent,node,label) %>%
  dplyr::mutate(mutList2 = paste0(coalesce(as.character(mutid3),""), collapse = ",")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(mutList2=gsub(",$","",gsub("^,","",gsub("(,)\\1+", "\\1", mutList2)))) %>%
  dplyr::mutate(branch.length.mapping=as.numeric(stringr::str_count(string = mutList, pattern = ">")))

mydatasummarizedperbranch <- mydata %>%
  # remove non-unique columns to obtain branches alone
  dplyr::select(-pval,-p_else_where,-mutid,-mutid2,-mutid3,-signif,-all) %>% unique()

newtree <- res$tree %>%
  tibble::as_tibble() %>%
  dplyr::left_join(mydatasummarizedperbranch, by=c("parent","node","branch.length","label")) %>%
  dplyr::mutate(branch.length=original.branch.length) %>%
  treeio::as.treedata()
tree_withzeros <- newtree


#################
## Save objects #
#################

newdir=paste0(mydir,"maxIt.",as.character(maxIt_pval))
dir.create(newdir)
newdir=paste0(mydir,"maxIt.",as.character(maxIt_pval),"/intervals/")
dir.create(newdir)
mydir <- newdir

saveRDS(object = res, file = paste0(mydir, "MutationMappingFunctionOutput.",interval,".rds"))
mutation_mapping_classification %>%
  write.table(file = paste0(mydir, "MutationMappingAnnotations.",interval,".bed"), sep = "\t", quote = F,
              append = F, row.names = F, col.names = T)
if (interval==1){
  write.table(profile_order[1:(length(profile_order)-1)], file = paste0(mydir, "/../PhylogenyMappingTipOrder.txt"), quote = F, col.names = F, row.names = F, eol = " ")
  saveRDS(object = profile_order, file = paste0(mydir, "/../MappingProfileOrder.rds"))
}
saveRDS(object = tree_withzeros, file = paste0(mydir, "TreeMutWithZeros.genome.",interval,".rds"))
