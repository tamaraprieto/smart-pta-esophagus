#!/usr/bin/env Rscript

library(data.table)
library(ape)
library(dplyr)
library(tidyr)
library(castor)
library(treeio)


####################################################
# Obtain the patient name and phylogeny directory #
####################################################

args = commandArgs(trailingOnly=TRUE)
mydir <- args[1]
treefile <- args[2]
donor_age <- as.numeric(args[3])
patient <- args[4]

print("Arguments: ")
print(mydir)
print(treefile)
print(donor_age)
print(patient)

main_dir <- gsub(paste0(patient,".*"),"",mydir)

#######################################
# load the tree and run the functions #
#######################################

treemut_RDS <- readRDS(file = paste0(mydir, treefile)) 

##########################################################################
# Update the RDS object with descendant celltypes for all internal nodes #
##########################################################################

# annotate the tree with immune/epithelial status  
phenotypes_rna_and_variablemarkers <- readRDS(paste0(main_dir, patient, "/results/phenotypes_rna_plusgeneexpression.rds"))

# Add all the metadata in the tree
treemut_RDS <- left_join(as_tibble(treemut_RDS),
                         phenotypes_rna_and_variablemarkers %>%
                           dplyr::rename(label=`sample name`) %>%
                           dplyr::filter(!is.na(label)), by = 'label') %>%
  treeio::as.treedata()



# get an array of states per node. Then define which ones have both immune and epithelial 
GetStateArrayFromNodes <- function(i,treedata=treemut_RDS){
  mynode <- i
  print(mynode)
  mysubtree <- treeio::tree_subset(treedata, node = mynode,levels_back = 0) %>%
    as_tibble() 
  out <- mysubtree %>%
    dplyr::select(parent,node,celltype) %>%
    dplyr::mutate(type=ifelse(is.na(celltype),"",ifelse(celltype=="EPCAM+","A","B"))) %>%
    dplyr::mutate(alldescendants_celltypes = paste(type, collapse = '')) %>%
    dplyr::select(alldescendants_celltypes) %>% unique() %>%
    dplyr::mutate(label=mynode)
  
  # get the names of the tips for the first split
  
  ancestral_node <- mysubtree %>% dplyr::filter(parent==node) %>% 
    dplyr::select(node) %>% unlist() %>% as.integer()
  two_subnodes <- mysubtree %>% dplyr::filter(parent==ancestral_node & parent!=node) %>%
    dplyr::select(label) %>% unlist()
  if (!grepl("^N",two_subnodes[1])){
    mytip1 <- two_subnodes[1]
  } else {
    mytip1 <- treeio::tree_subset(treedata, node = mynode,levels_back = 0) %>%
      as_tibble() %>% dplyr::filter(!grepl("^N",label)) %>%
      head(1) %>% dplyr::select(label) %>% unlist() %>% as.character()
  }
  if (!grepl("^N",two_subnodes[2])){
    mytip2 <- two_subnodes[2]
  } else {      
    mytip2 <- treeio::tree_subset(treedata, node = two_subnodes[2],levels_back = 0) %>%
      as_tibble() %>% dplyr::filter(!grepl("^N",label)) %>%
      head(1) %>% dplyr::select(label) %>% unlist() %>% as.character()
  }
  mytips=c(mytip1,mytip2)
  
  out <- cbind(out,taxonA=mytip1,taxonB=mytip2,MRCA=paste0(mytips,collapse = "-"))
  
  return(out)
}


# list the ancestral nodes  
mynodes <- treemut_RDS %>% as_tibble %>% dplyr::filter(grepl("^N",label)) %>%
  dplyr::select(label) %>% unlist() %>% as.character()
# obtain a array of descendant states for each ancestral node
my_node_characterization <- do.call("rbind", lapply(mynodes, GetStateArrayFromNodes, treedata=treemut_RDS)) %>%
  ## Nodes sharing both immune and epithelial should have diverged in utero
  dplyr::mutate(celltype_combination=ifelse((grepl("A",alldescendants_celltypes) & grepl("B",alldescendants_celltypes)),"two-cell-types","one-cell-type"))


treemut_RDS <- treemut_RDS %>%
  as_tibble() %>%
  dplyr::left_join(my_node_characterization) %>%
  treeio::as.treedata()
saveRDS(object = treemut_RDS, 
        file = paste0(mydir,"TreeMut.genome.embryolayerannotation.rds")) 


treemut_RDS <- readRDS(file = paste0(mydir,"TreeMut.genome.embryolayerannotation.rds")) 


##########################
# Create the newick file #
##########################

outdir=file.path(mydir, paste0("megacc_input_",patient))
dir.create(outdir)
ape::write.tree(treemut_RDS@phylo, file = paste0(outdir,"/phylogeny.nwk"))


################################# 
# create outgroup file (zeros?) #
#################################

write.table(paste0('"zeros"',"=outgroup"), file = paste0(outdir,"/outgroup.txt"), quote = F, row.names = F, col.names = F)

####################################
# CREATE THE TABLE FOR CALIBRATION #
####################################

# Root age is age+9 months max and min
fertilization_age <- donor_age + 0.75
birth_age <- donor_age

myroot <- castor::find_root(treemut_RDS@phylo)

  treemut_RDS %>% as_tibble() %>%
    dplyr::filter(celltype_combination=="two-cell-types") %>%
    dplyr::mutate(nodetype=ifelse(parent==node, "root","internal")) %>%
    dplyr::filter(!is.na(celltype_combination)) %>%
    dplyr::select(parent,node,branch.length, label, mut_count,nodetype,taxonA, taxonB,MRCA) %>%
    dplyr::filter(node!=myroot) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(calibrationtext=paste0("!MRCA=\'",gsub("_"," ",MRCA),"\' TaxonA=\'",gsub("_"," ",taxonA),"\' TaxonB=\'",gsub("_"," ",taxonB),"\' minTime=",birth_age," maxTime=",fertilization_age,";")) %>%
    dplyr::mutate(calibrationtext=ifelse(parent==myroot,paste0("!MRCA=\'",gsub("_"," ",MRCA),"\' TaxonA=\'",taxonA,"\' TaxonB=\'",taxonB,"\' time=",fertilization_age,";"),calibrationtext)) %>%
    dplyr::select(calibrationtext) %>%
    unique() %>%
    write.table(file = paste0(outdir,"/calibration.txt"), quote = F, row.names = F, col.names = F)
  
#}


