#!/usr/bin/env Rscript

library(data.table)
library(ape)
library(dplyr)
library(tidyr)
library(castor)
library(treeio)
source("/gpfs/commons/groups/landau_lab/ResolveOME/Scripts/rscripts/treefit_package/code.R")
source("/gpfs/commons/groups/landau_lab/ResolveOME/Scripts/rscripts/treefit_package/cna.R")
source("/gpfs/commons/groups/landau_lab/ResolveOME/Scripts/rscripts/treefit_package/plot_tree_annots_extra.R")
PD=readRDS("/gpfs/commons/groups/landau_lab/ResolveOME/Resources/calibration/PD6629.RDS")


####################################################
# Obtain the patient name and phylogeny directory #
####################################################

args = commandArgs(trailingOnly=TRUE)
#patient <- args[1]
mydir <- args[1]
treefile <- args[2]
donor_age <- as.numeric(args[3])
NITER=as.numeric(args[4])

print("Arguments: ")
print(mydir)
print(treefile)
print(donor_age)
print(NITER)

#######################################
# load the tree and run the functions #
#######################################

treemut_RDS <- readRDS(file = paste0(mydir, treefile)) 

new_data_table <- treemut_RDS@data
new_data_table <-  new_data_table %>%
  dplyr::select(-mutList, -mutList2)
empty_row <- which(is.na(new_data_table$original.branch.length))
new_data_table <- new_data_table[-c(empty_row),] # remove row with NA values
new_data_table$node[which(new_data_table$node >= empty_row)] <- new_data_table$node[which(new_data_table$node > empty_row)] - 1 # correcting node id's after removing zeros
treemut_RDS@data <- new_data_table

# create agedf from PD
# PD object to use as frame (loaded above)
# create a table with empty values for each tip to fill in
  agedf <- as.data.frame.matrix(matrix(NA, nrow = length(treemut_RDS@phylo$tip.label), ncol = length(colnames(PD$pdx$agedf))))
  colnames(agedf) <- colnames(PD$pdx$agedf)
  agedf$tip.label <- treemut_RDS@phylo$tip.label
  agedf$age_at_sample_pcy <- rep(donor_age, length(agedf$tip.label))
  agedf$telo_mean_length <- rnorm(length(agedf$tip.label), mean = 1416, sd = 100)
  agedf$age_at_sample_pcy[empty_row - 1] <- 0.000001 # substitute sample age of zeros
  agedf$telo_mean_length[empty_row - 1] <- NA 
  agedf$patient <- "patient"

agedf <- agedf %>%
  dplyr::mutate(per.sample.sensitivity.hybrid=1) %>%
  dplyr::mutate(per.sample.sensitivity.reg=1)

agedf$driver <- 'n/a'
agedf$driver3 <- 'n/a'

# make tree_ml object
tree_ml <- treemut_RDS@phylo
tree_ml$label <- treemut_RDS@data$profile.x 
tree_ml$el.snv.local.filtered <- treemut_RDS@phylo$edge.length # edited Sept 11 25
tree_ml$per.branch.sensitivity.hybrid.multi <- rep(1, length(tree_ml$edge.length))

# make nodes object 
# you probably will not need this bc we use the "null model" tree as our ultrametric tree
nodes <- as.data.frame.matrix(matrix(NA, nrow = 1, ncol = length(PD$nodes)))
colnames(nodes) <- colnames(PD$nodes)


myroot <- castor::find_root(tree_ml)
all_nodes <- seq(length(tree_ml$tip.label)+1,length(tree_ml$tip.label)+tree_ml$Nnode)
nodes_noroot <- setdiff(all_nodes, myroot)
mynode <- nodes_noroot[4]

nodes$node <- mynode 
nodes$driver <- 'n/a'
nodes$status <- 1
nodes$driver2 <- 'n/a'
nodes$driver3 <- 'n/a'
nodes$child_count <- ape::node.depth(treemut_RDS@phylo)[mynode]

# create PD object from our dataframes 
PD_ultra <- list()
PD_ultra$pdx <- list()
PD_ultra$pdx$agedf <- agedf
PD_ultra$pdx$tree_ml <- tree_ml
PD_ultra$nodes <- nodes
PD_ultra$localx.correction2 <- 1.01476   # check
PD_ultra$patient <- "patient"


# use rtreefit to get time tree
treemodel = "poisson_tree"
PD_ultra=wraptreefit(PD_ultra,niter=NITER,b.fit.null = TRUE, method=treemodel)
PD_ultra$fit$poisson_tree$nullmodel$ultratree
PD_ultra$fit$poisson_tree$altmodel$summary

write_rds(PD_ultra, file=paste0(mydir, "Treefit.RDS"))
write_rds(PD_ultra$fit$poisson_tree$nullmodel$ultratree, file=paste0(mydir, "TimeTree_NullModel.RDS"))
write_rds(PD_ultra$fit$poisson_tree$altmodel$ultratree, file=paste0(mydir, "TimeTree_AltModel.RDS"))
