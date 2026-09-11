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
patient <- as.character(args[3])
recall_object <- args[4]


print("Arguments: ")
print(mydir)
print(treefile)
print(patient)
print(recall_object)



treemut_RDS <- readRDS(file = paste0(mydir, treefile)) 

recall_percell <- readRDS(recall_object) %>%
  dplyr::filter(donor==patient) %>%
  dplyr::rename(label=cell) %>%
  dplyr::select(label, recall)

####################################################################################################
# correct external branches and internal branches using ADO product of descendant cells per branch #
####################################################################################################

cell_order <- readLines(file.path(mydir, "PhylogenyMappingTipOrder.txt"))
cell_order <- strsplit(trimws(cell_order), "\\s+")[[1]]
length(cell_order)
ado_vec <- 1 - recall_percell$recall[match(cell_order, recall_percell$label)]
sum(is.na(ado_vec))


compute_branch_ado_product <- function(profile, ado_vec) {
  if (is.na(profile)) return(NA_real_)
  profile <- substring(profile, 1, nchar(profile) - 1)   # drop the zeros output bit (end)
  bits <- as.integer(strsplit(profile, "")[[1]])
  if (length(bits) != length(ado_vec)) {
    stop(sprintf("profile length (%d) != number of cells (%d)", length(bits), length(ado_vec)))
  }
  if (sum(bits) == 0) return(NA_real_)   # mutation exclusive to the outgroup branch
  prod(ado_vec[bits == 1])
}

treemut_RDS_corrected <- treemut_RDS %>%
  treeio::as_tibble() %>%
  dplyr::left_join(recall_percell, by = "label") %>%
  dplyr::mutate(
    ado_product      = purrr::map_dbl(profile.x, compute_branch_ado_product, ado_vec = ado_vec),
    effective_recall = 1 - ado_product,
    branch.length = dplyr::case_when(
      is.na(effective_recall)  ~ branch.length,   # root or uncorrectable
      effective_recall <= 1e-6 ~ branch.length,   # degenerate recall: leave uncorrected
      TRUE                      ~ branch.length / effective_recall
    )
  ) %>%
  treeio::as.treedata()



saveRDS(object = treemut_RDS_corrected, file = paste0(mydir, "TreeMut.genome.recallcorrected.rds"))
