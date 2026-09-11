
library(data.table)
library(dplyr)
library(tidyr)
library(BSgenome.Hsapiens.UCSC.hg38)


##################
# Load arguments #
##################

args = commandArgs(trailingOnly=TRUE)
somatic_file <- args[1]
germline_file <- args[2]
mychrom <- args[3]
outdir=dirname(germline_file)


print("Arguments: ")
print(paste0("somatic: ",somatic_file))
print(paste0("germline: ",germline_file))
print(paste0("chrom: ",mychrom))


######################################################
# Create the required inputs. Swap the VAFs when 1|0 #
######################################################

# Load the somatic mutations
somatic_data_ori=data.table::fread(cmd = paste0("awk '{print $0}' ", somatic_file), dec = ".") %>%
  dplyr::filter(CHROM==mychrom)
somatic_data <- somatic_data_ori %>%
  dplyr::select(colnames(somatic_data_ori)[grepl("(CHROM|POS|^REF$|ALT$|.VAF$)",colnames(somatic_data_ori))])
somatic_data_dp <- somatic_data_ori %>%
  dplyr::select(colnames(somatic_data_ori)[grepl("(CHROM|POS|^REF$|ALT$|.AD$)",colnames(somatic_data_ori))])

# Load the phased phased germline SNPs
# Remove non-heterozygous sites
mytable_shapeit <- data.table::fread(germline_file, dec = ".") %>%
  dplyr::filter(AF>0.3 & AF<0.7 & Mean_Depth>3)
paste0("Number of SNPs: ",nrow(mytable_shapeit))
# remove variants with all missing data?

# obtain VAF
myvafs <- mytable_shapeit %>%
  dplyr::select(colnames(mytable_shapeit)[grepl(".VAF", colnames(mytable_shapeit))]) %>% as.data.frame()
# obtain Phased GT
myphased_gts <- mytable_shapeit %>%
  dplyr::select(colnames(mytable_shapeit)[grepl(".PhasedGT", colnames(mytable_shapeit))]) %>% as.data.frame()
# obtain DP
myphased_dps <- mytable_shapeit %>%
  dplyr::select(colnames(mytable_shapeit)[grepl(".DP", colnames(mytable_shapeit))]) %>% as.data.frame()
# Set up as NAs all VAFs supported by 3 reads or less 
myphased_gts[myphased_dps<=3] <- NA
myphased_gts[is.na(myphased_gts)] <- ".|."
myphased_gts2 <- myphased_gts
# swap VAF depending on phasing
rev_vafs <- 1-myvafs %>% as.data.frame()
myphased_gts2[myphased_gts == "1|0"] <- rev_vafs[myphased_gts == "1|0"]
myphased_gts2[myphased_gts != "1|0" & myphased_gts != ".|."] <- myvafs[myphased_gts != "1|0" & myphased_gts != ".|."]
myphased_gts2[myphased_gts == ".|."] <- NA
myphased_gts2 <- apply(myphased_gts2, 2, as.numeric)
rownames(myphased_gts2) <- mytable_shapeit$POS
myphased_dps <- data.matrix(myphased_dps)
rownames(myphased_dps) <- mytable_shapeit$POS
# Remove rows with all NAs
myphased_dps <- myphased_dps[!is.na(rowMeans(myphased_gts2, na.rm = T)),]
mytable_shapeit <- as.data.frame(mytable_shapeit)
mytable_shapeit <- as.data.frame(mytable_shapeit)[!is.na(rowMeans(myphased_gts2, na.rm = T)),]
myphased_gts2 <- myphased_gts2[!is.na(rowMeans(myphased_gts2, na.rm = T)),]


###################################################
# Create a list of 200kb intervals per chromosome #
###################################################

chrom_sizes <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
# Function to generate contiguous 200kb intervals for a chromosome
generate_contiguous_intervals <- function(chrom_length, interval_size = 200000) {
  num_intervals <- floor(chrom_length / interval_size)
  intervals <- data.frame(
    start = seq(1, by = interval_size, length.out = num_intervals),
    end = seq(interval_size, by = interval_size, length.out = num_intervals)
  )
  return(intervals)
}
intervals_list <- list()
for (chrom in names(chrom_sizes)) {
  chrom_length <- chrom_sizes[chrom]
  intervals_list[[chrom]] <- generate_contiguous_intervals(chrom_length)
}
myinterval_list <- intervals_list[[`mychrom`]]
print(paste0("Number of intervals: ",nrow(myinterval_list)))

####################################################################
# Fit a regression with the BAFs,                                  #
# and apply the model to predict the BAF of the somatic variant.   #
# Calculate values for ADO per region per cell                     #
####################################################################

variant_data <- cbind(start=as.numeric(mytable_shapeit$POS), 
                      vaf=as.numeric(mytable_shapeit$AF), 
                      total_ad=rowSums(myphased_dps, na.rm = T))  %>%
  as.data.frame()
# Convert to data.table for faster filtering
variant_data_dt <- as.data.table(variant_data)

# Ensure rownames are converted to a column before processing
myphased_gts2_dt <- as.data.table(myphased_gts2, keep.rownames = "start")
myphased_gts2_dt[, start := as.numeric(start)]  # Convert 'start' to numeric
myphased_dps_dt <- as.data.table(myphased_dps, keep.rownames = "start")
myphased_dps_dt[, start := as.numeric(start)]  # Convert 'start' to numeric


    
###################################################  
# Create the loess models, one per cell per chunk #
################################################### 
    
    #chunk_results <- list()
    # Loop through each cell in the myphased_gts2 matrix
    ADOandModelPerCell <- function(cell_idx, chunk_info, chunk_info_dp, pseudobulk) {
      cell_germ_vafs <- chunk_info[, cell_idx, with = FALSE] %>% as.data.frame()
      cell_germ_dps <- chunk_info_dp[, cell_idx, with = FALSE] %>% as.data.frame()
      # if all germline VAFs are zero, the loess will not work
      if (sum(!is.na(cell_germ_vafs))<=1){
        loess_cell <- NA
        # if there is more than 1 germline VAF:
      } else {
        cell_positions <- pseudobulk$start[!is.na(cell_germ_vafs)]
        # Clean data by removing NA values
        cell_germ_dps <- cell_germ_dps[!is.na(cell_germ_vafs),]
        cell_germ_vafs <- cell_germ_vafs[!is.na(cell_germ_vafs),]
        # Fit LOESS model for the cell and predict VAF and standard error
        loess_cell <- loess(cell_germ_vafs ~ cell_positions, weights = cell_germ_dps, degree = 2)
      }
      return(loess_cell)
    }


    #chunk_results <- list()
    # Loop through each cell in the myphased_gts2 matrix
    DiffPerCell <- function(cell_idx, chunk_info, chunk_info_dp, pseudobulk) {
      cell_germ_vafs <- chunk_info[, cell_idx, with = FALSE] %>% as.data.frame()
      cell_germ_dps <- chunk_info_dp[, cell_idx, with = FALSE] %>% as.data.frame()
      #print(cell_idx)
      # if all germline VAFs are zero, the loess will not work
      if (sum(!is.na(cell_germ_vafs))<=1){
        diff_median <- NA
        # if there is more than 1 germline VAF:
      } else {
        cell_positions <- pseudobulk$start[!is.na(cell_germ_vafs)]
        # Clean data by removing NA values
        cell_germ_dps <- cell_germ_dps[!is.na(cell_germ_vafs),]
        cell_germ_vafs <- cell_germ_vafs[!is.na(cell_germ_vafs),]
        # Calculate the absolute difference from 0.5
        # Compute the standard deviation of these differences
        differences <- abs(cell_germ_vafs - 0.5)
      }
      return(differences)
    }

    
    
    PerSiteCellResults <- function(cell_index, loess_model, 
                                   somatic_variant, somatic_vafs_vector,somatic_vafs_dp_vector ){
     # print(paste0("Cell index (get pbinom): ",cell_index))
      myloess_model <- loess_model[[cell_index]]
      if (is.na(myloess_model)){
        pval_binom <- NA
        return(pval_binom)
      } else {
      if (sum(is.na(myloess_model$fitted))>1 |  is.na(somatic_vafs_vector[cell_index]) | as.numeric(somatic_vafs_vector[cell_index])==0 ){
        pval_binom <- NA
      } else {
        pred <- stats::predict(myloess_model, somatic_variant, se = TRUE)
        # make sure predicted value is between zero and one
        if (is.na(pred$fit)){
          pval_binom <- NA
        } else {
          if ( pred$fit < 0 ) {
            pred$fit <- 0
          } else if ( pred$fit > 1 ) {
            pred$fit <- 1
          }
          # Switch vaf allele if that makes it closer to the predicted vaf.
          dist <- abs(as.numeric(somatic_vafs_vector[cell_index])-pred$fit)
          dist2 <- abs(abs(1-(as.numeric(somatic_vafs_vector[cell_index])))-pred$fit)
          dist3 <- min(dist,dist2)
          if (dist2<dist){
            reverse <- 1
          } else {
            reverse <- 0
          }
          # counts supporting the allele
          allele_counts <- strsplit(somatic_vafs_dp_vector[cell_index], ",")[[1]]
          a_count <- as.numeric(allele_counts[[1]])
          b_count <- as.numeric(allele_counts[[2]])
          ad_somatic <- sum(a_count,b_count, na.rm = T)
          if (reverse==1){
            alternative_count <- a_count
          } else {
            alternative_count <- b_count
          }
          # When 1 the predicted and observed are similar. The lowest, the most different are the VAFs
          pval_binom <- binom.test(alternative_count, n = ad_somatic, p = pred$fit, alternative = "two.sided")$p.value
          
        }
      }
        return(pval_binom)
    }
  }
    
    
    # Obtain expected VAF for each somatic position
    ObtainPredictionPerSite <- function(somatic_index, mysomatic, loess_model, somatic_data, somatic_data_dp, mysomatic_alt,mysomatic_ref){
      somatic_variant <- mysomatic[somatic_index]
      alternative_allele <- mysomatic_alt[somatic_index]
      reference_allele <- mysomatic_ref[somatic_index]
      #print(paste0("Somatic variant:",somatic_variant))
      # calculate the difference of the observed to the predicted  
      somatic_vafs_vector <- somatic_data %>%
        dplyr::filter(POS==somatic_variant) %>%
        unique() %>%
        dplyr::filter(REF==reference_allele) %>%
        dplyr::filter(ALT==alternative_allele) %>%
        dplyr::select(-CHROM,-POS,-ALT,-REF) %>%
        unlist()
      somatic_vafs_dp_vector  <- somatic_data_dp %>%
        dplyr::filter(POS==somatic_variant) %>%
        unique() %>%
        dplyr::filter(REF==reference_allele) %>%
        dplyr::filter(ALT==alternative_allele) %>%
        dplyr::select(-CHROM,-POS,-ALT,-REF) %>%
        unlist()
      if (length(somatic_vafs_vector[somatic_vafs_vector>0 & !is.na(somatic_vafs_vector)])==0){
        pvalues <- NA
      } else {
        index_cells_carrying_variants <- which(somatic_vafs_vector>0)
        pvalues <- sapply(as.vector(index_cells_carrying_variants), PerSiteCellResults, loess_model=loess_model, 
                          somatic_variant=somatic_variant, somatic_vafs_vector=somatic_vafs_vector, somatic_vafs_dp_vector=somatic_vafs_dp_vector)
      }
      if (length(pvalues)>1) {
        max_pval <- max(pvalues, na.rm = T)
        if (max_pval==-Inf){
          max_pval=NA
        }  
      } else {
        if (is.na(pvalues)) {
          max_pval=NA
        } else {
          max_pval <- pvalues 
        }
       }
      
      return(max_pval)
    }
    
     
    
    
PerChunk <- function(chunk_idx, myinterval_list, somatic_data, somatic_data_dp){
      
      print(paste0("Chunk ",chunk_idx))
      chunk_interval <- myinterval_list[chunk_idx,]
      a <- as.integer(chunk_interval$start)
      b <- as.integer(chunk_interval$end)

      ##################################################################
      # Fit germline SNPs on pseudobulk and obtain pseudobulk residual #
      ##################################################################
      pseudobulk <- variant_data_dt[start >= a & start <= b]      
      mysomatic <- somatic_data %>%
        dplyr::filter(POS>=a & POS<=b) %>%
        unique() %>%
        dplyr::select(POS) %>% unlist() %>% as.numeric()
      mysomatic_ref <- somatic_data %>%
        dplyr::filter(POS>=a & POS<=b) %>%
        unique() %>%
        dplyr::select(REF) %>% unlist() %>% as.character()
      mysomatic_alt <- somatic_data %>%
        dplyr::filter(POS>=a & POS<=b) %>%
        unique() %>%
        dplyr::select(ALT) %>% unlist() %>% as.character()


      if (nrow(pseudobulk)==0 | sum(!is.na(pseudobulk$vaf))<=1){
        return(NULL)
      } else {
        # Fit LOESS model for the pseudobulk and predict VAF and standard error
        loess_pseudobulk <- loess(pseudobulk$vaf ~ pseudobulk$start, weights = pseudobulk$total_ad, degree = 2)
        residual_standard_error_pseudobulk <- loess_pseudobulk$s
      
      ############################################################
      # Extract germline information for the chunk for all cells #
      ############################################################
      
      # Now filter the data based on 'start' and remove the 'start' column
      chunk_info <- myphased_gts2_dt[start >= a & start <= b, .SD]
      chunk_info_dp <- myphased_dps_dt[start >= a & start <= b, .SD]
      # Optionally remove the 'start' column if you no longer need it
      chunk_info[, start := NULL]
      chunk_info_dp[, start := NULL]
      
      #########################################################################
      # obtain annotations for the somatic mutations overlapping the interval #
      #########################################################################
      
      # Vectorized check without the need for sapply
      CellModels <- lapply(1:ncol(chunk_info), ADOandModelPerCell, chunk_info_dp=chunk_info_dp, chunk_info=chunk_info, pseudobulk=pseudobulk)
      #Differences <- lapply(1:ncol(chunk_info), DiffPerCell, chunk_info_dp=chunk_info_dp, chunk_info=chunk_info, pseudobulk=pseudobulk)
      #print(Differences)
      # Return NULL if the length is zero
      if (length(mysomatic) == 0) {
        return(NULL)
      } else {
        info_somatic <- lapply(1:length(mysomatic), ObtainPredictionPerSite, loess_model=CellModels, mysomatic=mysomatic, mysomatic_alt=mysomatic_alt, mysomatic_ref=mysomatic_ref, somatic_data=somatic_data, somatic_data_dp=somatic_data_dp)
        output <- cbind(chunk_idx,mysomatic,mysomatic_alt,loess_pseudobulk$s,unlist(info_somatic))
        return(output)
      }
     }
}

#PerChunk(chunk_idx = 282,somatic_data=somatic_data, somatic_data_dp=somatic_data_dp, myinterval_list=myinterval_list)
myresults <- lapply(seq(1,nrow(myinterval_list)), PerChunk, somatic_data=somatic_data, somatic_data_dp=somatic_data_dp, myinterval_list=myinterval_list)   
myADOannotations_forsomatic <- do.call("rbind",myresults) %>%
  as.data.frame() %>%
  dplyr::mutate(chrom=mychrom) %>%
  dplyr::select(chrom, mysomatic, mysomatic_alt, chunk_idx, V4, V5) %>%
  #dplyr::mutate(V4=ifelse(is.null(V4),NA, V4)) %>%
  dplyr::mutate(V5=ifelse(is.na(V5),NA, V5))
    
colnames(myADOannotations_forsomatic) <- c("CHROM","POS","ALT","PTA_ADO200kbinterval","PTA_ADOResidualLoessPseudobulk","PTA_ADOmaxbinompvalue")

write.table(x = myADOannotations_forsomatic, file=paste0(outdir, "/ADO_annotations.",mychrom,".txt"), quote = F, 
            row.names = F)

print("Finished!")





