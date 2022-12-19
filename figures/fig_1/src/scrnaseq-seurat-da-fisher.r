# # # # # # # # # # # # # # # # # # # # # # # #
#                                             #
#  Single cell RNA-sequencing data analysis   #
#                                             #
# # # # # # # # # # # # # # # # # # # # # # # #
#
# scrnaseq-seurat-da-fisher.r
#
# Differential cell type abundance analysis
#
# 0. Resources, variables and parameters ----
#
# Mandatory variables to be assigned before running the script:
#   - path_seurat_obj <- ASSIGN_PATH_SEURAT_OBJ
#     path to seurat object with raw counts and metadata (e.g. obtained from scraseq-get_scrna_obj.r script)
#   - path_results    <- ASSIGN_PATH_RESULTS
#     path to results folder
#   - metadata_feature_test <- ASSIGN_METADATA_FEATURE_TEST
#     character vector of cell type feature to test for differential abundance in cluster
#   - metadata_feature_group <- ASSIGN_METADATA_FEATURE_GROUP
#     character vector of cell type group to test for differential abundance in cluster
#
# 0. Resources ----
suppressWarnings(suppressMessages(library(Seurat)))

tab_to_write <- function(x, nm) 
{
  x <- cbind(rownames(x),x)
  x <- as.data.frame(x)
  colnames(x)[1] <- nm
  x
}

path_seurat_obj <- "seurat/results_seurat-cluster/seurat_obj.rds"
path_results    <- "seurat/results_seurat-cluster/da-fisher"

if(!dir.exists(path_results)) dir.create(path_results, recursive = T)

metadata_feature_test <- "Genotype" # Cell type feature to test for differential abundance in cluster
metadata_feature_group <- "seurat_clusters" # Cell type group to test for differential abundance in cluster

save_comparisons_obj = TRUE
# 1. Read data ----
seurat_obj <- readRDS(path_seurat_obj)

# 2. Create contingency tables ----
tab  <- table(seurat_obj[[metadata_feature_test, drop = T]],seurat_obj[[metadata_feature_group, drop = T]])
tot  <- rowSums(table(seurat_obj[[metadata_feature_test, drop = T]],seurat_obj[[metadata_feature_group, drop = T]]))
perc <- round(100*tab/tot,2)

# 3. Test pairwise differential abundance ----
my_contrasts <- as.data.frame(combn(rownames(tab),2))
colnames(my_contrasts) <- paste0("contrast-",apply(my_contrasts,2,paste0,collapse="vs"))
my_clusters  <- colnames(tab)

# re-order --
my_contrasts[[1]] <- rev(my_contrasts[[1]])
colnames(my_contrasts)[1] <- gsub("3AKOvs3BKO","3BKOvs3AKO",colnames(my_contrasts)[1])

comparisons       <- list()
comparisons_table <- list()
for(cl in my_clusters) {
  comparisons[[cl]] <- list()
  comparisons_table[[cl]] <- list()
  for(cn in colnames(my_contrasts)) {
    comparisons[[cl]][[cn]] <- list()
    comparisons_table[[cl]][[cn]] <- list()
    i <- my_contrasts[[cn]][1]
    j <- my_contrasts[[cn]][2]
    m <- matrix(c(tab[i,cl]
                  , tot[i]-tab[i,cl]
                  , tab[j,cl]
                  ,tot[j]-tab[j,cl])
                , nrow=2)
    comparisons[[cl]][[cn]][["m"]] <- m
    for(alternative in c("l","g","two.sided")) {
      comparisons[[cl]][[cn]][[alternative]] <- fisher.test(m, alternative=alternative)
      comparisons_table[[cl]][[cn]][[alternative]] <- data.frame("cluster"=cl
                                                                 ,"contrast"=cn
                                                                 ,"alternative"=alternative
                                                                 , "p.value"=comparisons[[cl]][[cn]][[alternative]]$p.value
                                                                 , "odds.ratio"=comparisons[[cl]][[cn]][[alternative]]$estimate)
    }
    comparisons_table[[cl]][[cn]] <- do.call(rbind.data.frame,comparisons_table[[cl]][[cn]])
  }
  comparisons_table[[cl]] <- do.call(rbind.data.frame,comparisons_table[[cl]])
}
comparisons_table <- do.call(rbind.data.frame,comparisons_table)

# 4. Save results ----
if(save_comparisons_obj) {
  outfile <- paste0(path_results,"/comparisons.rds")
  message(" -- saving to: ", outfile)
  saveRDS(comparisons, file = outfile)
}

outfile <- paste0(path_results,"/comparisons_table.txt")
write.table(comparisons_table, file = outfile, sep = "\t", quote = F, row.names = F)

tab_w   <- tab_to_write(tab, metadata_feature_test)
outfile <- paste0(path_results,"/tab.txt")
write.table(tab_w, file = outfile, sep = "\t", quote = F, row.names = F)

perc_w   <- tab_to_write(perc, metadata_feature_test)
outfile <- paste0(path_results,"/perc.txt")
write.table(perc_w, file = outfile, sep = "\t", quote = F, row.names = F)
