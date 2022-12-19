# # # # # # # # # # # # # # # # # # # # # # # #
#                                             #
#  Single cell RNA-sequencing data analysis   #
#                                             #
# # # # # # # # # # # # # # # # # # # # # # # #
#
# scrnaseq-pseudotime-monocle3.r
#
# Perform pseudotime analysis using Monocle (>=v3)
#
# 0. Resources, variables and parameters ----
#
# Mandatory variables to be assigned before running the script:
#   - path_scrna_obj <- ASSIGN_PATH_SCRNA_OBJ
#     path to seurat object with raw counts and metadata (e.g. obtained from scraseq-get_scrna_obj.r script)
#   - path_results    <- ASSIGN_PATH_RESULTS
#     path to results folder
#   - metadata_features <- ASSIGN_METADATA_FEATURES
#     character vector of features in metadata for grouping and visualization
# 
suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(SeuratWrappers)))
suppressWarnings(suppressMessages(library(monocle3)))
suppressWarnings(suppressMessages(library(scRNAseqRtools)))

# configuration ---
CORES = 12
# Set OpenMP threads
set_openmp_threads = TRUE
if(set_openmp_threads) {
  Sys.setenv("OMP_NUM_THREADS" = CORES)
}

# paths ---
path_scrna_obj <- "seurat/results_seurat-cluster/seurat_obj.rds"
path_results   <- "seurat/results_pseudotime-monocle3"

if(!dir.exists(path_results)) dir.create(path_results, recursive = T)

# global parameters ---
metadata_features <- c("TimeGenotype","seurat_clusters") # used for grouping in plots, character vector
seurat_assay = "RNA" # RNA | Default
# clustering ---
cluster_feature_name = "seurat_clusters" # Cluster
use_previous_clusters = F
reduction = "umap"
cluster_cells_resolution = NULL # NULL monocle3 default
cluster_cells_k = 5 # 20 monocle3 default

# visualization ---
p_cell_size = .65

# palette ---
get_palette_features <- function(metadata_features, pals = NULL)
{
  if(is.null(pals)) {
    # set defaults from ggsci
    pals <- ls('package:ggsci', pattern = 'pal')[1:length(metadata_features)] 
  }
  palette_features <- lapply(pals, function(p) {
    # fill missing items fromm RColorBrewer 
    rcb_pal <- sample(rownames(subset(RColorBrewer::brewer.pal.info,category=="div")),1)
    c(get(p)()(9),RColorBrewer::brewer.pal(9,rcb_pal))
  })
  names(palette_features) <- metadata_features
  return(palette_features)
}

set_feature_palettes = T

if(set_feature_palettes) {
  pal_default <- get_palette_features(metadata_features, pals = NULL)
} else {
  pal_default <- ggsci::pal_d3()(10)  
}
# alternatively, set pal_default to manually defined palettes
# pal_default <- ...
pal1 <- RColorBrewer::brewer.pal(6, "Greys")[c(3,6)]
names(pal1) <- c("EB_3D_WT","EB_9D_WT")
pal2 <- RColorBrewer::brewer.pal(6, "Oranges")[c(3,6)]
names(pal2) <- c("EB_3D_3BKO","EB_9D_3BKO")
pal3 <- RColorBrewer::brewer.pal(6, "Blues")[c(3,6)]
names(pal3) <- c("EB_3D_3AKO","EB_9D_3AKO")
pal_default[["TimeGenotype"]] <- c(pal1,pal2,pal3)

# read object ---
scrna_obj <- readRDS(path_scrna_obj)

save_object = FALSE # add pseudotime data to input object and save
save_cds_object = TRUE # save cds object
s33d = 1991
# 1. Preprocess data ----
# Convert object to cell_data_set ---
if(class(scrna_obj)=="SingleCellExperiment") {
  cds <- as.cell_data_set(as.Seurat(scrna_obj))
  rowData(cds) <- rowData(sce)
  rowData(cds)$gene_short_name <- rowData(cds)$gene_name
} else if(class(scrna_obj)=="Seurat") {
  if(seurat_assay=="Default") {
    cds <- as.cell_data_set(scrna_obj)
    rowData(cds) <- scrna_obj[[DefaultAssay(scrna_obj)]]@meta.features
    rowData(cds)$gene_short_name <- rownames(scrna_obj[[DefaultAssay(scrna_obj)]])
  } else{
    cds <- as.cell_data_set(scrna_obj, assay = seurat_assay)
    rowData(cds) <- scrna_obj[[seurat_assay]]@meta.features
    rowData(cds)$gene_short_name <- rownames(scrna_obj[[seurat_assay]])
  }
 
} else {
  cds <- scrna_obj
  rm(scrna_obj)
}


if(reduction!="umap") {
  reducedDim(cds,"UMAP") <- reducedDim(cds,toupper(reduction))
}
# 2. Pseudotime analysis ----
cds <- cluster_cells(cds = cds
			, reduction_method = "UMAP"
			, resolution=cluster_cells_resolution
			, k = cluster_cells_k)

if(use_previous_clusters) {
  # WORKAROUND to use existing clustering results ---
  to_clusters <- colData(cds)[[cluster_feature_name]]
  names(to_clusters) <- colnames(cds)
  cds@clusters[[toupper(reduction)]]$clusters <- to_clusters
  # ---
}

cds <- learn_graph(cds, use_partition = F)
cds <- order_cells(cds, root_pr_nodes = "Y_17") # Y_1 #Y_44

# Retreive pseudotime ---
colData(cds)$Pseudotime <- monocle3::pseudotime(cds)
colData(cds)$Pseudotime_scaled <- (colData(cds)$Pseudotime-min(colData(cds)$Pseudotime))/max(colData(cds)$Pseudotime)-min(colData(cds)$Pseudotime)

# 3. Visualization ----
pumap_trajectory_principal_points <- plot_cells(cds
                                                , color_cells_by = "pseudotime"
                                                , cell_size = p_cell_size
                                                , label_branch_points=T
                                                , label_leaves = T
                                                , label_roots = F
                                                , label_cell_groups = F
                                                , trajectory_graph_segment_size = .5
                                                , label_principal_points = T
                                                ) + theme_bw() + 
  my_theme_2 + theme(legend.key.size = unit(4,'mm'),legend.position = "right",aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

outfile <- paste0(path_results,"/dimplot.",reduction,".trajectory_principal_points.pdf")
pdf(file = outfile, paper = "a4r",w=unit(5,'cm'),h=unit(6,"cm"),useDingbats = F)
print(pumap_trajectory_principal_points)
dev.off()

pumap_trajectory_pseudotime <- plot_cells(cds
                                          , color_cells_by = "pseudotime"
                                          , cell_size = p_cell_size
                                          , label_branch_points=T
                                          , label_leaves = T
                                          , label_roots = F
                                          , label_cell_groups = F
                                          , trajectory_graph_segment_size = .5
                                          , label_principal_points = F
) + theme_bw() + 
  my_theme_2 + theme(legend.key.size = unit(4,'mm'),legend.position = "right",aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


pumap_trajectory <- list()
for(feature in metadata_features) {
  
  pumap_trajectory[[feature]] <- plot_cells(cds
                                            , color_cells_by = feature
                                            , cell_size = p_cell_size
                                            , label_branch_points=T
                                            , label_leaves = T
                                            , label_roots = F
                                            , label_cell_groups = F
                                            , trajectory_graph_segment_size = .5) + 
    theme_bw() + my_theme_2 + scale_color_manual(values = pal_default[[feature]]) + 
    theme(legend.key.size = unit(4,'mm'),legend.position = "right",aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  outfile <- paste0(path_results,"/dimplot.",reduction,".trajectory_",feature,".pdf")
  pdf(file = outfile, paper = "a4r",w=unit(5,'cm'),h=unit(6,"cm"),useDingbats = F)
  print(pumap_trajectory[[feature]])
  dev.off()
}

pumap_trajectory <- c(list("Pseudotime"=pumap_trajectory_pseudotime),pumap_trajectory)
pumap_trajectory_list <- ggpubr::ggarrange(plotlist = pumap_trajectory, nrow=1,align="hv")

outfile <- paste0(path_results,"/dimplot.",reduction,".trajectory.pdf")
pdf(file = outfile, paper = "a4r",w=unit(4.5*(length(metadata_features)+1),'cm'),h=unit(6.5,"cm"),useDingbats = F)
print(pumap_trajectory_list)
dev.off()

# 4. Save object ----
if(save_object) {
  scrna_obj@meta.data <- as.data.frame(colData(cds))
  outfile <- paste0(path_results,"/",basename(path_scrna_obj))
  message(" -- saving cds object to: ", outfile)
  saveRDS(cds, file = outfile)
}

if(save_cds_object) {
  outfile <- paste0(path_results,"/cds.rds")
  message(" -- saving cds object to: ", outfile)
  saveRDS(cds, file = outfile)
}
