# # # # # # # # # # # # # # # # # # # # # # # #
#                                             #
#  Single cell RNA-sequencing data analysis   #
#                                             #
# # # # # # # # # # # # # # # # # # # # # # # #
#
# scrnaseq-seurat-cluster.r
#
# Cluster cells with Seurat.
#
# 0. Resources, variables and parameters ----
#
# Mandatory variables to be assigned before running the script:
#   - path_seurat_obj <- ASSIGN_PATH_SEURAT_OBJ
#     path to seurat object with raw counts and metadata (e.g. obtained from scraseq-get_scrna_obj.r script)
#   - path_results    <- ASSIGN_PATH_RESULTS
#     path to results folder
#   - metadata_features <- ASSIGN_METADATA_FEATURES
#     character vector of features in metadata for grouping and visualization
# 
suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(scRNAseqRtools)))

# configuration ---
CORES = 12
# Set OpenMP threads
set_openmp_threads = TRUE
if(set_openmp_threads) {
  Sys.setenv("OMP_NUM_THREADS" = CORES)
}

# paths ---
path_seurat_obj <- "seurat/results_seurat-standard/seurat_obj.rds" # ASSIGN_PATH_SEURAT_OBJ
path_results    <- "seurat/results_seurat-cluster" # ASSIGN_PATH_RESULTS

if(!dir.exists(path_results)) dir.create(path_results, recursive = T)

# global parameters ---
metadata_features <- c("Library_PreP_Batch","TimeGenotype") # used for grouping in plots, character vector

# Dimensionality reduction
max_pcs = 100
npcs = 8
npcs_method = "elbow" # elbow | jackstraw
npc_loadings = 4
npcs_jackstraw = 20
perplexity = 10
tsne_reduction = "pca"
umap_reduction = "pca"

# Clustering
clustering_reduction = "pca"
clustering_resolution = .8 # .8 seurat default

# Cluster markers 
find_all_markers_assay = "RNA"
find_all_markers_latent_vars = c("Library_PreP_Batch","percent.mt")
find_all_markers_test_use = "LR" # LR, bimod, t, negbinom, poisson, LR, MAST, DESeq2
find_all_markers_only_pos = TRUE
find_all_markers_min_pct  = 0.05
find_all_markers_logfc_th = 0.1
markers_do_heatmap = FALSE
markers_do_vln = FALSE
markers_heatmap_n_top = 20
markers_vln_n_top = 20
markers_vln_ncol = 3
markers_reduction_n_top = 5

# Visualization
reductions = c("tsne","umap")
extended_plots = F # make additional plots for variable features, PCA component scores, markers
reductions_markers = "umap"
reductions_technical_features = "umap"

# read object ---
seurat_obj <- readRDS(path_seurat_obj)

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

save_object = TRUE # save obj 
rename_seurat_obj = FALSE # rename obj when saving if exists
s33d = 1991
# 1. Clustering ----
seurat_obj <- FindNeighbors(seurat_obj
                            , reduction = clustering_reduction
                            , dims = 1:npcs
                            , graph.name = paste0(clustering_reduction,c("_nn","_snn")))
seurat_obj <- FindClusters(seurat_obj
                           , resolution = clustering_resolution
                           , graph.name = paste0(clustering_reduction,"_snn"))

# 2. Find cluster markers ----
markers <- FindAllMarkers(seurat_obj
                          , assay = find_all_markers_assay
                          , latent.vars = find_all_markers_latent_vars
                          , only.pos = find_all_markers_only_pos
                          , min.pct = find_all_markers_min_pct
                          , logfc.threshold = find_all_markers_logfc_th
                          , test.use = find_all_markers_test_use)

outfile <- paste0(path_results,"/seurat_clusters_markers.txt.gz")
write.table(markers, file = gzfile(outfile), row.names = F, col.names = T, sep = "\t", quote = F)

# 3. Visualization ----
dimplot_wrapper <- function(seurat_obj
                            , metadata_features
                            , pal = pal_default
                            , reduction = "pca") 
{
  pdimred_list <- list()
  if(!is.list(pal)) {
    pal_dimred <- list("seurat_clusters"=pal)
  } else {
    pal_dimred <- pal
  }
  pdimred_list[["seurat_clusters"]] <- DimPlot(seurat_obj, reduction = reduction, group.by = "seurat_clusters", pt.size = .5) +
    theme_bw() + my_theme + theme(panel.grid.major = element_blank(), legend.key.size = unit(4,'mm'), aspect.ratio = 1) +
    xlab(paste0(gsub("PCA","PC",toupper(reduction))," 1")) + ylab(paste0(gsub("PCA","PC",toupper(reduction))," 2")) + ggtitle("seurat_clusters") + scale_color_manual(values = pal_dimred[["seurat_clusters"]])
  
  for(feature in metadata_features) {
    if(!is.list(pal)) pal_dimred[[feature]] <- pal
    pdimred_list[[feature]] <- DimPlot(seurat_obj, reduction = reduction, group.by = feature, pt.size = .5) +
      theme_bw() + my_theme + theme(panel.grid.major = element_blank(), legend.key.size = unit(4,'mm'), aspect.ratio = 1) +
      xlab(paste0(gsub("PCA","PC",toupper(reduction))," 1")) + ylab(paste0(gsub("PCA","PC",toupper(reduction))," 2")) + ggtitle(feature) + scale_color_manual(values = pal_dimred[[feature]])
  }
  
  pdimred <- ggpubr::ggarrange(plotlist = pdimred_list, align = "hv",nrow=1)
  pdimred <- ggpubr::annotate_figure(pdimred, top=text_grob(paste0("Dimensionality reduction - ",toupper(reduction)), hjust=0.5, size=8))
  
  return(pdimred)
}

# palettes
if(set_feature_palettes) {
  pal_cluster <- get_palette_features(metadata_features = "seurat_clusters"
                                      ,pals = ls('package:ggsci', pattern = 'pal')[length(metadata_features)+1])
  pal_default <- c(pal_default, pal_cluster)
}

pdimred <- list()
for(reduction in reductions) {
  pdimred[[reduction]] <- dimplot_wrapper(seurat_obj, metadata_features, reduction = reduction, pal = pal_default)
  
  outfile <- paste0(path_results, "/dimplot.",reduction,".pdf")
  message(" -- saving to: ",outfile)
  pdf(file = outfile, paper = "a4r", w = unit(10, 'cm'), h = unit(6, 'cm'), useDingbats = F)
  print(pdimred[[reduction]])
  dev.off()  
}

# markers
if(markers_do_vln) {
  
  top_n_genes <- lapply(split(markers, markers$cluster), function(x) head(x, n = markers_vln_n_top))
  top_n_genes <- do.call(rbind.data.frame, top_n_genes)
  
  if(length(top_n_genes$gene)>markers_vln_ncol*5) {
    top_n_genes_list <- split(top_n_genes$gene, ceiling(seq_along(top_n_genes$gene)/as.integer(markers_vln_ncol*5)))
  } else {
    top_n_genes_list <- list(top_n_genes$gene)
  }
  pmarkers_vln <- list()
  for(i in 1:length(top_n_genes_list)) {
    pmarkers_vln[[i]] <- VlnPlot(seurat_obj, features = top_n_genes_list[[i]]
                                 , group.by = "seurat_clusters"
                                 , cols = pal_default[["seurat_clusters"]]
                                 , pt.size = 0.05,ncol = markers_vln_ncol) &
      theme_bw() + my_theme_2 + theme(panel.grid.major = element_blank(), legend.key.size = unit(4,'mm'))
    
  }
  
  outfile <- paste0(path_results, "/seurat_clusters_markers.vln.top_",markers_vln_n_top,".pdf")
  message(" -- writing to: ", outfile)
  pdf(file = outfile, paper = "a4", w = unit(20, 'cm'), h = unit(20, 'cm'), useDingbats = F)
  print(pmarkers_vln)
  dev.off()
  
}
if(markers_do_heatmap) {
  top_n_genes <- lapply(split(markers, markers$cluster), function(x) head(x, n = markers_heatmap_n_top))
  top_n_genes <- do.call(rbind.data.frame, top_n_genes)
  
  outfile <- paste0(path_results, "/seurat_clusters_markers.heatmap.top_",markers_heatmap_n_top,".pdf")
  pdf(file = outfile, paper = "a4", w = unit(8, 'cm'), h = unit(10, 'cm'), useDingbats = F)
  DoHeatmap(seurat_obj, features = top_n_genes$gene) + NoLegend()
  dev.off()
}

if(extended_plots) {
  # technical feature plots
  pf1 <- FeaturePlot(seurat_obj, features = "nFeature_RNA",reduction = reductions_technical_features, pt.size = .5) +
    theme_bw() + my_theme_2 + theme(panel.grid.major = element_blank(), legend.key.size = unit(4,'mm')) +
    xlab(paste0(gsub("PCA","PC",toupper(reductions_technical_features))," 1")) + ylab(paste0(gsub("PCA","PC",toupper(reductions_technical_features))," 2")) 
  pf2 <- FeaturePlot(seurat_obj, features = "nCount_RNA",reduction = reductions_technical_features, pt.size = .5) +
    theme_bw() + my_theme_2 + theme(panel.grid.major = element_blank(), legend.key.size = unit(4,'mm')) +
    xlab(paste0(gsub("PCA","PC",toupper(reductions_technical_features))," 1")) + ylab(paste0(gsub("PCA","PC",toupper(reductions_technical_features))," 2")) 
  pf3 <- FeaturePlot(seurat_obj, features = "percent.mt",reduction = reductions_technical_features, pt.size = .5) +
    theme_bw() + my_theme_2 + theme(panel.grid.major = element_blank(), legend.key.size = unit(4,'mm')) +
    xlab(paste0(gsub("PCA","PC",toupper(reductions_technical_features))," 1")) + ylab(paste0(gsub("PCA","PC",toupper(reductions_technical_features))," 2")) 
  
  pdimred_technical_features <- ggpubr::ggarrange(pf1,pf2,pf3,align = "hv", nrow = 1, ncol = 3)
  outfile <- paste0(path_results, "/dimplot.",reductions_technical_features,"_technical_features.pdf")
  pdf(file = outfile, paper = "a4r", w = unit(10, 'cm'), h = unit(4, 'cm'), useDingbats = F)
  print(pdimred_technical_features)
  dev.off()
  
  top_n_genes <- lapply(split(markers, markers$cluster), function(x) head(x, n = markers_reduction_n_top))
  top_n_genes <- do.call(rbind.data.frame, top_n_genes)
  
  pdimred_markers <- list()
  for(reduction in reductions_markers) {
    pdimred_markers[[reduction]] <- FeaturePlot(seurat_obj, features = top_n_genes$gene, reduction = reduction, pt.size = .1, ncol = 5) &
      theme_bw() + my_theme_2 + theme(panel.grid.major = element_blank(), legend.key.size = unit(4,'mm'))
    
    outfile <- paste0(path_results,"/seurat_clusters_markers.dimred_",reduction,".pdf")
    message(" -- output file: ", outfile)
    pdf(file = outfile, paper = "a4", useDingbats = F, w=unit(30,'cm'), height = unit(30,'cm'))
    print(pdimred_markers[[reduction]])
    dev.off()
  }
}
# 4. Save object ----
if(save_object) {
  outfile <- paste0(path_results,"/seurat_obj.rds")
  if(file.exists(outfile) && rename_seurat_obj) {
    message("[!] Object exists, renaming")
    outfile <- gsub("\\.rds$",".proc.rds",outfile)
  }
  message(" -- saving seurat object to: ", outfile)
  saveRDS(seurat_obj, file = outfile)
  
}
