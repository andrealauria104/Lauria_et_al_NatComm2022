#!/usr/bin/env Rscript
# # # # # # # # # # # # # # # # # #
#                                 #
#  RNA-sequencing data analysis   #
#                                 #
# # # # # # # # # # # # # # # # # #
#
# rnaseqrtools-dge-edger.r
#
# Perfom Differential Expression Analysis with edgeR
# 
# 0. Resources, variables and parameters ----
#
# Mandatory variables to be assigned before running the script:
#   - path_dgelist_obj <- ASSIGN_PATH_DGELIST_OBJ
#     path to dgelist object with raw counts and metadata (e.g. obtained from rnaseq-get_dgelist_obj.r script)
#   - path_results    <- ASSIGN_PATH_RESULTS
#     path to results folder
#   - metadata_features_col <- ASSIGN_METADATA_FEATURES_COL
#     character vector of features in metadata for grouping and visualization
# 
suppressWarnings(suppressMessages(library(RNAseqRtools)))

# paths ---
path_dgelist_obj <- "data/dgelist_obj.rds" # ASSIGN_PATH_DGELIST_OBJ
path_results <- "results/dge-edger" # ASSIGN_PATH_RESULTS

if(!dir.exists(path_results)) dir.create(path_results, recursive = T)

# global parameters ---
metadata_features_col <- c("StageGenotype") # ASSIGN_METADATA_FEATURES_COL # used for grouping in plots, character vector

# Differential gene expression analysis ---
de_method = "qlf" # qlf, lrt, exact
de_robust.dispersion = TRUE
de_return.y = FALSE
de_formula = "~0+StageGenotype"
de_cf = NULL
de_contrast  = list("ESC_3BKO_vs_WT"=c("StageGenotype","ESC_3BKO","ESC_WT")
                    ,"EpiLC_3BKO_vs_WT"=c("StageGenotype","EpiLC_3BKO","EpiLC_WT")
                    ,"ME24h_3BKO_vs_WT"=c("StageGenotype","ME24h_3BKO","ME24h_WT")
                    ,"ME48h_3BKO_vs_WT"=c("StageGenotype","ME48h_3BKO","ME48h_WT"))

de_anovalike = FALSE
de_analysis  = "timecourse_contrasts" # dge

fdrTh  = 0.05
fcTh   = 1
save_table = TRUE
save_excel = FALSE

# Clustering ---
expression.unit = "logRPKM"
de_clustering_method = "kmeans" # hclust, kmeans, pam, manual. none = no clustering
n_clust = 4
de_clustering_scaledata = TRUE
save_clustering_res = TRUE

# hclust
get_hclust_distance = "euclidean"
get_hclust_method   = "complete"
get_hclust_cut_tree = "static" # static, dynamic
get_hclust_minClusterSize = 50
# kmeans
kmeans_nclust_method = "manual" # silhouette
# pam
pam_metric = "euclidean"

# Visualization ---
# heatmaps
hm_annotate_columns = TRUE
hm_use_raster = FALSE
hm_show_row_dend = FALSE
hm_show_column_dend = TRUE
hm_cluster_rows = TRUE
hm_cluster_columns = FALSE
hm_show_row_names = FALSE

do_full_heatmap = TRUE # do heatmap with all samples 
do_contrast_heatmap = TRUE # do heatmap with contrast-specific samples

# violin plots
do_contrast_vln = TRUE # do violin plots with contrast-specific samples
save_contrast_vln = FALSE # save violin plots with contrast-specific samples
save_contrast_vln_arranged = TRUE # violin plots arranged
contrast_vln_arranged_ncol = NULL 

contrast_vln_top = 5 # top genes to label
contrast_vln_label_selection = NULL # genes to label in violin plot

# cluster expression plots 
do_cluster_expression_plots   = FALSE
cluster_expression_plot_type  = "summary_median"
cluster_expression_facet_ncol = 2
cluster_expression_scale = TRUE

# palette ---
set_feature_palettes = T

if(set_feature_palettes) {
  pal_default <- utilsRtools::get_palette_features(metadata_features_col
                                                   , pals = NULL)
} else {
  pal_default <- ggsci::pal_d3()(10)  
}
# alternatively, set pal_default to manually defined palettes
# pal_default <- ...
pal_1 <- RColorBrewer::brewer.pal(9, "Oranges")[c(5:8)]
names(pal_1) <-c("ESC_3BKO","EpiLC_3BKO","ME24h_3BKO","ME48h_3BKO")
pal_2 <- RColorBrewer::brewer.pal(9, "Greys")[c(5:8)]
names(pal_2) <- c("ESC_WT","EpiLC_WT","ME24h_WT","ME48h_WT")

pal_default[["StageGenotype"]] <- c(pal_1,pal_2)

# read object ---
dgelist_obj <- readRDS(path_dgelist_obj)

s33d = 1991
# 1. Differential expression analysis ----
if(!is.null(de_contrast) && grepl("~0+",de_formula)) {
  # test defined contrasts (parametrization: ~ 0 + group ...) ---
  if(!is.list(de_contrast)) de_contrast <- list("de_contrast"=de_contrast)
  de <- lapply(de_contrast, function(x) {
    RNAseqRtools::calculateDiffExprEdgeR(y = dgelist_obj
                                         , design = de_formula
                                         , cf = de_cf
                                         , method = de_method
                                         , contrast = x
                                         , robust.dispersion = de_robust.dispersion
                                         , return.y = de_return.y
                                         , anovalike = de_anovalike)
  })
} else {
  # test model coefficients (parametrization: ~ group ... )---
  if(!is.list(de_cf)) de_cf <- list("de_cf"=de_cf)
  de <- lapply(de_cf, function(x) {
    RNAseqRtools::calculateDiffExprEdgeR(y = dgelist_obj
                                         , design = de_formula
                                         , cf = x
                                         , method = de_method
                                         , contrast = de_contrast
                                         , robust.dispersion = de_robust.dispersion
                                         , return.y = de_return.y
                                         , anovalike = de_anovalike)
  })
}

# Get differentially expressed genes ---
de <- lapply(de, function(x) {
  x$sig <- tryCatch(RNAseqRtools::getDEgenes(x
                                             , fdrTh = fdrTh
                                             , fcTh = fcTh)
                    , error = function(e) {message(e);return(NA)})
  return(x)
})

# Save results ---
if(save_table) {
  for(i in names(de)) {
    outfile <- paste0(path_results,"/dge-edger.",de_method,".",i,".table.txt.gz")
    message(" -- writing to: ", outfile)
    write.table(de[[i]]$table
                , file = gzfile(outfile)
                , row.names = F
                , col.names = T
                , sep = "\t"
                , quote = F)
  }
}

na_idx <- lapply(de, function(x) !is.data.frame(x$sig) && is.na(x$sig))
na_idx <- which(unlist(na_idx))
if(length(na_idx)!=0) de <- de[-na_idx] # filter non significant results

for(i in names(de)) {
  outfile <- paste0(path_results,"/dge-edger.",de_method,".",i,".fcTh_",fcTh,".fdrTh_",fdrTh,".txt.gz")
  message(" -- writing to: ", outfile)
  write.table(de[[i]]$sig
              , file = gzfile(outfile)
              , row.names = F
              , col.names = T
              , sep = "\t"
              , quote = F)
}

if(save_excel) {
  outfile <- paste0(path_results,"/dge-edger.",de_method,".",de_analysis,".fcTh_",fcTh,".fdrTh_",fdrTh,".xlsx")
  RNAseqRtools::saveXLSresEdgeR(res = lapply(de,"[[","sig")
                                , outfile = outfile
                                , force = T)
  
}

# 2. Clustering ----
de_genes <- unique(unlist(lapply(de, function(x) x$sig$gene_name)))
m_de <- dgelist_obj[[expression.unit]][de_genes,]

if(de_clustering_method=="hclust" && !is.null(n_clust)) {
  # perform hierarchical clustering 
  m_de_hclust <- RNAseqRtools::get_hclust(m_de
                                          , cut_tree_k = n_clust
                                          , distance   = get_hclust_distance
                                          , method     = get_hclust_method
                                          , cut_tree   = get_hclust_cut_tree
                                          , minClusterSize = get_hclust_minClusterSize)
  hm_split <- m_de_hclust$cl
  
} else if(de_clustering_method=="kmeans") {
  # perform kmeans clustering 
  if(kmeans_nclust_method!="manual" || is.null(n_clust)) {
    n_clust <- RNAseqRtools::setKMClusters(m_de
                                           , method=kmeans_nclust_method
                                           , ret=T)
  }
  m_de_kmeans <- RNAseqRtools::getKMeans(m_de
                                         , centers = n_clust
                                         , scaledata = de_clustering_scaledata)
  hm_split <- m_de_kmeans$cluster
} else if(de_clustering_method=="pam" && !is.null(n_clust)) {
  # perform pam clustering
  if(is.null(n_clust)) {
    n_clust <- RNAseqRtools::setPAMCnumber(m_de
                                           , metric = pam_metric
                                           , scaledata = de_clustering_scaledata
                                           , ret=T)
  }
  set.seed(s33d)
  m_de_pam <- RNAseqRtools::getPAM(m_de
                                   , k = n_clust
                                   , metric = pam_metric
                                   , scaledata = de_clustering_scaledata)
  hm_split <- m_de_pam$cluster$clustering
} else if(de_clustering_method=="manual") {
  # fill section with rules for manual clustering
  # hm_split <- ...
  
} else {
  # no clustering
  hm_split <- NULL
}

if(de_clustering_method!="none" && save_clustering_res) {
  clustering_res <- reshape2::melt(hm_split, value.name="cluster")
  clustering_res$gene_name <- rownames(clustering_res)
  clustering_res <- clustering_res[order(clustering_res$cluster),]
  
  outfile <- paste0(path_results,"/dge-edger.",de_method,".",de_analysis,".fcTh_",fcTh,".fdrTh_",fdrTh,".clustering_",de_clustering_method,".txt.gz")
  message(" -- writing to: ", outfile)
  write.table(clustering_res
              , file = gzfile(outfile)
              , row.names = F
              , col.names = T
              , sep = "\t"
              , quote = F)
}
# 3. Visualization ----
# heatmaps ---
if(do_full_heatmap) {
  if(hm_annotate_columns & !is.null(metadata_features_col)) {
    hm_annotDF  <- dgelist_obj$samples[,metadata_features_col, drop=F]
    hm_annotCol <- list(pal_default[[metadata_features_col]])
    names(hm_annotCol) <- metadata_features_col
    hm_annotDF <- hm_annotDF[order(hm_annotDF[[metadata_features_col]]),,drop=F]
    cidx <- rownames(hm_annotDF)
  } else {
    hm_annotDF  <- NULL
    hm_annotCol <- NULL
    cidx        <- colnames(m_de)
  }
  if(!is.null(hm_split)) {
    gidx <- names(hm_split)
  } else {
    gidx <- rownames(m_de)
  }
  set.seed(s33d)
  hm_de_genes <- RNAseqRtools::get_heatmap4(m = m_de[gidx,cidx]
                                            , show_row_dend=hm_show_row_dend
                                            , show_column_dend=hm_show_column_dend
                                            , show_row_names=hm_show_row_names
                                            , cluster_columns=hm_cluster_columns
                                            , cluster_rows=hm_cluster_rows
                                            , myZscale = c(-1.5,0,1.5)
                                            , row_title_gp = gpar(fontsize=8)
                                            , retHm    = T
                                            , split = hm_split
                                            , use_raster = hm_use_raster
                                            , annotDF  = hm_annotDF
                                            , annotCol = hm_annotCol)
  
  outfile <- paste0(path_results,"/dge-edger.",de_method,".",de_analysis,".fcTh_",fcTh,".fdrTh_",fdrTh,".heatmap.pdf")
  message(" -- writing to: ", outfile)
  pdf(file = outfile, paper = "a4", width = unit(5,'cm'),height = unit(8,'cm'))
  draw(hm_de_genes, heatmap_legend_side = "bottom")
  dev.off()
}

if(do_contrast_heatmap && (length(de_contrast) > 1 || length(de_cf) > 1)) {
  hm_de <- list()
  for(i in names(de)) {
    if(length(de_contrast) > 1) {
      idx <- which(dgelist_obj$samples[,de_contrast[[i]][1]]%in%de_contrast[[i]][-1])
      m_de_contrast <- dgelist_obj[[expression.unit]][de[[i]]$sig$gene_name,idx]  
    } else if(length(de_cf) > 1) {
      # TO DO ...
    }
    
    if(hm_annotate_columns & !is.null(metadata_features_col)) {
      hm_annotDF  <- dgelist_obj$samples[idx,metadata_features_col, drop=F]
      hm_annotCol <- list(pal_default[[metadata_features_col]])
      names(hm_annotCol) <- metadata_features_col
      hm_annotDF <- hm_annotDF[order(hm_annotDF[[metadata_features_col]]),,drop=F]
      cidx <- rownames(hm_annotDF)
    } else {
      hm_annotDF  <- NULL
      hm_annotCol <- NULL
      cidx        <- colnames(m_de_contrast)
    }
    
    de_genes_i <- list("all"=rownames(de[[i]]$sig)
                       ,"up-regulated"=rownames(subset(de[[i]]$sig, logFC>0))
                       ,"down-regulated"=rownames(subset(de[[i]]$sig, logFC<0)))
    
    hm_split_i <- gsub("\\d+","",names(unlist(de_genes_i[2:3])))
    names(hm_split_i) <- unlist(de_genes_i[2:3])
    
    set.seed(s33d)
    hm_de[[i]] <-  RNAseqRtools::get_heatmap4(m = m_de_contrast[names(hm_split_i),cidx]
                                              , show_row_dend=hm_show_row_dend
                                              , show_column_dend=hm_show_column_dend
                                              , show_row_names=hm_show_row_names
                                              , cluster_columns=hm_cluster_columns
                                              , cluster_rows=hm_cluster_rows
					      , cluster_row_slices = F
                                              , myZscale = c(-1.5,0,1.5)
                                              , row_title_gp = gpar(fontsize=8)
                                              , retHm    = T
                                              , split = hm_split_i
                                              , use_raster = hm_use_raster
                                              , annotDF  = hm_annotDF
                                              , annotCol = hm_annotCol)
    
    outfile <- paste0(path_results,"/dge-edger.",de_method,".",i,".fcTh_",fcTh,".fdrTh_",fdrTh,".heatmap.pdf")
    message(" -- writing to: ", outfile)
    pdf(file = outfile, paper = "a4", width = unit(2.6,'cm'),height = unit(4.8,'cm'))
    draw(hm_de[[i]], heatmap_legend_side = "bottom")
    dev.off()
  }
}

# violin plots 
if(do_contrast_vln) {
  vln_de <- list()
  for(i in names(de)) {
    vln_de[[i]] <- RNAseqRtools::plotDiffExprRes(de[[i]]$table
                                                 , type = "volcano"
                                                 , lfcTh = fcTh
                                                 , pvTh = fdrTh
                                                 , top = contrast_vln_top
                                                 , label_selection = contrast_vln_label_selection) + ggtitle(i)
    if(save_contrast_vln) {
      outfile <- paste0(path_results,"/dge-edger.",de_method,".",i,".fcTh_",fcTh,".fdrTh_",fdrTh,".vln.pdf")
      message(" -- writing to: ", outfile)
      pdf(file = outfile, paper = "a4", width = unit(3.5,'cm'),height = unit(4,'cm'))
      print(vln_de[[i]])
      dev.off()
    }
  }
  if(save_contrast_vln_arranged && length(de) > 1) {
    vln_de_arranged <- ggpubr::ggarrange(plotlist = vln_de, align = "hv", common.legend = T, ncol = contrast_vln_arranged_ncol)
    
    outfile <- paste0(path_results,"/dge-edger.",de_method,".",de_analysis,".fcTh_",fcTh,".fdrTh_",fdrTh,".vln_arranged.pdf")
    message(" -- writing to: ", outfile)
    pdf(file = outfile, paper = "a4", width = unit(6,'cm'),height = unit(8,'cm'))
    print(vln_de_arranged)
    dev.off()
  }
}

# cluster expression 
if(de_clustering_method!="none" && do_cluster_expression_plots) {
  cl_expr <- RNAseqRtools::plot_cluster_expression2(y = dgelist_obj
                                                    , cl = hm_split
                                                    , group_by = metadata_features_col
                                                    , pal = pal_default[[metadata_features_col]]
                                                    , expression.unit = expression.unit
                                                    , plot.type = cluster_expression_plot_type
                                                    , facet_ncol = cluster_expression_facet_ncol
                                                    , scale = cluster_expression_scale)
  
  outfile <- paste0(path_results,"/dge-edger.",de_method,".",de_analysis,".fcTh_",fcTh,".fdrTh_",fdrTh,".clustering_",de_clustering_method,".cluster_expression.pdf")
  message(" -- writing to: ", outfile)
  pdf(file = outfile, paper = "a4", width = unit(3,'cm'),height = unit(6.5,'cm'))
  print(cl_expr)
  dev.off()
}
