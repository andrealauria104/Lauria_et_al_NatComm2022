# 
# Generate DMR heatmap
#
# 0. Resources ----
suppressWarnings(suppressMessages(library(BSseqRtools)))

path_mdata_list <- c("DMR/dssDMR_WT-time-course_3BKOvsWT.mregion.ratio.rds")
path_metadata <- "metadata.txt"

metadata_features_summarize = "StageGenotypeFix"
metadata_features_profile = c("StageFix","Genotype")
metadata_features_profile_precedence <- list(c(2,1,3,4),c(2,1))
names(metadata_features_profile_precedence) <- metadata_features_profile

nm <- gsub("DMR/dssDMR_|.mregion.ratio.rds","",path_mdata_list)
nm <- gsub("_"," ",nm)
names(nm) <- path_mdata_list

exclude_pattern <- list(NULL,NULL,NULL,"B126|B77","WT")
names(exclude_pattern) <- path_mdata_list
pdf_widths <- c(rep(2.8,3),1.8,1.8)
names(pdf_widths)  <- path_mdata_list
pdf_heights <- c(rep(3.95,5))
names(pdf_heights)  <- path_mdata_list
nclust <- list("DMR/dssDMR_WT-time-course_3BKOvsWT.coverage_filt.WT_dmrclust_denovo.labels")
names(nclust) <- path_mdata_list

# palette
bPalette <- c("black", "red4","red3", "red2")
ramp <- circlize::colorRamp2(c(0, 60, 80, 100), bPalette)

pal_default <- list()
pal_1 <- RColorBrewer::brewer.pal(9, "Oranges")[c(5:8)]
names(pal_1) <- c("ESC-3BKO","EpiLC-3BKO","ME24h-3BKO","ME48h-3BKO")
pal_2 <- RColorBrewer::brewer.pal(9, "Greys")[c(5:8)]
names(pal_2) <- c("ESC-WT","EpiLC-WT","ME24h-WT","ME48h-WT")

pal_default[["StageGenotypeFix"]] <- c(pal_1,pal_2)

dmrclust <- list()
for(path_mdata in path_mdata_list) {
  # 1. Read data ----
  mdata <- readRDS(path_mdata)
  colnames(mdata) <- gsub("_","-",colnames(mdata)) 
  colnames(mdata) <- gsub("-Rep"," #",colnames(mdata)) 
  colnames(mdata) <- gsub("WGBS-","",colnames(mdata))
  colnames(mdata) <- gsub("MESO","ME",colnames(mdata))
  colnames(mdata) <- gsub("EpiSC","EpiLC",colnames(mdata))
  
  if(!is.null(exclude_pattern[[path_mdata]])) {
    mdata <- mdata[,grep(exclude_pattern[[path_mdata]],colnames(mdata),invert = T)] 
  }
  
  # 2. Generate heatmap ----
  if(is.integer(nclust[[path_mdata]])) {
    set.seed(91)
    dmrclust[[nm[path_mdata]]] <- kmeans(x = t(scale(t(mdata))), centers = nclust[[path_mdata]])
  } else if(is.character(nclust[[path_mdata]])) {
    dmrlabels <- read.delim(nclust[[path_mdata]],header = F)
    dmrclust[[nm[path_mdata]]] <- list()
    dmrclust[[nm[path_mdata]]][["cluster"]] <- dmrlabels[,2]
    names(dmrclust[[nm[path_mdata]]][["cluster"]]) <- dmrlabels[,1]
    mdata <- mdata[which(rownames(mdata)%in%names(dmrclust[[nm[path_mdata]]][["cluster"]])),]
  } else{
    dmrclust[[nm[path_mdata]]] <- NULL
  }
  col_idx <- c(grep("WT",colnames(mdata)),grep("\\-B",colnames(mdata)))
  hm <- ComplexHeatmap::Heatmap(mdata[names(dmrclust[[nm[path_mdata]]]$cluster),col_idx]
                                , cluster_rows = F
                                , cluster_columns = F
                                , split = dmrclust[[nm[path_mdata]]]$cluster
                                , col = ramp
                                , row_names_gp = gpar(fontsize = 8)
                                , row_title_gp = gpar(fontsize = 8)
                                , row_title_rot = 0
                                , column_names_gp = gpar(fontsize = 8)
                                , column_title_gp = gpar(fontsize = 8, fontface = "plain")
                                , heatmap_legend_param = list(title = "Methylation [%]", 
                                                              title_gp = gpar(fontsize = 8), title_position = "topcenter", 
                                                              values_gp = gpar(fontsize = 8), legend_direction = "horizontal")
                                , show_row_names = F
                                , show_row_dend = F
                                , use_raster = T
                                , raster_device = "tiff"
                                , raster_quality = 4
                                , border = T
                                , na_col = "black"
                                , row_title = paste0(as.roman(names(table(dmrclust[[nm[path_mdata]]]$cluster))),"\n [n=",table(dmrclust[[nm[path_mdata]]]$cluster),"]")
                                , column_title = paste0(nm[path_mdata],"\n",nrow(mdata)," DMRs")
  )
  
  # Save ---
  outfile <- gsub(".rds$",".heatmap.col_idx.pdf",path_mdata)
  message("-- saving to: ",outfile)
  pdf(file=outfile, paper = "a4",w=unit(pdf_widths[path_mdata],'cm'), height=unit(pdf_heights[path_mdata],'cm'))
  draw(hm, heatmap_legend_side="bottom")
  dev.off()
  
  outfile <- gsub(".rds$",".dmrclust.col_idx.rds",path_mdata)
  message("-- saving to: ",outfile)
  saveRDS(dmrclust[[nm[path_mdata]]], file = outfile)
}

