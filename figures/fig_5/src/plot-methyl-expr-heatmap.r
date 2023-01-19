#
# plot mid/late heatmap
#
#
suppressPackageStartupMessages(library(utilsRtools))

gene_clusters_split <- readRDS("data/gene_clusters_split.rds")
gene_clusters_split <- gene_clusters_split[-which(gene_clusters_split=="early")]
# rnaseq ----
m_rnaseq <- readRDS("data/m_rnaseq.rds")
m_rnaseq <- m_rnaseq[names(gene_clusters_split),]

set.seed(1991)
hm_rnaseq_1 <- RNAseqRtools::get_heatmap4(m = m_rnaseq
                                        , show_row_names = F
                                        , show_row_dend = F
                                        , cluster_columns=F
                                        , cluster_rows=T
                                        , cluster_row_slices=F
                                        , myLegend = "logRPKM"
                                        # , myZscale = c(-1.5,0,1.5)
                                        , split = gene_clusters_split
                                        , use_raster = T
                                        , row_title_gp = gpar(fontsize=8)
                                        , retHm = T
                                        , bm = F
                                        , column_title = "RNA-seq")

gene_set <- c("Gli2","Tubb3","Nnat","Sox2","Id2","Sox1","Olig3")
gene_set_idx <- intersect(gene_set,rownames(m_rnaseq))
r_idx <- ComplexHeatmap::row_order(draw(hm_rnaseq_1))
r_idx_nm <- lapply(r_idx,function(i) rownames(m_rnaseq)[i])
r_idx_nm_at_i <- lapply(r_idx_nm,function(i) which(i%in%gene_set_idx))
r_idx_nm_at <- mapply(function(a,b) a[b], a=r_idx,b=r_idx_nm_at_i)
r_idx_nm_lab <- lapply(r_idx_nm,function(i) i[which(i%in%gene_set_idx)])

anno <- anno_mark(at = unlist(r_idx_nm_at)
                  , labels = unlist(r_idx_nm_lab)
                  , which = "row"
                  , side="left"
                  , labels_gp=gpar(fontsize=8)
                  , labels_rot = 0
                  , link_width = unit(2,'mm'))
set.seed(1991)
hm_rnaseq <- RNAseqRtools::get_heatmap4(m = m_rnaseq
                                        , show_row_names = F
                                        , show_row_dend = F
                                        , cluster_columns=F
                                        , cluster_rows=T
                                        , cluster_row_slices=F
                                        , myLegend = "logRPKM"
                                        , myZscale = c(-1.5,0,1.5)
                                        , split = gene_clusters_split
                                        , use_raster = T
                                        , row_title_gp = gpar(fontsize=8)
                                        , retHm = T
                                        , bm = F
                                        , column_title = "RNA-seq"
                                        , left_annotation = rowAnnotation(mark=anno)
                                        )
# wgbs ----
m_wgbs <- readRDS("data/m_wgbs.rds")
m_wgbs <- m_wgbs[names(gene_clusters_split),]

bPalette <- c("black", "red4","red3", "red2")
ramp <- circlize::colorRamp2(c(0, 60, 80, 100), bPalette)
hm_wgbs <- ComplexHeatmap::Heatmap(m_wgbs
                                   , cluster_rows = F
                                   , cluster_columns = F
                                   , col = ramp
                                   , row_names_gp = gpar(fontsize = 6)
                                   , row_title_gp = gpar(fontsize = 8)
                                   , row_title_rot = 0
                                   , column_names_gp = gpar(fontsize = 8)
                                   , column_title_gp = gpar(fontsize = 8, fontface = "plain")
                                   , heatmap_legend_param = list(title = "Methylation [%]", 
                                                                 title_gp = gpar(fontsize = 8), title_position = "topcenter", 
                                                                 values_gp = gpar(fontsize = 8), legend_direction = "horizontal")
                                   , show_row_names = F
                                   , show_column_names = T
                                   , show_row_dend = F
                                   , show_column_dend = F
                                   , use_raster = T
                                   , raster_device = "tiff"
                                   , raster_quality = 4
                                   , border = T
                                   , na_col = "black"
                                   , split = gene_clusters_split
                                   , column_title = "WGBS")

# regulatory region ----
dmr_dge_common <- read.delim("data/dmr_dge_common.denovo_3b_dmr.txt.gz")
m_reg <- unique(dmr_dge_common[which(dmr_dge_common$gene%in%names(gene_clusters_split) & !is.na(dmr_dge_common$region)),c("gene","region")])
m_reg <- ddply(m_reg, .(gene), summarise, region_av=paste0(region, collapse = ","))
m_reg$region_av[grep("Super",m_reg$region_av)] <- "SuperEnhancers"
m_reg$region_av[grep("Typical",m_reg$region_av)] <- "TypicalEnhancers"
m_reg$region_av[grep("Promoters",m_reg$region_av)] <- "Promoters"
rownames(m_reg) <- m_reg$gene
m_reg <- as.matrix(m_reg[names(gene_clusters_split),-1,drop=F])
hm_reg <- ComplexHeatmap::Heatmap(m_reg
                                  , cluster_rows = F
                                  , cluster_columns = F
                                  , col = pal_d3()(3)
                                  , heatmap_legend_param = list(title = "Regulatory Region", 
                                                                 title_gp = gpar(fontsize = 8), title_position = "topcenter", 
                                                                 values_gp = gpar(fontsize = 8), legend_direction = "horizontal")
                                   , show_row_names = F
                                   , show_column_names = F
                                   , show_row_dend = F
                                   , show_column_dend = F
                                   , use_raster = T
                                   , raster_device = "tiff"
                                   , raster_quality = 4
                                   , border = T
                                   , na_col = "black"
                                   , split = gene_clusters_split)
# save ----
outfile <- "results/methyl_expr_dynamics.heatmap.pdf"
message(" -- saving to: ", outfile)
pdf(file = outfile, paper = "a4", width = unit(5.5,'cm'), height = unit(4.4,'cm'))
draw(hm_rnaseq+hm_wgbs+hm_reg, heatmap_legend_side="bottom")
dev.off()

