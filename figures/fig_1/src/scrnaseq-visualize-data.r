#
# Data visualization ---
# 
# 0. Resources ----
suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(scRNAseqRtools)))

# paths ---
path_seurat_obj <- "seurat/results_seurat-cluster/seurat_obj.rds"
path_cds <- "seurat/results_pseudotime-monocle3/cds.rds"
path_results <- "seurat/results_visualize-data"

if(!dir.exists(path_results)) dir.create(path_results, recursive = T)

metadata_features <- c("seurat_clusters","Genotype","TimeGenotype")
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

pal_default[["seurat_clusters"]] <- ggsci::pal_lancet()(9)[c(1:4,6)]

rename_seurat_obj = FALSE # rename obj when saving if exists
s33d = 1991

seurat_obj <- readRDS(path_seurat_obj)

seurat_obj$Genotype <- factor(seurat_obj$Genotype, levels = c("WT","3AKO","3BKO"))
seurat_obj$GenotypeLabel <- paste0(seurat_obj$Genotype,"\ [n=",table(seurat_obj$Genotype)[seurat_obj$Genotype],"]")
seurat_obj$GenotypeLabel <- factor(seurat_obj$GenotypeLabel, levels = unique(seurat_obj$GenotypeLabel))

cds <- readRDS(path_cds)

cds$Genotype <- factor(cds$Genotype, levels = c("WT","3AKO","3BKO"))
cds$GenotypeLabel <- paste0(cds$Genotype,"\ [n=",table(cds$Genotype)[cds$Genotype],"]")
cds$GenotypeLabel <- factor(cds$GenotypeLabel, levels = unique(cds$GenotypeLabel))

p_cell_size = 0.65
# 1. Visualize ----
pal_dimred = pal_default
feature = "TimeGenotype"
split.by = "GenotypeLabel"
shape.by = "GenotypeFull"
# reductions = c("umap","umap_mnn")
reductions = "umap"

pdimred <- list()
for(reduction in reductions) {
  pdimred[[reduction]] <- DimPlot(seurat_obj, reduction = reduction, group.by = feature, pt.size = .65, split.by = split.by, shape.by = shape.by) +
    theme_bw() + my_theme_2 + theme(panel.grid.major = element_blank(), legend.key.size = unit(4,'mm'), aspect.ratio = 1) +
    guides(shape=guide_legend(nrow=5)) +
    xlab(paste0(gsub("PCA","PC",toupper(reduction))," 1")) + ylab(paste0(gsub("PCA","PC",toupper(reduction))," 2")) + ggtitle(feature) + 
    scale_color_manual(values = pal_dimred[[feature]])
  pdimred[[reduction]] <- ggrastr::rasterise(pdimred[[reduction]],dpi=600)
  outfile <- paste0(path_results, "/dimplot.",reduction,".",feature,".split.pdf")
  message(" -- saving to: ",outfile)
  pdf(file = outfile, paper = "a4r", w = unit(4.5, 'cm'), h = unit(3.5, 'cm'), useDingbats = F)
  print(pdimred[[reduction]])
  dev.off()
}

# clusters ---
feature = "seurat_clusters"
split.by = NULL
shape.by = "GenotypeFull"
# reductions = c("umap","umap_mnn")
reductions = "umap"

pdimred <- list()
for(reduction in reductions) {
  pdimred[[reduction]] <- DimPlot(seurat_obj, reduction = reduction, group.by = feature, pt.size = .65, split.by = split.by, shape.by = shape.by) +
    theme_bw() + my_theme_2 + theme(panel.grid.major = element_blank(), legend.key.size = unit(4,'mm'), aspect.ratio = 1) +
    guides(shape=guide_legend(nrow=5)) +
    xlab(paste0(gsub("PCA","PC",toupper(reduction))," 1")) + ylab(paste0(gsub("PCA","PC",toupper(reduction))," 2")) + ggtitle(feature) + 
    scale_color_manual(values = pal_dimred[[feature]])
  pdimred[[reduction]] <- ggrastr::rasterise(pdimred[[reduction]],dpi=600)
  outfile <- paste0(path_results, "/dimplot.",reduction,".",feature,".pdf")
  message(" -- saving to: ",outfile)
  pdf(file = outfile, paper = "a4r", w = unit(3, 'cm'), h = unit(3, 'cm'), useDingbats = F)
  print(pdimred[[reduction]])
  dev.off()
}

feature = "TimeGenotype"
split.by = NULL
shape.by = "GenotypeFull"
# reductions = c("umap","umap_mnn")
reductions = "umap"

pdimred <- list()
for(reduction in reductions) {
  pdimred[[reduction]] <- DimPlot(seurat_obj, reduction = reduction, group.by = feature, pt.size = .65, split.by = split.by, shape.by = shape.by) +
    theme_bw() + my_theme_2 + theme(panel.grid.major = element_blank(), legend.key.size = unit(4,'mm'), aspect.ratio = 1) +
    guides(shape=guide_legend(nrow=5)) +
    xlab(paste0(gsub("PCA","PC",toupper(reduction))," 1")) + ylab(paste0(gsub("PCA","PC",toupper(reduction))," 2")) + ggtitle(feature) + 
    scale_color_manual(values = pal_dimred[[feature]])
  pdimred[[reduction]] <- ggrastr::rasterise(pdimred[[reduction]],dpi=600)
  outfile <- paste0(path_results, "/dimplot.",reduction,".",feature,".pdf")
  message(" -- saving to: ",outfile)
  pdf(file = outfile, paper = "a4r", w = unit(3, 'cm'), h = unit(3, 'cm'), useDingbats = F)
  print(pdimred[[reduction]])
  dev.off()
}

split.by = "GenotypeFull"
pdimred <- list()
for(reduction in reductions) {
  pdimred[[reduction]] <- DimPlot(seurat_obj, reduction = reduction, group.by = feature, pt.size = .5, split.by = split.by, shape.by = shape.by) +
    theme_bw() + my_theme_2 + theme(panel.grid.major = element_blank(), legend.key.size = unit(4,'mm'), aspect.ratio = 1) +
    guides(shape=guide_legend(nrow=5)) +
    xlab(paste0(gsub("PCA","PC",toupper(reduction))," 1")) + ylab(paste0(gsub("PCA","PC",toupper(reduction))," 2")) + ggtitle(feature) + 
    scale_color_manual(values = pal_dimred[[feature]])
  
  outfile <- paste0(path_results, "/dimplot.split_genotypefull.",reduction,".pdf")
  message(" -- saving to: ",outfile)
  pdf(file = outfile, paper = "a4r", w = unit(8, 'cm'), h = unit(5, 'cm'), useDingbats = F)
  print(pdimred[[reduction]])
  dev.off()
}

pal_dimred = pal_default
feature = "TimeGenotype"
split.by = "TimeGenotype"
shape.by = "GenotypeFull"
# reductions = c("umap","umap_mnn")
reductions = "umap"

pdimred <- list()
for(reduction in reductions) {
  pdimred[[reduction]] <- DimPlot(seurat_obj, reduction = reduction, group.by = feature, pt.size = .5, split.by = split.by, shape.by = shape.by, ncol=3) +
    theme_bw() + my_theme_2 + theme(panel.grid.major = element_blank(), legend.key.size = unit(4,'mm'), aspect.ratio = 1) +
    guides(shape=guide_legend(nrow=5)) +
    xlab(paste0(gsub("PCA","PC",toupper(reduction))," 1")) + ylab(paste0(gsub("PCA","PC",toupper(reduction))," 2")) + ggtitle(feature) + 
    scale_color_manual(values = pal_dimred[[feature]])
  
  outfile <- paste0(path_results, "/dimplot.split_all.",reduction,".pdf")
  message(" -- saving to: ",outfile)
  pdf(file = outfile, paper = "a4", w = unit(8, 'cm'), h = unit(12, 'cm'), useDingbats = F)
  print(pdimred[[reduction]])
  dev.off()
}

# 2. Visualize monocle3 ----
pumap_trajectory_pseudotime <- plot_cells(cds
                                          , color_cells_by = "pseudotime"
                                          , cell_size = p_cell_size
                                          , label_branch_points=T
                                          , label_leaves = F
                                          , label_roots = F
                                          , label_cell_groups = F
                                          , trajectory_graph_segment_size = .5
                                          , label_principal_points = F
) + theme_bw() + 
  my_theme_2 + theme(legend.key.size = unit(4,'mm'),legend.position = "right",aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pumap_trajectory_pseudotime <- ggrastr::rasterise(pumap_trajectory_pseudotime,dpi=600)
outfile <- paste0(path_results,"/dimplot.",reduction,".trajectory_pseudotime.pdf")
message(" -- saving to: ",outfile)
pdf(file = outfile, paper = "a4r", w = unit(2.8, 'cm'), h = unit(2.8, 'cm'), useDingbats = F)
print(pumap_trajectory_pseudotime)
dev.off()

colData(cds)$Branch_1 <- "other"
colData(cds)$Branch_1[colData(cds)$seurat_clusters%in%c(0,3)] <- "Epi_to_ME"
colData(cds)$Branch_2 <- "other"
colData(cds)$Branch_2[colData(cds)$seurat_clusters%in%c(0,1,2)] <- "Epi_to_NE"
p_branch_1 <- plot_cells(cds
           , color_cells_by = "Branch_1"
           , cell_size = 0.3
           , label_branch_points=F
           , label_leaves = F
           , label_roots = F
           , label_cell_groups = F
           , trajectory_graph_segment_size = .3
           , rasterize = T) + 
  theme_bw() + my_theme_2 + 
  # scale_color_manual(values = c("other"="lightgrey","Epi_to_ME"="#3F007D")) + 
  scale_color_manual(values = c("other"="lightgrey","Epi_to_ME"="#003333")) + 
  # facet_grid(~Genotype)+
  theme(legend.key.size = unit(4,'mm'),legend.position = NULL
        ,aspect.ratio = 1
        , panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        , axis.text = element_blank()
        , axis.title = element_blank()
        , axis.ticks.length = unit(0,'mm'))
p_branch_2 <- plot_cells(cds
                         , color_cells_by = "Branch_2"
                         , cell_size = 0.3
                         , label_branch_points=F
                         , label_leaves = F
                         , label_roots = F
                         , label_cell_groups = F
                         , trajectory_graph_segment_size = .3
                         , rasterize = T) + 
  theme_bw() + my_theme_2 + 
  # scale_color_manual(values = c("other"="lightgrey","Epi_to_NE"="#3F007D")) + 
  scale_color_manual(values = c("other"="lightgrey","Epi_to_NE"="#003333")) + 
  # facet_grid(~Genotype)+
  theme(legend.key.size = unit(4,'mm'),legend.position = NULL
        ,aspect.ratio = 1
        , panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        , axis.text = element_blank()
        , axis.title = element_blank()
        , axis.ticks.length = unit(0,'mm'))
p_branch <- ggpubr::ggarrange(p_branch_1,p_branch_2, align = "hv",nrow = 2,common.legend = T,legend = "bottom")

outfile <- paste0(path_results,"/dimplot.",reduction,".trajectory_branch.pdf")
message(" -- saving to: ",outfile)
pdf(file = outfile, paper = "a4r", w = unit(1.5, 'cm'), h = unit(2, 'cm'), useDingbats = F)
print(p_branch)
dev.off()

colData(cds)$Branch_1_pt <- NA
colData(cds)$Branch_1_pt[colData(cds)$seurat_clusters%in%c(0,3)] <- colData(cds)$Pseudotime[colData(cds)$seurat_clusters%in%c(0,3)]
colData(cds)$Branch_2_pt <- NA
colData(cds)$Branch_2_pt[colData(cds)$seurat_clusters%in%c(0,1,2)] <- colData(cds)$Pseudotime[colData(cds)$seurat_clusters%in%c(0,1,2)]
p_branch_1_pt <- plot_cells(cds
                         , color_cells_by = "Branch_1_pt"
                         , cell_size = 0.3
                         , label_branch_points=F
                         , label_leaves = F
                         , label_roots = F
                         , label_cell_groups = F
                         , trajectory_graph_segment_size = .3
                         , rasterize = T) + 
  theme_bw() + my_theme_2 + 
  scale_color_viridis_c(option = "plasma",limits=c(0,9.661588), na.value="lightgrey") +
  guides(col=FALSE) + 
  # scale_color_manual(values = c("other"="lightgrey","Epi_to_ME"="#3F007D")) + 
  # facet_grid(~Genotype)+
  theme(legend.key.size = unit(4,'mm'),legend.position = NULL
        ,aspect.ratio = 1
        , panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        , axis.text = element_blank()
        , axis.title = element_blank()
        , axis.ticks.length = unit(0,'mm'))
p_branch_2_pt <- plot_cells(cds
                         , color_cells_by = "Branch_2_pt"
                         , cell_size = 0.3
                         , label_branch_points=F
                         , label_leaves = F
                         , label_roots = F
                         , label_cell_groups = F
                         , trajectory_graph_segment_size = .3
                         , rasterize = T) + 
  theme_bw() + my_theme_2 +
  scale_color_viridis_c(option = "plasma",limits=c(0,9.661588), na.value="lightgrey") +
  guides(col=FALSE) + 
  # scale_color_manual(values = c("other"="lightgrey","Epi_to_NE"="#3F007D")) + 
  # facet_grid(~Genotype)+
  theme(legend.key.size = unit(4,'mm'),legend.position = NULL
        ,aspect.ratio = 1
        , panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        , axis.text = element_blank()
        , axis.title = element_blank()
        , axis.ticks.length = unit(0,'mm'))
p_branch_pt <- ggpubr::ggarrange(p_branch_1_pt,p_branch_2_pt, align = "hv",nrow = 2,common.legend = T,legend = "bottom")

outfile <- paste0(path_results,"/dimplot.",reduction,".trajectory_branch_pseudotime.pdf")
message(" -- saving to: ",outfile)
pdf(file = outfile, paper = "a4r", w = unit(1.5, 'cm'), h = unit(1.7, 'cm'), useDingbats = F)
print(p_branch_pt)
dev.off()

pumap_trajectory_pseudotime_genotype <- plot_cells(cds
           , color_cells_by = feature
           , cell_size = p_cell_size
           , label_branch_points=F
           , label_leaves = F
           , label_roots = F
           , label_cell_groups = F
           , trajectory_graph_segment_size = .3) + 
  theme_bw() + my_theme_2 + scale_color_manual(values = pal_default[[feature]]) + 
  facet_grid(~Genotype)+
  theme(legend.key.size = unit(4,'mm'),aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pumap_trajectory_pseudotime_genotype <- ggrastr::rasterise(pumap_trajectory_pseudotime_genotype,dpi=600)
outfile <- paste0(path_results,"/dimplot.",reduction,".trajectory_pseudotime_genotype.pdf")
message(" -- saving to: ",outfile)
pdf(file = outfile, paper = "a4r", w = unit(4.5, 'cm'), h = unit(3.5, 'cm'), useDingbats = F)
print(pumap_trajectory_pseudotime_genotype)
dev.off()

# 3. Markers ----
plot_genes_violin_wrapper <- function(seurat_obj, pal, genes
                                      , assay="RNA"
                                      , point_size=0.1
                                      , jitter_width=0.4
                                      , summary_point_size=1.5
                                      , ident="seurat_clusters")
{
  to_plot <- as.data.frame(seurat_obj@assays[[assay]][genes,])
  to_plot <- cbind("gene_name"=rownames(to_plot),to_plot)
  to_plot <- melt(to_plot, id.vars = "gene_name")
  colnames(to_plot)[2:3] <- c("Sample","expression")
  to_plot <- merge(to_plot,seurat_obj@meta.data, by = "Sample")
  
  p <- ggplot(to_plot, aes_string(x=ident,y="expression",fill=ident,col=ident)) +
    facet_wrap(~gene_name, scales = "free") + 
    geom_violin(trim=T,show.legend = F,scale="width") + 
    geom_jitter(size=point_size,width = jitter_width, show.legend = F) +
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) +
    stat_summary(fun = "mean"
                 ,geom = "point"
                 ,color = "black"
                 ,fill = "white"
                 ,show.legend = F
                 ,size=summary_point_size) +
    theme_bw() + my_theme_2 + 
    theme(axis.text = element_text(size=6)
          , axis.title = element_text(size=6)
          , strip.text = element_text(size=6))
  
  return(p)
}

p_1 <- plot_genes_violin_wrapper(seurat_obj
                                 , pal = pal_default[["seurat_clusters"]]
                                 , genes = c("Pou5f1","Fgf5","Lefty1")
                                 , jitter_width=0.3
                                 , summary_point_size=1)
p_2 <- plot_genes_violin_wrapper(seurat_obj
                                 , pal = pal_default[["seurat_clusters"]]
                                 , genes = c("Msx1","T","Gata4")
                                 , jitter_width=0.3
                                 , summary_point_size=1)

p_3 <- plot_genes_violin_wrapper(seurat_obj
                                 , pal = pal_default[["seurat_clusters"]]
                                 , genes = c("Sox2","Sox1","Tubb3")
                                 , jitter_width=0.3
                                 , summary_point_size=1)
p_markers <- ggpubr::ggarrange(p_1,p_2,p_3,ncol = 1 ,align = "hv",common.legend = T)

outfile <- paste0(path_results,"/seurat_clusters_markers.violin.pdf")
message(" -- output file: ", outfile)
pdf(file = outfile, paper = "a4", useDingbats = F, w=unit(3,'cm'), height = unit(3.5,'cm'))
print(p_markers)
dev.off()

