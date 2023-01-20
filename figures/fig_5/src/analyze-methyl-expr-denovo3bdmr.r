#
# Methylation/Expression analysis - Dnmt3b target denovo DMRs
#
# 0. Resources ----
suppressPackageStartupMessages(library(BSseqRtools))
suppressPackageStartupMessages(library(RNAseqRtools))

# paths ---
# rnaseq 
path_dge  = "data/rnaseq/dea_edgeR.anovaqlf_pairwiseqlf.ESC_EpiSC_Meso.fcTh_1.fdrTh_0.05.rds"
path_expr = "data/rnaseq/processed_counts.rds"
# wgbs 
path_denovo_3b_dmr_genes <- "data/wgbs/dssDMR_WT-time-course_3BKOvsWT.coverage_filt.WT_dmrclust_denovo_rGREAT.region_gene_map_default.gz"
path_denovo_3b_dmr_signal <- "data/wgbs/dssDMR_WT-time-course_3BKOvsWT.mregion.ratio.average.rds"
path_denovo_3b_dmr_regulatory_regions <- "data/wgbs/dssDMR_WT-time-course_3BKOvsWT.coverage_filt.WT_dmrclust_denovo.nm.Regulatory_Regions.bed"

# 1. Read data ----
# rnaseq
expr <- readRDS(path_expr)
dge  <- readRDS(path_dge)
dge  <- lapply(dge, function(x) {
  x[,-1] <- lapply(x[,-1], as.numeric)
  x
})

# wgbs
denovo_3b_dmr_genes <- read.delim(path_denovo_3b_dmr_genes)
denovo_3b_dmr_signal <- readRDS(path_denovo_3b_dmr_signal)

chr <- gsub("^(.*)\\.(\\d+)\\.(\\d+)$", "\\1", rownames(denovo_3b_dmr_signal))
start <- as.numeric(gsub("^(.*)\\.(\\d+)\\.(\\d+)$", "\\2", rownames(denovo_3b_dmr_signal)))
end <- as.numeric(gsub("^(.*)\\.(\\d+)\\.(\\d+)$", "\\3", rownames(denovo_3b_dmr_signal)))
denovo_3b_dmr_signal <- as.data.frame(denovo_3b_dmr_signal)
denovo_3b_dmr_signal$idx <- paste0(chr,":",start,"-",end)
denovo_3b_dmr_signal_genes <- merge(denovo_3b_dmr_signal, denovo_3b_dmr_genes, by ="idx")

# regulatory regions ---
denovo_3b_dmr_regulatory_regions <- read.delim(path_denovo_3b_dmr_regulatory_regions, header = F)
denovo_3b_dmr_regulatory_regions <- unique(denovo_3b_dmr_regulatory_regions[,4:5])
colnames(denovo_3b_dmr_regulatory_regions) <- c("idx_tmp","region")

chr <- gsub("^(.*)\\.(\\d+)\\.(\\d+)$", "\\1", denovo_3b_dmr_regulatory_regions$idx_tmp)
start <- as.numeric(gsub("^(.*)\\.(\\d+)\\.(\\d+)$", "\\2", denovo_3b_dmr_regulatory_regions$idx_tmp))
end <- as.numeric(gsub("^(.*)\\.(\\d+)\\.(\\d+)$", "\\3", denovo_3b_dmr_regulatory_regions$idx_tmp))
denovo_3b_dmr_regulatory_regions$idx<- paste0(chr,":",start,"-",end)
denovo_3b_dmr_regulatory_regions <- denovo_3b_dmr_regulatory_regions[,-1]

# 2. Map DMR/gene expression ----
contrast_list <- list("3BKOvsWT_ESC"=c("ESC-3BKO","ESC-WT")
                      ,"3BKOvsWT_EpiLC"=c("EpiLC-3BKO","EpiLC-WT")
                      ,"3BKOvsWT_ME24h"=c("ME24h-3BKO","ME24h-WT")
                      ,"3BKOvsWT_ME48h"=c("ME48h-3BKO","ME48h-WT"))
contrast_name_map <- names(contrast_list)
names(contrast_name_map) <- names(dge)
names(dge) <- contrast_name_map[names(dge)]

# retain dmr with de gene
denovo_3b_dmr_signal_genes_list <- lapply(names(contrast_list), function(i) 
  {
  message(" -- processing: ", i)
  x <- unique(denovo_3b_dmr_signal_genes[,c("idx","gene", "distTSS",contrast_list[[i]])])
  x$average_meth_diff_x <- x[,contrast_list[[i]][1]]-x[,contrast_list[[i]][2]]
  x <- merge(x, denovo_3b_dmr_regulatory_regions, by = "idx", all.x = T)
  x <- merge(x, dge[[i]], by.x="gene",by.y="genes", all = F)
  x <- ddply(x, .(gene), mutate
             , average_meth_diff=mean(average_meth_diff_x)
             , average_idx=paste0(idx,collapse=";")
             , average_dist=paste0(distTSS,collapse=";")
             , average_region=paste0(region,collapse=";"))
  return(x)
})
names(denovo_3b_dmr_signal_genes_list) <- names(contrast_list)

# retain all dmr 
denovo_3b_dmr_signal_genes_list_full <- lapply(names(contrast_list), function(i) 
{
  message(" -- processing: ", i)
  x <- unique(denovo_3b_dmr_signal_genes[,c("idx","gene", "distTSS",contrast_list[[i]])])
  x$average_meth_diff_x <- x[,contrast_list[[i]][1]]-x[,contrast_list[[i]][2]]
  x <- merge(x, denovo_3b_dmr_regulatory_regions, by = "idx", all.x = T)
  x <- merge(x, dge[[i]], by.x="gene",by.y="genes", all.x = T)
  x <- ddply(x, .(gene), mutate
             , average_meth_diff=mean(average_meth_diff_x)
             , average_idx=paste0(idx,collapse=";")
             , average_dist=paste0(distTSS,collapse=";")
             , average_region=paste0(region,collapse=";"))
  return(x)
})
names(denovo_3b_dmr_signal_genes_list_full) <- names(contrast_list)

outfile <- "results/dmr_dge.denovo_3b_dmr.rds"
saveRDS(denovo_3b_dmr_signal_genes_list, file = outfile)

# common gene/dmr list
# expression
denovo_3b_dmr_dge_common_l <- lapply(denovo_3b_dmr_signal_genes_list, function(x) 
{
  unique(x[x[,"average_meth_diff"]<=(-20),"idx"])
})
dmr_dge_common <- unique(unlist(denovo_3b_dmr_dge_common_l))
dmr_dge_common <- lapply(denovo_3b_dmr_signal_genes_list, function(x) 
{
  unique(x[x[,"logFC"]>=1 & x[,"idx"]%in%dmr_dge_common,c("gene","idx")])
})
dmr_dge_common <- melt(dmr_dge_common)
colnames(dmr_dge_common)[3] <- "contrast_expression"
dmr_dge_common <- ddply(dmr_dge_common, .(idx,gene), summarize, contrast_expression=paste0(contrast_expression,collapse=";"))

# methylation
denovo_3b_dmr_dge_common_full_l <- lapply(denovo_3b_dmr_signal_genes_list_full, function(x) 
{
  unique(x[x[,"average_meth_diff"]<=(-20),"idx"])
})
dmr_common <- melt(denovo_3b_dmr_dge_common_full_l)
colnames(dmr_common) <- c("idx","contrast_methylation")
dmr_common <- ddply(dmr_common, .(idx), summarize, contrast_methylation=paste0(contrast_methylation,collapse=";"))

# regulatory regions
denovo_3b_dmr_regulatory_regions <- read.delim(path_denovo_3b_dmr_regulatory_regions, header = F)
denovo_3b_dmr_regulatory_regions <- unique(denovo_3b_dmr_regulatory_regions[,4:5])
colnames(denovo_3b_dmr_regulatory_regions) <- c("idx_tmp","region")

chr <- gsub("^(.*)\\.(\\d+)\\.(\\d+)$", "\\1", denovo_3b_dmr_regulatory_regions$idx_tmp)
start <- as.numeric(gsub("^(.*)\\.(\\d+)\\.(\\d+)$", "\\2", denovo_3b_dmr_regulatory_regions$idx_tmp))
end <- as.numeric(gsub("^(.*)\\.(\\d+)\\.(\\d+)$", "\\3", denovo_3b_dmr_regulatory_regions$idx_tmp))
denovo_3b_dmr_regulatory_regions$idx<- paste0(chr,":",start,"-",end)
denovo_3b_dmr_regulatory_regions <- denovo_3b_dmr_regulatory_regions[,-1]

dmr_dge_common <- merge(dmr_dge_common, denovo_3b_dmr_regulatory_regions, by="idx",all.x=T)
dmr_dge_common <- merge(dmr_dge_common, dmr_common, by="idx")

outfile <- "results/dmr_dge_common.denovo_3b_dmr.txt.gz"
write.table(dmr_dge_common
            , file = gzfile(outfile)
            , sep = "\t"
            , quote = F
            , row.names = F
            , col.names = T)
