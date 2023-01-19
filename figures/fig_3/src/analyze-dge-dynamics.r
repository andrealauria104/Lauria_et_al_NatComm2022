#
# Analyze DGE dynamics
#
# 0. Resources ----
suppressPackageStartupMessages(library(RNAseqRtools))
suppressPackageStartupMessages(library(ggrastr))
# paths 
path_results_dge <- list.files("results/dge-edger"
                               , pattern=".fcTh_0.5.fdrTh_0.05.txt.gz"
                               , full.names = T
                               , recursive = T)
path_dgelist_obj <- "data/dgelist_obj.rds"

# 1. Read data ----
dgelist_obj <- readRDS(path_dgelist_obj)
results_dge <- lapply(path_results_dge, read.delim)
names(results_dge) <- gsub(".*dge-edger.qlf.|.fc.*","",path_results_dge)

# 2. Classify genes ----
degenes <- unique(unlist(lapply(results_dge,"[[","gene_name")))
degenes_stage_fc <- lapply(names(results_dge),function(x) {
 y <- results_dge[[x]][,c("gene_name","logFC")]
 y$contrast <- x
 y
})
degenes_stage_fc <- do.call(rbind,degenes_stage_fc)
degenes_stage_fc <- ddply(degenes_stage_fc
                          , .(gene_name)
                          , mutate
                          , max_contrast = contrast[which.max(abs(logFC))]
                          , max_contrast_direction = sign(logFC[which.max(abs(logFC))]))
degenes_stage_fc_clust <- unique(degenes_stage_fc[,c("gene_name","max_contrast","max_contrast_direction")])
degenes_stage_fc_clust$clust <- NA
degenes_stage_fc_clust$clust[which(degenes_stage_fc_clust$max_contrast=="EpiLC_3B_vs_none" & degenes_stage_fc_clust$max_contrast_direction==1)] <- "early_up"
degenes_stage_fc_clust$clust[which(degenes_stage_fc_clust$max_contrast=="EpiLC_3B_vs_none" & degenes_stage_fc_clust$max_contrast_direction==(-1))] <- "early_down"
degenes_stage_fc_clust$clust[which(degenes_stage_fc_clust$max_contrast=="ME24h_3B_vs_none" & degenes_stage_fc_clust$max_contrast_direction==1)] <- "mid_up"
degenes_stage_fc_clust$clust[which(degenes_stage_fc_clust$max_contrast=="ME24h_3B_vs_none" & degenes_stage_fc_clust$max_contrast_direction==(-1))] <- "mid_down"
degenes_stage_fc_clust$clust[which(degenes_stage_fc_clust$max_contrast=="ME48h_3B_vs_none" & degenes_stage_fc_clust$max_contrast_direction==1)] <- "late_up"
degenes_stage_fc_clust$clust[which(degenes_stage_fc_clust$max_contrast=="ME48h_3B_vs_none" & degenes_stage_fc_clust$max_contrast_direction==(-1))] <- "late_down"

gene_clust <- degenes_stage_fc_clust$clust
names(gene_clust) <- degenes_stage_fc_clust$gene_name
outfile <- "results/gene_clust.rds"
saveRDS(gene_clust, file = outfile)

# 3. Clustering/Visualization ----
# lines
m_rnaseq <- dgelist_obj$logRPKM[names(gene_clust),]
colnames(m_rnaseq) <- gsub("B(126|77)_rep","B\\1_none_rep",colnames(m_rnaseq))
m_rnaseq_scaled <- t(scale(t(m_rnaseq)))
m_rnaseq_m <- melt(m_rnaseq_scaled)
colnames(m_rnaseq_m) <- c("gene","sample","expression")
m_rnaseq_m$cluster_direction <- gene_clust[as.character(m_rnaseq_m$gene)]
m_rnaseq_m$cluster_direction <- gene_clust[as.character(m_rnaseq_m$gene)]
m_rnaseq_m$stage <- sapply(strsplit(as.character(m_rnaseq_m$sample),"\\_"),"[[",1)
m_rnaseq_m$stage <- gsub("MesoT(24|48)","ME\\1h",m_rnaseq_m$stage)
m_rnaseq_m$overexpression <- paste0(sapply(strsplit(as.character(m_rnaseq_m$sample),"\\_"),"[[",2)
                                    ,"-"
                                    ,sapply(strsplit(as.character(m_rnaseq_m$sample),"\\_"),"[[",3))
m_rnaseq_m$overexpression <- gsub("(B126|B77)\\-","",m_rnaseq_m$overexpression)
m_rnaseq_m <- ddply(m_rnaseq_m, .(gene,stage,overexpression,cluster_direction),summarize
                    , average_expression=mean(expression))
m_rnaseq_m$cluster <- gsub("^(\\w+)\\_(\\w+)$","\\1",m_rnaseq_m$cluster_direction)
m_rnaseq_m$cluster <- factor(m_rnaseq_m$cluster, levels=c("early","mid","late"))
m_rnaseq_m$direction <- gsub("^(\\w+)\\_(\\w+)$","\\2",m_rnaseq_m$cluster_direction)
m_rnaseq_m$direction <- factor(m_rnaseq_m$direction, levels=c("up","down"))
levels(m_rnaseq_m$direction) <- c("activated","repressed")
m_rnaseq_m$stage <- factor(m_rnaseq_m$stage, levels=c("EpiLC", "ME24h", "ME48h"))

plines_rnaseq <- ggplot(m_rnaseq_m, aes(x=stage,y=average_expression, group=paste0(overexpression,"_",gene),col=overexpression,alpha=cluster)) +
  facet_grid(direction~cluster)+
  rasterise(geom_line(),dpi=300) + theme_bw() + my_theme_2 + xlab(NULL) + ylab("Average expression [scaled RPKM]") +
  geom_vline(data=subset(m_rnaseq_m, cluster=="early"),aes(xintercept = "EpiLC"), col="black", linetype=2)+
  geom_vline(data=subset(m_rnaseq_m, cluster=="mid"),aes(xintercept = "ME24h"), col="black", linetype=2)+
  geom_vline(data=subset(m_rnaseq_m, cluster=="late"),aes(xintercept = "ME48h"), col="black", linetype=2)+
  scale_color_manual(values = c("none"="#F16913","3B"="#807DBA")) + scale_alpha_manual(values=c(.1,.1,.1)) + ggtitle("[3BKO] DNMT3B-overexpression regulated genes") +
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))

outfile <- "results/expr_dynamics.lines.pdf"
message(" -- saving to: ", outfile)
pdf(file = outfile, paper = 'a4', w=unit(3.8,'cm'), h=unit(2.8,'cm'), useDingbats = F)
print(plines_rnaseq)
dev.off()

m_rnaseq_m_trend <- ddply(m_rnaseq_m, .(stage,overexpression,cluster_direction),summarize
                          , trend_median_expression=median(average_expression)
                          , q1=quantile(average_expression)["25%"]
                          , q3=quantile(average_expression)["75%"])
m_rnaseq_m_trend$cluster <- gsub("_.*","",m_rnaseq_m_trend$cluster_direction)
m_rnaseq_m_trend$direction <- gsub(".*_","",m_rnaseq_m_trend$cluster_direction)
m_rnaseq_m_trend$direction <- gsub("up","Activated",m_rnaseq_m_trend$direction)
m_rnaseq_m_trend$direction <- gsub("down","Repressed",m_rnaseq_m_trend$direction)
m_rnaseq_m_trend$cluster <- factor(m_rnaseq_m_trend$cluster, levels = unique(m_rnaseq_m_trend$cluster)[c(1,3,2)])

plines_rnaseq_trend <- ggplot(m_rnaseq_m_trend, aes(x=stage,y=trend_median_expression, group=overexpression,col=overexpression)) +
  facet_grid(direction~cluster)+
  geom_line(alpha=1,lwd=1) + 
  geom_ribbon(aes(ymin=q1, ymax=q3,fill=overexpression), col="white", alpha=0.3,lwd=0.2) +
  theme_bw() + my_theme_2 + xlab(NULL) + ylab("Average expression [scaled RPKM]") +
  geom_vline(data=subset(m_rnaseq_m_trend, cluster=="early"),aes(xintercept = "EpiLC"), col="black", linetype=2)+
  geom_vline(data=subset(m_rnaseq_m_trend, cluster=="mid"),aes(xintercept = "ME24h"), col="black", linetype=2)+
  geom_vline(data=subset(m_rnaseq_m_trend, cluster=="late"),aes(xintercept = "ME48h"), col="black", linetype=2)+
  scale_color_manual(values = c("none"="#F16913","3B"="#807DBA")) + 
  scale_fill_manual(values = c("none"="#F16913","3B"="#807DBA")) + 
  ggtitle("[3BKO] DNMT3B-overexpression regulated genes") +
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1), plot.title = element_text(face = "plain"))

outfile <- "results/expr_dynamics.lines.trend.pdf"
message(" -- saving to: ", outfile)
pdf(file = outfile, paper = 'a4', w=unit(3.8,'cm'), h=unit(2.8,'cm'), useDingbats = F)
print(plines_rnaseq_trend)
dev.off()
