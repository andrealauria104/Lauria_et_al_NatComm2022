#
# plot dynamics
#
#
suppressPackageStartupMessages(library(utilsRtools))
suppressPackageStartupMessages(library(ggrastr))

gene_clusters_split <- readRDS("data/gene_clusters_split.rds")
# rnaseq
m_rnaseq <- readRDS("data/m_rnaseq.rds")
m_rnaseq_scaled <- t(scale(t(m_rnaseq)))
m_rnaseq_m <- melt(m_rnaseq_scaled)
colnames(m_rnaseq_m) <- c("gene","sample","expression")
m_rnaseq_m$cluster <- gene_clusters_split[as.character(m_rnaseq_m$gene)]
m_rnaseq_m$stage <- gsub("^(\\w+)\\-(.*)","\\1",m_rnaseq_m$sample)
m_rnaseq_m$genotype <- gsub("^(\\w+)\\-(\\w+)(\\s|\\-).*","\\2",m_rnaseq_m$sample)
m_rnaseq_m <- ddply(m_rnaseq_m, .(gene,stage,genotype,cluster),summarize
                    , average_expression=mean(expression))
m_rnaseq_m$stage <- factor(m_rnaseq_m$stage,levels=unique(m_rnaseq_m$stage)[c(2,1,3,4)])

plines_rnaseq <- ggplot(m_rnaseq_m, aes(x=stage,y=average_expression, group=paste0(genotype,"_",gene),col=genotype,alpha=cluster)) +
  facet_grid(~cluster)+
  geom_line() + theme_bw() + my_theme_2 + xlab(NULL) + ylab("Average expression [scaled RPKM]") +
  geom_vline(data=subset(m_rnaseq_m, cluster=="early"),aes(xintercept = "EpiLC"), col="black", linetype=2)+
  geom_vline(data=subset(m_rnaseq_m, cluster=="mid"),aes(xintercept = "ME24h"), col="black", linetype=2)+
  geom_vline(data=subset(m_rnaseq_m, cluster=="late"),aes(xintercept = "ME48h"), col="black", linetype=2)+
  scale_color_manual(values = c("#F16913","#737373")) + scale_alpha_manual(values=c(.3,.3,.3)) + ggtitle("RNA-seq") +
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))

# wgbs
m_wgbs <- readRDS("data/m_wgbs.rds")
m_wgbs_m <- melt(m_wgbs)
colnames(m_wgbs_m) <- c("gene","sample","methylation")
m_wgbs_m$cluster <- gene_clusters_split[as.character(m_wgbs_m$gene)]
m_wgbs_m$stage <- gsub("^(\\w+)\\-(.*)","\\1",m_wgbs_m$sample)
m_wgbs_m$genotype <- gsub("^(\\w+)\\-(\\w+)(\\s|\\-).*","\\2",m_wgbs_m$sample)
m_wgbs_m <- ddply(m_wgbs_m, .(gene,stage,genotype,cluster),summarize
                    , average_methylation=mean(methylation))
m_wgbs_m$stage <- factor(m_wgbs_m$stage,levels=unique(m_wgbs_m$stage)[c(2,1,3,4)])

plines_wgbs <- ggplot(m_wgbs_m, aes(x=stage,y=average_methylation, group=paste0(genotype,"_",gene),col=genotype,alpha=cluster)) +
  facet_grid(~cluster)+
  geom_line() + theme_bw() + my_theme_2 + xlab(NULL) + ylab("Average methylation [%]") +
  geom_vline(data=subset(m_rnaseq_m, cluster=="early"),aes(xintercept = "EpiLC"), col="black", linetype=2)+
  geom_vline(data=subset(m_rnaseq_m, cluster=="mid"),aes(xintercept = "EpiLC"), col="black", linetype=2)+
  geom_vline(data=subset(m_rnaseq_m, cluster=="late"),aes(xintercept = "EpiLC"), col="black", linetype=2)+
  scale_color_manual(values = c("#F16913","#737373")) + scale_alpha_manual(values=c(.3,.3,.3)) + ggtitle("WGBS") +
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))

plines_arranged <- ggpubr::ggarrange(plines_rnaseq,plines_wgbs, nrow = 2, ncol = 1, align = "hv", common.legend = T, legend = "bottom")

outfile <- "results/methyl_expr_dynamics.lines.pdf"
message(" -- saving to: ", outfile)
pdf(file = outfile, paper = 'a4', w=unit(3.8,'cm'), h=unit(4,'cm'), useDingbats = F)
print(plines_arranged)
dev.off()

# trend ---
m_rnaseq_m_trend <- ddply(m_rnaseq_m, .(stage,genotype,cluster),summarize
                          , trend_median_expression=median(average_expression)
                          , q1=quantile(average_expression)["25%"]
                          , q3=quantile(average_expression)["75%"])
m_rnaseq_m_trend$cluster <- gsub("_.*","",m_rnaseq_m_trend$cluster)
m_rnaseq_m_trend$cluster <- factor(m_rnaseq_m_trend$cluster, levels = unique(m_rnaseq_m_trend$cluster))

plines_rnaseq_trend <- ggplot(m_rnaseq_m_trend, aes(x=stage,y=trend_median_expression, group=genotype,col=genotype)) +
  facet_grid(~cluster)+
  geom_line(alpha=1,lwd=1) + 
  geom_ribbon(aes(ymin=q1, ymax=q3,fill=genotype), col="white", alpha=0.3,lwd=0.2) +
  theme_bw() + my_theme_2 + xlab(NULL) + ylab("Average expression [scaled RPKM]") +
  geom_vline(data=subset(m_rnaseq_m_trend, cluster=="early"),aes(xintercept = "EpiLC"), col="black", linetype=2)+
  geom_vline(data=subset(m_rnaseq_m_trend, cluster=="mid"),aes(xintercept = "ME24h"), col="black", linetype=2)+
  geom_vline(data=subset(m_rnaseq_m_trend, cluster=="late"),aes(xintercept = "ME48h"), col="black", linetype=2)+
  scale_color_manual(values = c("#F16913","#737373")) + 
  scale_fill_manual(values = c("#F16913","#737373")) + 
  ggtitle("RNA-seq") +
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1), plot.title = element_text(face = "plain"))

m_wgbs_m_trend <- ddply(m_wgbs_m, .(stage,genotype,cluster),summarize
                          , trend_median_methylation=median(average_methylation)
                          , q1=quantile(average_methylation)["25%"]
                          , q3=quantile(average_methylation)["75%"])
m_wgbs_m_trend$cluster <- gsub("_.*","",m_wgbs_m_trend$cluster)
m_wgbs_m_trend$cluster <- factor(m_wgbs_m_trend$cluster, levels = unique(m_wgbs_m_trend$cluster))

plines_wgbs_trend <- ggplot(m_wgbs_m_trend, aes(x=stage,y=trend_median_methylation, group=genotype,col=genotype)) +
  facet_grid(~cluster)+
  geom_line(alpha=1,lwd=1) + 
  geom_ribbon(aes(ymin=q1, ymax=q3,fill=genotype), col="white", alpha=0.3,lwd=0.2) +
  theme_bw() + my_theme_2 + xlab(NULL) + ylab("Average methylation [%]") +
  geom_vline(data=subset(m_wgbs_m_trend, cluster=="early"),aes(xintercept = "EpiLC"), col="black", linetype=2)+
  geom_vline(data=subset(m_wgbs_m_trend, cluster=="mid"),aes(xintercept = "ME24h"), col="black", linetype=2)+
  geom_vline(data=subset(m_wgbs_m_trend, cluster=="late"),aes(xintercept = "ME48h"), col="black", linetype=2)+
  scale_color_manual(values = c("#F16913","#737373")) + 
  scale_fill_manual(values = c("#F16913","#737373")) + 
  ggtitle("WGBS") +
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1), plot.title = element_text(face = "plain"))

plines_trend_arranged <- ggpubr::ggarrange(plines_rnaseq_trend,plines_wgbs_trend, nrow = 2, ncol = 1, align = "hv", common.legend = T, legend = "bottom")

outfile <- "results/methyl_expr_dynamics.lines.trend.pdf"
message(" -- saving to: ", outfile)
pdf(file = outfile, paper = 'a4', w=unit(3.8,'cm'), h=unit(4,'cm'), useDingbats = F)
print(plines_trend_arranged)
dev.off()
