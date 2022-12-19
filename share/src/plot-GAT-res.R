#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

if(length(args)<2) {
  message("\n[!] Invalid arguments")
  message("\nUsage: plot-GAT-res.R <path_gat_tsv> <pvalue_th> <plot_title>")
  quit(save = "no", status = 0, runLast = TRUE)
}

suppressWarnings(suppressMessages(require(utilsRtools)))

path_gat <- as.character(args[1])
pvth <- as.numeric(args[2])
plot_title <- as.character(args[3])

if(is.na(plot_title)) plot_title <- NULL

message(" -- reading: ",path_gat)
gat <- read.delim(path_gat,header = T, stringsAsFactors = F)
gat$signif <- "n.s."
gat$signif[which(gat$qvalue<pvth & gat$l2fold>0)] <- paste0("Enriched [Q < ",pvth,"]")
gat$signif[which(gat$qvalue<pvth & gat$l2fold<0)] <- paste0("Depleted [Q < ",pvth,"]")
gat$signif <- factor(gat$signif, levels = c(paste0("Enriched [Q < ",pvth,"]"),paste0("Depleted [Q < ",pvth,"]"),"n.s."))
pal <- c("#FF3333","#3399FF","#E0E0E0")

#gat$annotation <- gsub("\\-.*","",gat$annotation)
if(length(unique(gat$track))==1) {
  ord_idx <- with(gat, order(l2fold, decreasing = T))
  gat$annotation <- factor(gat$annotation, levels=gat$annotation[ord_idx])
  pdf_w <- unit(4,'cm')
  pdf_h <- unit(2,'cm')
} else {
  pdf_w <- unit(length(unique(gat$track))*2,'cm')
  pdf_h <- unit(2,'cm')
}
pgat <- ggplot(gat, aes(x=annotation,y=l2fold,fill=signif)) +
  geom_col(col="black",lwd=0.25) + ggtitle(plot_title) + theme_bw() + my_theme_2 +
  theme(legend.position = "right", legend.key.size = unit(4,'mm'),axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,size = 6)
        ,plot.title = element_text(face = "plain",size=8)) +
  ylab("log2[Obs/Exp]") + guides(fill=guide_legend(title=NULL)) + scale_fill_manual(values = pal) +  xlab(NULL)

if(length(unique(gat$track))>1) pgat <- pgat + facet_wrap(~track,nrow=1)

outfile <- gsub("\\.tsv$",paste0(".enrichment_plot.q_",pvth,".pdf"),path_gat)
message(" -- writing to: ",outfile)
pdf(file = outfile, paper = "a4r",useDingbats = F, width = pdf_w,height = pdf_h)
pgat
dev.off()
