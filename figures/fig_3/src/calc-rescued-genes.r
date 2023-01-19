# 
# compare rescued
# 
# 0. Resources ---
suppressPackageStartupMessages(library(RNAseqRtools))
suppressPackageStartupMessages(library(GeneOverlap))

path_tables_wt_ko <- list.files("results/dge-edger"
                         , pattern = "_3BKO_vs_WT.table.txt.gz"
                         , full.names = T)
path_tables_ko_re <- list.files("results/dge-edger"
                         , pattern = "_3B_vs_none.table.txt.gz"
                         , full.names = T)

# 1. Read data ----
path_tables_wt_ko <- path_tables_wt_ko[c(2,1,3:4)]
tables_wt_ko <- lapply(path_tables_wt_ko, read.delim)
names(tables_wt_ko) <- gsub("dge-edger.qlf.|_3BKO_vs_WT.table.txt.gz","",basename(path_tables_wt_ko))

genes_wt_ko <- lapply(tables_wt_ko, function(x) {
  list("up"=subset(x, logFC>=.5 & FDR<=0.05)$gene_name
       ,"down"=subset(x, logFC<(-.5) & FDR<=0.05)$gene_name)
})

tables_ko_re <- lapply(path_tables_ko_re, read.delim)
names(tables_ko_re) <- gsub("dge-edger.qlf.|_3B_vs_none.table.txt.gz","",basename(path_tables_ko_re))

genes_ko_re <- lapply(tables_ko_re, function(x) {
  list("up"=subset(x, logFC>=.5 & FDR<=0.05)$gene_name
       ,"down"=subset(x, logFC<(-.5) & FDR<=0.05)$gene_name)
})

# 2. Calculate percentage overlap ----
rescue_overlaps <- mapply(function(a,b) {
  res <- data.frame("tot_up"=length(a$up), 
                    "n_rescued_up"=length(intersect(a$up,b$down)),
                    "perc_rescued_up"=100*length(intersect(a$up,b$down))/length(a$up),
                    "tot_down"=length(a$down),
                    "n_rescued_down"=length(intersect(a$down,b$up)),
                    "perc_rescued_down"=100*length(intersect(a$down,b$up))/length(a$down))
  res$n_non_rescued_up <- res$tot_up-res$n_rescued_up
  res$n_non_rescued_down <- res$tot_down-res$n_rescued_down
  return(res)
}, a = genes_wt_ko[2:4], b = genes_ko_re, SIMPLIFY = F)
rescue_overlaps <- do.call(rbind.data.frame,rescue_overlaps)
rescue_overlaps$Stage <- rownames(rescue_overlaps)
rescue_overlaps_m <- rescue_overlaps

rescue_overlaps <- reshape2::melt(rescue_overlaps, id.vars=c("Stage","perc_rescued_down","perc_rescued_up","tot_down","tot_up"))
rescue_overlaps$direction <- gsub(".*(up|down)","\\1",rescue_overlaps$variable)
rescue_overlaps$direction <- factor(rescue_overlaps$direction, levels = unique(rescue_overlaps$direction))
rescue_overlaps$variable <- gsub("^n_|_up|_down","",rescue_overlaps$variable)
rescue_overlaps$variable <- factor(rescue_overlaps$variable, levels = unique(rescue_overlaps$variable)[c(2,1)])
rescue_overlaps$perc <- NA
rescue_overlaps$perc[which(rescue_overlaps$direction=="up")] <- rescue_overlaps$perc_rescued_up[which(rescue_overlaps$direction=="up")]
rescue_overlaps$perc[which(rescue_overlaps$direction=="down")] <- rescue_overlaps$perc_rescued_down[which(rescue_overlaps$direction=="down")]
rescue_overlaps$perc[which(rescue_overlaps$variable=="non_rescued")] <- NA

prescue_overlaps <- ggplot(rescue_overlaps, aes(x=Stage,y=value,fill=variable)) + 
  geom_col(width=0.8,col="black") + ggtitle("3BKO vs WT rescued genes") +
  theme_bw() + my_theme + facet_wrap(~direction) +
  ylab("Number of genes") + 
  scale_fill_manual(values = c("white","purple"), labels=c("n.s.","rescued [FDR < 0.05 & |logFC| > 0.5]")) +
  theme(legend.key.size = unit(4,'mm')) +
  geom_text(aes(label=round(perc,2)),show.legend = F,size=2, nudge_y = 90) +
  guides(fill=guide_legend(title = NULL))
  
outfile <- "results/rescue_overlaps.pdf"
message(" -- saving to: ", outfile)
pdf(file = outfile, paper = "a4r",width = unit(3.2,'cm'), height = unit(2.4,'cm'))
print(prescue_overlaps)
dev.off()

# all stages ---
rescued_all_stages <- data.frame("tot_up"=length(unique(unlist(lapply(genes_wt_ko[2:4], "[[","up"))))
                                 ,"n_rescued_up"=length(intersect(unique(unlist(lapply(genes_wt_ko[2:4], "[[","up")))
                                                                ,unique(unlist(lapply(genes_ko_re, "[[","down")))))
                                 ,"tot_down"=length(unique(unlist(lapply(genes_wt_ko[2:4], "[[","down"))))
                                 ,"n_rescued_down"=length(intersect(unique(unlist(lapply(genes_wt_ko[2:4], "[[","down")))
                                                                  ,unique(unlist(lapply(genes_ko_re, "[[","up"))))))
rescued_all_stages$n_non_rescued_up <- rescued_all_stages$tot_up-rescued_all_stages$n_rescued_up
rescued_all_stages$n_non_rescued_down <- rescued_all_stages$tot_down-rescued_all_stages$n_rescued_down
rescued_all_stages$perc_rescued_up <- round(100*rescued_all_stages$n_rescued_up/rescued_all_stages$tot_up,2)
rescued_all_stages$perc_rescued_down <- round(100*rescued_all_stages$n_rescued_down/rescued_all_stages$tot_down,2)

rescued_all_stages_m <- rescued_all_stages
rescued_all_stages_m$Stage <- "All"

rescued_all_stages <- reshape2::melt(rescued_all_stages, id.vars=c("perc_rescued_down","perc_rescued_up","tot_down","tot_up"))
rescued_all_stages$direction <- gsub(".*(up|down)","\\1",rescued_all_stages$variable)
rescued_all_stages$direction <- factor(rescued_all_stages$direction, levels = unique(rescued_all_stages$direction))
rescued_all_stages$variable <- gsub("^n_|_up|_down","",rescued_all_stages$variable)
rescued_all_stages$variable <- factor(rescued_all_stages$variable, levels = unique(rescued_all_stages$variable)[c(2,1)])
rescued_all_stages$perc <- NA
rescued_all_stages$perc[which(rescued_all_stages$direction=="up")] <- rescued_all_stages$perc_rescued_up[which(rescued_all_stages$direction=="up")]
rescued_all_stages$perc[which(rescued_all_stages$direction=="down")] <- rescued_all_stages$perc_rescued_down[which(rescued_all_stages$direction=="down")]
rescued_all_stages$perc[which(rescued_all_stages$variable=="non_rescued")] <- NA
prescued_all_stages <- ggplot(rescued_all_stages, aes(x=direction,y=value,fill=variable)) + 
  geom_col(width=0.8,col="black") + ggtitle("3BKO vs WT rescued genes") +
  theme_bw() + my_theme + 
  ylab("Number of genes") + 
  scale_fill_manual(values = c("white","purple"), labels=c("n.s.","rescued [FDR < 0.05 & |logFC| > 0.5]")) +
  theme(legend.key.size = unit(4,'mm'), legend.position = "right") +
  geom_text(aes(label=round(perc,2)),show.legend = F,size=2, nudge_y = 90) +
  guides(fill=guide_legend(title = NULL))

outfile <- "results/rescued_all_stages.pdf"
message(" -- saving to: ", outfile)
pdf(file = outfile, paper = "a4r",width = unit(4,'cm'), height = unit(2.4,'cm'))
print(prescued_all_stages)
dev.off()

# full ---
rescue_overlaps_full <- rbind.data.frame(rescue_overlaps_m,rescued_all_stages_m[,colnames(rescue_overlaps_m)])
rescue_overlaps_full <- reshape2::melt(rescue_overlaps_full, id.vars=c("Stage","perc_rescued_down","perc_rescued_up","tot_down","tot_up"))
rescue_overlaps_full$direction <- gsub(".*(up|down)","\\1",rescue_overlaps_full$variable)
rescue_overlaps_full$direction <- factor(rescue_overlaps_full$direction, levels = unique(rescue_overlaps_full$direction))
rescue_overlaps_full$variable <- gsub("^n_|_up|_down","",rescue_overlaps_full$variable)
rescue_overlaps_full$variable <- factor(rescue_overlaps_full$variable, levels = unique(rescue_overlaps_full$variable)[c(2,1)])
rescue_overlaps_full$perc <- NA
rescue_overlaps_full$perc[which(rescue_overlaps_full$direction=="up")] <- rescue_overlaps_full$perc_rescued_up[which(rescue_overlaps_full$direction=="up")]
rescue_overlaps_full$perc[which(rescue_overlaps_full$direction=="down")] <- rescue_overlaps_full$perc_rescued_down[which(rescue_overlaps_full$direction=="down")]
rescue_overlaps_full$perc[which(rescue_overlaps_full$variable=="non_rescued")] <- NA

prescue_overlaps_full <- ggplot(rescue_overlaps_full, aes(x=Stage,y=value,fill=variable)) + 
  geom_col(width=0.8,col="black") + ggtitle("3BKO vs WT rescued genes") +
  theme_bw() + my_theme + facet_wrap(~direction) +
  ylab("Number of genes") + 
  scale_fill_manual(values = c("white","purple"), labels=c("n.s.","rescued [FDR < 0.05 & |logFC| > 0.5]")) +
  theme(legend.key.size = unit(4,'mm'), axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
  geom_text(aes(label=round(perc,2)),show.legend = F,size=2, nudge_y = 90) +
  guides(fill=guide_legend(title = NULL))

outfile <- "results/rescue_overlaps_full.pdf"
message(" -- saving to: ", outfile)
pdf(file = outfile, paper = "a4r",width = unit(3.2,'cm'), height = unit(2.8,'cm'))
print(prescue_overlaps_full)
dev.off()

# 3. Gene Overlaps ----
universe <- unique(c(unlist(lapply(tables_wt_ko, "[[","gene_name"))
              ,unlist(lapply(tables_ko_re, "[[","gene_name"))))
l1 <- unlist(genes_wt_ko,recursive = F)
names(l1) <- paste0("wt_ko.",names(l1))
l1$wt_ko.all_up <- unique(unlist(l1[grep("(EpiLC.*|ME.*)up",names(l1),value = T)]))
l1$wt_ko.all_down <- unique(unlist(l1[grep("(EpiLC.*|ME.*)down",names(l1),value = T)]))
l2 <- unlist(genes_ko_re,recursive = F)
names(l2) <- paste0("ko_re.",names(l2))
l2$ko_re.all_down <- unique(unlist(l2[grep("(EpiLC.*|ME.*)down",names(l2),value = T)]))
l2$ko_re.all_up <- unique(unlist(l2[grep("(EpiLC.*|ME.*)up",names(l2),value = T)]))
gom.obj <- GeneOverlap::newGOM(gsetA         = l2
                               , gsetB       = l1
                               , genome.size = length(universe))

outfile <- "results/stage_rescue.overlaps.mat.rds"
message(" -- saving to: ", outfile)
saveRDS(gom.obj, file = outfile)

# 4. Pie chart ----
rescued_up <- subset(rescued_all_stages,direction=="up")$value
names(rescued_up) <- c("rescued [FDR < 0.05 & |logFC| > 0.5]","n.s.")
ppie_up <- plot_ggpie(data = rescued_up, variable = "rescued", pal = c("purple","white"),nudge_x_val = 0.6, title = "up in 3BKO vs WT") 
rescued_down <- subset(rescued_all_stages,direction=="down")$value
names(rescued_down) <- c("rescued [FDR < 0.05 & |logFC| > 0.5]","n.s.")
ppie_down <- plot_ggpie(data = rescued_down, variable = "rescued", pal = c("purple","white"),nudge_x_val = 0.6, title = "down in 3BKO vs WT") 

outfile <- "results/rescued_all_stages.pie.pdf"
message(" -- saving to: ", outfile)
pdf(file = outfile, paper = "a4",width = unit(4,'cm'), height = unit(4,'cm'))
print(ppie_up)
print(ppie_down)
dev.off()
