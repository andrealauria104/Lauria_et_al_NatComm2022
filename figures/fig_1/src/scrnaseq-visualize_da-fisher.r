# Visualize differential abundance results
# 0. Resources ----
library(utilsRtools)

path_comparison_table <- "seurat/results_seurat-cluster/da-fisher/comparisons_table.txt"
path_perc <- "seurat/results_seurat-cluster/da-fisher/perc.txt"
path_results <- "seurat/results_seurat-cluster/da-fisher"

pal_default <- list()
pal_default[["Genotype"]] <- c("WT"="#969696","3AKO"="#6BAED6","3BKO"="#FD8D3C")

# 1. Read data ----
comparison_table <- read.delim(path_comparison_table)

perc <- read.delim(path_perc)
colnames(perc) <- gsub("^X","cluster_",colnames(perc))
perc <- reshape2::melt(perc)
colnames(perc)[2] <- "cluster"
perc$cluster <- as.factor(as.integer(gsub("cluster_","",perc$cluster)))

perc_contrasts <- list()
for(i in unique(comparison_table$contrast)) {
  print(i)
  i_1 <- sapply(strsplit(gsub("contrast-","",i),"vs"),"[[",1)
  i_2 <- sapply(strsplit(gsub("contrast-","",i),"vs"),"[[",2)
  perc_contrasts[[i]] <- subset(perc, Genotype%in%c(i_1,i_2))
  perc_contrasts[[i]]$contrast <- i
}
perc_contrasts <- do.call(rbind.data.frame,perc_contrasts)
rownames(perc_contrasts) <- NULL
perc_contrasts$Genotype <- factor(perc_contrasts$Genotype, levels = c("WT","3AKO","3BKO"))
perc_contrasts$contrast <- factor(perc_contrasts$contrast, levels = unique(perc_contrasts$contrast)[c(2,3,1)])
# 2. Plot results ----
p_perc <- ggplot(perc_contrasts, aes(x=cluster, y=value, fill=Genotype)) +geom_col(col="black",position = "dodge") +
  facet_wrap(~contrast, scales = "fixed") + theme_bw() + ylab("% of cells") +
  my_theme_2 + xlab("cluster") + scale_fill_manual(values = pal_default[["Genotype"]]) +
  theme(legend.key.size = unit(4,'mm'), legend.position = "top", axis.title.x = element_blank()) 

p_perc_lineage <- ggplot(subset(perc_contrasts,cluster%in%c(4,5,6,9)), aes(x=cluster, y=value, fill=Genotype)) +geom_col(col="black",position = "dodge") +
  facet_wrap(~contrast, scales = "free") + theme_bw() + ylab("% of cells") +
  my_theme_2 + xlab("cluster") + scale_fill_manual(values = pal_default[["Genotype"]]) +
  theme(legend.key.size = unit(4,'mm'), legend.position = "top", axis.title.x = element_blank()) 

comparison_table_sel_1 <- subset(comparison_table, p.value<=0.05 & alternative!="two.sided")
comparison_table_sel_1$pval_mask <- ""
comparison_table_sel_1$pval_mask[which(comparison_table_sel_1$p.value <= 0.05)] <- "*"
comparison_table_sel_1$pval_mask[which(comparison_table_sel_1$p.value <= 0.01)] <- "**"
comparison_table_sel_1$pval_mask[which(comparison_table_sel_1$p.value <= 0.001)] <- "***"
comparison_table_sel_1$pval_mask[which(comparison_table_sel_1$p.value <= 0.0001)] <- "****"

comparison_table_sel <- subset(comparison_table, alternative=="two.sided")
comparison_table_sel$alternative_sig <- "none"
comparison_table_sel$pval_mask <- ""
comparison_table_sel_1$idx <- paste0(comparison_table_sel_1[,2],"_",comparison_table_sel_1[,1])
comparison_table_sel$idx <- paste0(comparison_table_sel[,2],"_",comparison_table_sel[,1])
match_idx <- comparison_table_sel$idx%in%comparison_table_sel_1$idx
comparison_table_sel$alternative_sig[match_idx] <-comparison_table_sel_1$alternative
comparison_table_sel$pval_mask[match_idx] <-comparison_table_sel_1$pval_mask
comparison_table_sel$alternative_sig <- gsub("g","greater",comparison_table_sel$alternative_sig)
comparison_table_sel$alternative_sig <- gsub("l","less",comparison_table_sel$alternative_sig)
comparison_table_sel$contrast <- factor(comparison_table_sel$contrast, levels = unique(comparison_table_sel$contrast)[c(2,3,1)])

p_table <- ggplot(comparison_table_sel, aes(x=as.factor(cluster),y=odds.ratio, fill=alternative_sig, label=pval_mask))+
  geom_col(col="black",position = "dodge") +
  facet_wrap(~contrast, scales = "fixed") + geom_hline(yintercept = 1, linetype = "dashed",col="#FF33FF",lwd=0.5) +
  theme_bw() + my_theme_2 + theme(legend.key.size = unit(4,'mm')) + geom_text(nudge_y = .01, size=3) +
  xlab("cluster") + scale_fill_manual(values = c("red","blue","white")) + theme(strip.text = element_blank()) +
  guides(fill=guide_legend(title = "Fisher's exact test \n significance")) + scale_y_log10()

# p_table_lineage <- ggplot(subset(comparison_table_sel, cluster%in%c(4,5,6,9)), aes(x=as.factor(cluster),y=odds.ratio, fill=alternative_sig, label=pval_mask))+geom_col(col="black",position = "dodge") +
#   facet_wrap(~contrast, scales = "free_x") + geom_hline(yintercept = 1, linetype = "dashed",col="#FF33FF",lwd=0.5) +
#   theme_bw() + my_theme_2 + theme(legend.key.size = unit(4,'mm')) + geom_text(nudge_y = .01, size=3) +
#   xlab("cluster") + scale_fill_manual(values = c("red","blue","white")) + theme(strip.text = element_blank()) + 
#   guides(fill=guide_legend(title = "Fisher's exact test \n significance"))

# arranged ---
p_perc_table <- ggpubr::ggarrange(p_perc,p_table, nrow=2, align = "v")

outfile <- paste0(path_results,"/p_perc_table.pdf")
message(" -- output file: ", outfile)
pdf(file = outfile, paper = "a4", useDingbats = F, w=unit(3.8,'cm'), height = unit(3.5,'cm'))
print(p_perc_table)
dev.off()

# p_perc_table_lineage <- ggpubr::ggarrange(p_perc_lineage,p_table_lineage, nrow=2, align = "v")
# 
# outfile <- paste0(path_results,"/p_perc_table_lineage.pdf")
# message(" -- output file: ", outfile)
# pdf(file = outfile, paper = "a4", useDingbats = F, w=unit(5,'cm'), height = unit(4,'cm'))
# print(p_perc_table_lineage)
# dev.off()
