# 
# Plot DMR gene CpG
#
# 0. Resources ----
suppressPackageStartupMessages(library(BSseqRtools))

# paths ---
path_dmr_gene <- list.files("."
                            , pattern = "DMR_hyper_CRISPRoff_none.bedGraph.gz"
                            , full.names = T
                            , recursive = T)
metadata <- read.delim("metadata.txt")
# palette ---
pal_default <- list()
pal_2 <- rep(RColorBrewer::brewer.pal(9, "Purples")[4],4)
names(pal_2) <- c("CRISPRoff-gRNAs-E1","CRISPRoff-gRNAs-E2","CRISPRoff-gRNAs-E3","CRISPRoff-gRNAs-E4")
pal_3 <- rep(RColorBrewer::brewer.pal(9, "Oranges")[4],2)
names(pal_3) <- c("none","CRISPRoff")

pal_default[["Condition"]] <- c(pal_2,pal_3)

# 1. Read data ----
dmr_gene <- lapply(path_dmr_gene, read.delim, header=F)
names(dmr_gene) <- gsub("\\.\\/|\\.trimmed_bismark_.*","",path_dmr_gene)
dmr_gene <- do.call(rbind,dmr_gene)
colnames(dmr_gene) <- c("chr","start","end","name","score")
dmr_gene$sample <- gsub("\\..*","", rownames(dmr_gene))
rownames(dmr_gene) <- NULL
dmr_gene <- merge(dmr_gene, metadata, by.x="sample",by.y="Sample")
dmr_gene$sample <- factor(dmr_gene$sample, levels=unique(dmr_gene$sample)[c(7,1,8,2,9,3,10,4,11,5,12,6)])
pal_default$Condition[grep("CRISPRoff-gRNAs",names(pal_default$Condition))] <- "#F0037F"

p_dmr_gene_row <- ggplot(dmr_gene, aes(x=sample, y=score, col=Condition)) +
  geom_boxplot(outlier.shape = NA, show.legend = F, fatten=1) + geom_jitter(size=.4, width = .1, show.legend = F,alpha=0.7) + 
  facet_grid(name~1, scales = "free_x") +
  theme_bw() + my_theme_2 + 
  scale_color_manual(values = pal_default[["Condition"]]) +
  ylab("CpG methylation [%]") + scale_y_continuous(limits = c(0,100), breaks = c(0,50,100)) +
  # stat_compare_means(comparisons = my_comparisons, paired = T, label = "p.signif",size=2, step.increase = 0) +
  theme(axis.title.x = element_blank()
        , axis.text.x = element_text(angle=90, hjust = 1, vjust = .5, size=6)
        , axis.text.y = element_text(size=6)
        , axis.title.y = element_text(size=6)
        , strip.text = element_text(size=6))

outfile <- "DMR_hyper_CRISPRoff_none.boxplot.row_col.pdf"
pdf(file=outfile, paper = "a4", width = unit(1.9,'cm'), height = unit(4,'cm'))
print(p_dmr_gene_row)
dev.off()
