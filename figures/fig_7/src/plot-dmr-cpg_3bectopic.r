# 
# Plot DMR gene CpG
#
# 0. Resources ----
suppressPackageStartupMessages(library(BSseqRtools))

# paths ---
dmr_gene <- "Sox2"
path_dmr_gene <- list.files("."
                            , pattern = paste0("DMR_",dmr_gene,".bedGraph.gz")
                            , full.names = T
                            , recursive = T)
metadata <- read.delim("metadata.txt")
# palette ---
pal_default <- list()
pal_1 <- RColorBrewer::brewer.pal(9, "Greys")[5]
names(pal_1) <- c("WT-none")
pal_2 <- rep(RColorBrewer::brewer.pal(9, "Purples")[5],2)
names(pal_2) <- c("3BKO-B77-3B","3BKO-B126-3B")
pal_3 <- rep(RColorBrewer::brewer.pal(9, "Oranges")[5],2)
names(pal_3) <- c("3BKO-B77-none","3BKO-B126-none")

pal_default[["GenotypeOverexpression"]] <- c(pal_1,pal_2,pal_3)

# 1. Read data ----
dmr_gene <- lapply(path_dmr_gene, read.delim, header=F)
names(dmr_gene) <- gsub("\\.\\/|\\.trimmed_bismark_.*","",path_dmr_gene)
dmr_gene <- do.call(rbind,dmr_gene)
colnames(dmr_gene) <- c("chr","start","end","name","score")
dmr_gene$sample <- gsub("\\..*","", rownames(dmr_gene))
rownames(dmr_gene) <- NULL
dmr_gene$sample <- gsub("^WT$","WT-none",dmr_gene$sample)
dmr_gene$sample <- gsub("(B77|B126)$","3BKO-\\1-none",dmr_gene$sample)
dmr_gene$sample <- gsub("(B77|B126)_3b_WT","3BKO-\\1-3B",dmr_gene$sample)
dmr_gene$sample <- factor(dmr_gene$sample, levels = c("WT-none"
                                                      ,"3BKO-B77-none"
                                                      ,"3BKO-B126-none"
                                                      ,"3BKO-B77-3B"
                                                      ,"3BKO-B126-3B"))

p_dmr_gene_row <- ggplot(dmr_gene, aes(x=sample, y=score, col=sample)) +
  geom_boxplot(outlier.shape = NA, show.legend = F, fatten=1) + geom_jitter(size=.4, width = .1, show.legend = F, alpha=0.6) + 
  facet_grid(name~1, scales = "free_x") +
  theme_bw() + my_theme_2 + 
  scale_color_manual(values = pal_default[["GenotypeOverexpression"]]) +
  ylab("CpG methylation [%]") + scale_y_continuous(limits = c(0,100), breaks = c(0,50,100)) +
  # stat_compare_means(comparisons = my_comparisons, paired = T, label = "p.signif",size=2, step.increase = 0) +
  theme(axis.title.x = element_blank()
        , axis.text.x = element_text(angle=90, hjust = 1, vjust = .5, size=6)
        , axis.text.y = element_text(size=6)
        , axis.title.y = element_text(size=6)
        , strip.text = element_text(size=6))

outfile <- "dmr.gene_Sox2.boxplot.row.pdf"
pdf(file=outfile, paper = "a4", width = unit(1.4,'cm'), height = unit(4,'cm'))
print(p_dmr_gene_row)
dev.off()