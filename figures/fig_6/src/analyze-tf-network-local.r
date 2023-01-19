#
# Analyze TF network locally: TF target enrichment
#
# 0. Resources ----
suppressPackageStartupMessages(library(utilsRtools))

# 1. Read data ----
# full network
tf_target <- read.delim("results/tf_target.gz")

# expressed genes
y <- readRDS("data/rnaseq/processed_counts.rds")
expressed_genes <- rownames(y$counts)

# 3B targets
reg_genes_with_rescue <- read.delim("results/reg_genes_with_rescue.gz", header = F)[,1]

# degenes
path_dge  <- "data/rnaseq/dea_edgeR.anovaqlf_pairwiseqlf.ESC_EpiSC_Meso.fcTh_1.fdrTh_0.05.rds"
dge <- readRDS(path_dge)
degenes <- unique(unlist(lapply(dge, "[[","genes")))

# 2. Test TF target enrichment ----
tf_direct <- unique(subset(tf_target, from %in% reg_genes_with_rescue)$from)
tf_direct_fisher <- list()
for(i in tf_direct) {
  i_1 <- subset(tf_target, from%in%i & to%in%degenes)$to # target & diff expressed
  i_2 <- setdiff(degenes,subset(tf_target, from%in%i & to%in%degenes)$to) # non target & diff expressed
  i_3 <- subset(tf_target, from%in%i & !(to%in%degenes))$to # target & non diff expressed
  i_4 <- setdiff(expressed_genes,unique(c(degenes,subset(tf_target, from%in%i)$to))) # non target & non diff expressed
  m <- matrix(c(length(i_1),length(i_2)
                ,length(i_3),length(i_4)), nrow = 2
              , dimnames = list(c("target","non target"),c("de","non de")))
  tf_direct_fisher[[i]] <- list("test"=fisher.test(m, alternative = "g"),"m"=m)
}

# save table 
tf_comparisons_table <- list()
for(i in names(tf_direct_fisher)) {
  tf_comparisons_table[[i]] <- data.frame("TF"=i
                                       ,"alternative"="g"
                                       , "p.value"=tf_direct_fisher[[i]]$test$p.value
                                       , "odds.ratio"=tf_direct_fisher[[i]]$test$estimate
                                       , row.names = NULL)
}
tf_comparisons_table <- do.call(rbind,tf_comparisons_table)
tf_comparisons_table <- tf_comparisons_table[order(tf_comparisons_table$p.value,decreasing = F),]
write.table(tf_comparisons_table
            , file = gzfile("tf_comparisons_table.gz")
            , sep = "\t"
            , quote = F
            , row.names = F
            , col.names = T)
# 3. Visualize results ----
tf_comparisons_table$TF <- factor(tf_comparisons_table$TF, levels = rev(tf_comparisons_table$TF))

p_tf_comparisons_table <- ggplot(tf_comparisons_table, aes(y=TF,x=-log10(p.value),size=odds.ratio)) + 
  geom_vline(xintercept = -log10(0.05), linetype="dashed", col="red") +
  geom_point() + theme_bw() + my_theme + xlab("-log10 [Pvalue]") + 
  ylab(NULL) + ggtitle("DNMT3B-direct TFs activity\n on target genes") +
  scale_size(range = c(0,2)) + 
  theme(axis.text = element_text(size=6),axis.title = element_text(size=8), plot.title = element_text(size=7),
        legend.position = "right", legend.text = element_text(size = 6))

outfile <-"results/p_tf_comparisons_table.pdf"
pdf(file=outfile, paper = "a4", width = unit(2.2,'cm'),height = unit(2.4,'cm'))
print(p_tf_comparisons_table)
dev.off()

p_tf_comparisons_table_2 <- ggplot(tf_comparisons_table, aes(y=TF,x=-log10(p.value),size=odds.ratio,fill=p.value<=0.05)) + 
  geom_point(shape=21) + theme_bw() + my_theme + xlab("-log10 [Pvalue]") + 
  geom_vline(xintercept = -log10(0.05), linetype="dashed", col="red") +
  ylab(NULL) + ggtitle("DNMT3B-direct TFs activity\n on target genes") +
  scale_size(range = c(0,3)) + scale_fill_manual(values = c("white","black")) +
  theme(axis.text = element_text(size=6),axis.title = element_text(size=6), plot.title = element_text(size=7),
        legend.position = "right", legend.text = element_text(size = 6))

outfile <-"results/p_tf_comparisons_table_2.pdf"
pdf(file=outfile, paper = "a4", width = unit(2.4,'cm'),height = unit(2.5,'cm'))
print(p_tf_comparisons_table_2)
dev.off()

p_tf_comparisons_table_top <- ggplot(head(tf_comparisons_table,n=10), aes(y=TF,x=-log10(p.value),size=odds.ratio,fill=p.value<=0.05)) + 
  geom_point(shape=21) + theme_bw() + my_theme + xlab("-log10 [Pvalue]") + 
  geom_vline(xintercept = -log10(0.05), linetype="dashed", col="red") +
  ylab(NULL) + ggtitle("Top 10 DNMT3B-direct TFs \nactivity on target genes") +
  scale_size(range = c(0,3.5)) + scale_fill_manual(values = c("white","black")) +
  theme(axis.text = element_text(size=7),axis.title = element_text(size=7), plot.title = element_text(size=7),
        legend.position = "right", legend.text = element_text(size = 6))

outfile <-"results/p_tf_comparisons_table_top.pdf"
pdf(file=outfile, paper = "a4", width = unit(2.6,'cm'),height = unit(2.2,'cm'))
print(p_tf_comparisons_table_top)
dev.off()
