#
# Methylation/Expression analysis:
# 
# Multiple Factor Analysis (MFA)
#
#
# 0. Resources ----
suppressPackageStartupMessages(library(FactoMineR))
suppressPackageStartupMessages(library(utilsRtools))

# paths ---
path_expr = "data/rnaseq/processed_counts.rds"
path_dge  = "data/rnaseq/dea_edgeR.anovaqlf_pairwiseqlf.ESC_EpiSC_Meso.fcTh_1.fdrTh_0.05.rds"
path_meth = "data/wgbs/dssDMR_WT-time-course_3BKOvsWT.mregion.ratio.rds"

# 1. Read data ----
# rnaseq
expr   <- readRDS(path_expr)
expr$samples$condition_fix <- gsub("(B126|B77)_Rep.*","3BKO_\\1",expr$samples$sample)
m_expr <- as.matrix(edgeR::cpmByGroup(expr, group=expr$samples$condition_fix), normalized.lib.sizes = TRUE)
colnames(m_expr) <- gsub("EpiSC","EpiLC",colnames(m_expr))
colnames(m_expr) <- gsub("Meso_","ME",colnames(m_expr))
colnames(m_expr) <- gsub("T(24|48)","\\1h",colnames(m_expr))
colnames(m_expr) <- gsub("\\_","\\-",colnames(m_expr))

# wgbs
m_meth <- as.matrix(readRDS(path_meth))
colnames(m_meth) <- gsub("WGBS_","",colnames(m_meth))
colnames(m_meth) <- gsub("EpiSC","EpiLC",colnames(m_meth))
colnames(m_meth) <- gsub("MESO","ME",colnames(m_meth))
colnames(m_meth) <-  gsub("(B126|B77)_Rep.*","3BKO_\\1",colnames(m_meth))
colnames(m_meth) <- gsub("T(24|48)","\\1h",colnames(m_meth))
colnames(m_meth) <- gsub("\\_","\\-",colnames(m_meth))

m_expr <- m_expr[,colnames(m_meth)]

dge  <- readRDS(path_dge)
dge  <- lapply(dge, function(x) {
  x[,-1] <- lapply(x[,-1], as.numeric)
  x
})
dge_genes <- unique(unlist(lapply(dge, "[[","genes")))
m_expr <- m_expr[dge_genes,]

# metadata 
metadata <- data.frame("sample"=colnames(m_expr)
                       ,"genotype"=sapply(strsplit(colnames(m_expr),"\\-"),"[[",2)
                       ,"clone"=sapply(strsplit(colnames(m_expr),"\\-"),"[[",3)
                       ,"stage_genotype"=paste0(sapply(strsplit(colnames(m_expr),"\\-"),"[[",1),"-",sapply(strsplit(colnames(m_expr),"\\-"),"[[",2)))
metadata$clone[metadata$genotype=="WT"] <- "WT"
metadata$stage_genotype <- factor(metadata$stage_genotype, levels=unique(metadata$stage_genotype))

# palette 
pal_default <- list()
pal_1 <- RColorBrewer::brewer.pal(9, "Oranges")[c(5:8)]
names(pal_1) <- c("ESC-3BKO","EpiLC-3BKO","ME24h-3BKO","ME48h-3BKO")
pal_2 <- RColorBrewer::brewer.pal(9, "Greys")[c(5:8)]
names(pal_2) <- c("ESC-WT","EpiLC-WT","ME24h-WT","ME48h-WT")

pal_default[["stage_genotype"]] <- c(pal_1,pal_2)

# 2. MFA ----
r.mfa <- FactoMineR::MFA(base = t(rbind(m_expr,m_meth))
                       , group = c(nrow(m_expr), nrow(m_meth))
                       , graph=FALSE)
outfile <- "results/r.mfa.rds"
message(" -- saving to: ", outfile)
saveRDS(r.mfa, file=outfile)
# extract the H and W matrices from the MFA run result
mfa.h <- r.mfa$global.pca$ind$coord
mfa.w <- r.mfa$quanti.var$coord

# 3. Visualize results ----
# create a dataframe with the H matrix and the group label
mfa_df <- as.data.frame(mfa.h)
mfa_df$sample <- rownames(mfa_df)
mfa_df <- merge(mfa_df, metadata, by="sample")

# create the plot
p_mfa <- ggplot(mfa_df, aes(x=Dim.1, y=Dim.2, color=stage_genotype, shape=clone)) +
  geom_point() + ggtitle("MFA [ RNAseq + WGBS ]") +
  scale_color_manual(values = pal_default$stage_genotype) + 
  theme_bw() + my_theme + 
  theme(legend.position = "right", aspect.ratio = 1, legend.key.size = unit(4,'mm')) +
  guides(col=guide_legend(title = "StageGenotype"))

outfile <- "results/mfa.pdf"
message(" -- saving to: ", outfile)
pdf(file = outfile, paper = 'a4', w=unit(3.6,'cm'), h=unit(3.6,'cm'), useDingbats = F)
print(p_mfa)
dev.off()
