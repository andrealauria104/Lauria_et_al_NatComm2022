#
# Build Dnmt3b-dependent TF network
#
# 0. Resources ----
library(igraph)
library(utilsRtools)
# Functions ---
process_chea <- function(path_gmt, human_mouse_map)
{
  chea <- read.gmt.file(path_gmt)
  names(chea) <- gsub("_.*","",names(chea))
  chea <- melt(chea)
  chea <- chea[,c(2,1)]
  colnames(chea) <- c("from","to")
  human_mouse_map <- human_mouse_map[which(human_mouse_map$human_gene_name%in%unique(c(chea$from,chea$to))),]
  chea <- chea[which(chea$from%in%c(human_mouse_map$human_gene_name)),]
  chea <- chea[which(chea$to%in%c(human_mouse_map$human_gene_name)),]
  chea$from <- human_mouse_map[match(chea$from, human_mouse_map$human_gene_name),"mouse_gene_name"]
  chea$to <- human_mouse_map[match(chea$to, human_mouse_map$human_gene_name),"mouse_gene_name"]
  return(chea)
}
read.gmt.file <- function(pathMsigDbFile) 
{
  inputFile <- pathMsigDbFile
  con  <- file(inputFile, open = "r")
  
  c = 1
  pathway.list <- vector(mode="list",length=0)
  while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0)
  {
    myVector <- do.call("rbind",strsplit(oneLine, "\t"))
    t = vector(mode="list",length=1)
    t[[1]] = myVector[3:length(myVector)]
    names(t) = myVector[1]
    pathway.list = c(pathway.list,t)
    c = c+1
  }
  
  close(con)
  return(pathway.list)
}

# Paths ---
path_trust <- "data/trrust_rawdata.mouse.tsv"
path_chea <- list.files("data", pattern=".gmt", recursive = T, full.names = T)
path_chea <- path_chea[-3]
path_orthologs <- "data/ensembl_mart.hsapiens.homolog.one2one_through_gene_name.header_added.gz"

# 1. Read TF-target data ----
orthologs <- read.delim(path_orthologs)
human_mouse_map <- unique(orthologs[,c("external_gene_name","homolog_associated_gene_name")])
colnames(human_mouse_map) <- c("mouse_gene_name","human_gene_name")

# chea ---
tf_target_chea <- lapply(path_chea, process_chea, human_mouse_map=human_mouse_map)
names(tf_target_chea) <- gsub(".*\\/|\\.gmt","",path_chea)
tf_target_chea_coexpr <- do.call(rbind,tf_target_chea[c(1,3)])
rownames(tf_target_chea_coexpr) <- NULL
tf_target_chea_coexpr <- unique(tf_target_chea_coexpr)
tf_target_chea_bind <- do.call(rbind,tf_target_chea[c(2,4,5)])
rownames(tf_target_chea_bind) <- NULL
tf_target_chea_bind <- unique(tf_target_chea_bind)

tf_idx <- intersect(tf_target_chea_coexpr$from,tf_target_chea_bind$from)
target_idx <- list()
for(i in tf_idx) {
  target_idx[[i]] <- intersect(subset(tf_target_chea_coexpr, from==i)$to,subset(tf_target_chea_bind, from==i)$to)
}
tf_target_chea_df <- melt(target_idx)
tf_target_chea_df <- tf_target_chea_df[,c(2,1)]
colnames(tf_target_chea_df) <- c("from","to")


# trust ---
tf_target_trust <- read.delim2(path_trust,header = F)[,1:2]
colnames(tf_target_trust) <- c("from","to")

# unite ---
tf_target <- rbind.data.frame(tf_target_chea_df, tf_target_trust)
tf_target <- tf_target[order(tf_target$from, tf_target$to),]
tf_target <- unique(tf_target)

# save ---
write.table(tf_target
            , file = gzfile("results/tf_target.gz")
            , sep = "\t"
            , quote = F
            , row.names = F
            , col.names = T)
# 2. Read Dnmt3b regulatory data - WGBS, ChIP-seq ----
gene_clusters_split <- readRDS("data/gene_clusters_split.rds")
reg_genes <- names(gene_clusters_split)

# rescued genes 
rescued_genes <- readRDS("data/gene_clust.rds")
rescued_genes_down <- names(rescued_genes[grepl("_down",rescued_genes)])

reg_genes <- intersect(reg_genes, rescued_genes_down)
reg_genes_with_rescue <- intersect(reg_genes, rescued_genes_down)

write.table(reg_genes_with_rescue
            , file = gzfile("results/reg_genes_with_rescue.gz")
            , sep = "\t"
            , quote = F
            , row.names = F
            , col.names = F)
# 3. Read Dnmt3b regulatory data - RNA-seq ----
path_dge  <- "data/rnaseq/dea_edgeR.anovaqlf_pairwiseqlf.ESC_EpiSC_Meso.fcTh_1.fdrTh_0.05.rds"
dge <- readRDS(path_dge)
degenes <- unique(unlist(lapply(dge, "[[","genes")))

# 4. Build regulatory network ----
detf <- intersect(degenes, unique(tf_target$from))
reg_network_tf <- subset(tf_target, from%in%detf & to%in%degenes)
nodes <- data.frame("id"=unique(unlist(reg_network_tf)))
links <- reg_network_tf
net   <- igraph::graph_from_data_frame(d=links, vertices=nodes, directed=T) 

d <- igraph::degree(net,mode="out") 
nm <- names(d)
nm[-which(d>=quantile(d,seq(0,1,.01))["99%"])] <- NA
nm[which(!is.na(nm))]<-NA
e <- igraph::get.edgelist(net)
nm_col <- c(6,3)
names(nm_col) <- c("TRUE","FALSE")
V(net)$color <- nm_col[as.character(V(net)$name%in%detf)]

coul <- igraph::categorical_pal(2)
set.seed(1991)
igraph::plot.igraph(simplify(net)
                    , vertex.size= 2
                    , edge.arrow.size=.1
                    , vertex.label.cex = d/100
                    , vertex.label.color = "#FF8000"
                    , vertex.frame.color=nm_col[as.character(V(net)$name%in%detf)]
                    , edge.color=adjustcolor("darkgrey", alpha.f = .1)
                    , display.isolates=F
                    , vertex.label=nm
                    , layout=layout_in_circle(simplify(net))
                    )

outfile <-"results/Dnmt3b_tf_network_2.full.pdf"
pdf(file=outfile, paper = "a4", width = unit(4,'cm'),height = unit(5,'cm'))
set.seed(1991)
igraph::plot.igraph(simplify(net)
                    , vertex.size= 2
                    , edge.arrow.size=.1
                    , vertex.label.cex = d/100
                    , vertex.label.color = "#FF8000"
                    , vertex.frame.color=nm_col[as.character(V(net)$name%in%detf)]
                    , edge.color=adjustcolor("darkgrey", alpha.f = .1)
                    , display.isolates=F
                    , vertex.label=nm
                    , layout=layout_in_circle(simplify(net))
)
dev.off()

# 5. Network metrics - Centrality ----
net_centr_degree <- igraph::centr_degree(net, mode="out", loops = T, normalized = T)
to_plot <- melt(net_centr_degree$res)
to_plot$TF <- V(net)$name
to_plot <- subset(to_plot, TF%in%detf)
idx <- to_plot$TF[order(to_plot$value,decreasing = F)]
to_plot$TF<- factor(to_plot$TF, levels = idx)
to_plot$label <- as.character(to_plot$TF)
top_idx <- which(to_plot$value>=quantile(to_plot$value,seq(0,1,.01))["98%"])
lab_top_idx <- which(to_plot$TF%in%head(to_plot$label[order(to_plot$value,decreasing = T)],n=5))
to_plot$label[-top_idx] <- ""
to_plot$col <- ""
to_plot$col[top_idx] <- "top 2%"
to_plot$col <- factor(to_plot$col, levels=c("top 2%",""))

pdegree <- ggplot(to_plot, aes(x=TF,y=value,label=label,col=col)) + geom_point(size=.5) + theme_bw() + my_theme_2 +
  theme(panel.grid.major = element_blank(),axis.text.x = element_blank(),axis.line.x = element_blank(),axis.ticks.x = element_blank()) + 
  ylab("Out degree centrality") +xlab("DE TFs") + geom_text_repel(size=3, show.legend = F, force_pull = .8, force = 1.2, segment.linetype=3,max.overlaps = 100) + 
  scale_color_manual(values = c("black","darkgrey")) +
  guides(col=guide_legend(title = NULL))

outfile <-"results/Dnmt3b_tf_network_degree_2.pdf"
pdf(file=outfile, paper = "a4",width=unit(1.8,'cm'),height=unit(2,'cm'),useDingbats = F)
pdegree
dev.off()

# 6. Network metrics - Others ----
# Nodes
N <- length(V(net))
# Degree
k     <- igraph::degree(net,mode="all")
k_in  <- igraph::degree(net,mode="in")
k_out <- igraph::degree(net,mode="out")
# Tot links
L <- gsize(net) # directed: sum(k_in) | sum(k_out) # undirected: sum(k)/2
# Average degree
average_k <- L/N # directed: mean(k_in) | mean(k_out) # undirected: 2*L/N

# Degree distribution
p_k <- degree_distribution(net, cumulative=F, mode="all")
# p_k <- p_k[which(p_k!=0)]

pwl_range <- 100:(length(p_k)-600)
p_k_pwl <- p_k[order(p_k,decreasing = F)]
p_k_fit <- fit_power_law(p_k_pwl)
p_k_df <- data.frame(k=1:length(p_k),p_k=p_k)
p_k_df$p_k_fitted <- NA
p_k_df$p_k_fitted<- p_k_df$k^(-p_k_fit$alpha)

plot_p_k <- ggplot(p_k_df,aes(x=k,y=p_k)) + geom_point(size=1) +
  geom_line(aes(x=k,y=p_k_fitted),linetype="dashed") +
  # geom_vline(xintercept = average_k, linetype="dashed")+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() + my_theme_2 + xlab("k") + ylab("p(k)")
plot_p_k <- plot_p_k + annotation_logticks()

# Degree distribution - CCDF
dd <- igraph::degree_distribution(net, cumulative=T, mode="all")
plot_dd_cdf <- ggplot(data.frame(x=0:max(k),y=1-dd),aes(x=x,y=y)) + geom_point(size=1,col="orange") +
  theme_bw() + my_theme_2 + xlab("Degree") + ylab("CDF") + 
  theme(plot.title = element_text(face = "plain"))

outfile <-"results/Dnmt3b_tf_network_degree.cdf_2.pdf"
pdf(file=outfile, paper = "a4",width=unit(2.5,'cm'),height=unit(2,'cm'),useDingbats = F)
print(plot_dd_cdf)
dev.off()

plot_dd_ccdf <- ggplot(data.frame(x=0:max(k),y=dd),aes(x=x,y=y)) + geom_point(size=1,col="orange") +
  theme_bw() + my_theme_2 + xlab("Degree") + ylab("CCDF") + 
  scale_y_log10() + scale_x_log10() + theme(plot.title = element_text(face = "plain"))

outfile <-"results/Dnmt3b_tf_network_degree.ccdf_2.pdf"
pdf(file=outfile, paper = "a4",width=unit(2.5,'cm'),height=unit(2,'cm'),useDingbats = F)
print(plot_dd)
dev.off()
