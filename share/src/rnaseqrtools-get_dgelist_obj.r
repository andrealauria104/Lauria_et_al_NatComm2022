#!/usr/bin/env Rscript
# # # # # # # # # # # # # # # # # #
#                                 #
#  RNA-sequencing data analysis   #
#                                 #
# # # # # # # # # # # # # # # # # #
#
# rnaseqrtools-get_dgelist_obj.r
#
# Create data structure object for downstream analysis with edgeR
# 
#
# 0. Resources ----
suppressWarnings(suppressMessages(library(docopt)))

'RNA-sequencing data analysis
 
 rnaseq-get_dgelist_obj.r
 
 Create data structure (DGElist) object for downstream analysis with edgeR

Usage:
   rnaseqrtools-get_dgelist_obj.r [-c <counts> -m <metadata> -g <gencode> -o <outdir> -u <unit> -n <norm> -t <tlen> -e <vexpr> -s <vsample> -x <extended> -i <group> -r <reference>]
              
Options:
   -c, --counts Path to raw counts (comma separated if > 1) (required). 
   -m, --metadata Path to experiment metadata (comma separated if > 1) (required). 
   -g, --gencode Path to gene annotation info (gencode) (optional).
   -o, --outdir Path to output directory [default: .].
   -u, --unit Normalized expression units [default: cpm].   
   -n, --norm Normalization factors method [default: TMM].
   -t, --tlen Path to transcript lenghts (optional).
   -e, --vexpr Filter expression value [default: 1].
   -s, --vsample Filter sample value [default: 3].
   -x, --extended Save extended data [default: TRUE].
   -i, --group Group by condition (optional).
   -r, --reference Reference group by condition (optional).

Author:
   Andrea Lauria' -> doc

opts <- docopt(doc)
# required arguments ---
required_args <- opts[1:2]
if(any(sapply(required_args, is.null))) {
  missing_idx <- sapply(required_args, is.null)
  missing_args <- gsub("--"," ",names(required_args)[missing_idx])
  message("\n[!] Missing required arguments: ",missing_args,"\n")
  message(doc)
  quit(save = "no", status = 0, runLast = TRUE)
}

# check arguments ---
available_unit <- c("cpm","rpkm")
if(!opts$unit%in%available_unit) {
  message("\n[!] Invalid expression units\n")
  message(available_unit)
  message(doc)
  quit(save = "no", status = 0, runLast = TRUE)
}

available_norm <- c("TMM","TMMwsp","RLE","upperquartile","none")
if(!opts$norm%in%available_norm) {
  message("\n[!] Invalid normalization method\n")
  message(available_norm)
  message(doc)
  quit(save = "no", status = 0, runLast = TRUE)
}
# resources ----
suppressWarnings(suppressMessages(library(RNAseqRtools)))

# 1. Parse arguments ----
opts_idx <- c("counts","metadata")
if(!is.null(opts$gencode))  opts_idx <- c(opts_idx,"gencode")
if(!is.null(opts$tlen)) opts_idx <- c(opts_idx,"tlen")

opts[opts_idx] <- lapply(opts[opts_idx], function(x) 
  {
  x <- as.character(x)
  if(any(grepl("\\,",x))) {
    x <- unlist(strsplit(x,"\\,"))
  }
  return(x)
  }
)

# 2. Read data ----
read_data_list <- list()
for(i in opts_idx) {
  if(length(opts[[i]])>1) {
    # integrate multiple runs/datasets
    read_data_list_tmp <- lapply(opts[[i]], function(x) {
      message(" -- reading: ",x)
      read.delim(x
                 , stringsAsFactors = F
                 , header = T
                 , sep = ifelse(grepl("\\.csv$",x),",","\t"))
      
    })
    if(i=="counts" || (i=="gencode" && length(unique(opts[[i]]))!=1)) {
      read_data_list[[i]] <- Reduce(
        function(x, y, ...) merge(x, y, all = TRUE, by = 1),
        read_data_list_tmp
      )
    } else if(i=="metadata" || i=="qcmatrix") {
      read_data_list[[i]] <- do.call(rbind.data.frame, read_data_list_tmp)
    } else if(i=="gencode" && length(unique(opts[[i]]))==1){
      # N.B.: gencode annotation should be the same!
      read_data_list[[i]] <- read_data_list_tmp
    } else if(i=="tlen" && length(unique(opts[[i]]))==1){
      # N.B.: gencode annotation should be the same!
      read_data_list[[i]] <- read_data_list_tmp
    }
  } else {
    message(" -- reading: ",opts[[i]])
    read_data_list[[i]] <- read.delim(opts[[i]]
                                      , stringsAsFactors = F
                                      , header = T
                                      , sep = ifelse(grepl("\\.csv$",opts[[i]]),",","\t"))  
  }
}

rownames(read_data_list$counts) <- read_data_list$counts[,1]
read_data_list$counts <- read_data_list$counts[,-1]
rownames(read_data_list$metadata) <- read_data_list$metadata[,grep("Sample",colnames(read_data_list$metadata),ignore.case = T)]
read_data_list$metadata <- read_data_list$metadata[colnames(read_data_list$counts),]

if(!is.null(opts$gencode)) {
  read_data_list$gene_info <- read.delim(as.character(opts$gencode), stringsAsFactors = F, header = T)
  read_data_list$gene_info <- ddply(read_data_list$gene_info, .(gene_name),summarize
                                    , gene_id=paste0(gene_id, collapse = ";")
                                    , gene_type=paste0(gene_type, collapse = ";"))
  rownames(read_data_list$gene_info) <- read_data_list$gene_info$gene_name
  read_data_list$gene_info <- read_data_list$gene_info[rownames(read_data_list$counts),]
}

if(!is.null(opts$tlen)) {
  tlen_nm <- read_data_list$tlen[,1]
  read_data_list$tlen <- read_data_list$tlen[,2]
  names(read_data_list$tlen) <- tlen_nm
}
# checks ---
if(!all(read_data_list$metadata[,grep("sample",colnames(read_data_list$metadata),ignore.case = T)]%in%colnames(read_data_list$counts))) {
  message("\n[!] Colnames in count matrix do not match metadata samples!\n")
  message(doc)
  quit(save = "no", status = 0, runLast = TRUE)
}
# 3. Prepare dgelist (y) object ---- 
factor_ord_idx <- grep("_factor_ord$",colnames(read_data_list$metadata),value = T)
if(length(factor_ord_idx)!=0) {
  # set factor orders in metadata fields
  for(i in factor_ord_idx) {
    levels_ord <- unique(read_data_list$metadata[,gsub("_factor_ord$","",i)][order(read_data_list$metadata[,i])])
    read_data_list$metadata[,gsub("_factor_ord$","",i)] <- factor(read_data_list$metadata[,gsub("_factor_ord$","",i)]
                                                                  , levels = levels_ord)
  }
}

y <- processRNAseqEdgeR(m = read_data_list$counts
                        , experimental_info = read_data_list$metadata
                        , gene_info = read_data_list$gene_info
                        , group = opts$group
                        , reference = opts$reference
                        , filter.expr.th = as.integer(as.character(opts$vexpr))
                        , filter.sample.th = as.integer(as.character(opts$vsample))
                        , norm.fact.method = as.character(opts$norm)
                        , expression.unit = as.character(opts$unit)
                        , tlen = read_data_list$tlen
                        )

outfile <- paste0(as.character(opts$outdir),"/dgelist_obj.rds")
message(" -- saving object to: ", outfile)
saveRDS(y, file = outfile)

if(as.logical(as.character(opts$extended))) {
  outfile <- paste0(as.character(opts$outdir),"/norm_counts.",tolower(as.character(opts$unit)),".",tolower(as.character(opts$norm)),".txt.gz")
  message(" -- saving normalized expression to: ", outfile)
  write.table(cbind.data.frame("gene_name"=rownames(y[[toupper(as.character(opts$unit))]]),y[[toupper(as.character(opts$unit))]])
              , file = gzfile(outfile)
              , sep = "\t"
              , quote = F
              , row.names = F
              , col.names = T)
}

