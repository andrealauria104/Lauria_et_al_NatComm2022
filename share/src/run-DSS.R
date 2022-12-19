#!/usr/bin/env Rscript
# 
#  Differential Methylation Analysis with DSS 
#                                             
# 0. Resources ----
suppressWarnings(suppressMessages(library(docopt)))

'Differential Methylation Analysis with DSS 

Usage:
   run-DSS.R [-m <metadata>] [-c <contrast>] [-d <datadir>] [-n <ncores>] [-p <pipeline>] [-a <analysis>] [-o <outdir>] [--min_cov <min_cov>] [--smoothing <smoothing>] [--equal_disp <equal_disp>] [--smoothing_span <smoothing_span>] [--dml_delta <dml_delta>] [--dml_p_threshold <dml_p_threshold>] [--dmr_delta <dmr_delta>] [--dmr_p_threshold <dmr_p_threshold>] [--dmr_minlen <dmr_minlen>] [--dmr_minCG <dmr_minCG>] [--dmr_dis_merge <dmr_dis_merge>] [--dmr_pct_sig <dmr_pct_sig>] [--show_one_dmr <show_one_dmr>]

Options:
   -m, --metadata Path to metadata.txt (required)
   -c, --contrast String specifing contrast from metadata. Format is METADATAFIELD_COND1_COND2,METADATAFIELD_COND3_COND4 (required)
   -d, --datadir Path to directory containing bismark/bsmap methylation data or to BSobj (as .rds file) [default: .]
   -o, --outdir Path to output directory [default: .]
   -n, --ncores Number of cores [default: 32]
   -p, --pipeline Processing pipeline (bismark/bsmap) [default: bismark]
   -a, --analysis Analysis contrast name [default: NULL]
   --min_cov Minimum cytosine coverage [default: 1]
   --smoothing DSS parameter smoothing [default: TRUE]
   --equal_disp DSS parameter equal_disp [default: FALSE]
   --smoothing_span DSS parameter smoothing_span [default: 500]
   --dml_delta DSS parameter dml_delta [default: 0.1]
   --dml_p_threshold DSS parameter dml_p_threshold [default: 0.001]
   --dmr_delta DSS parameter dmr_delta [default: 0.2]
   --dmr_p_threshold DSS parameter dmr_p_threshold [default: 0.05]
   --dmr_minlen DSS parameter dmr_minlen [default: 50]
   --dmr_minCG DSS parameter dmr_minCG [default: 3]
   --dmr_dis_merge DSS parameter dmr_dis_merge [default: 100]
   --dmr_pct_sig DSS parameter dmr_pct_sig [default: 0.5]
   --show_one_dmr Show randomly selected DMR [default: FALSE]
   --save_all Save all files [default: FALSE]
' -> doc

opts <- docopt(doc)
print(opts)
# required arguments ---
if(is.null(opts$metadata) || is.null(opts$contrast)) {
  message("\n [!] Missing required arguments \n ")
  message(doc)
  quit(save = "no", status = 0, runLast = TRUE)
}

# resources ---
suppressWarnings(suppressMessages(library(BSseqRtools)))
suppressWarnings(suppressMessages(library(DSS)))

# 1. Load data ----
metadata       <- read.delim(as.character(opts$metadata), stringsAsFactors=F)

if(grepl("BSobj.rds",opts$datadir)) {
  BSobj <- readRDS(opts$datadir)
} else {
  BSobj <- prepareDSS(path_files = as.character(opts$datadir)
                    , pipeline   = as.character(opts$pipeline)
                    , min.cov    = as.character(opts$min_cov)
                    , samples_file_pattern = metadata[,1])
}

if(as.logical(opts$save_all)) saveRDS(BSobj, file=paste0(as.character(opts$outdir),"/BSobj.rds"))

# parse contrasts 
my_comparisons <- list()

if(grepl("\\,",opts$contrast)) {
	
	contrast_vec <- unlist(strsplit(as.character(opts$contrast),"\\,"))
	contrast_nm_vec <- unlist(strsplit(as.character(opts$analysis),"\\,"))
	
	if(length(contrast_vec)!=length(contrast_nm_vec)) {
		stop(message("[!] Invalid contrast names, lengths do not match."))
	}

	for(i in 1:length(contrast_vec)) {
		contrast_c <- unlist(strsplit(contrast_vec[i],"\\_"))
    		if(is.null(opts$analysis)) {
      			contrast_nm <- contrast_vec[i]
     		} else {
      			contrast_nm <- contrast_nm_vec[i]
     		}
    	idx1 <- which(metadata[,contrast_c[1]]==contrast_c[2])
    	idx2 <- which(metadata[,contrast_c[1]]==contrast_c[3])

    	my_comparisons[[contrast_nm]] <- list(metadata[idx1,1],metadata[idx2,1])
    }
} else {
	contrast_c <- unlist(strsplit(as.character(opts$contrast),"\\_"))

  	if(is.null(opts$analysis)) {
    	contrast_nm <- as.character(opts$contrast)
    } else {
    	contrast_nm <- as.character(opts$analysis)
  	}

  	idx1 <- which(metadata[,contrast_c[1]]==contrast_c[2])
  	idx2 <- which(metadata[,contrast_c[1]]==contrast_c[3])

  	my_comparisons[[contrast_nm]] <- list(metadata[idx1,1],metadata[idx2,1])
}

# 2. Run DSS ----
dss_analysis <- list()
for(i in names(my_comparisons)) {
  dss_analysis[[i]] <- runTwoGroupDSS(BSobj=BSobj
                                    , group1=my_comparisons[[i]][[1]]
                                    , group2=my_comparisons[[i]][[2]]
                                    , smoothing=as.logical(opts$smoothing)
                                    , equal.disp=as.logical(opts$equal_disp)
                                    , smoothing.span=as.integer(opts$smoothing_span)
                                    , dml.delta=as.numeric(opts$dml_delta)
                                    , dml.p.threshold=as.numeric(opts$dml_p_threshold)
                                    , dmr.delta=as.numeric(opts$dmr_delta)
                                    , dmr.p.threshold=as.numeric(opts$dmr_p_threshold)
                                    , dmr.minlen=as.integer(opts$dmr_minlen)
                                    , dmr.minCG=as.integer(opts$dmr_minCG)
                                    , dmr.dis.merge=as.integer(opts$dmr_dis_merge)
                                    , dmr.pct.sig=as.numeric(opts$dmr_pct_sig)
                                    , n.cores=as.integer(opts$ncores)
                                    , outdir=as.character(opts$outdir)
                                    , analysis=i)
}

outfile <- paste0(as.character(opts$outdir),"/dss_analysis.rds")
message("\n -- saving complete run to: ", outfile)
saveRDS(dss_analysis, file=outfile)

# 3. Show DMR ----
if(as.logical(opts$show_one_dmr)) {
for(nm in names(dss_analysis)) {
  ridx <- sample(2:nrow(dss_analysis[[nm]]$dmr),4, replace=F)
  outfile <- paste0(as.character(opts$outdir),"/dssOneDMR_",nm,".pdf")
  pdf(file=outfile, paper = "a4", w=unit(3.5,'cm'),h=unit(8,'cm'))
  showOneDMR(dss_analysis[[nm]]$dmr[ridx[1],], BSobj) 
  showOneDMR(dss_analysis[[nm]]$dmr[ridx[2],], BSobj)
  showOneDMR(dss_analysis[[nm]]$dmr[ridx[3],], BSobj)
  showOneDMR(dss_analysis[[nm]]$dmr[ridx[4],], BSobj)
  dev.off()
  }
}
