#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

if(length(args)<2) {
	message("\n[!] Invalid arguments")
	message("\nUsage: run-rGREAT.R <path_bed> <species> <rule> <distance>")
	quit(save = "no", status = 0, runLast = TRUE)
}

suppressWarnings(suppressMessages(library(rGREAT)))
suppressWarnings(suppressMessages(library(utilsRtools)))

path_bed <- as.character(args[1])
species  <- as.character(args[2])

if(length(args)<3) {
	rule <- "basalPlusExt"
} else {
	rule <- as.character(args[3])
}

if(length(args)<4) {
	adv_span <- 1000.0
    adv_twoDistance <- 1000.0
    adv_oneDistance <- 1000.0
} else {
	adv_span <- as.numeric(args[4])
    adv_twoDistance <- as.numeric(args[4])
    adv_oneDistance <- as.numeric(args[4])
}

outdir <- gsub(".bed","_rGREAT",path_bed)
if(!dir.exists(outdir)) dir.create(outdir,recursive=T)

message("[+] Reading bed file: ", path_bed)
bed <- read.delim(path_bed,header=F)
if(ncol(bed)<4) {
   bed$id <- paste0(bed[,1],".",bed[,2],".",bed[,3])
}
colnames(bed) <- c("chr","start","end","id")

message("[+] Running submitGreatJob ...")
message(" -- species: ", species)
job <- submitGreatJob(bed
	, species=species
	, rule=rule
	, adv_span=adv_span
	, adv_twoDistance=adv_twoDistance
	, adv_oneDistance=adv_oneDistance)

message("[+] Running getEnrichmentTables ...")
tb <- list()
tb_genes <- list()
for(ont in availableOntologies(job)) {
	message(" -- retrieving: ",ont)
	ont_idx <- gsub("[[:space:]]","_",ont)
	tb[[ont_idx]] <- getEnrichmentTables(job, ontology = ont)[[1]]
	outfile <- paste0(outdir,"/",ont_idx,".txt.gz")
	message(" -- writing to: ",outfile)
	write.table(tb[[ont_idx]], file=gzfile(outfile), sep="\t",quote=F, row.names=F,col.names=T)
	
	tb_genes[[ont_idx]] <- getEnrichmentTables(job, download_by = "tsv", ontology = ont)[[1]]
	outfile <- paste0(outdir,"/",ont_idx,".genes.txt.gz")
	message(" -- writing to: ",outfile)
	write.table(tb_genes[[ont_idx]], file=gzfile(outfile), sep="\t",quote=F, row.names=F,col.names=T)
}

message("[+] Get RegionGeneAssociation ...")
res <- plotRegionGeneAssociationGraphs(job)
res <- as.data.frame(res)
colnames(res)[1] <- "chr" 
res$idx <- with(res, paste0(chr,":",start,"-",end))
bed$idx <- with(bed, paste0(chr,":",start,"-",end))
res <- merge(res, bed[,c("id","idx")], by="idx")

res <- res[,c(2:4,1,8,5:7)]
outfile <- paste0(outdir,"/RegionGeneAssociation.txt.gz")
message(" -- writing to: ",outfile)
write.table(res, file=gzfile(outfile), sep="\t",quote=F, row.names=F,col.names=T)


outfile <- paste0(outdir,"/RegionGeneAssociationGraphs.pdf")
pdf(file=outfile, paper="a4r", width=unit(8,'cm'),height=unit(3,'cm'))
plotRegionGeneAssociationGraphs(job)
dev.off()
