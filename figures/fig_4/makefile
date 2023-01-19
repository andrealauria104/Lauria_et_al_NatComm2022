################
# DMR analysis # ----
################
# Contrasts from metadata, space separated
CONTRAST=StageGenotype_EpiSC-WT_ESC-WT StageGenotype_MESO24h-WT_ESC-WT StageGenotype_MESO48h-WT_ESC-WT
CONTRAST+=StageGenotype_MESO24h-WT_EpiSC-WT StageGenotype_MESO48h-WT_EpiSC-WT StageGenotype_MESO48h-WT_MESO24h-WT
CONTRAST+=StageGenotype_EpiSC-3BKO_ESC-3BKO StageGenotype_MESO24h-3BKO_ESC-3BKO StageGenotype_MESO48h-3BKO_ESC-3BKO
CONTRAST+=StageGenotype_MESO24h-3BKO_EpiSC-3BKO StageGenotype_MESO48h-3BKO_EpiSC-3BKO StageGenotype_MESO48h-3BKO_MESO24h-3BKO
CONTRAST+=StageGenotype_ESC-3BKO_ESC-WT StageGenotype_EpiSC-3BKO_EpiSC-WT StageGenotype_MESO24h-3BKO_MESO24h-WT StageGenotype_MESO48h-3BKO_MESO48h-WT 
# Contrasts in DMR filenames, space separated
CONTRAST_NAME=EpiSCvsESC_WT MESO24hvsESC_WT MESO48hvsESC_WT MESO24hvsEpiSC_WT MESO48hvsEpiSC_WT MESO48hvsMESO24h_WT
CONTRAST_NAME+=EpiSCvsESC_3BKO MESO24hvsESC_3BKO MESO48hvsESC_3BKO MESO24hvsEpiSC_3BKO MESO48hvsEpiSC_3BKO MESO48hvsMESO24h_3BKO
CONTRAST_NAME+=3BKOvsWT_ESC 3BKOvsWT_EpiSC 3BKOvsWT_MESO24h 3BKOvsWT_MESO48h

GENOME=mm10

# GAT
GAT_ANNOTATIONS=data/genomic_regions/Full_Annotations.bed.gz
GAT_WORKSPACE=data/genomic_regions/contigs.no_blacklist_gaps.bed
# ------------- #
# Configuration #
# ------------- #
DATADIR?=data/

CONTRAST?=ASSIGN_IN_MAKEFILE
CONTRAST_NAME?=ASSIGN_IN_MAKEFILE
DESIGN?=
DSS_CORES?=4

.SECONDARY:

.PHONY: all_dss

# DSS parameters ---
SMOOTHING?=TRUE
EQUAL_DISP?=FALSE

DML_DELTA?=0.1
DML_PVALUE?=0.001
DMR_DELTA?=0.2
DMR_PVALUE?=0.05
DMR_MINCG?=5
DMR_MINLEN?=50
DMR_DIS_MERGE?=100
DMR_PCT_SIG?=0.5

GENOME?=mm10

ALL+=$(ALL_DSS)
ALL+=$(ALL_SCREEN)
ALL+=$(ALL_GREAT)

# utilities for make 
comma := ,
empty :=
space := $(empty) $(empty)
# --------------------------------------- #
# Differential Methylation Analysis - DSS #
# --------------------------------------- #
ALL_DSS=dss_analysis.rds
ALL_DSS_TXT=$(addprefix dssDMR_, $(addsuffix _delta_$(DMR_DELTA)_p_$(DMR_PVALUE).$(DMR_MINLEN)_$(DMR_MINCG)_$(DMR_DIS_MERGE)_$(DMR_PCT_SIG).txt, $(CONTRAST_NAME)))
ALL_DSS_BED=$(addprefix dssDMR_, $(addsuffix _delta_$(DMR_DELTA)_p_$(DMR_PVALUE).$(DMR_MINLEN)_$(DMR_MINCG)_$(DMR_DIS_MERGE)_$(DMR_PCT_SIG).bed, $(CONTRAST_NAME)))
ALL_DSS_TRACKS=$(addprefix dssDMR_, $(addsuffix _delta_$(DMR_DELTA)_p_$(DMR_PVALUE).$(DMR_MINLEN)_$(DMR_MINCG)_$(DMR_DIS_MERGE)_$(DMR_PCT_SIG).track.bed, $(CONTRAST_NAME)))
ALL_DSS+=$(ALL_DSS_TXT) $(ALL_DSS_BED) $(ALL_DSS_TRACKS)

all_dss: $(ALL_DSS)
	@echo done

dss_analysis.rds: data/metadata.txt
	conda activate r-environment;\
	bin/run-DSS.R -m $< -d $(DATADIR) \
	-c $(subst ${space},${comma},${CONTRAST}) \
	-a $(subst ${space},${comma},${CONTRAST_NAME}) \
	--smoothing $(SMOOTHING) --equal_disp $(EQUAL_DISP) \
	--ncores $(DSS_CORES) -o .

dssDMR_%.track.bed: dssDMR_%.txt
	awk 'BEGIN{OFS="\t";print "track useScore=1"}(NR>1 && $$1 ~ "chr"){print $$1,$$2,$$3,$$8,sqrt($$8^2)*1000}' $< > $@

dssDMR_%.bed: dssDMR_%.txt
	awk 'BEGIN{OFS="\t"}(NR>1 && $$1 ~ "chr"){print $$1,$$2,$$3,$$8}' $< | sort -k1,1 -k2,2n > $@

dssDMR_%.hypo.bed: dssDMR_%.bed
	awk '$$4<0' $< > $@
dssDMR_%.hyper.bed: dssDMR_%.bed
	awk '$$4>0' $< > $@
dssDMR_%.hypo.names: dssDMR_%.bed
	awk '$$4<0{print $$1":"$$2"-"$$3}' $< > $@
dssDMR_%.hyper.names: dssDMR_%.bed
	awk '$$4>0{print $$1":"$$2"-"$$3}' $< > $@

dssDMR_%_delta_$(DMR_DELTA)_p_$(DMR_PVALUE).$(DMR_MINLEN)_$(DMR_MINCG)_$(DMR_DIS_MERGE)_$(DMR_PCT_SIG).txt: dss_analysis.rds
	conda activate r-environment;\
	bin/call-DMR.R -d $< -o . \
	--dmr_delta $(DMR_DELTA) --dmr_p_threshold $(DMR_PVALUE) \
	--dmr_minlen $(DMR_MINLEN) --dmr_minCG $(DMR_MINCG) \
	--dmr_dis_merge $(DMR_DIS_MERGE) --dmr_pct_sig $(DMR_PCT_SIG)

# -------------- #
# DMR annotation #
# -------------- #

# GAT
GAT_ANNOTATIONS?=Full_Annotations.bed.gz
GAT_NM?=$(shell basename $(GAT_ANNOTATIONS) | sed -e 's/\.bed.*//g')
GAT_WORKSPACE?=contigs.no_blacklist_gaps.bed
GAT_NSAMPLES?=10000
GAT_CORES?=8
GAT_PLOT_TITLE?="DMRs annotation"

%.$(GAT_NM).nsamples_$(GAT_NSAMPLES).gat.tsv: %.bed $(GAT_ANNOTATIONS) $(GAT_WORKSPACE)
	conda activate gat;\
	gat-run.py --segments=$< \
	--annotations=$(word 2, $^) \
	--workspace=$(word 3, $^) \
	--ignore-segment-tracks \
	--num-samples=$(GAT_NSAMPLES) \
	--num-threads=$(GAT_CORES) \
	--log=$*.$(GAT_NM).nsamples_$(GAT_NSAMPLES).gat.log --stdout=$@

%.gat.enrichment_plot.q_0.05.pdf: %.gat.tsv
	conda activate r-environment;\
	bin/plot-GAT-res $< 0.05 $(GAT_PLOT_TITLE)

# ------------------------------ #
# ENCODE SCREEN ccREs annotation #
# ------------------------------ #
ALL_SCREEN=$(addprefix dssDMR_, $(addsuffix _delta_$(DMR_DELTA)_p_$(DMR_PVALUE).$(DMR_MINLEN)_$(DMR_MINCG)_$(DMR_DIS_MERGE)_$(DMR_PCT_SIG).ccREs.bed, $(CONTRAST_NAME)))
ALL_SCREEN+=$(addprefix dssDMR_, $(addsuffix _delta_$(DMR_DELTA)_p_$(DMR_PVALUE).$(DMR_MINLEN)_$(DMR_MINCG)_$(DMR_DIS_MERGE)_$(DMR_PCT_SIG).ccREs_map.gz, $(CONTRAST_NAME)))

all_screen: $(ALL_SCREEN)

ifeq ($(strip $(GENOME)),)
else
ENCODE_SCREEN?=/sto1/ref/bioinfotree/task/encode-screen/dataset/v13/$(GENOME)-ccREs.bed
endif

%.ccREs.bed: %.bed $(ENCODE_SCREEN) 
	intersectBed -a $< -b $(word 2, $^) -wa -wb | cut -f1-4,10 \
	| bedtools sort -i stdin | bedtools merge -i stdin -c 4,5 -o distinct \
	> $@

%.ccREs_map.gz: %.ccREs.bed
	awk 'BEGIN{print "idx\tccRE"}{print $$1":"$$2"-"$$3"\t"$$5}' $< | gzip -c - >  $@

# ------------------------------------------ #
# GREAT annotation and functional enrichment #
# ------------------------------------------ #
GREAT_MAX_DIST_CUT?=default
# default = 100kb for basalPlusExt rule

ALL_GREAT=$(addprefix dssDMR_, $(addsuffix _delta_$(DMR_DELTA)_p_$(DMR_PVALUE).$(DMR_MINLEN)_$(DMR_MINCG)_$(DMR_DIS_MERGE)_$(DMR_PCT_SIG)_rGREAT/RegionGeneAssociation.txt.gz, $(CONTRAST_NAME)))
ALL_GREAT+=$(addprefix dssDMR_, $(addsuffix _delta_$(DMR_DELTA)_p_$(DMR_PVALUE).$(DMR_MINLEN)_$(DMR_MINCG)_$(DMR_DIS_MERGE)_$(DMR_PCT_SIG)_rGREAT/region_gene_map_$(GREAT_MAX_DIST_CUT).gz, $(CONTRAST_NAME)))

all_great: $(ALL_GREAT)

%_rGREAT/RegionGeneAssociation.txt.gz: %.bed
	conda activate rstudio_Rv4.0.3; export LD_LIBRARY_PATH=/sto1/ref/miniconda2/envs/rstudio_Rv4.0.3/jre/lib/amd64/server;\
	bin/run-rGREAT.R $< $(GENOME) "basalPlusExt" 100

ifeq ($(GREAT_MAX_DIST_CUT),default)
%_rGREAT/region_gene_map_$(GREAT_MAX_DIST_CUT).gz: %_rGREAT/RegionGeneAssociation.txt.gz
	gzip -cd $< | awk 'BEGIN{OFS="\t"};(NR==1 || $$5!="NA"){print $$4,$$8,$$5}' | gzip -c - > $@
else
%_rGREAT/region_gene_map_$(GREAT_MAX_DIST_CUT).gz: %_rGREAT/RegionGeneAssociation.txt.gz
	gzip -cd $< | awk 'BEGIN{OFS="\t"};(NR==1 || ($$5!="NA" && $$5<=$(GREAT_MAX_DIST_CUT) && $$5>=(-$(GREAT_MAX_DIST_CUT)))){print $$4,$$8,$$5}' | gzip -c - > $@
endif

# ----- #
# Utils #
# ----- #
dss%txt.gz: dss%txt
	gzip $<

dssDMR_%.dmr_info.txt.gz: dssDMR_%.bed dssDMR_%.ccREs_map.gz dssDMR_%_rGREAT/region_gene_map_default.gz dssDMR_%.hyper.names dssDMR_%.hypo.names
	bin/collect-dmr-info.r dssDMR_$* 

# ---- #
DMR/dssDMR_WT-time-course_3BKOvsWT.mregion.ratio.heatmap.col_idx.pdf: src/get-dmr-heatmap_WT-time-course_3BKO-vs-WT.r
	conda activate rstudio_Rv4.0.3; export LD_LIBRARY_PATH=/sto1/ref/miniconda2/envs/rstudio_Rv4.0.3/jre/lib/amd64/server;\
	Rscript --vanilla $<

