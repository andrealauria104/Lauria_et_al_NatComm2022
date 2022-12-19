#################
# Configuration # ----
#################
CONDA_ROOT?=/opt/conda
CONDA_VERSION?=miniconda3
CONDA_ACTIVATE=source $(CONDA_ROOT)/$(CONDA_VERSION)/etc/profile.d/conda.sh; conda activate
ENV_PREFIX?=$(CONDA_ROOT)/$(CONDA_VERSION)
GENCODE_BIOTYPES?=$(GENCODE_DIR)/$(GENCODE_ANNOTATION).ensg2gene_symbol2biotype.map.header_added

#################
# Preprocessing # ----
#################
# create single-cell object ---
data/seurat_obj.rds: data/GEP.count.gz \
		     data/metadata.txt \
		     $(GENCODE_BIOTYPES)
	$(CONDA_ACTIVATE) rstudio_Rv4.0.3;export LD_LIBRARY_PATH=$(ENV_PREFIX)/envs/rstudio_Rv4.0.3/jre/lib/amd64/server;\
	bin/scrnaseq-get_scrna_obj.r -c $< \
	-m $^2 -g $^3 -w seurat -a $(ANALYSIS)

#####################
# Analysis - Seurat # ----
#####################
INPUT_SEURAT_STANDARD?=data/seurat_obj.rds
INPUT_SEURAT_CLUSTER?=seurat/results_seurat-standard/seurat_obj.rds
INPUT_SEURAT_PSEUDOTIME_MONOCLE3?=seurat/results_seurat-cluster/seurat_obj.rds

seurat/results_seurat-standard/seurat_obj.rds: src/scrnaseq-seurat-standard.r $(INPUT_SEURAT_STANDARD)
	$(CONDA_ACTIVATE) rstudio_Rv4.0.3;export LD_LIBRARY_PATH=$(ENV_PREFIX)/envs/rstudio_Rv4.0.3/jre/lib/amd64/server;\
	OMP_NUM_THREADS=$(CORES) Rscript --vanilla $<
seurat/results_seurat-cluster/seurat_obj.rds: src/scrnaseq-seurat-cluster.r $(INPUT_SEURAT_CLUSTER)
	$(CONDA_ACTIVATE) rstudio_Rv4.0.3;export LD_LIBRARY_PATH=$(ENV_PREFIX)/envs/rstudio_Rv4.0.3/jre/lib/amd64/server;\
	OMP_NUM_THREADS=$(CORES) Rscript --vanilla $<
seurat/results_pseudotime-monocle3/cds.rds: src/scrnaseq-pseudotime-monocle3.r $(INPUT_SEURAT_PSEUDOTIME_MONOCLE3)
	$(CONDA_ACTIVATE) rstudio_Rv4.0.3;export LD_LIBRARY_PATH=$(ENV_PREFIX)/envs/rstudio_Rv4.0.3/jre/lib/amd64/server;\
	OMP_NUM_THREADS=$(CORES) Rscript --vanilla $<

seurat/results_visualize-data/dimplot.umap.pdf: src/scrnaseq-visualize-data.r src/results_seurat-cluster/seurat_obj.rds src/results_pseudotime-monocle3/cds.rds
	$(CONDA_ACTIVATE) rstudio_Rv4.0.3;export LD_LIBRARY_PATH=$(ENV_PREFIX)/envs/rstudio_Rv4.0.3/jre/lib/amd64/server;\
	OMP_NUM_THREADS=$(CORES) Rscript --vanilla $<	

seurat/results_seurat-transfer/done: src/scrnaseq-seurat-transfer.r seurat/results_seurat-cluster/seurat_obj.rds
	$(CONDA_ACTIVATE) rstudio_Rv4.0.3;export LD_LIBRARY_PATH=$(ENV_PREFIX)/envs/rstudio_Rv4.0.3/jre/lib/amd64/server;\
	OMP_NUM_THREADS=$(CORES) Rscript --vanilla $<
	touch $@

seurat/results_seurat-cluster/da-fisher/comparisons.rds: src/scrnaseq-seurat-da-fisher.r seurat/results_seurat-cluster/seurat_obj.rds
	$(CONDA_ACTIVATE) rstudio_Rv4.0.3;export LD_LIBRARY_PATH=$(ENV_PREFIX)/envs/rstudio_Rv4.0.3/jre/lib/amd64/server;\
	OMP_NUM_THREADS=$(CORES) Rscript --vanilla $<
	touch $@
seurat/results_seurat-cluster/da-fisher/p_perc_table_lineage.pdf: src/scrnaseq-visualize_da-fisher.r seurat/results_seurat-cluster/da-fisher/comparisons.rds
	$(CONDA_ACTIVATE) rstudio_Rv4.0.3;export LD_LIBRARY_PATH=$(ENV_PREFIX)/envs/rstudio_Rv4.0.3/jre/lib/amd64/server;\
	OMP_NUM_THREADS=$(CORES) Rscript --vanilla $<

PVAL_TH?=0.01
seurat/results_seurat-cluster/seurat_clusters_markers.pval_$(PVAL_TH).txt.gz: seurat/results_seurat-cluster/seurat_clusters_markers.txt.gz seurat/results_seurat-cluster/seurat_obj.rds
	zcat $< | awk 'NR==1 || $$1<$(PVAL_TH)' | gzip -c - > $@

