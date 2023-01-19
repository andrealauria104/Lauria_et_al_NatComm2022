# ----- #
# Setup #
# ----- #
# Conda ----
CONDA_ROOT=/opt/conda
CONDA_VERSION=miniconda3
CONDA_ACTIVATE=source $(CONDA_ROOT)/$(CONDA_VERSION)/etc/profile.d/conda.sh; conda activate
ENV_PREFIX?=$(CONDA_ROOT)/$(CONDA_VERSION)

GENOME?=mm10
PATH_ANN_INFO?=$(REFERENCE_ROOT)/indexes/hisat2/grcm38_tran/Mus_musculus.GRCm38.84.gene_info.fix_header.tsv.gz
PATH_TLEN?=$(REFERENCE_ROOT)/indexes/hisat2/grcm38_tran/Mus_musculus.GRCm38.84.gene_length.tsv

# ---------------------------------- #
# Downstream analysis - RNAseqRtools #
# ---------------------------------- #
data/dgelist_obj.rds: data/raw_counts.csv \
                      data/metadata.txt $(PATH_ANN_INFO) $(PATH_TLEN) 
	$(CONDA_ACTIVATE) rstudio_Rv4.0.3; export LD_LIBRARY_PATH=$(ENV_PREFIX)/envs/rstudio_Rv4.0.3/jre/lib/amd64/server;\
	bin/rnaseqrtools-get_dgelist_obj.r -c $< -m $(word 2,$^) -g $(word 3,$^) -t $(word 4,$^) -o `dirname $@` -u rpkm -s 2
data/norm_counts.rpkm.tmm.txt.gz: rnaseqrtools/dgelist_obj.rds
	@echo depends on $<
	@echo done

results/dge-edger/done: src/rnaseqrtools-dge-edger.r data/dgelist_obj.rds
	$(CONDA_ACTIVATE) rstudio_Rv4.0.3; export LD_LIBRARY_PATH=$(ENV_PREFIX)/envs/rstudio_Rv4.0.3/jre/lib/amd64/server;\
	Rscript --vanilla $<
	date > $@

results/expr_dynamics.lines.trend.pdf: src/analyze-dge-dynamics.r data/dgelist_obj.rds results/dge-edger/done
	$(CONDA_ACTIVATE) rstudio_Rv4.0.3; export LD_LIBRARY_PATH=$(ENV_PREFIX)/envs/rstudio_Rv4.0.3/jre/lib/amd64/server;\
	Rscript --vanilla $<

results/rescue_overlaps_full.pdf: src/calc-rescued-genes.r results/dge-edger/done
	$(CONDA_ACTIVATE) rstudio_Rv4.0.3; export LD_LIBRARY_PATH=$(ENV_PREFIX)/envs/rstudio_Rv4.0.3/jre/lib/amd64/server;\
	Rscript --vanilla $<
