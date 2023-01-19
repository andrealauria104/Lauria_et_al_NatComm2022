# ----- #
# Setup #
# ----- #
# Conda ----
CONDA_ROOT=/opt/conda
CONDA_VERSION=miniconda3
CONDA_ACTIVATE=source $(CONDA_ROOT)/$(CONDA_VERSION)/etc/profile.d/conda.sh; conda activate
ENV_PREFIX?=$(CONDA_ROOT)/$(CONDA_VERSION)

# ------------ #
# DMR Gene CpG #
# ------------ #
dmr.gene_Sox2.boxplot.row.pdf: src/plot-dmr-cpg_3bectopic.r
	$(CONDA_ACTIVATE) rstudio_Rv4.0.3; export LD_LIBRARY_PATH=$(ENV_PREFIX)/envs/rstudio_Rv4.0.3/jre/lib/amd64/server;\
	Rscript --vanilla $<

DMR_hyper_CRISPRoff_none.boxplot.row_col.pdf: src/plot-dmr-cpg_crisproff.r
	$(CONDA_ACTIVATE) rstudio_Rv4.0.3; export LD_LIBRARY_PATH=$(ENV_PREFIX)/envs/rstudio_Rv4.0.3/jre/lib/amd64/server;\
	Rscript --vanilla $<
