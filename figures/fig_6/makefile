# ----- #
# Setup #
# ----- #
# Conda ----
CONDA_ROOT=/opt/conda
CONDA_VERSION=miniconda3
CONDA_ACTIVATE=source $(CONDA_ROOT)/$(CONDA_VERSION)/etc/profile.d/conda.sh; conda activate
ENV_PREFIX?=$(CONDA_ROOT)/$(CONDA_VERSION)

# -------- #
# Analysis #
# -------- #
results/tf_target.gz: src/build-Dnmt3b-network-2.r
	$(CONDA_ACTIVATE) rstudio_Rv4.0.3; export LD_LIBRARY_PATH=$(ENV_PREFIX)/envs/rstudio_Rv4.0.3/jre/lib/amd64/server;\
	Rscript --vanilla $<

results/Dnmt3b_tf_network_2.full.pdf: results/tf_target.gz
	@echo done
results/Dnmt3b_tf_network_degree_2.pdf: results/tf_target.gz
	@echo done
results/reg_genes_with_rescue.gz: results/tf_target.gz
	@echo done

results/p_tf_comparisons_table_top.pdf: src/analyze-tf-network-local.r results/tf_target.gz results/reg_genes_with_rescue.gz
	$(CONDA_ACTIVATE) rstudio_Rv4.0.3; export LD_LIBRARY_PATH=$(ENV_PREFIX)/envs/rstudio_Rv4.0.3/jre/lib/amd64/server;\
	Rscript --vanilla $<
