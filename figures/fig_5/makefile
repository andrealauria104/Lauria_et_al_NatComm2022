# ----- #
# Setup #
# ----- #
# Conda ----
CONDA_ROOT=/opt/conda
CONDA_VERSION=miniconda3
CONDA_ACTIVATE=source $(CONDA_ROOT)/$(CONDA_VERSION)/etc/profile.d/conda.sh; conda activate
ENV_PREFIX?=$(CONDA_ROOT)/$(CONDA_VERSION)

GENOME?=mm10

# -------- #
# Analysis #
# -------- #
results/r.mfa.rds: src/mfa-methyl-expr.r
	$(CONDA_ACTIVATE) rstudio_Rv4.0.3; export LD_LIBRARY_PATH=$(ENV_PREFIX)/envs/rstudio_Rv4.0.3/jre/lib/amd64/server;\
	Rscript --vanilla $<
results/mfa.pdf: results/r.mfa.rds
	@echo done

results/methyl_expr_dynamics.heatmap.pdf: src/plot-methyl-expr-heatmap.r
	$(CONDA_ACTIVATE) rstudio_Rv4.0.3; export LD_LIBRARY_PATH=$(ENV_PREFIX)/envs/rstudio_Rv4.0.3/jre/lib/amd64/server;\
	Rscript --vanilla $<
results/methyl_expr_dynamics.lines.trend.pdf: src/plot-methyl-expr-dynamics.r
	$(CONDA_ACTIVATE) rstudio_Rv4.0.3; export LD_LIBRARY_PATH=$(ENV_PREFIX)/envs/rstudio_Rv4.0.3/jre/lib/amd64/server;\
	Rscript --vanilla $<

