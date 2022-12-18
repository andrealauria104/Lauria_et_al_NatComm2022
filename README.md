# DNMT3B supports meso-endoderm differentiation from mouse embryonic stem cells

Code repository for the paper

-- <cite>Andrea Lauria*, Guohua Meng*, Valentina Proserpio*, Stefania Rapelli*, Mara Maldotti,
Isabelle Laurence Polignano, Francesca Anselmi, Chiara Levra Levron, Giacomo Donati,
Danny Incarnato, Anna Krepelova, Francesco Neri, Ivan Molineris, Salvatore Oliviero, Lack
of DNMT3B in mouse epiblasts impairs meso-endoderm differentiation.</cite>

## Abstract

## Repository structure

The folder structure of the repository is organized in the following way:

```bash
.
├── figures
│   ├── fig_1
│   ├── fig_2
│   ├── fig_3
│   ├── fig_4
│   ├── fig_5
│   ├── fig_6
│   └── fig_7
├── README.md
└── share
    ├── bin
    ├── env
    └── src
```

where:

- `figures`: code to reproduce the results presented in the figures.
- `share`: code and conda environment files shared among the analyses.

Each figure folder is organized as follows, e.g. for `fig_1`:

```bash

figures/fig_1
├── bin
│   └── scrnaseq-get_scrna_obj.r -> ../src/scrnaseq-get_scrna_obj.r
├── data
├── env
│   └── rstudio_Rv4.0.3.yml -> ../../../share/env/rstudio_Rv4.0.3.yml
├── makefile
└── src
    ├── scrnaseq-get_scrna_obj.r
    ├── scrnaseq-pseudotime-monocle3.r
    ├── scrnaseq-seurat-da-fisher.r
    ├── scrnaseq-seurat-standard.r
    ├── scrnaseq-test_branch_expression.r
    ├── scrnaseq-visualize_da-fisher.r
    └── scrnaseq-visualize-data.r
```

where:

- `bin`: executable scripts/tools.
- `data`: processed data files.
- `env`: conda environment files.
- `makefile`: makefile with rules to run the analysis.
- `src`: analysis scripts.

## Data availability

Raw and processed sequencing data are available from the Gene Expression Omnibus (GEO) database under accession code GSEXXXX

## Additional sofware required for the analysis:

- `ngsRtools`: https://github.com/andrealauria104/ngsRtools
- `epigeNGS-pipe`: https://github.com/andrealauria104/epigeNGS-pipe

