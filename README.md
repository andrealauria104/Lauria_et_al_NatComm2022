# DNMT3B supports meso-endoderm differentiation from mouse embryonic stem cells

Code repository for the paper:

<cite>Andrea Lauria*, Guohua Meng*, Valentina Proserpio*, Stefania Rapelli*, Mara Maldotti,
Isabelle Laurence Polignano, Francesca Anselmi, Chiara Levra Levron, Giacomo Donati,
Danny Incarnato, Anna Krepelova, Francesco Neri, Ivan Molineris, Salvatore Oliviero, "DNMT3B supports meso-endoderm differentiation from mouse embryonic stem cells." Nat Commun 14, 367 (2023). https://doi.org/10.1038/s41467-023-35938-x </cite>

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
    ├── scrnaseq-seurat-cluster.r
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

Raw and processed sequencing data are available from the Gene Expression Omnibus (GEO) database under accession code [GSE168415](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168415).

## Used software

- STAR 2.7.1a: https://github.com/alexdobin/STAR
- HISAT2 v2.2.1: http://daehwankimlab.github.io/hisat2/ 
- featureCounts v1.6.4: https://subread.sourceforge.net/ 
- Bismark v0.22.3: https://www.bioinformatics.babraham.ac.uk/projects/bismark/
- Bowtie2 v2.3.5: https://bowtie-bio.sourceforge.net/bowtie2/index.shtml
- bedtools v2.29.2: https://bedtools.readthedocs.io/en/latest/
- deepTools v3.3.1: https://deeptools.readthedocs.io/en/develop/
- MACS2 v2.2.5: https://pypi.org/project/MACS2/ 
- R v3.5.3, v3.6.1, v4.0.3: https://www.r-project.org/
- Seurat v4.0.0: https://satijalab.org/seurat/index.html
- Monocle3 v1.0.0: https://cole-trapnell-lab.github.io/monocle3/
- edgeR v3.32.1: https://bioconductor.org/packages/release/bioc/html/edgeR.html
- DSS v2.34.0: https://bioconductor.org/packages/release/bioc/html/DSS.html
- MethylKit v1.16.1: https://www.bioconductor.org/packages/release/bioc/html/methylKit.html

## Additional sofware required for the analysis

- `ngsRtools`: https://github.com/andrealauria104/ngsRtools
- `epigeNGS-pipe`: https://github.com/andrealauria104/epigeNGS-pipe

