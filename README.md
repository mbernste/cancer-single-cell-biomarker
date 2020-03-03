# CACTUS: CharACTerizing TUmor Subpopulations

CACTUS is a R/Shiny application with associated data processing pipeline for exploring tumor subpopulations in single-cell RNA-seq data.

CACTUS implements the following features: 
* Visualize tumor scRNA-seq datasets with [PHATE](https://github.com/KrishnaswamyLab/PHATE) 
* Cluster each tumor 
* Perform differential expression (DE) analysis on each cluster for each tumor 
* Pefrom gene set enrichment analysis on each cluster's DE genes using [Enrichr](https://amp.pharm.mssm.edu/Enrichr/) 
* Align tumor datasets with [Seurat](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867419305598%3Fshowall%3Dtrue) 
* Classify cell types using [CellAssign](https://shahlab.ca/projects/cellassign/)

These analyses are all automated using a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline. 

## Dependencies

CACTUS has a number of R and Python dependencies. 

Required Python packages include:
* [NumPy](https://numpy.org) 
* [Pandas](https://pandas.pydata.org) 
* [Scanpy](https://icb-scanpy.readthedocs-hosted.com/en/stable/) 
* [Louvain](https://louvain-igraph.readthedocs.io/en/latest/) 
* [GSEApy](https://gseapy.readthedocs.io/en/latest/)
* [PHATE](https://github.com/KrishnaswamyLab/PHATE) 
* [Snakemake]((https://snakemake.readthedocs.io/en/stable/)

Required R packages include:
* [Shiny](https://www.google.com/search?client=safari&rls=en&q=R+Shiny&ie=UTF-8&oe=UTF-8)
* [heatmaply](https://cran.r-project.org/web/packages/heatmaply/index.html) 
* [viridis](https://cran.r-project.org/web/packages/viridis/index.html) 
* [CellAssign](https://shahlab.ca/projects/cellassign/) 
* [dplyr](https://www.rdocumentation.org/packages/dplyr/versions/0.7.8) 
* [Seurat](https://satijalab.org/seurat/) 

## Running CACTUS

The Snakemake pipeline:
![DAG](https://github.com/mbernste/cancer-single-cell-biomarker/blob/master/dag.png)
