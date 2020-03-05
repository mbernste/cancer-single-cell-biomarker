# CACTUS: CharACTerizing TUmor Subpopulations

CACTUS is a R/Shiny application with associated data processing pipeline for exploring tumor subpopulations in single-cell RNA-seq data.

Currently, CACTUS supports exploring high-grade glioma scRNA-seq data from [Yuan et al.](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-018-0567-9). Future additions to CACTUS will enable exploration of other datasets and cancer types.

CACTUS implements the following features: 
* Visualize tumor scRNA-seq datasets with [PHATE](https://github.com/KrishnaswamyLab/PHATE) 
* Cluster each tumor 
* Perform differential expression (DE) analysis on each cluster for each tumor 
* Perform [gene set enrichment analysis](https://www.pnas.org/content/102/43/15545) on each cluster's DE genes using [Enrichr](https://amp.pharm.mssm.edu/Enrichr/) 
* Align tumor datasets with [Seurat](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867419305598%3Fshowall%3Dtrue) 
* Classify cell types using [CellAssign](https://shahlab.ca/projects/cellassign/) 

These analyses are all automated using a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline. 

A screenshot of the CACTUS homepage: 

![screenshot](https://github.com/mbernste/cancer-single-cell-biomarker/blob/master/screenshot.png) 

## Dependencies

CACTUS has a number of R and Python dependencies. 

Required Python packages include:
* [NumPy](https://numpy.org) 
* [Pandas](https://pandas.pydata.org) 
* [h5py](https://pypi.org/project/h5py/) 
* [Scanpy](https://icb-scanpy.readthedocs-hosted.com/en/stable/) 
* [Louvain](https://louvain-igraph.readthedocs.io/en/latest/) 
* [GSEApy](https://gseapy.readthedocs.io/en/latest/)
* [PHATE](https://github.com/KrishnaswamyLab/PHATE) 
* [Snakemake](https://snakemake.readthedocs.io/en/stable/)

Required R packages include:
* [Shiny](https://www.google.com/search?client=safari&rls=en&q=R+Shiny&ie=UTF-8&oe=UTF-8)
* [heatmaply](https://cran.r-project.org/web/packages/heatmaply/index.html) 
* [viridis](https://cran.r-project.org/web/packages/viridis/index.html) 
* [CellAssign](https://shahlab.ca/projects/cellassign/) 
* [dplyr](https://www.rdocumentation.org/packages/dplyr/versions/0.7.8) 
* [Seurat](https://satijalab.org/seurat/) 
* [ggplot2](https://ggplot2.tidyverse.org) 
* [rhdf5](http://bioconductor.org/packages/release/bioc/html/rhdf5.html) 
* [shinycssloaders](https://cran.r-project.org/web/packages/shinycssloaders/index.html) 
* [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html) 
* [rjson](https://cran.r-project.org/web/packages/rjson/index.html) 
* [DT](https://cran.r-project.org/web/packages/DT/index.html) 

## Run the Shiny application with pre-computed data

The R/Shiny application only requires the R dependencies outlined above. Once these dependencies are installed, the Shiny application can be run using pre-computed data from [Yuan et al.](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-018-0567-9) by following these steps:
1. Download the [zipped-data](https://uwmadison.box.com/s/8gxeyb7ropvi0up1ydaiy2jud0amlo0f)
2. Unzip the folder:
``unzip cactus_data.zip``
3. Run Shiny application with the command:
``Rscript run_app.R <absolute path to cactus_data>``

## Running CACTUS

There exist two Snakemake workflows for processing the RNA-seq data in [Yuan et al.](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-018-0567-9). The first downloads the raw data from [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103224) (GEO). The second computes all of the analyses that are viewable in the R/Shiny application. 

### Downloading the data

The data download and preparation procedure is implemented by the following [Snakemake workflow](https://github.com/mbernste/cancer-single-cell-biomarker/blob/master/src/preprocess/Snakefile):

![DAG](https://github.com/mbernste/cancer-single-cell-biomarker/blob/master/prep_data_dag.png)

To run this Snakemake workflow, perform the following:
1. Modify the [config.json](https://github.com/mbernste/cancer-single-cell-biomarker/blob/master/src/config.json) file to specify the paths for which CACTUS should store the raw data from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103224)
2. Change the working directory to the [preprocessing directory](https://github.com/mbernste/cancer-single-cell-biomarker/tree/master/src/preprocess): ``cd src/preprocess``
3. Run Snakemake: ``snakemake``

### Running the data processing pipeline

The [second Snakemake pipeline](https://github.com/mbernste/cancer-single-cell-biomarker/blob/master/src/analysis/Snakefile) implements the follwoing:

![DAG](https://github.com/mbernste/cancer-single-cell-biomarker/blob/master/dag.png)
