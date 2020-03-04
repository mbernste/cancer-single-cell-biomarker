#############
# Modified by Zijian Ni
#############

# Input: Raw data from multiple tumors.
# Output: Normalized and integrated data.
# Usage: 
# Rscript integration.R <Input h5 file> <Output h5 prefix> <# of high-variant features. Default=2000>
# E.g. 
# Rscript integration.R In.h5 Out


if (!requireNamespace("rhdf5", quietly = TRUE)) install.packages("rhdf5")
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")

library(Seurat)
library(Matrix)
library(rhdf5)

# Command line parameters
Args = commandArgs(trailingOnly=TRUE)

DATA_DIR = Args[1]
OUT_DIR = Args[2]
if(length(Args)==3) Nfeature=Args[3] else Nfeature=2000

# Load raw data

fname <- h5ls(DATA_DIR)
Raw_Dat <- list()
for(i in fname$name){
    Raw_Dat[[i]] <- h5read(DATA_DIR,i)
}
h5closeAll()

Raw_Dat$gene_name <- make.unique(Raw_Dat$gene_name)
Raw_Dat$gene_name <- gsub("_","-",Raw_Dat$gene_name)
Raw_Dat$count <- as(Raw_Dat$count,"dgCMatrix")
rownames(Raw_Dat$count) <- Raw_Dat$gene_name
colnames(Raw_Dat$count) <- Raw_Dat$cell


##Run Seurat integration pipeline

datasets <- list()
for(tumor in unique(Raw_Dat$tumor)){
#for(tumor in c('PJ025', 'PJ048', 'TCGA_GBM')){
    datasets[[tumor]] <- Raw_Dat$count[,grep(tumor,Raw_Dat$tumor)]
}
#print(dim(datasets[['TCGA_GBM']]))

for (i in seq_along(names(datasets))) {
    datasets[[i]] <- CreateSeuratObject(counts = datasets[[i]], project = "NYGC", min.cells = 1)
    datasets[[i]] <- NormalizeData(datasets[[i]], scale.factor = 10e6)
    datasets[[i]] <- FindVariableFeatures(datasets[[i]],selection.method = "vst", nfeatures = Nfeature)
}

# Integrate together
tumor.anchors <- FindIntegrationAnchors(object.list = datasets, dims = 1:50)
tumors.integrated <- IntegrateData(anchorset = tumor.anchors, dims = 1:50)


str(tumors.integrated@assays$integrated@data)
# Write output H5 file
Raw_Dat$count <- as.matrix(tumors.integrated@assays$integrated@data)
gene_all <- Raw_Dat$gene_name
id_all <- Raw_Dat$gene_id
names(gene_all) <- gene_all
names(id_all) <- gene_all

Raw_Dat$gene_name <- unname(gene_all[rownames(Raw_Dat$count)])
Raw_Dat$gene_id <- unname(id_all[rownames(Raw_Dat$count)])

h5createFile(OUT_DIR)
for(i in seq_along(Raw_Dat)){
    h5write(Raw_Dat[[i]], file = OUT_DIR, name = names(Raw_Dat)[i])
}
h5closeAll()
