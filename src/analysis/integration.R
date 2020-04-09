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
library(rjson)

# Command line parameters
Args = commandArgs(trailingOnly=TRUE)

DATA_F = Args[1]
METADATA_F = Args[2]
OUT_DIR = Args[3]
#if(length(Args)==3) Nfeature=Args[3] else Nfeature=2000
Nfeature=2000

# Load metadata
meta_df <- fromJSON(file = METADATA_F)

# Load raw data

fname <- h5ls(DATA_F)
Raw_Dat <- list()
for(i in fname$name){
  Raw_Dat[[i]] <- h5read(DATA_F,i)
}
h5closeAll()

Raw_Dat$gene_name <- make.unique(Raw_Dat$gene_name)
Raw_Dat$gene_name <- gsub("_","-",Raw_Dat$gene_name)
Raw_Dat$count <- as(Raw_Dat$count,"dgCMatrix")
rownames(Raw_Dat$count) <- Raw_Dat$gene_name
colnames(Raw_Dat$count) <- Raw_Dat$cell


##Run Seurat integration pipeline

datasets <- list()
for(ds in unique(Raw_Dat$dataset)){
  datasets[[ds]] <- Raw_Dat$count[,grep(ds, Raw_Dat$dataset)]
}

for (i in seq_along(names(datasets))) {
    datasets[[i]] <- CreateSeuratObject(counts = datasets[[i]], project = "NYGC", min.cells = 1)
    #datasets[[i]] <- NormalizeData(datasets[[i]], scale.factor = 10e6)
    datasets[[i]] <- FindVariableFeatures(datasets[[i]],selection.method = "vst", nfeatures = Nfeature)
}

# Integrate together
for (integ_set_name in names(meta_df$integration_sets)) {
    # Get relavent datasets
    integ_set <- meta_df$integration_sets[[integ_set_name]]
    print(paste('Integrating datasets: ', paste(unlist(integ_set), collapse=' ')))
    integ_datasets <- datasets[integ_set]

    datasets_l <- c()
    for (ds in integ_set) {
        datasets_l <- c(datasets_l, rep(ds, length(grep(ds, Raw_Dat$cell))))
    }

    # Set integration parameters
    min_len <- min(unlist(lapply(integ_datasets, function(x) length(colnames(x)))))
    filt <- as.integer(min_len * 0.1)
    weight <- as.integer(filt / 2)
    print(paste('The filter is ', filt))
  
    tumor.anchors <- FindIntegrationAnchors(
        object.list = integ_datasets, 
        k.filter = filt, 
        dims = 1:50
    )
    tumors.integrated <- IntegrateData(
        anchorset = tumor.anchors, 
        k.weight = weight, 
        dims = 1:50
    )
 
    integ_data <- list()
    integ_data$count <- as.matrix(tumors.integrated@assays$integrated@data)
    curr_cells <- c(unlist(lapply(integ_datasets, function(x) unlist(colnames(x)))))
    integ_data$cell <- curr_cells
    integ_data$dataset <- datasets_l
  
    out_f <- paste0(OUT_DIR, '/', integ_set_name, '_aligned.h5')
    h5createFile(out_f)
    for(i in seq_along(integ_data)){
        h5write(integ_data[[i]], file = out_f, name = names(integ_data)[i])
    }
}
h5closeAll()
