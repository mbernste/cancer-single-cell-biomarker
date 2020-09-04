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
tum_1 = Args[2]
tum_2 = Args[3]
#DATA_F = "/Users/matthewbernstein/Development/single-cell-hackathon/charts.h5"
OUT_F = Args[4]
#if(length(Args)==3) Nfeature=Args[3] else Nfeature=2000
Nfeature=2000

# Load raw data
fname <- h5ls(DATA_F)
Raw_Dat <- list()
for(i in fname$name){
  Raw_Dat[[i]] <- h5read(DATA_F,i)
}
h5closeAll()


gene_names_tum_1 = paste0(tum_1, "_gene_name")
cells_tum_1 = paste0(tum_1, "_cell")
expr_tum_1 = paste0(tum_1, "_log1_tpm")
gene_names_tum_2 = paste0(tum_2, "_gene_name")
cells_tum_2 = paste0(tum_2, "_cell")
expr_tum_2 = paste0(tum_2, "_log1_tpm")


Raw_Dat[[gene_names_tum_1]] <- make.unique(Raw_Dat[[gene_names_tum_1]])
Raw_Dat[[gene_names_tum_1]] <- gsub("_","-",Raw_Dat[[gene_names_tum_1]])
Raw_Dat[[expr_tum_1]] <- as(Raw_Dat[[expr_tum_1]],"dgCMatrix")


Raw_Dat[[gene_names_tum_2]] <- make.unique(as(Raw_Dat[[gene_names_tum_2]], "vector"))
Raw_Dat[[gene_names_tum_2]] <- gsub("_","-",Raw_Dat[[gene_names_tum_2]])
Raw_Dat[[expr_tum_2]] <- as(Raw_Dat[[expr_tum_2]],"dgCMatrix")

##Run Seurat integration pipeline

datasets <- list()

ds_1 <- Raw_Dat[[expr_tum_1]]
ds_2 <- Raw_Dat[[expr_tum_2]]

rownames(ds_1) <-  Raw_Dat[[gene_names_tum_1]]
colnames(ds_1) <- Raw_Dat[[cells_tum_1]]

rownames(ds_2) <-  Raw_Dat[[gene_names_tum_2]]
colnames(ds_2) <- Raw_Dat[[cells_tum_2]]

print("HERE")
print(Raw_Dat[[cells_tum_1]])
print(colnames(ds_1))

ds_1 <- CreateSeuratObject(counts = ds_1, project = "CHARTS", min.cells = 1)
ds_1 <- FindVariableFeatures(ds_1, selection.method = "vst", nfeatures = Nfeature)
ds_2 <- CreateSeuratObject(counts = ds_2, project = "CHARTS", min.cells = 1)
ds_2 <- FindVariableFeatures(ds_2, selection.method = "vst", nfeatures = Nfeature)

integ_datasets <- list(ds_1, ds_2)

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

combined_expr <- paste0(tum_1, '.', tum_2, '_log1_tpm')
combined_cells <- paste0(tum_1, '.', tum_2, '_cell')

integ_data <- list()
integ_data[[combined_expr]] <- as.matrix(tumors.integrated@assays$integrated@data)
curr_cells <- c(unlist(lapply(integ_datasets, function(x) unlist(colnames(x)))))
integ_data[[combined_cells]] <- curr_cells

out_f <- paste0('./', OUT_F)
h5createFile(out_f)
for(i in seq_along(integ_data)){
  h5write(integ_data[[i]], file = out_f, name = names(integ_data)[i])
}
h5closeAll()

