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
cells_df <- list()
for(tumor in unique(Raw_Dat$tumor)){
    cells <- Raw_Dat$cell
    datasets[[tumor]] <- Raw_Dat$count[,grep(tumor,cells)]
    cells_df[[tumor]] <- cells[grep(tumor,cells)]
    
    datasets[[tumor]] <- CreateSeuratObject(counts = datasets[[tumor]], project = "NYGC", min.cells = 1)
    datasets[[tumor]] <- NormalizeData(datasets[[tumor]], scale.factor = 10e6)
    datasets[[tumor]] <- FindVariableFeatures(datasets[[tumor]],selection.method = "vst", nfeatures = Nfeature)
}

tumor_ids <- sort(names(datasets))
print(tumor_ids)
for (t_i in tumor_ids) {
    for (t_j in tumor_ids) {
        print(paste("Integrating tumor ", t_i, "and", t_j))
        curr_datasets = list(datasets[[t_i]], datasets[[t_j]])
        
        # Integrate together
        tumor.anchors <- FindIntegrationAnchors(object.list = curr_datasets, dims = 1:50)
        tumors.integrated <- IntegrateData(anchorset = tumor.anchors, dims = 1:50)
     
        str(tumors.integrated@assays$integrated@data)
        pair_dat <- list()
        # Write output H5 file
        pair_dat$count <- as.matrix(tumors.integrated@assays$integrated@data)
        gene_all <- Raw_Dat$gene_name
        id_all <- Raw_Dat$gene_id
        names(gene_all) <- gene_all
        names(id_all) <- gene_all

        pair_dat$gene_name <- unname(gene_all[rownames(Raw_Dat$count)])
        pair_dat$gene_id <- unname(id_all[rownames(Raw_Dat$count)])
        pair_dat$cell <- c(as.character(cells_df[[t_i]]), as.character(cells_df[[t_j]]))
        print(length(pair_dat$cell))
        tumor_i_ids <- as.character(rep(t_i, length(cells_df[[t_i]])))
        tumor_j_ids <- as.character(rep(t_j, length(cells_df[[t_j]])))
        tumor_all_ids <- c(tumor_i_ids, tumor_j_ids)
        print(unique(tumor_all_ids))

        h5createFile(paste0(Args[2], "/", t_i, "_", t_j, ".h5"))
        rownames(pair_dat$cell) <- NULL
        for(i in seq_along(pair_dat)){
            h5write(as.character(pair_dat$cell), file = paste0(Args[2], "/", t_i, "_", t_j, ".h5"), name = "cell")
            h5write(as.character(tumor_all_ids), file = paste0(Args[2], "/", t_i, "_", t_j, ".h5"), name = "tumor")
            h5write(unname(pair_dat$count), file = paste0(Args[2], "/", t_i, "_", t_j, ".h5"), name = "count")
        }
    }
}
h5closeAll()
