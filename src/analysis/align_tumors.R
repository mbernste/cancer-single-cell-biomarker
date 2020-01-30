library(jsonlite)
library(sets)
library("rhdf5")
library("SingleCellExperiment")
library("scran")
library("cellassign")
library("readr")
library("Seurat")

args = commandArgs(trailingOnly=TRUE)

DATA_DIR = args[1]

TUMORS = c('PJ016', 'PJ018', 'PJ025', 'PJ048', 'PJ030', 'PJ035', 'PJ017', 'PJ032')
#TUMORS = c('PJ016', 'PJ018')

print('Loading counts data...')
counts = h5read(paste0(DATA_DIR, "/GSE103224.h5"),"count")
cells = h5read(paste0(DATA_DIR, "/GSE103224.h5"), "cell")
tumor_ids = h5read(paste0(DATA_DIR, "/GSE103224.h5"), "tumor")
genes = h5read(paste0(DATA_DIR, "/GSE103224.h5"), "gene_name")
print('done.')

sos <- list()
for (tumor in TUMORS) {
    print(paste('Processing tumor', tumor))
    tumor_indices = c()
    for (i in 1:length(cells)) {
        if (tumor_ids[[i]] == tumor) {
            tumor_indices <- append(tumor_indices, i)
        }
    }
    tumor_cells <- cells[tumor_indices]
    tumor_counts <- counts[,tumor_indices]
    rownames(tumor_counts) <- genes
    colnames(tumor_counts) <- tumor_cells                            

    so <- CreateSeuratObject(
        tumor_counts,
        project = tumor,
        assay = "RNA",
        min.cells = 1,
        min.features = 0,
        names.field = NULL,
        names.delim = NULL,
        meta.data = NULL
    )
    so <- NormalizeData(so)
    sos[[tumor]] <- so
}

anchors <- FindIntegrationAnchors(object.list = sos, dims = 1:50)
integ <- IntegrateData(anchorset = anchors, dims = 1:50)

DefaultAssay(integ) <- "integrated"


#print(integ)
#print(integ@assays)
#print(integ@assays$integrated)

#integrated_df <- as.matrix(integ@assays$data)
write.csv(integ@assays$integrated@data, "integrated_scaled.tsv", row.names = TRUE, sep = "\t")
  
