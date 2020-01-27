library(jsonlite)
library(sets)
library("rhdf5")
library("SingleCellExperiment")
library("scran")
library("cellassign")
library("readr")

args = commandArgs(trailingOnly=TRUE)

DATA_DIR = args[1]
OUT_DIR = args[2]
JSON_F = args[3]

TUMORS = c('PJ016', 'PJ018', 'PJ025', 'PJ048', 'PJ030', 'PJ035', 'PJ017', 'PJ032')

json_txt <- read_file(JSON_F)
j <- fromJSON(json_txt)

all_genes <- c()
all_groups <- c()

for (x in names(j)) {
  all_genes <- c(all_genes, j[[x]])
  all_groups <- c(all_groups, x)
  #all_genes <- set(union(set(all_genes), set(j[x])))
}
all_genes <- unique(all_genes)
all_groups <- unique(all_groups)

mat <- c()
for (gene in all_genes) {
  row <- c()
  for (group in all_groups) {
    if (gene %in% j[[group]]) {
      row <- append(row, 1)
    }
    else {
      row <- append(row, 0)
    }
  }
  mat <- append(mat, row)
  print(mat)
}
marker_mat <- matrix(mat, nrow = length(all_genes), ncol = length(all_groups))
rownames(marker_mat) <- all_genes
colnames(marker_mat) <- all_groups

print('Loading counts data...')
counts = h5read(paste0(DATA_DIR, "/GSE103224.h5"),"count")
cells = h5read(paste0(DATA_DIR, "/GSE103224.h5"), "cell")
tumor_ids = h5read(paste0(DATA_DIR, "/GSE103224.h5"), "tumor")
genes = h5read(paste0(DATA_DIR, "/GSE103224.h5"), "gene_name")
print('done.')
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
  
  sce <- SingleCellExperiment(assays = list(counts=tumor_counts))
  rownames(sce) <- genes
  colnames(sce) <- tumor_cells                            
  
  print('Computing sum factors...')
  sce <- computeSumFactors(sce)
  s <- sizeFactors(sce)
  print('done.')
  
  sce_small <- sce[rownames(marker_mat),]
  print('Running CellAssign...')
  fit <- cellassign(
    exprs_obj = sce_small, 
    marker_gene_info = marker_mat, 
    s = s, 
    learning_rate = 1e-2, 
    shrinkage = TRUE,
    verbose = FALSE
  )
  print('done.')
  celltypes_mat <- cbind(celltypes(fit))
  rownames(celltypes_mat) <- tumor_cells
  colnames(celltypes_mat) <- c('predicted_cell_type')
  
  cellprobs_mat <- cbind(cellprobs(fit))
  rownames(cellprobs_mat) <- tumor_cells
  colnames(cellprobs_mat) <- colnames(marker_mat)

  write.table(celltypes_mat, file = paste0(OUT_DIR, '/', tumor, '_predicted_cell_types.tsv'), sep = '\t', quote = FALSE)
  write.table(cellprobs_mat, file = paste0(OUT_DIR, '/', tumor, '_predicted_cell_type_probs.tsv'), sep = '\t', quote = FALSE)
}
  
