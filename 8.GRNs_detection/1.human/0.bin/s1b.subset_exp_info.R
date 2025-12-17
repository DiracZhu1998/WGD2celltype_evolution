library(Seurat)
library(dplyr)
library(Matrix)

args <- commandArgs(T)
input <- args[1]
output <- args[2]

set.seed(42)
obj <- readRDS(input)
exprMat <- GetAssayData(obj, layer = 'counts')
#LogexprMat <- readRDS("LogexprMat.rds")
cellInfo <- obj@meta.data
cellInfo$CellID <- rownames(cellInfo)
# Sample 500 cells for each unique integrated clusters
chosen_cells <- cellInfo %>% filter(`Refined family` != 'Mixture') %>% group_by(`Refined family`) %>% slice_sample(n = 2000, replace = FALSE) %>% ungroup()
chosen_cells <- chosen_cells$CellID

cellInfo <- cellInfo[match(chosen_cells, rownames(cellInfo)),]
exprMat <- as.matrix(exprMat[, chosen_cells])
#LogexprMat <- as.matrix(LogexprMat[, chosen_cells])

write.table(exprMat, file = paste0(output, "/0.exprMat.subset.csv"), sep = ",")
#write.table(LogexprMat, file = paste0(output, "/0.LogexprMat.subset.csv"), sep = ",")
write.table(cellInfo, file = paste0(output, "/0.cellInfo.subset.csv"), sep = ",")
