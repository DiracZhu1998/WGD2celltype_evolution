suppressMessages(suppressWarnings({
    library(Seurat)
}))
source('Trinarization_score.R')


parser = argparse::ArgumentParser(description="")
parser$add_argument('-I','--input', help='input atlas')
parser$add_argument('-S','--species', help='input species name')
args = parser$parse_args()

obj <- readRDS(args$input)
label = 'Refined family'
obj@meta.data[,label] <- as.character(obj@meta.data[,label])
Idents(obj) <- label

# Extract count matrix from Seurat object (assumes RNA assay)
counts <- GetAssayData(obj, layer = "counts")
cluster_ids <- obj@meta.data[, label]
rm(obj)
counts <- as.matrix(counts)


result <- trinarize(counts, cluster_ids, 0.1)

result[result < 0.05] = 0
result[result != 0] = 1

result <- reshape2::melt(result)
result <- result[result$value != 0, ]
result$value <- NULL
colnames(result) <- c('gene', 'cluster')
result$species <- args$species
result$species_gene <- paste(args$species, result$gene, sep = '_')
write.table(result, paste0("exp_domain/", args$species, '.exp_tri_score_0.1_0.05.txt'), quote = F, sep = '\t', row.names = F, col.names = T)
