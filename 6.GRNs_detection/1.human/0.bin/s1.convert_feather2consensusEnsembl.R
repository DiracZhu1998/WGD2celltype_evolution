
library(arrow)
library(Seurat)
library(dplyr)

gene_info <- read.delim("0.human_gene_info_from_Ensembl.mod", header =T, sep = "\t")
feather1 <- read_feather("cisTargets/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
feather2 <- read_feather("cisTargets/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")

gene_info[gene_info$Gene.name == NULL, "Gene.name"] <- gene_info[gene_info$Gene.name == NULL, "Gene.stable.ID"]
gene_info[gene_info$Gene.name == NULL, "Gene.Synonym"] <- gene_info[gene_info$Gene.name == NULL, "Gene.stable.ID"]

Ensembl <- list()
for (i in 1:nrow(gene_info)){
    tryCatch({Ensembl[[gene_info[i, 2]]] <- unique(c( Ensembl[[gene_info[i, 2]]], gene_info[i, 1]))},
        error = function(e){Ensembl[[gene_info[i, 2]]] <- gene_info[i, 1] })
    tryCatch({Ensembl[[gene_info[i, 3]]] <- unique(c( Ensembl[[gene_info[i, 3]]], gene_info[i, 1]))},
        error = function(e){Ensembl[[gene_info[i, 3]]] <- gene_info[i, 1] })
}

# some genes are not exist in Ensembl biomart, keep it as original name, e.g., 0610010B08Rik
# some lncRNA has relationships with two gene IDs, like 1700030C10Rik corresponding to ENSMUSG00000099759 and ENSMUSG00000091071
# I also wouldn't change these genes ID, it doesn't matter actually since for our scRNA analysis we only focus on coding genes.
colnames(feather1)[2:ncol(feather1)] <- setNames(vapply(colnames(feather1)[2:ncol(feather1)], 
                                                        FUN = function(g){ if (length(Ensembl[[g]]) == 1){
                                                            Ensembl[[g]] } else {
                                                                g }}, FUN.VALUE=character(1)), NULL)
colnames(feather2)[2:ncol(feather2)] <- setNames(vapply(colnames(feather2)[2:ncol(feather2)], 
                                                        FUN = function(g){ if (length(Ensembl[[g]]) == 1){
                                                            Ensembl[[g]] } else {
                                                                g }}, FUN.VALUE=character(1)), NULL)
# A730008H23Rik Hjurp is same gene in latest version Ensembl, but feather file contain both column, there are several other few cases like this.
# only keep first occurence, this is an essential step
duplicated1 <- table(colnames(feather1))
duplicated1 <- names(duplicated1)[duplicated1 >= 2]
for (i in duplicated1){
    tmp <- which(colnames(feather1) == i)
    feather1[, tmp[1:length(tmp)]] <- NULL # remove all occurence!
}
duplicated2 <- table(colnames(feather2))
duplicated2 <- names(duplicated2)[duplicated2 >= 2]
for (i in duplicated2){
    tmp <- which(colnames(feather2) == i)
    feather2[, tmp[1:length(tmp)]] <- NULL # remove all occurence!
}

#write_feather(feather1, sink = "cisTargets/hg38-500bp.feather", version = 1)
#write_feather(feather2, sink = "cisTargets/hg38-tss-centered-10kb.feather", version = 1)

atlas <- readRDS("/mnt/data01/yuanzhen/01.Vertebrate_cell_evo/01.data/02.atlas/human/human_atlas_test_subset_100k.rds")
cell2cluster <- read.table("/mnt/data01/yuanzhen/01.Vertebrate_cell_evo/01.data/02.atlas/human_mouse_lizard/s4.harmony_final/cell_clusters.integrated.txt", header = T)
atlas@meta.data$species <- "Human"
atlas@meta.data$integrated_cluster <- vapply(colnames(atlas), FUN = function(x){
                                                 tmp <- cell2cluster[cell2cluster$cell_id == x, "clusters"]
                                                 if (length(tmp) == 0){
                                                     tmp <- 0
                                                 }
                                                 return(tmp)
        }, FUN.VALUE = double(1))

exprMat <- GetAssayData(atlas, layer = "counts")
rownames(exprMat) <- setNames(vapply(rownames(exprMat), FUN = function(g) {
                                  if (length(Ensembl[[g]]) == 1){
                                      Ensembl[[g]] } else { g }}, FUN.VALUE=character(1)), NULL)

LogexprMat <- GetAssayData(atlas, layer = "data")
rownames(LogexprMat) <- setNames(vapply(rownames(LogexprMat), FUN = function(g) {
                                    if (length(Ensembl[[g]]) == 1){
                                        Ensembl[[g]] } else { g }}, FUN.VALUE=character(1)), NULL)

cellInfo <- atlas@meta.data
cellInfo$CellType <- cellInfo$supercluster

saveRDS(exprMat, file = "exprMat.rds")
saveRDS(LogexprMat, file = "LogexprMat.rds")
saveRDS(cellInfo, file = "cellInfo.rds")
