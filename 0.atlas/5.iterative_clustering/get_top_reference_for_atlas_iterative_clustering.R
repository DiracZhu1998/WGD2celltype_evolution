# script for getting ratio of major cell type and subtype labelling in our iterative clusters
suppressMessages({
    library(dplyr)
    library(Seurat)
})
source('ClusterPlotFunctions.R')

process_table <- function(obj, major, sub){
    table1 <- calculate_ratios(obj, order = levels(obj@meta.data$iter.clustering), major)
    table2 <- calculate_ratios(obj, order = levels(obj@meta.data$iter.clustering), sub)
    table1 <- table1 %>% group_by(c) %>% filter(ratio == max(ratio)) %>% slice(1) %>% ungroup()
    table2 <- table2 %>% group_by(c) %>% filter(ratio == max(ratio)) %>% slice(1) %>% ungroup()
    # in one case, max has two lines, same ratio
    table <- cbind(table1, table2[,-c(1,2)])
    colnames(table) <- c('Cluster','Cell number','Major labelling of reference',
                         '%Major labelling of reference', 'subtype labelling of reference', '%subtype labelling of reference')
    return(table)
}

Pmar <- readRDS("Pmar/Pmar.count_mean.dataset.rds")
Pvit <- readRDS("Pvit/Pvit.count_mean.dataset.rds")
Mmus <- readRDS("Mmus/Mmus.count_mean.dataset.rds")
Hsap <- readRDS("Hsap/Hsap.count_mean.dataset.rds")

Pmar_table <- process_table(Pmar, "cell_type", "cell_types_def")
Pvit_table <- process_table(Pvit, "cell_type", "Taxonomy3")
Mmus_table <- process_table(Mmus, "cell_type", "ClusterName")
Hsap_table <- process_table(Hsap, "cell_type", "cellname")

write.table(x = Pmar_table, file = "Pmar/Pmar.iter_clustering.pct_of_reference_annotation.count_mean.txt", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(x = Pvit_table, file = "Pvit/Pvit.iter_clustering.pct_of_reference_annotation.count_mean.txt", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(x = Mmus_table, file = "Mmus/Mmus.iter_clustering.pct_of_reference_annotation.count_mean.txt", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(x = Hsap_table, file = "Hsap/Hsap.iter_clustering.pct_of_reference_annotation.count_mean.txt", row.names = F, col.names = T, sep = "\t", quote = F)
