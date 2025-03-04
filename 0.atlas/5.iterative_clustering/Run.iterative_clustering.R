suppressMessages({
    library(scrattch.hicat)
    library(dendextend)
    library(dplyr)
    library(matrixStats)
    library(Matrix)
    library(Seurat)
    library(RColorBrewer)
    library(scrattch.vis)
    library(ggdendro)
    library(tidyr)
    library(ggExtra)
})

args = commandArgs(trailingOnly=TRUE)

obj <- readRDS(args[1])
label <- args[3]
#obj <- subset(obj, cells = sample(colnames(obj), size = 5000))

dgeMatrix_count <- as.matrix(obj@assays$RNA$counts)
dgeMatrix_cpm <- cpm(dgeMatrix_count)
norm.dat <- log2(dgeMatrix_cpm + 1)
gene.counts <- log2(colSums(norm.dat > 0))
nUMI <- log2(colSums(dgeMatrix_count))
if (args[2] == 'NA'){
    rm.eigen <- as.matrix(cbind(gene.counts,nUMI)) #for human
} else {
    perctMito <- obj@meta.data[, args[2]]
    rm.eigen <- as.matrix(cbind(gene.counts,nUMI,perctMito))
}

de.param <- de_param(padj.th     = 0.05, 
                     low.th      = 1, 
                     q1.th       = 0.4, 
                     q2.th       = NULL,
                     q.diff.th   = 0.5,
                     de.score.th = 100,
                     min.cells = 100,
                     min.genes = 6)

iter.result <- iter_clust(norm.dat, 
                          counts = norm.dat,
                          dim.method = "pca",
                          max.dim = 80,
                          de.param = de.param,
                          rm.eigen = rm.eigen,
                          method = "louvain",
                          prefix = label,
                          verbose = TRUE,
                          split.size = 350,
                          #sampleSize = 20000,
                          #max.cl.size = 500
                          )
warnings()

rd.dat <- t(dgeMatrix_count[iter.result$markers,])
merge.result <- scrattch.bigcat::merge_cl(norm.dat, 
                         cl = iter.result$cl,
                         rd.dat = rd.dat,
                         de.param = de.param, verbose = 1)

cat(length(unique(merge.result$markers)), " DE genes")
cat(length(unique(merge.result$cl))," Clusters\n")

pdf(paste0(label, '.pdf'))
obj@meta.data$iter.clustering <-as.character(merge.result$cl)
#colors <-  c("#ebcb2e", "#9ec22f", "#a9961b", "#cc3a1b", "#cc8778" , "#d14c8d", "#4cabdc", "#5ab793", "#e7823a","#e6bb9b", "#046c9a", "#4784a2" , "#4990c9")
DimPlot(obj, group.by = 'iter.clustering',
                label.size = 4) + NoLegend()
dev.off()

save.image(file = paste0(label,".RData"))
