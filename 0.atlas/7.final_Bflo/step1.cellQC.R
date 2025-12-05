
## ---------------------------
##
## Script name:  
##
## Purpose of script: filter low quality cell and generate report graph and statistics
##
## Author: Yuanzhen Zhu
##
## Date Created: 2020-11-11, Data Modified: 2021-06-01, 2025-06-30
##
## Version: V2.0
##
## ---------------------------
##
## Notes: This script is universal and general used!
##
## ---------------------------

### Get the parameters

parser = argparse::ArgumentParser(description="Script to QC and Cluster scRNA data")
parser$add_argument('-I','--input', help = 'input raw matrix or 10X-like dir')
parser$add_argument('-F','--nf',help = 'HVG number set')
parser$add_argument('-D','--dim',help = 'PCA dim usage')
parser$add_argument('-P','--percentage',help = 'doublets percentage')
parser$add_argument('-R','--res',help = 'Map resolution usage')
parser$add_argument('-O','--out',help = 'out directory')
parser$add_argument('-S','--sample',help = 'sample sample name')
parser$add_argument('-G','--gene',help = 'file containing retain gene names')
args = parser$parse_args()

### Load library
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(DoubletFinder)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(patchwork)))

#EC <- readRDS(args$input)
EC.data <- Read10X(data.dir = args$input, gene.column = 1)

nf.usage <- as.numeric(if(!is.null(args$nf)) args$nf else 3000)
dim.usage <- as.numeric(if(!is.null(args$dim)) args$dim else 30)
doublets.percentage <- if(!is.null(args$percentage)) args$percentage else 0.05
doublets.percentage <- as.numeric(doublets.percentage)
res.usage <- as.numeric(if(!is.null(args$res)) args$res else 0.6)

### Creat Seurat object, basic filtering and statistics
EC <- CreateSeuratObject(EC.data, project = args$sample, min.cells = 3, min.features = 200)
# usually retain only protein_coding genes
if (!is.null(args$gene)){
    gene_ <- read.delim(args$gene, header = F)$V1
    EC <- subset(EC, features = gene_)
}
MTgenes <- rownames(EC)[!grepl("^bf|^Bflo", rownames(EC))]
#MTgenes <- rownames(EC)[grepl('^MT\\-', rownames(EC))]
MTgenes
EC[["percent.mt"]] <- PercentageFeatureSet(EC, features = MTgenes)

meanT <- sum(EC@meta.data$nCount_RNA)/nrow(EC@meta.data)
meanG <- sum(EC@meta.data$nFeature_RNA)/nrow(EC@meta.data)
meanMT <- sum(EC[["percent.mt"]])/nrow(EC@meta.data)
info <- c(meanT, meanG, meanMT)

pdf(paste0(args$out, "/count_mt.cor.pdf"))
FeatureScatter(EC, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()
pdf(paste0(args$out, "/count_gene.cor.pdf"))
FeatureScatter(EC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

p1 <- VlnPlot(EC, features = "nFeature_RNA",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#1B9E77") + NoLegend() + xlab("Gene") + labs(title="")+ theme(axis.text.x = element_blank())
p2 <- VlnPlot(EC, features = "nCount_RNA",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#D95F02") + NoLegend() + xlab("Transcript") + labs(title="")+ theme(axis.text.x = element_blank()) #+ ylim(0,15000)
p3 <- VlnPlot(EC, features = "percent.mt",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#7570B3") + NoLegend() + xlab("percent.mt") + labs(title="")+ theme(axis.text.x = element_blank())
p <-p1|p2|p3
pdf(paste0(args$out,"/QC_before.dis.pdf"))
print(p)
dev.off()

#ECmeta <- EC@meta.data[order(-EC@meta.data$nFeature_RNA),]
#n95 <- as.numeric(as.integer(nrow(ECmeta) * 0.05))
#n95features <- as.numeric(ECmeta[n95, "nFeature_RNA"])
# you need to manually tune below filtering
EC <- subset(EC, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10 & nCount_RNA < 10000  )
#if (type == 'sn'){
#    EC <- subset(EC, subset = nFeature_RNA > 400 & nFeature_RNA < 3500 & percent.mt < 5)
#} else if (type == 'sc'){
#    EC <- subset(EC, subset = nFeature_RNA > 400 & nFeature_RNA < 6500 & percent.mt < 5)
#}
p1 <- VlnPlot(EC, features = "nFeature_RNA",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#1B9E77") + NoLegend() + xlab("Gene") + labs(title="")+ theme(axis.text.x = element_blank())
p2 <- VlnPlot(EC, features = "nCount_RNA",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#D95F02") + NoLegend() + xlab("Transcript") + labs(title="")+ theme(axis.text.x = element_blank()) #+ ylim(0,15000)
p3 <- VlnPlot(EC, features = "percent.mt",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#7570B3") + NoLegend() + xlab("percent.mt") + labs(title="")+ theme(axis.text.x = element_blank())
p <-p1|p2|p3
pdf(paste0(args$out,"/QC_after.dis.pdf"))
print(p)
dev.off()


meanT <- sum(EC@meta.data$nCount_RNA)/nrow(EC@meta.data)
meanG <- sum(EC@meta.data$nFeature_RNA)/nrow(EC@meta.data)
meanMT <- sum(EC[["percent.mt"]])/nrow(EC@meta.data)
info <- rbind(info, c(meanT, meanG, meanMT))
rownames(info) <- c("before", "after")
colnames(info) <- c("Transcripts", "Genes", "pct.mt")
write.table(info, paste0(args$out, "/QC.mean.txt"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

EC <- NormalizeData(EC, normalization.method = "LogNormalize", scale.factor = 10000)
EC <- FindVariableFeatures(EC, selection.method = "vst", nfeatures = nf.usage)
EC <- ScaleData(EC)

pdf(paste0(args$out, "/VariableFeaturePlot.pdf"))
VariableFeaturePlot(EC)
dev.off()

### PCA, statistics
EC <- RunPCA(EC, npcs = dim.usage)
pdf(paste0(args$out,"/PCA_DimHeatmap.pdf"), height = 20, width = 20)
DimHeatmap(EC, dims = 1:dim.usage, cells = 100, balanced = TRUE)
dev.off()
EC <- JackStraw(EC, num.replicate = 100)
EC <- ScoreJackStraw(EC, dims = 1:20)
pdf(paste0(args$out, "/PCA_JackStrawPlot.pdf"))
JackStrawPlot(EC, dims = 1:20)
dev.off()
pdf(paste0(args$out, "/PCA_ElbowPlot.pdf"))
ElbowPlot(EC)
dev.off()

### Choose suitable dim for UMAP & T-SNE: dim reduce unfavored signal, determine percent of variation associated with each PC; calculate cumulative percents for each PCs; Determine the difference between variation of PC and subsequent PC;
pct <- EC[["pca"]]@stdev / sum(EC[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 80 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
pdf(paste0(args$out, "/PCA_choose.pdf"))
plot(EC[["pca"]]@stdev, xlab="PC", ylab = "SD explaineded (%)")
abline(v=pcs, col="red")
dev.off()
print( paste0("PCA dim: ", as.character(dim.usage), "; Favored choose dim 1:", as.character(pcs)))

### Find neighbors, UMAP and T-SNE
EC <- FindNeighbors(EC, dims = 1:(dim.usage-5))
EC <- FindClusters(EC, resolution = res.usage)
EC <- RunUMAP(EC, dims = 1:(dim.usage-5))
EC <- RunTSNE(EC, dims = 1:(dim.usage-5))

p1 <- DimPlot(EC, reduction = "umap")
p2 <- DimPlot(EC, reduction = "tsne")
p <- p1|p2
ggsave(filename = paste0(args$out, "/Umap_TSNE.pdf"), plot = p, device = "pdf", width = 15, height = 7)

saveRDS(EC,paste(args$out,"/",args$sample,"_QCed.rds",sep=""))


