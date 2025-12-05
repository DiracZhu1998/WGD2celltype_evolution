suppressMessages(suppressWarnings({
    library(Seurat)
    library(ggplot2)
    library(dplyr)
    library(cowplot)
    library(plyr)
}))
    
    


suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(cowplot)))
suppressMessages(suppressWarnings(library(plyr)))
require(harmony)

options(Seurat.object.assay.version = "v5")

## ---------------------------
##
## Script name:
##
## Purpose of script: Integration (logN/SCT + harmony correction + Downstream Seurat) and clustering based on given parameters
##
## Author: Yuanzhen Zhu
##
## Date Created: 2022-07-15, Data Modified: 2022-07-21, 2023-11-25
##
## Version: V1.1
##
## ---------------------------
##
## Notes: This script can be universal and general used for Seurat v5

parser = argparse::ArgumentParser(description = "Script to integrate scRNA data")
parser$add_argument('-I', '--input', help = 'Input file containing path of RDS file needing integration')
parser$add_argument('-F', '--nf', help = 'Number of HVGs')
parser$add_argument('-H', '--hvg', help = 'scale and PCA on HVG or all genes, HVG is T, all genes F')
parser$add_argument('-N', '--normalize', help = 'RNA or SCT')
parser$add_argument('-D', '--dim', help = 'PCA dim usage and harmony used PCA dim')
parser$add_argument('-V', '--var', help = 'Input confounding variable, separate by ;(no spaces), default is orig.ident')
parser$add_argument('-T', '--theta', help = 'default = 2,  theta=0 does not encourage any diversity. Larger values of theta result in more diverse clusters')
parser$add_argument('-L', '--lambda', help = 'default = 1, Bigger values protect against over correction')
parser$add_argument('-C', '--nclust', help = 'nclust for harmony running, default is 50, if you assume your data contain more than 50 cell types, please set another value')
parser$add_argument('-R', '--res', help = 'Map resolution usage')
parser$add_argument('-O', '--out', help = 'Out directory')
parser$add_argument('-M', '--marker', help = 'input marker list, containing marker gene   annotation')
args = parser$parse_args()



nf.usage <- as.numeric(if(!is.null(args$nf)) args$nf else 3000)
scale_on <- as.character(if(!is.null(args$hvg)) args$hvg else "T")
dim.usage <- as.numeric(if(!is.null(args$dim)) args$dim else 30)
res.usage <- as.numeric(if(!is.null(args$res)) args$res else 0.6)
var.usage <- as.character(if(!is.null(args$var)) args$var else "orig.ident")
theta.usage <- as.double(if(!is.null(args$theta)) args$theta else 2)
lambda.usage <- as.numeric(if(!is.null(args$lambda)) args$lambda else 1)

var.usage <- strsplit(var.usage, split = ";")[[1]]
nclust.usage <- as.numeric(if(!is.null(args$nclust)) args$nclust else 50)
n.usage <- if (!is.null(args$normalize)){
    if (args$normalize %in% c('RNA', "SCT")){
        args$normalize
    } else {'RNA'}
} else { 'RNA' }

# read input RDS, and for loop it
path.list <- read.table(args$input)
id.list <- c()
for (i in 1:nrow(path.list)){
    abc = strsplit(path.list[i,1], split = "/")
    id <- abc[[1]][lengths(abc)]
    del_id <- path.list[i,2]
    assign(del_id, readRDS(path.list[i,1]))
    id.list <- c(id.list, del_id)
}

# build list containing objects
rds.list <- lapply(id.list, get)
print(id.list)
names(rds.list) <- id.list

Combined <- merge(rds.list[[1]], y = as.vector(rds.list)[-1])
Combined <- subset(Combined, subset = nFeature_RNA >= 200)
dim(Combined)

DefaultAssay(Combined) <- "RNA"
Combined <- Combined %>% NormalizeData(verbose = FALSE)  %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = nf.usage, assay = n.usage)
if (scale_on == "T"){
    Combined <- Combined %>% 
        ScaleData(assay = n.usage, verbose = FALSE) %>%
        RunPCA(npcs = dim.usage, verbose = FALSE, assay = n.usage)
} else if (scale_on == "F"){
    Combined <- Combined %>% 
        ScaleData(features = rownames(Combined), assay = n.usage, verbose = FALSE) %>%
        RunPCA(features = rownames(Combined), npcs = dim.usage, verbose = FALSE, assay = n.usage)
}

Combined <- Combined %>% IntegrateLayers(method = HarmonyIntegration,
    key = var.usage, npcs = dim.usage, nclust = nclust.usage, 
    max.iter.harmony = 10, theta = theta.usage, lamdba = lambda.usage,
    .options = harmony_options(max.iter.cluster = 30, epsilon.harmony = -Inf, epsilon.cluster = -Inf), 
    orig.reduction = "pca", new.reduction = "harmony", verbose = TRUE)

Combined[["RNA"]] <- JoinLayers(Combined[["RNA"]])

Combined <- RunUMAP(object = Combined, assay = n.usage, reduction = "harmony", dims = 1:(dim.usage - 5))
Combined <- FindNeighbors(object = Combined, assay = n.usage, reduction = "harmony", dims = 1:(dim.usage - 5))
Combined <- FindClusters(object = Combined, resolution = res.usage, graph.name = paste0(n.usage, "_snn"))

# pca then heatmap 
DefaultAssay(Combined) <- "RNA"
pdf(paste0(args$out,"/DimHeatmap.pdf"), width = 30, height = 20)
DimHeatmap(Combined, dims = 1:dim.usage, cells = 100, balanced = TRUE, assays = n.usage)
dev.off()

p1 <- DimPlot(Combined, reduction = "umap", group.by = var.usage)  
p2 <- DimPlot(Combined, reduction = "umap", label = T, pt.size=0.5)
#p3 <- DimPlot(Combined, reduction = "tsne", group.by = var.usage) + NoLegend()
#p4 <- DimPlot(Combined, reduction = "tsne", label = T, pt.size=0.5)



pdf(paste0(args$out, "/", n.usage, "_harmony_cell_cluster.umap.pdf"), height = 10, width = 20)
plot_grid(p1, p2, ncol = 2)
dev.off()

pdf(paste0(args$out, "/", n.usage, "_harmony_cell_cluster_umi.pdf"), height=10, width=10)
p5 <- ggplot(as.data.frame(Combined@reductions$umap@cell.embeddings), aes(x=umap_1, y=umap_2, color=Combined@meta.data$nCount_RNA)) + geom_point(size=0.5) + theme(legend.position = "none") + theme_bw() + labs(color="UMIs") + scale_color_viridis_c(direction = -1)
#p6 <- ggplot(as.data.frame(Combined@reductions$tsne@cell.embeddings), aes(x=tSNE_1, y=tSNE_2, color=Combined@meta.data$nCount_RNA)) + geom_point(size=0.5) + theme(legend.position = "none") + theme_bw() + labs(color="UMIs") + scale_color_viridis_c(direction = -1)
#plot_grid(p5, p6)
p5
dev.off()

pdf(paste0(args$out,"/", n.usage, "_harmony_cell_cluster_gene.pdf"),height=10,width=10)
p7 <- ggplot(as.data.frame(Combined@reductions$umap@cell.embeddings), aes(x=umap_1, y=umap_2, color=Combined@meta.data$nFeature_RNA)) + geom_point(size=0.5) + theme(legend.position = "none") + theme_bw() + labs(color="genes") + scale_color_viridis_c(direction = -1)
#p8 <- ggplot(as.data.frame(Combined@reductions$tsne@cell.embeddings), aes(x=tSNE_1, y=tSNE_2, color=Combined@meta.data$nFeature_RNA)) + geom_point(size=0.5) + theme(legend.position = "none") + theme_bw() + labs(color="genes") + scale_color_viridis_c(direction = -1)
#plot_grid(p7, p8)
p7
dev.off()

#markers <- FindAllMarkers(Combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#write.table(markers,paste0(args$out,"/", n.usage, "_harmony_markers.txt"),sep = '\t',row.names =FALSE, col.names =TRUE,quote =FALSE) 
#top <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
#tmp <- WhichCells(Combined, idents = levels(Combined), downsample = 1000)
#tmp <- subset(Combined, cells = tmp)
#pdf(paste0(args$out,"/top_marker_heatmap.pdf"), width = 40, height = 20)
#DoHeatmap(
#  tmp,
#  features = top$gene,
#  cells = NULL,
#  group.by = "ident",
#  group.bar = TRUE,
#  group.colors = NULL,
#  disp.min = -2.5,
#  disp.max = NULL,
#  slot = "scale.data",
#  assay = n.usage,
#  label = TRUE,
#  size = 5.5,
#  hjust = 0,
#  angle = 45,
#  raster = TRUE,
#  draw.lines = TRUE,
#  lines.width = NULL,
#  group.bar.height = 0.02,
#  combine = TRUE
#)
#dev.off()

#### customu set ####

if(!is.null(args$marker)){
    sample_marker <- read.table(args$marker,header=F)
    colnames(sample_marker) = c("marker", "anno")
    plot_theme<-theme(panel.background = element_blank(),axis.line = element_line(size=0.1),axis.ticks = element_line(size=0.1),axis.text = element_text(size=7),axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),axis.ticks.length = unit(1, "pt"),plot.title = element_text(size = 7, face = "bold",margin=margin(0,0,4,0)),plot.margin = unit(c(0,0,0,0),"cm"),legend.key.size = unit(0.1, 'cm'),legend.key.height = unit(0.1, 'cm'),legend.key.width = unit(0.1, 'cm'),legend.title = element_text(size=6),legend.text = element_text(size=6))
    dot_p <- DotPlot(Combined, assay="RNA", features=sample_marker$marker, cols = c("grey","red"), dot.scale=3) + labs(x="",y="") + plot_theme
    ggsave(filename = paste0(args$out, "/Marker.DotPlot.pdf"), plot = dot_p, device = "pdf", width = round(nrow(sample_marker)/8+2), height = 4)
}

saveRDS(Combined, paste0(args$out,"/", n.usage, "_harmony_integrated.RDS"))

