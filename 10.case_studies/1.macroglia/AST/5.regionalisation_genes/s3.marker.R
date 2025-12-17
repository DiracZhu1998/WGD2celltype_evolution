suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(future)))

parser = argparse::ArgumentParser(description="Find markers")

parser$add_argument('-I','--input', help='input Rdata')
parser$add_argument('-O','--out', help='input output filename prefix')
args = parser$parse_args()


obj <- readRDS(args$input)
Idents(obj) <- 'Refined family'
obj@meta.data$type <- vapply(obj@meta.data$`Refined family`, FUN = function(x){
                                                         paste0(strsplit(x, split = '\\ ')[[1]][-1], collapse = "_")
                                                             }, FUN.VALUE = character(1))

AST <- subset(obj, cells = colnames(obj)[which(obj@meta.data$type == 'Astrocytes')])
GABA <- subset(obj, cells = colnames(obj)[which(obj@meta.data$type == 'GABAergic_neurons')])
Glut <- subset(obj, cells = colnames(obj)[which(obj@meta.data$type == 'glutamatergic_neurons')])

plan("multisession", workers = 8)
Markers1 <- FindAllMarkers(AST, test.use = "roc", only.pos = T)
Markers2 <- FindAllMarkers(GABA, test.use = "roc", only.pos = T)
Markers3 <- FindAllMarkers(Glut, test.use = "roc", only.pos = T)
Markers <- rbind(Markers1, Markers2, Markers3)
#Markers <- Markers[Markers$p_val_adj < 0.01, ]
write.table(x = Markers, file = paste0(args$out, '.roc.txt'), sep = "\t", quote = F, row.names = F, col.names = T)
