suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(future)))
suppressMessages(suppressWarnings(library(dplyr)))
### Get the parameters
parser = argparse::ArgumentParser(description="Find markers with Seurat")
parser$add_argument('-I','--input', help='input Rdata')
parser$add_argument('-L','--label', help='input Seurat Idents label within metadata')
parser$add_argument('-O','--out', help='input output filename prefix')
args = parser$parse_args()

set.seed(42)

obj <- readRDS(args$input)
#obj <- SCTransform(object = obj, verbose = FALSE, return.only.var.genes = FALSE)
# let's try SAM normalisation value first
sampled_df <- obj@meta.data
sampled_df$cellname <- rownames(sampled_df)
sampled_df <- sampled_df %>% group_by(args$label) %>% slice_sample(n = 3000) %>% ungroup()
obj <- subset(obj, cells = as.character(sampled_df$cellname))

Idents(obj) <- args$label
plan("multisession", workers = 8)
Markers1 <- FindAllMarkers(obj, test.use = "roc", only.pos = T)
Markers2 <- FindAllMarkers(obj, test.use = "wilcox", only.pos = T)
Markers2 <- Markers2[Markers2$p_val_adj < 0.01, ]
write.table(x = Markers1, file = paste0(args$out, '.roc.txt'), sep = "\t", quote = F, row.names = F, col.names = T)
write.table(x = Markers2, file = paste0(args$out, '.wilcox.FDR001.txt'), sep = "\t", quote = F, row.names = F, col.names = T)

