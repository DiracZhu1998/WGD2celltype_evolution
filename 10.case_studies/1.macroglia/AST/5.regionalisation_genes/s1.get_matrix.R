suppressMessages({
    library(scrattch.hicat)
    library(dendextend)
    library(dplyr)
    library(Seurat)
    library(DESeq2)
})

parser = argparse::ArgumentParser(description="get DESeq2 normalised matrix")
parser$add_argument('-I','--input', help='input seurat object amniote/vertebrate and ortholog/metagene')
parser$add_argument('-L','--label', help='input label')
parser$add_argument('-T','--table', help='input region conversion table')
parser$add_argument('-R','--region', help='input region column name')
args = parser$parse_args()

label <- args$label
obj <- readRDS(args$input)
obj@meta.data$`Refined family` <- as.character(obj@meta.data$`Refined family`)
selected <- c('Astrocytes', 'Diencephalon GABAergic neurons', 'Mesencephalon GABAergic neurons',
              'Telencephalon GABAergic neurons','Rhombencephalon GABAergic neurons',
              'Rhombencephalon glutamatergic neurons', 'Telencephalon glutamatergic neurons',
              'Mesencephalon glutamatergic neurons','Diencephalon glutamatergic neurons')
obj <- subset(obj, cells = colnames(obj)[which(obj@meta.data$`Refined family` %in% selected)])

table <- read.delim(args$table, header = F)
location <- table$V1
names(location) <- table$V2
region <- args$region
obj@meta.data$location <- vapply(as.character(obj@meta.data[,region]), FUN = function(x){
                        names(location)[location == x]}, FUN.VALUE = character(1))
obj@meta.data[which(obj@meta.data$`Refined family` == 'Astrocytes'), 'Refined family'] <- 
        paste(as.character(obj@meta.data[which(obj@meta.data$`Refined family` == 'Astrocytes'), 'location']), 'Astrocytes', sep = ' ')

tmp <- obj@meta.data
tmp$cellID <- rownames(obj@meta.data)
tmp <- tmp %>% group_by(`Refined family`) %>% sample_n(min(2000, n()), replace = FALSE) %>% ungroup()
obj <- subset(obj, cells = tmp$cellID)
saveRDS(obj, paste0(label, '.subsetted.rds'))

cl.fact <- setNames(factor(obj@meta.data[,'Refined family']), rownames(obj@meta.data))
dgeMatrix_count <- as.matrix(obj@assays$RNA$counts)
cl.sum <- get_cl_sums(dgeMatrix_count, cl.fact)
metadata <- obj@meta.data
experiment <- unique(data.frame(Refined.family = metadata$`Refined family`))
experiment$location <- vapply(experiment$`Refined.family`, FUN = function(x){
                                  strsplit(x, split = ' ')[[1]][1]
                        }, FUN.VALUE = character(1))
experiment$celltype <- vapply(experiment$`Refined.family`, FUN = function(x){
                                  paste0(strsplit(x, split = ' ')[[1]][-1], collapse  = '_')
                        }, FUN.VALUE = character(1))


dataset <- DESeqDataSetFromMatrix(countData = cl.sum, colData = experiment, design = ~ location + celltype)
dds <- DESeq(dataset)
rld <- vst(dds)
mat <- assay(rld)
write.table(mat, paste0(label,'.DEseq2.txt'), col.names = T, row.names = T, sep = '\t', quote = F)

