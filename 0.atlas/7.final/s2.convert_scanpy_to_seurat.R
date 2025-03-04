library(sceasy)
library(Seurat)
library(reticulate)
use_python("/home/zoo/ball6395/software/miniconda3/envs/scpy39/bin/python3.9")
loompy <- reticulate::import('loompy')


parser = argparse::ArgumentParser(description="Test paralog importance on cell marker genes")
parser$add_argument('-I','--input', help='input h5ad')
parser$add_argument('-O','--out', help='output rds')
args = parser$parse_args()

sceasy::convertFormat(args$input, from="anndata", to="seurat",
                      outFile=args$out, main_layer = 'data')

obj <- readRDS(args$out)
obj <- UpdateSeuratObject(obj)

saveRDS(obj, file = args$out)
