library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(igraph)
require(VennDiagram)


parser = argparse::ArgumentParser(description="Test paralog importance on cell marker genes")
parser$add_argument('-I','--input', help='input Seurat rds')
parser$add_argument('-P','--paralog', help='input Seurat paralog file retrieved from Ensembl')
parser$add_argument('-M','--marker', help='input marker gene list')
parser$add_argument('-S','--species', help='input species name')
parser$add_argument('-O','--out', help='input output path')
args = parser$parse_args()

species <- args$species
if (!dir.exists(args$out)) {
  dir.create(args$out, recursive = TRUE)
  print(paste("Path:", args$out, " created."))
}

########## Paralogs #############
sobj <- readRDS(args$input)
DefaultAssay(sobj) <- "RNA"
# read paralog file
paralogs <- read.delim(args$paralog, header = T, sep = ",")
# remove species name for convenience
colnames(paralogs) <- vapply(colnames(paralogs), FUN = function(x){
  gsub(paste0("\\.?",species,"\\.?"), "", x) 
}, FUN.VALUE = character(1))

paralogs <- paralogs[which(paralogs$Gene.stable.ID != "" & paralogs$paralogue.gene.stable.ID != ""), ]
paralogs[paralogs$Gene.name == "", "Gene.name"] <- paralogs[paralogs$Gene.name == "", "Gene.stable.ID"]
paralogs[paralogs$paralogue.associated.gene.name == "", "paralogue.associated.gene.name"] <- 
  paralogs[paralogs$paralogue.associated.gene.name == "", "paralogue.gene.stable.ID"] 

# get family level paralogs
g <- graph_from_data_frame(paralogs[,c(2,4)], directed = FALSE)
components <- clusters(g)$membership
paralog_pairs <- split(names(components), components)

p0 <- ggplot(data = NULL, aes(x = sapply(paralog_pairs, function(x) length(x) ))) +
  geom_histogram(binwidth = 1, fill = "grey", color = "black") +
  labs(title = "Distribution of the number of paralogs \n in paralog families") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(args$out, '/', species,".whole_paralog_family_size.pdf"), p0, width = 8, height = 5)

##### paralog ratio within DEGs #########
subcluster_markers <- read.delim(args$marker, header = T, sep = "\t")
stats <- subcluster_markers %>% as_tibble() %>% group_by(cluster) %>% 
  summarise(total = length(cluster), 
    n_paralogs = sum(gene %in% unique(c(paralogs$Gene.name, paralogs$paralogue.associated.gene.name))), # calculate number of paralogs
    n_paralog_pairs = length(unique(vapply(gene[gene %in% unique(c(paralogs$Gene.name, paralogs$paralogue.associated.gene.name))], 
        FUN =  function(y) {which(sapply(paralog_pairs, function(x) y %in% x))}, FUN.VALUE = double(1)))),
    "paralogs%" = n_paralogs/total,
    "families_divided_by_paralogs%" = n_paralog_pairs/n_paralogs
    ) %>%
  ungroup()

paralog_number_in_ourdata <- sum(unique(c(paralogs$Gene.name, paralogs$paralogue.associated.gene.name)) %in%  rownames(sobj))
paralog_ratio_in_ourdata <- round(paralog_number_in_ourdata/length(rownames(sobj)),2)
non_paralog_number_in_ourdata <- sum(!rownames(sobj) %in% unique(c(paralogs$Gene.name, paralogs$paralogue.associated.gene.name)))
  

p1 <- stats %>% as.data.frame() %>% ggplot(aes(`paralogs%`)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.6) + 
  geom_vline(aes(xintercept=mean(`paralogs%`)),
             color="blue", linetype="dashed", linewidth=1) + 
  geom_vline(aes(xintercept=paralog_ratio_in_ourdata),
             color="red", linetype="dashed", linewidth=1) + xlim(0,1)

ggsave(paste0(args$out, '/', species,".paralog_ratio_inDEGs.pdf"), p1, width = 6, height = 2)

### get stats of paralogs of one-one-one orthologs ###
paralogs_in_our_dataset <- paralogs[paralogs$Gene.name %in% rownames(sobj) & paralogs$paralogue.associated.gene.name %in% rownames(sobj), ]
g <- graph_from_data_frame(paralogs_in_our_dataset[,c(2,4)], directed = FALSE)
components <- clusters(g)$membership
paralog_pairs_in_our_dataset <- split(names(components), components)
######
p2 <- stats %>% as.data.frame() %>% ggplot(aes(`families_divided_by_paralogs%`)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.6) +
  geom_vline(aes(xintercept=mean(`families_divided_by_paralogs%`)),
             color="blue", linetype="dashed", linewidth=1) +
  geom_vline(aes(xintercept=length(paralog_pairs_in_our_dataset)/length(unlist(paralog_pairs_in_our_dataset))), # you should use
             color="red", linetype="dashed", linewidth=1) + xlim(0,1)
ggsave(paste0(args$out, '/', species,".paralog_families_divided_paralogs_inDEGs.pdf"), p2, width = 6, height = 2)


######### Label all genes with different DEGs and paralogs as T/F #########
stats_for_chi2 <- data.frame(gene = rownames(sobj), DEG = NA, paralog = NA)
stats_for_chi2$DEG <- stats_for_chi2$gene %in% subcluster_markers$gene
stats_for_chi2$paralog <- stats_for_chi2$gene %in% unique(c(paralogs$Gene.name, paralogs$paralogue.associated.gene.name))

######## Venn plot show DEGs and paralogs #########
venn.plot <- venn.diagram(
  x = list(set1 = stats_for_chi2[stats_for_chi2$DEG, "gene"], 
           set2 = stats_for_chi2[stats_for_chi2$paralog, "gene"]
           ),
  category.names = c("DEG", "paralogs"),
  filename = paste0(args$out, '/', species,".DEG_paralogs.Venn.png"),
  col=c("#440154ff", '#21908dff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  output = TRUE
)
############## Chi squared test of paralogs and DEGs ##################
chi <- data.frame(matrix(data = c(sum(stats_for_chi2$DEG & stats_for_chi2$paralog),
                                        sum(stats_for_chi2$DEG & !stats_for_chi2$paralog),
                                        sum(!stats_for_chi2$DEG & stats_for_chi2$paralog),
                                        sum(!stats_for_chi2$DEG & !stats_for_chi2$paralog) 
                                        ), ncol = 2, byrow = T))

rownames(chi) <- c("DEGs", "non-DEGs")
colnames(chi) <- c("paralogs", "non-paralogs")
write.table(x = chi, file = paste0(args$out, '/', species, ".paralog_DEGs.chi2.txt"), sep = "\t", quote = F,
            col.names = T, row.names = T)
write.table(x = capture.output(print(chisq.test(chi, correct = F)))[5], 
      file = paste0(args$out, '/', species, ".paralog_DEGs.chi2.txt"), sep = "\t", quote = F,
      col.names = F, row.names = F, append = T)
write(paste0("OR = ",(chi[1,1]*chi[2,2])/(chi[1,2]*chi[2,1])), file = paste0(args$out, '/', species, ".paralog_DEGs.chi2.txt"), sep = "\t", append = T)
