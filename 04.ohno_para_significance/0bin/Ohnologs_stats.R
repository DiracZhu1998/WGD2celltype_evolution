library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(igraph)
require(VennDiagram)

parser = argparse::ArgumentParser(description="Test paralog importance on cell marker genes")
parser$add_argument('-I','--input', help='input Seurat rds')
parser$add_argument('-P','--ohnolog', help='input Seurat ohnolog file retrieved from Ohnolog v2')
parser$add_argument('-M','--marker', help='input marker gene list')
parser$add_argument('-S','--species', help='input species name')
parser$add_argument('-O','--out', help='input output path')
args = parser$parse_args()

species <- args$species
if (!dir.exists(args$out)) {
  dir.create(args$out, recursive = TRUE)
  print(paste("Path:", args$out, " created."))
}

########## Ohnologs #############
sobj <- readRDS(args$input)
DefaultAssay(sobj) <- "RNA"
# read paralog file
ohnologs <- read.delim(args$ohnolog, header = T)
### get family level ohnologs ###
g <- graph_from_data_frame(ohnologs[,1:2], directed = FALSE)
components <- clusters(g)$membership
ohnolog_pairs <- split(names(components), components)

p0 <- ggplot(data = NULL, aes(x = sapply(ohnolog_pairs, function(x) length(x) ))) +
  geom_histogram(binwidth = 1, fill = "grey", color = "black") +
  labs(title = "Distribution of the number of ohnologs \n in ohnolog families") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(args$out, '/', species,".whole_ohnolog_family_size.pdf"), p0, width = 8, height = 5)

##### ohnolog ratio within DEGs #########
subcluster_markers <- read.delim(args$marker, header = T, sep = "\t")
stats <- subcluster_markers %>% as_tibble() %>% group_by(cluster) %>%
  summarise(total = length(cluster),
            n_ohnologs = sum(gene %in% unique(c(ohnologs$Symbol1, ohnologs$Symbol2))), # calculate number of ohnologs
            n_ohnolog_pairs = length(unique(vapply(gene[gene %in% unique(c(ohnologs$Symbol1, ohnologs$Symbol2))],
                                                   FUN =  function(y) {which(sapply(ohnolog_pairs, function(x) y %in% x))}, FUN.VALUE = double(1)))),
            "ohnologs%" = n_ohnologs/total,
            "families_divided_by_ohnologs%" = n_ohnolog_pairs/n_ohnologs
            ) %>%
  ungroup()

ohnolog_number_in_ourdata <- sum(unique(c(ohnologs$Symbol1, ohnologs$Symbol2)) %in%  rownames(sobj))
ohnolog_ratio_in_ourdata <- round(ohnolog_number_in_ourdata/length(rownames(sobj)),2)
non_ohnolog_number_in_ourdata <- sum(!rownames(sobj) %in% unique(c(ohnologs$Symbol1, ohnologs$Symbol2)))

p1 <- stats %>% as.data.frame() %>% ggplot(aes(`ohnologs%`)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.6) +
  geom_vline(aes(xintercept=mean(`ohnologs%`)),
             color="blue", linetype="dashed", linewidth=1) +
  geom_vline(aes(xintercept=ohnolog_ratio_in_ourdata),
             color="red", linetype="dashed", linewidth=1) + xlim(0,1)

ggsave(paste0(args$out, '/', species,".ohnolog_ratio_inDEGs.pdf"), p1, width = 6, height = 2)

### get stats of ohnologs ###
ohnologs_in_our_dataset <- ohnologs[ohnologs$Symbol1 %in% rownames(sobj) & ohnologs$Symbol2 %in% rownames(sobj), ]
g <- graph_from_data_frame(ohnologs_in_our_dataset[,c(3,4)], directed = FALSE)
components <- clusters(g)$membership
ohnolog_pairs_in_our_dataset <- split(names(components), components)

######
p2 <- stats %>% as.data.frame() %>% ggplot(aes(`families_divided_by_ohnologs%`)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.6) +
  geom_vline(aes(xintercept=mean(`families_divided_by_ohnologs%`)),
             color="blue", linetype="dashed", linewidth=1) +
  geom_vline(aes(xintercept=length(ohnolog_pairs_in_our_dataset)/ohnolog_number_in_ourdata),
             color="red", linetype="dashed", linewidth=1) + xlim(0,1)

ggsave(paste0(args$out, '/', species,".ohnolog_families_divided_ohnologs_inDEGs.pdf"), p2, width = 6, height = 2)

######### Label all genes with different DEGs and paralogs as T/F #########
stats_for_chi2 <- data.frame(gene = rownames(sobj), DEG = NA, ohnolog = NA)
stats_for_chi2$DEG <- stats_for_chi2$gene %in% subcluster_markers$gene
stats_for_chi2$ohnolog <- stats_for_chi2$gene %in% unique(c(ohnologs$Symbol1, ohnologs$Symbol2))

######## Venn plot show DEGs and ohnologs #########
venn.plot <- venn.diagram(
  x = list(set1 = stats_for_chi2[stats_for_chi2$DEG, "gene"],
           set2 = stats_for_chi2[stats_for_chi2$ohnolog, "gene"]
  ),
  category.names = c("DEG", "ohnologs"),
  filename = paste0(args$out, '/', species,".DEG_ohnologs.Venn.png"),
  col=c("#440154ff", '#21908dff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  output = TRUE
)
############## Chi squared test of ohnologs and DEGs ##################
chi <- data.frame(matrix(data = c(sum(stats_for_chi2$DEG & stats_for_chi2$ohnolog),
                                        sum(stats_for_chi2$DEG & !stats_for_chi2$ohnolog),
                                        sum(!stats_for_chi2$DEG & stats_for_chi2$ohnolog),
                                        sum(!stats_for_chi2$DEG & !stats_for_chi2$ohnolog)
                                        ), ncol = 2, byrow = T))


rownames(chi) <- c("DEGs", "non-DEGs")
colnames(chi) <- c("ohnologs", "non-ohnologs")

write.table(x = chi, file = paste0(args$out, '/', species, ".ohnolog_DEGs.chi2.txt"), sep = "\t", quote = F,
            col.names = T, row.names = T)
write.table(x = capture.output(print(chisq.test(chi, correct = F)))[5], 
            file = paste0(args$out, '/', species, ".ohnolog_DEGs.chi2.txt"), sep = "\t", quote = F,
            col.names = F, row.names = F, append = T)
write(paste0("OR = ",(chi[1,1]*chi[2,2])/(chi[1,2]*chi[2,1])), file = paste0(args$out, '/', species, ".ohnolog_DEGs.chi2.txt"), sep = "\t", append = T)







