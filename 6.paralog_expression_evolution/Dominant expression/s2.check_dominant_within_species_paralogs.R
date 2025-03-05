suppressMessages(suppressWarnings({
    library(Seurat)
    library(ggplot2)
    library(igraph)
    library(reshape2)
    library(tidyverse)
    library(ggpubr)
    library(rstatix)
}))


parser = argparse::ArgumentParser(description="Test paralog/ohnolog shift")
parser$add_argument('-I','--input', help='input processed pseudobulk exp/pct rds')
parser$add_argument('-P1','--paralog', help='input paralog pairs')
parser$add_argument('-P2','--ohnologs', help='input ohnolog pairs')
parser$add_argument('-L','--label', help='output file label')
parser$add_argument('-O','--out', help='input output path')
args = parser$parse_args()

### load atlas ###
output <- args$out
tmp <- paste0(output, '/boxplot')
if (!dir.exists(tmp)) {
    dir.create(tmp, recursive = TRUE)
    print(paste("Path:", output, " created."))
}
label <- args$label
mtx <- readRDS(args$input) 
paralogs <- read.table(args$paralog, header = T, sep = "\t")
paralogs <- paralogs[paralogs$Dup1 %in% rownames(mtx) & paralogs$Dup2 %in% rownames(mtx), ] # paralogs in our dataset
### get family level paralogs ###
g <- graph_from_data_frame(paralogs, directed = FALSE)
components <- components(g)$membership
paralog_family <- split(names(components), components)

ohnologs <- read.delim(args$ohnolog, header = T)
### get family level ohnologs ###
g <- graph_from_data_frame(ohnologs[,1:2], directed = FALSE)
ohnologs <- ohnologs[ohnologs$Ohno1 %in% rownames(mtx) & ohnologs$Ohno2 %in% rownames(mtx), ] # ohnologs in our dataset
components <- clusters(g)$membership
ohnolog_family <- split(names(components), components)

# filter and to get non-ohnolog paralogs
SSDparalogs <- paralogs[!paralogs$Dup1 %in% unique(c(ohnologs$Ohno1, ohnologs$Ohno2)),]
SSDparalogs <- SSDparalogs[!SSDparalogs$Dup2 %in% unique(c(ohnologs$Ohno1, ohnologs$Ohno2)),]
### get family level paralogs ###
g <- graph_from_data_frame(SSDparalogs, directed = FALSE)
components <- clusters(g)$membership
SSDparalog_family <- split(names(components), components)


title1 <- c('group1','group2','n1','n2','statistic','p','p.adj','p.adj.signif')
title2 <- c('Var1','variable','n','min','max','median','iqr','mean','sd','se','ci')
title3 <- c('family','n','statistic','df','p','method')
get_ <- function(paralog_family, type){
    lapply(seq_along(paralog_family), FUN = function(s){
           paralog_family[[s]] <- paralog_family[[s]][paralog_family[[s]] %in% rownames(mtx)]
            
           label = paste0(label, '.', type)
           # add title
           if (s == 1){
               write.table(data.frame(t(title1)), file = paste0(output,'/', label, ".pairwise_wilcox_test.txt"), 
                           sep = "\t", row.names = FALSE, col.names = FALSE, quote = F, append = F)
               write.table(data.frame(t(title2)), file = paste0(output,'/', label, ".mean.txt"),
                           sep = "\t", row.names = FALSE, col.names = FALSE, quote = F, append = F)
               write.table(data.frame(t(title3)), file = paste0(output,'/',label,'.friedman_test.txt'), 
                           sep = "\t", row.names = FALSE, col.names = FALSE, quote = F, append = F)

           }

           if (length(paralog_family[[s]]) >=2){
               tmp = melt(as.matrix(mtx[paralog_family[[s]], ]))
               tmp$Var1 <- factor(tmp$Var1) # gene
               tmp$Var2 <- factor(tmp$Var2) # cluster
               res.fried <- tmp %>% friedman_test(value ~ Var1 |Var2)
               res.fried$`.y.` <- paste0(paralog_family[[s]], collapse = ";")


               pwc <- tmp %>% wilcox_test(value ~ Var1, paired = TRUE, p.adjust.method = "bonferroni")
               mean <- tmp %>% group_by(Var1) %>% get_summary_stats(value, type = "common")
               
               write.table(data.frame(pwc)[,-1], file = paste0(output,'/', label, ".pairwise_wilcox_test.txt"), 
                           sep = "\t", row.names = FALSE, col.names = FALSE, quote = F, append = T)
               write.table(data.frame(mean), file = paste0(output,'/', label, ".mean.txt"),
                           sep = "\t", row.names = FALSE, col.names = FALSE, quote = F, append = T)
               write.table(data.frame(res.fried), file = paste0(output,'/',label,'.friedman_test.txt'), 
                           sep = "\t", row.names = FALSE, col.names = FALSE, quote = F, append = T)
               
               
               pwc <- pwc %>% add_xy_position(x = "Var1")
               
               if ( length(unique(tmp$Var1)) <= 2){ #only plot boxplot for paralogs <= 4
                   p <- ggpaired(tmp, x = "Var1", y = "value", add = "point", color = "Var1", id = 'Var2', line.color = "gray", line.size = 0.4, palette = "jco") + 
                       #ggboxplot(tmp, x = "Var1", y = "value", add = "point", color = "Var1", line.color = "gray", line.size = 0.4) + 
                       stat_pvalue_manual(pwc, hide.ns = TRUE) + 
                       labs(
                            subtitle = get_test_label(res.fried,  detailed = TRUE),
                            caption = get_pwc_label(pwc)
                            )  
                   ggsave(filename = paste0(output, "/boxplot/boxplot.", type, ".", paste( paralog_family[[s]] , collapse ="-"), ".pdf"), 
                          plot = p,  width = 5+round(length(paralog_family[[s]])/3), height = 8)
               }
           }
})
}

get_(ohnolog_family, type = "WGD")
get_(SSDparalog_family, type = "SSD")


