suppressPackageStartupMessages({
        require(igraph)
            library(dplyr)
            library(ggvenn)
                library(tidyr)
                library(ggpubr)
                    library(rstatix)
})

species <- c('Hsap','Mmus','Pvit','Pmar')

duplicate_pairs <- readRDS("../2.paralog_switching/Combined.SSD_WGD.rds")

Hsap2Mmus <- read.table("/mnt/data01/yuanzhen/01.Vertebrate_cell_evo/02.gene_relationships/run4/results/Ortho_pipeline/OrthoFinder/RBH_results/RBH.Hsap.Mmus.tsv", header = T)
Mmus2Hsap <- read.table("/mnt/data01/yuanzhen/01.Vertebrate_cell_evo/02.gene_relationships/run4/results/Ortho_pipeline/OrthoFinder/RBH_results/RBH.Mmus.Hsap.tsv", header = T)
Hsap2Pvit <- read.table("/mnt/data01/yuanzhen/01.Vertebrate_cell_evo/02.gene_relationships/run4/results/Ortho_pipeline/OrthoFinder/RBH_results/RBH.Hsap.Pvit.tsv", header = T)
Pvit2Hsap <- read.table("/mnt/data01/yuanzhen/01.Vertebrate_cell_evo/02.gene_relationships/run4/results/Ortho_pipeline/OrthoFinder/RBH_results/RBH.Pvit.Hsap.tsv", header = T)
Mmus2Pvit <- read.table("/mnt/data01/yuanzhen/01.Vertebrate_cell_evo/02.gene_relationships/run4/results/Ortho_pipeline/OrthoFinder/RBH_results/RBH.Mmus.Pvit.tsv", header = T)
Pvit2Mmus <- read.table("/mnt/data01/yuanzhen/01.Vertebrate_cell_evo/02.gene_relationships/run4/results/Ortho_pipeline/OrthoFinder/RBH_results/RBH.Pvit.Mmus.tsv", header = T)

Hsap2Pmar <- read.table("/mnt/data01/yuanzhen/01.Vertebrate_cell_evo/02.gene_relationships/run4/results/Ortho_pipeline/OrthoFinder/RBH_results/RBH.Hsap.Pmar.tsv", header = T)
Pmar2Hsap <- read.table("/mnt/data01/yuanzhen/01.Vertebrate_cell_evo/02.gene_relationships/run4/results/Ortho_pipeline/OrthoFinder/RBH_results/RBH.Pmar.Hsap.tsv", header = T)
Mmus2Pmar <- read.table("/mnt/data01/yuanzhen/01.Vertebrate_cell_evo/02.gene_relationships/run4/results/Ortho_pipeline/OrthoFinder/RBH_results/RBH.Mmus.Pmar.tsv", header = T)
Pmar2Mmus <- read.table("/mnt/data01/yuanzhen/01.Vertebrate_cell_evo/02.gene_relationships/run4/results/Ortho_pipeline/OrthoFinder/RBH_results/RBH.Pmar.Mmus.tsv", header = T)
Pvit2Pmar <- read.table("/mnt/data01/yuanzhen/01.Vertebrate_cell_evo/02.gene_relationships/run4/results/Ortho_pipeline/OrthoFinder/RBH_results/RBH.Pvit.Pmar.tsv", header = T)
Pmar2Pvit <- read.table("/mnt/data01/yuanzhen/01.Vertebrate_cell_evo/02.gene_relationships/run4/results/Ortho_pipeline/OrthoFinder/RBH_results/RBH.Pmar.Pvit.tsv", header = T)



# pay attention, hsap2mmus just species linkage orthologs, not dup1 not dup2 
colnames(Hsap2Mmus) <- c("dup1","dup2")
colnames(Mmus2Hsap) <- c("dup1","dup2")
colnames(Hsap2Pvit) <- c("dup1","dup2")
colnames(Pvit2Hsap) <- c("dup1","dup2")
colnames(Mmus2Pvit) <- c("dup1","dup2")
colnames(Pvit2Mmus) <- c("dup1","dup2")
colnames(Hsap2Pmar) <- c("dup1","dup2")
colnames(Pmar2Hsap) <- c("dup1","dup2")
colnames(Mmus2Pmar) <- c("dup1","dup2")
colnames(Pmar2Mmus) <- c("dup1","dup2")
colnames(Pvit2Pmar) <- c("dup1","dup2")
colnames(Pmar2Pvit) <- c("dup1","dup2")


orthologs <- rbind(Hsap2Mmus,Mmus2Hsap,Hsap2Pvit,Pvit2Hsap,Mmus2Pvit,Pvit2Mmus,
                                      Hsap2Pmar,Pmar2Hsap,Mmus2Pmar,Pmar2Mmus,Pvit2Pmar,Pmar2Pvit)

# combined paralogs and orthologs file together, and build graph
Combined_pairs_ssd <- duplicate_pairs %>% filter(species %in% species & type == "SSD") %>% select(c("dup1", "dup2"))
Combined_pairs_wgd <- duplicate_pairs %>% filter(species %in% species & type == "WGD") %>% select(c("dup1", "dup2"))

Combined_pairs_ssd <- rbind(orthologs, Combined_pairs_ssd)
Combined_pairs_wgd <- rbind(orthologs, Combined_pairs_wgd)

# load gene expressed in our sc dataset
background <- readRDS("/mnt/data01/yuanzhen/01.Vertebrate_cell_evo/04.ohno_para_significance/0bin/gene_bg.rds")

# retain genes only in my dataset
Combined_pairs_ssd <- Combined_pairs_ssd[Combined_pairs_ssd$dup1 %in% background$gene & Combined_pairs_ssd$dup2 %in% background$gene, ]
Combined_pairs_wgd <- Combined_pairs_wgd[Combined_pairs_wgd$dup1 %in% background$gene & Combined_pairs_wgd$dup2 %in% background$gene, ]

# get SSD/WGD families
g <- graph_from_data_frame(Combined_pairs_ssd, directed = FALSE)
components <- components(g)$membership
SSD_pairs <- split(names(components), components)
g <- graph_from_data_frame(Combined_pairs_wgd, directed = FALSE)
components <- components(g)$membership
WGD_pairs <- split(names(components), components)

# load wilcox test results and plot barplot to show celltype non-speicfic dominate expression across species
species <- c("Hsap", "Mmus", "Pvit", "Pmar")
type <- c("WGD", "SSD")
mtx <- c("exp", "pct")

pairwise_res <- Reduce(rbind, lapply(species, FUN = function(s){
                     Reduce(rbind, lapply(type, FUN = function(t){
                        Reduce(rbind, lapply(mtx, FUN = function(m){
                            info <- read.delim(paste0(s,"/",m,"/",s, ".", m, ".", t, ".pairwise_wilcox_test.txt"), header = T)
                            info$species <- s
                            info$type <- t
                            info$mtx <- m
                            return(info)
                        }))
    }))
}))

exp_pct_res <- Reduce(rbind, lapply(species, FUN = function(s){
                    Reduce(rbind, lapply(type, FUN = function(t){
                        Reduce(rbind, lapply(mtx, FUN = function(m){
                            info <- read.delim(paste0(s,"/",m,"/",s, ".", m, ".", t, ".mean.txt"), header = T)
                            info$species <- s
                            info$type <- t
                            info$mtx <- m
                            return(info)
                        }))
                    }))
}))

# for each 
get_res <- function(diff_, mean_, pairs, label){
    # the following condition is based on pairwise wilcox test comparison
    if(length(pairs) == 2){
        n_sign <- diff_[diff_$group1 %in% pairs, "p"]
    } else if (length(pairs) > 2) {
        n_sign <- diff_[diff_$group1 %in% pairs, "p.adj.signif"]
    } else {
        n_sign <- NA
    }
    n_sign <- n_sign[!is.na(n_sign) & n_sign < 0.01]

    # load expression data and use mean expression to represent expression among cell types
    exp_mean_tmp <- setNames(mean_[mean_$Var1 %in% pairs, "mean"],  mean_[mean_$Var1 %in% pairs, "Var1"])
        
    if (length(n_sign >= 1)){
        info <- c(length(pairs), length(n_sign), names(exp_mean_tmp)[which.max(exp_mean_tmp)])
    } else {
        if (length(pairs) == 1){
            info <- c(length(pairs), length(n_sign), NA)
        } else {
            info <- c(length(pairs), length(n_sign), NA)
        }   
    }
    info <- c(info, label)
    return(info)
}

dominant_copy_species <- function(family, type_, mtx_){
    info <- Reduce(rbind, lapply(seq_along(family), FUN = function(x){
                                     pairs <- family[[x]]
                                     pairs_Hsap <- pairs[pairs %in% background[background$species == "Hsap", "gene"]]
                                     pairs_Mmus <- pairs[pairs %in% background[background$species == "Mmus", "gene"]]
                                     pairs_Pvit <- pairs[pairs %in% background[background$species == "Pvit", "gene"]]
                                     pairs_Pmar <- pairs[pairs %in% background[background$species == "Pmar", "gene"]]
                                     # pay attention, type == type wouldn't filter any rows, you can use type == !!type  to disambiguate second type as a variable
                                     diff_Hsap <- pairwise_res %>% filter(species == "Hsap" & type == type_ & mtx == mtx_)
                                     diff_Mmus <- pairwise_res %>% filter(species == "Mmus" & type == type_ & mtx == mtx_)
                                     diff_Pvit <- pairwise_res %>% filter(species == "Pvit" & type == type_ & mtx == mtx_)
                                     diff_Pmar <- pairwise_res %>% filter(species == "Pmar" & type == type_ & mtx == mtx_)
                                     
                                     mean_Hsap <- exp_pct_res %>% filter(species == "Hsap" & type == type_ & mtx == mtx_)
                                     mean_Mmus <- exp_pct_res %>% filter(species == "Mmus" & type == type_ & mtx == mtx_)
                                     mean_Pvit <- exp_pct_res %>% filter(species == "Pvit" & type == type_ & mtx == mtx_)
                                     mean_Pmar <- exp_pct_res %>% filter(species == "Pmar" & type == type_ & mtx == mtx_)
                                                                         
                                     res <- rbind(get_res(diff_Hsap, mean_Hsap, pairs_Hsap, "Hsap"),
                                                  get_res(diff_Mmus, mean_Mmus, pairs_Mmus, "Mmus"),
                                                  get_res(diff_Pvit, mean_Pvit, pairs_Pvit, "Pvit"),
                                                  get_res(diff_Pmar, mean_Pmar, pairs_Pmar, "Pmar")
                                                )
                                     res <- data.frame(res)
                                     res$type = type_
                                     res$mtx = mtx_
                                     res$family_index <- x
                                     return(res)
                                    }))
        info <- data.frame(info)
        colnames(info) <- c("number_of_paralogs", "number_of_sign", "Max", "species", "type", "mtx", "family_index")
        return(info)
}

WGD_exp <- dominant_copy_species(family = WGD_pairs, type_ = "WGD", mtx_ = "exp")
WGD_pct <- dominant_copy_species(family = WGD_pairs, type_ = "WGD", mtx_ = "pct")
SSD_exp <- dominant_copy_species(family = SSD_pairs, type_ = "SSD", mtx_ = "exp")
SSD_pct <- dominant_copy_species(family = SSD_pairs, type_ = "SSD", mtx_ = "pct")

# get dominant copy and several other informations
get_dominant_copy_across_species <- function(results, species1, species2, orthologs){
    # filter to certain condition
    results <- results %>% filter(species %in% c(species1, species2)) %>% group_by(family_index) %>% 
        filter(all(number_of_paralogs == 2) & all(number_of_sign > 0) ) %>%
        # make sure there are at least 2 copies for both species in that family
        #filter(any(number_of_paralogs >= 2) & all(number_of_paralogs >= 1 & number_of_paralogs <= 4))%>% 
        ungroup()
    res <- Reduce(rbind, lapply(unique(results$family_index), FUN = function(i){
        tmp = results %>% filter(family_index == i)
        tmp$orthologs <- ifelse(tmp$Max %in% orthologs$dup1, orthologs[match(tmp$Max, orthologs$dup1), 2], 
                                ifelse(tmp$Max %in% orthologs$dup2, orthologs[match(tmp$Max, orthologs$dup2), 1], NA))
        tmp$orthologs_is_dom <- vapply(tmp$orthologs, FUN = function(x){
                                           if (is.na(x)){ return(FALSE) } else { x %in% tmp$Max }
                                }, FUN.VALUE = logical(1))
        return(tmp)
                                    }))
    res$number_of_paralogs <- as.numeric(res$number_of_paralogs)
    res$number_of_sign <- as.numeric(res$number_of_sign)
    return(res)
}



dominant_usage_wgd_exp_Hsap2Mmus <- get_dominant_copy_across_species(WGD_exp, "Hsap", "Mmus", orthologs = Hsap2Mmus)
dominant_usage_wgd_exp_Hsap2Pvit <- get_dominant_copy_across_species(WGD_exp, "Hsap", "Pvit", orthologs = Hsap2Pvit)
dominant_usage_wgd_exp_Mmus2Pvit <- get_dominant_copy_across_species(WGD_exp, "Mmus", "Pvit", orthologs = Mmus2Pvit)
dominant_usage_wgd_exp_Hsap2Pmar <- get_dominant_copy_across_species(WGD_exp, "Hsap", "Pmar", orthologs = Hsap2Pmar)
dominant_usage_wgd_exp_Mmus2Pmar <- get_dominant_copy_across_species(WGD_exp, "Mmus", "Pmar", orthologs = Mmus2Pmar)
dominant_usage_wgd_exp_Pvit2Pmar <- get_dominant_copy_across_species(WGD_exp, "Pvit", "Pmar", orthologs = Pvit2Pmar)

dominant_usage_ssd_exp_Hsap2Mmus <- get_dominant_copy_across_species(SSD_exp, "Hsap", "Mmus", orthologs = Hsap2Mmus)
dominant_usage_ssd_exp_Hsap2Pvit <- get_dominant_copy_across_species(SSD_exp, "Hsap", "Pvit", orthologs = Hsap2Pvit)
dominant_usage_ssd_exp_Mmus2Pvit <- get_dominant_copy_across_species(SSD_exp, "Mmus", "Pvit", orthologs = Mmus2Pvit)
dominant_usage_ssd_exp_Hsap2Pmar <- get_dominant_copy_across_species(SSD_exp, "Hsap", "Pmar", orthologs = Hsap2Pmar)
dominant_usage_ssd_exp_Mmus2Pmar <- get_dominant_copy_across_species(SSD_exp, "Mmus", "Pmar", orthologs = Mmus2Pmar)
dominant_usage_ssd_exp_Pvit2Pmar <- get_dominant_copy_across_species(SSD_exp, "Pvit", "Pmar", orthologs = Pvit2Pmar)

dominant_usage_wgd_pct_Hsap2Mmus <- get_dominant_copy_across_species(WGD_pct, "Hsap", "Mmus", orthologs = Hsap2Mmus)
dominant_usage_wgd_pct_Hsap2Pvit <- get_dominant_copy_across_species(WGD_pct, "Hsap", "Pvit", orthologs = Hsap2Pvit)
dominant_usage_wgd_pct_Mmus2Pvit <- get_dominant_copy_across_species(WGD_pct, "Mmus", "Pvit", orthologs = Mmus2Pvit)
dominant_usage_wgd_pct_Hsap2Pmar <- get_dominant_copy_across_species(WGD_pct, "Hsap", "Pmar", orthologs = Hsap2Pmar)
dominant_usage_wgd_pct_Mmus2Pmar <- get_dominant_copy_across_species(WGD_pct, "Mmus", "Pmar", orthologs = Mmus2Pmar)
dominant_usage_wgd_pct_Pvit2Pmar <- get_dominant_copy_across_species(WGD_pct, "Pvit", "Pmar", orthologs = Pvit2Pmar)

dominant_usage_ssd_pct_Hsap2Mmus <- get_dominant_copy_across_species(SSD_pct, "Hsap", "Mmus", orthologs = Hsap2Mmus)
dominant_usage_ssd_pct_Hsap2Pvit <- get_dominant_copy_across_species(SSD_pct, "Hsap", "Pvit", orthologs = Hsap2Pvit)
dominant_usage_ssd_pct_Mmus2Pvit <- get_dominant_copy_across_species(SSD_pct, "Mmus", "Pvit", orthologs = Mmus2Pvit)
dominant_usage_ssd_pct_Hsap2Pmar <- get_dominant_copy_across_species(SSD_pct, "Hsap", "Pmar", orthologs = Hsap2Pmar)
dominant_usage_ssd_pct_Mmus2Pmar <- get_dominant_copy_across_species(SSD_pct, "Mmus", "Pmar", orthologs = Mmus2Pmar)
dominant_usage_ssd_pct_Pvit2Pmar <- get_dominant_copy_across_species(SSD_pct, "Pvit", "Pmar", orthologs = Pvit2Pmar)

# summary of results
summary_exp <- rbind(dominant_usage_wgd_exp_Hsap2Mmus %>% mutate(species_pair = "Hsap_vs_Mmus"), 
                     dominant_usage_wgd_exp_Hsap2Pvit %>% mutate(species_pair = "Hsap_vs_Pvit"),
                     dominant_usage_wgd_exp_Mmus2Pvit %>% mutate(species_pair = "Mmus_vs_Pvit"),
                     dominant_usage_wgd_exp_Hsap2Pmar %>% mutate(species_pair = "Hsap_vs_Pmar"), 
                     dominant_usage_wgd_exp_Mmus2Pmar %>% mutate(species_pair = "Mmus_vs_Pmar"),
                     dominant_usage_wgd_exp_Pvit2Pmar %>% mutate(species_pair = "Pvit_vs_Pmar"),
                     dominant_usage_ssd_exp_Hsap2Mmus %>% mutate(species_pair = "Hsap_vs_Mmus"),
                     dominant_usage_ssd_exp_Hsap2Pvit %>% mutate(species_pair = "Hsap_vs_Pvit"),
                     dominant_usage_ssd_exp_Mmus2Pvit %>% mutate(species_pair = "Mmus_vs_Pvit"),
                     dominant_usage_ssd_exp_Hsap2Pmar %>% mutate(species_pair = "Hsap_vs_Pmar"), 
                     dominant_usage_ssd_exp_Mmus2Pmar %>% mutate(species_pair = "Mmus_vs_Pmar"),
                     dominant_usage_ssd_exp_Pvit2Pmar %>% mutate(species_pair = "Pvit_vs_Pmar"))

summary_pct <- rbind(dominant_usage_wgd_pct_Hsap2Mmus %>% mutate(species_pair = "Hsap_vs_Mmus"), 
                     dominant_usage_wgd_pct_Hsap2Pvit %>% mutate(species_pair = "Hsap_vs_Pvit"),
                     dominant_usage_wgd_pct_Mmus2Pvit %>% mutate(species_pair = "Mmus_vs_Pvit"),
                     dominant_usage_wgd_pct_Hsap2Pmar %>% mutate(species_pair = "Hsap_vs_Pmar"), 
                     dominant_usage_wgd_pct_Mmus2Pmar %>% mutate(species_pair = "Mmus_vs_Pmar"),
                     dominant_usage_wgd_pct_Pvit2Pmar %>% mutate(species_pair = "Pvit_vs_Pmar"),
                     dominant_usage_ssd_pct_Hsap2Mmus %>% mutate(species_pair = "Hsap_vs_Mmus"),
                     dominant_usage_ssd_pct_Hsap2Pvit %>% mutate(species_pair = "Hsap_vs_Pvit"),
                     dominant_usage_ssd_pct_Mmus2Pvit %>% mutate(species_pair = "Mmus_vs_Pvit"),
                     dominant_usage_ssd_pct_Hsap2Pmar %>% mutate(species_pair = "Hsap_vs_Pmar"), 
                     dominant_usage_ssd_pct_Mmus2Pmar %>% mutate(species_pair = "Mmus_vs_Pmar"),
                     dominant_usage_ssd_pct_Pvit2Pmar %>% mutate(species_pair = "Pvit_vs_Pmar"))

# investigate conservation and divergence at family level
exp_details <- summary_exp %>% group_by(type, species_pair, family_index) %>% 
        summarise(
                  conserved = ifelse(all(number_of_sign > 0) & sum(orthologs_is_dom) > 0,"YES", "NO"),
                  divergent = ifelse(all(number_of_sign > 0) & sum(!is.na(Max) & !is.na(orthologs) & orthologs_is_dom == FALSE) > 0, "YES", "NO")
                  ) %>% ungroup()
pct_details <- summary_pct %>% group_by(type, species_pair, family_index) %>% 
        summarise(
                  conserved = ifelse(all(number_of_sign > 0) & sum(orthologs_is_dom) > 0,"YES", "NO"),
                  divergent = ifelse(all(number_of_sign > 0) & sum(!is.na(Max) & !is.na(orthologs) & orthologs_is_dom == FALSE) > 0, "YES", "NO")
                  ) %>% ungroup()

# calculate the number of conserved pattern and divergent pattern across species
dominance_pattern_exp <- exp_details %>% group_by(species_pair, type) %>% summarise(
                                                                                    total_number = length(family_index),
                                                                                    conserved = sum(conserved == "YES", na.rm = TRUE),
                                                                                    divergent = sum(divergent == "YES", na.rm = TRUE)
                                                                                    ) %>% ungroup()

dominance_pattern_pct <- pct_details %>% group_by(species_pair, type) %>% summarise(
                                                                                    total_number = length(family_index),
                                                                                    conserved = sum(conserved == "YES", na.rm = TRUE),
                                                                                    divergent = sum(divergent == "YES", na.rm = TRUE)
                                                                                    ) %>% ungroup()

# prepare data for plotting
dominance_pattern_exp <- dominance_pattern_exp %>% pivot_longer(cols = c(conserved, divergent), names_to = "variable",  values_to = "value")
dominance_pattern_exp$type <- factor(dominance_pattern_exp$type, levels = c("WGD","SSD"))
dominance_pattern_exp$species_pair <- factor(dominance_pattern_exp$species_pair, 
                                             levels = c("Hsap_vs_Mmus","Hsap_vs_Pvit","Mmus_vs_Pvit","Hsap_vs_Pmar","Mmus_vs_Pmar","Pvit_vs_Pmar"))

dominance_pattern_pct <- dominance_pattern_pct %>% pivot_longer(cols = c(conserved, divergent), names_to = "variable", values_to = "value")
dominance_pattern_pct$type <- factor(dominance_pattern_pct$type, levels = c("WGD","SSD"))
dominance_pattern_pct$species_pair <- factor(dominance_pattern_pct$species_pair, 
                                             levels = c("Hsap_vs_Mmus","Hsap_vs_Pvit","Mmus_vs_Pvit","Hsap_vs_Pmar","Mmus_vs_Pmar","Pvit_vs_Pmar"))

pdf("2.Fig3_plots/dominant_expression_pattern_across_species.exp.2copies.pdf", width = 6, height = 4)
dominance_pattern_exp %>% ggbarplot(x = "variable", y = "value", color = "variable", fill = "variable") + 
    scale_fill_manual(values=c("#99990066","#66666666"))+
    scale_color_manual(values=c("#999900","#666666")) + 
    facet_grid(vars(type), vars(species_pair))
dev.off()

pdf("2.Fig3_plots/dominant_expression_pattern_across_species.pct.2copies.pdf", width = 6, height = 4)
dominance_pattern_pct %>% ggbarplot(x = "variable", y = "value", color = "variable", fill = "variable") + 
    scale_fill_manual(values=c("#99990066","#66666666"))+
    scale_color_manual(values=c("#999900","#666666")) +
    facet_grid(vars(type), vars(species_pair))
dev.off()

save.image(file = "plot_dominant_expression_across_species.RData")
