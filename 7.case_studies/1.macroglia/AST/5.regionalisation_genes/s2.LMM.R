suppressMessages(suppressWarnings({
    library(lme4)
    library(ggplot2)
    library(dplyr)
}))

parser = argparse::ArgumentParser(description="Test LMM on different datasets")
parser$add_argument('-I','--input', help='input normalised expression matrix')
parser$add_argument('-L','--label', help='input label')
args = parser$parse_args()

label <- args$label
mat <- read.delim(args$input, header = T, sep = '\t')
mat$Gene <- rownames(mat)
exp <- reshape2::melt(mat)

exp$location <- vapply(as.character(exp$variable), FUN = function(x) {
                           strsplit(x, split = '\\.')[[1]][1]
}, FUN.VALUE = character(1))
exp$celltype <- vapply(as.character(exp$variable), FUN = function(x) {
                           paste0(strsplit(x, split = '\\.')[[1]][-1], collapse = '_')
}, FUN.VALUE = character(1))

exp_AST_GABA <- exp[exp$celltype %in% c('GABAergic_neurons', 'Astrocytes'), ]
tmp <- unique(exp_AST_GABA$location)[rowSums(table(unique(exp_AST_GABA[,c('location', 'celltype')]))) !=2] # 2 cell types
exp_AST_GABA <- exp_AST_GABA[!exp_AST_GABA$location %in% tmp, ]

exp_AST_Glut <- exp[exp$celltype %in% c('glutamatergic_neurons', 'Astrocytes'), ]
tmp <- unique(exp_AST_Glut$location)[rowSums(table(unique(exp_AST_Glut[,c('location', 'celltype')]))) !=2] # 2 cell types
exp_AST_Glut <- exp_AST_Glut[!exp_AST_Glut$location %in% tmp, ]

# variance decomp
fit_model_for_gene <- function(data_subset) {
    # Fit the mixed-effects model
    model <- lmer(value ~ (1 | celltype) + (1 | location), data = data_subset)
    variance_components <- as.data.frame(VarCorr(model))
    total_variance <- sum(variance_components$vcov)
    variance_proportion <- variance_components$vcov / total_variance
    result <- c(unique(data_subset$Gene),variance_proportion[variance_components$grp == "celltype"],
                variance_proportion[variance_components$grp == "location"],
                variance_proportion[variance_components$grp == "Residual"]
                )
    return(result)
}

get_res <- function(exp, type){
    results <- do.call(rbind, lapply(split(exp, exp$Gene), fit_model_for_gene))
    results <- as.data.frame(results)
    colnames(results) <- c('Gene', 'celltype', 'location', 'Residual')

    results$celltype <- as.double(results$celltype)
    results$location <- as.double(results$location)
    results$Residual <- as.double(results$Residual)
    results$Category <- with(results, ifelse(
                                         location > 0.5 , "High across location",
                                         ifelse(celltype > 0.5, "High across cell type family", "None")
                                         ))
    results <- results[complete.cases(results), ]
    
    p <- ggplot(results, aes(x = celltype, y = location, color = Category)) +
        geom_point(alpha = 0.6, size = 2) + 
        scale_color_manual(values = c("High across location" = "#A4CB9E",
                                    "High across cell type family" = "#E17327",
                                    "None" = "#A3A5A6"
                                    )) +  
        labs(x = "Fraction of variance across cell type family", y = "Fraction of variance across location", color = "Category") +
        theme_minimal() + theme( legend.position = "top", text = element_text(size = 12)) + coord_fixed(ratio = 1)
    ggsave(paste0(label,'.',type,".LMM_variance_decomp.cell_vs_location.pdf"), p, width = 6, height = 6)
    write.table(results, paste0(label,'.', type ,'.LMM_variance_decomp.cell_vs_location.txt'), col.names = T, row.names = F, quote = F, sep = '\t')
}

get_res(exp_AST_GABA, 'AST_GABA')
get_res(exp_AST_Glut, 'AST_Glut')


