# Load necessary libraries
library(Seurat)
library(matrixStats)  # For efficient row/column operations

# Define the trinarize function
trinarize <- function(counts, cluster_ids, f = 0.2) {
    # Extract count matrix from Seurat object (assumes RNA assay)
    #counts <- GetAssayData(seurat_obj, layer = "counts")
    #counts <- as.matrix(counts)
    # Get the cluster assignments and gene expressions
    #cluster_ids <- seurat_obj@meta.data[, label]
      
    # Calculate the number of non-zero entries (expressed genes) per cluster per gene
    Nonzeros <- t(aggregate(t(counts > 0), by = list(cluster_ids), FUN = sum, drop = TRUE))   
    # Calculate the number of cells per cluster
    NCells <- table(cluster_ids)    
    # Function to calculate the trinarization probability for one gene in one cluster
    p_half <- function(k, n, f) {
        a <- 1.5
        b <- 2
        incb <- pbeta(f, a + k, b - k + n)      
        if (incb == 0) {
            return(1.0)
        } else {
            log_incb <- log(incb)
            betaln_val <- lbeta(a + k, b - k + n)
            gamma_diff <- lgamma(a + b + n) - (lgamma(a + k) + lgamma(b - k + n))
            p <- 1.0 - exp(log_incb + betaln_val + gamma_diff)
            return(p)
        }
    }
    # Vectorize the function to apply to all clusters and genes
    p_half_vectorized <- Vectorize(p_half)
    # Compute trinarization probabilities for all clusters and genes
    result <- t(as.data.frame(sapply(2:nrow(Nonzeros), function(i) {
                                         p_half_vectorized(as.integer(Nonzeros[i, ]), NCells, f)
    })))
    rownames(result) <- rownames(counts)
    colnames(result) <- names(NCells)
    return(result)
}
