
#############################################################################
##### This function takes a list of markers and a Seurat v2.3.4 object ######
##### and return a dotplot of the markers expression across clusters   ######
#############################################################################

Cluster.dotplot <- function(Dataset,
                            Marker.genes,
                            min.expression,
                            percent.min,
                            maxdot.size){
  #Extract and transpose the matrix of marker genes expression
  data.to.plot <- data.frame(t(as.matrix(Dataset[["RNA"]]@data[Marker.genes, ])))
  
  data.to.plot$Cell <- rownames(data.to.plot)
  data.to.plot$id <- Idents(Dataset)
  
  #Reshape the dataframe
  data.to.plot <- data.to.plot %>% tidyr::gather(key = Marker.genes, value = expression, -c(Cell, id)) 
  
  #For each genes in each cluster calculate mean expression and percent cell with norm expression > 0
  data.to.plot <- data.to.plot %>%
		  group_by(id, Marker.genes) %>% 
    		  summarize(avg.exp = mean(expm1(x = expression)), pct.exp = length(x = expression[expression > min.expression]) / length(x = expression) ) 
  
  data.to.plot <- data.to.plot %>% ungroup() %>%
    group_by(Marker.genes) %>% 
    mutate(avg.exp.scale = scale(x = avg.exp)) %>%
    mutate(avg.exp.scale = MinMax(data = avg.exp.scale, max = 2, min = -2)) # add column with scaled expression values from -2 to 2
  
  data.to.plot$genes.plot <- factor(x = data.to.plot$Marker.genes, levels = rev(x = Marker.genes)) #Put gene names as factor 
  
  data.to.plot$pct.exp[data.to.plot$pct.exp < percent.min] <- NA #Set to Na if less than percent.min of cells expresse the gene
  data.to.plot$pct.exp <- data.to.plot$pct.exp * 100
  
  p <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot, y = id)) +
	    geom_point(mapping = aes(size = pct.exp, fill = avg.exp.scale), shape = 21, color = "black", stroke = 0.5) + # modify the colors by if want to color by domains or cluster ident
	    scale_size_area(max_size= maxdot.size) + # Scale the radius of the dot from 0 to 6 
	    scale_y_discrete(position = "right") +
	    theme(axis.text.x = element_text(angle = 90, vjust = 1), axis.title.y = element_blank()) +
	    xlab("") + ylab("") +
	    scale_fill_gradientn(colours = brewer.pal(9,"Reds")) + coord_flip()
} 





#############################################################################
##### This function takes a list of markers and a Seurat v2.3.4 object ######
##### and return a dotplot of the markers expression across clusters   ######
#############################################################################

Cluster.dotplot.ann <- function(Dataset,
                            Marker.genes,
                            min.expression,
                            percent.min,
                            maxdot.size, annotation){
  #Extract and transpose the matrix of marker genes expression
  data.to.plot <- data.frame(t(as.matrix(Dataset[["RNA"]]@data[Marker.genes, ])))
  
  data.to.plot$Cell <- rownames(data.to.plot)
  data.to.plot$id <- Idents(Dataset)
  data.to.plot$ann <- Dataset@meta.data[, annotation]
  
  #Reshape the dataframe
  data.to.plot <- data.to.plot %>% tidyr::gather(key = Marker.genes, value = expression, -c(Cell, id, ann)) 
  
  #For each genes in each cluster calculate mean expression and percent cell with norm expression > 0
  data.to.plot <- data.to.plot %>%
		  group_by(ann, id, Marker.genes) %>% 
    		  summarize(avg.exp = mean(expm1(x = expression)), pct.exp = length(x = expression[expression > min.expression]) / length(x = expression) ) 
  
  data.to.plot <- data.to.plot %>% ungroup() %>%
    group_by(Marker.genes) %>% 
    mutate(avg.exp.scale = scale(x = avg.exp)) %>%
    mutate(avg.exp.scale = MinMax(data = avg.exp.scale, max = 2, min = -2)) # add column with scaled expression values from -2 to 2
  
  data.to.plot$genes.plot <- factor(x = data.to.plot$Marker.genes, levels = rev(x = Marker.genes)) #Put gene names as factor 
  
  data.to.plot$pct.exp[data.to.plot$pct.exp < percent.min] <- NA #Set to Na if less than percent.min of cells expresse the gene
  data.to.plot$pct.exp <- data.to.plot$pct.exp * 100
  
  p <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot, y = id)) +
	    geom_point(mapping = aes(size = pct.exp, fill = ann), shape = 21, color = "black", stroke = 0.5) + # modify the colors by if want to color by domains or cluster ident
	    scale_size_area(max_size= maxdot.size) + # Scale the radius of the dot from 0 to 6 
	    scale_y_discrete(position = "right") +
	    theme(axis.text.x = element_text(angle = 90, vjust = 1), axis.title.y = element_blank()) +
	    xlab("") + ylab("") +
	    scale_fill_manual(values= c(brewer.pal(12,'Set3'), brewer.pal(6, 'Set1'), brewer.pal(12,'Paired'), brewer.pal(8, 'Dark2') )) + coord_flip()
} 







###################################################################################
############ This function returns a list of all nodes marker as defined ##########
############ using the roc test implemented in Seurat FindMarkersNode function ####
###################################################################################

Find.nodes.markers <- function(Datasrt,
                               ngenes) {
  Branch.markers <- list()
  Node <- list()
  tree <- Datasrt@cluster.tree[[1]]
  
  for (i in  unique(tree$edge[,1])) {
    
    Node[[i]] <- FindMarkersNode(object = Datasrt,
                                 test.use = "roc",
                                 node = i,
                                 min.pct = 0,
                                 logfc.threshold = 0.6,
                                 print.bar = F,
                                 only.pos = F)
    
    Node[[i]]$gene <- rownames(Node[[i]])
    Node[[i]]$q.diff <- abs(Node[[i]]$pct.1 - Node[[i]]$pct.2)/max(Node[[i]]$pct.1, Node[[i]]$pct.2)
    Node[[i]]$Specificity.index <- Node[[i]]$avg_logFC/ (1-Node[[i]]$q.diff)
    
    pos.genes  <- Node[[i]] %>%
      filter(Specificity.index > 1 & pct.1 > 0.65) %>% #New
      top_n(ngenes, Specificity.index) %>%
      arrange(Specificity.index) %>%
      pull(gene)
    
    neg.genes  <- Node[[i]] %>%
      filter(Specificity.index < -1 & pct.2 > 0.65) %>% #New
      top_n(-ngenes, Specificity.index) %>%
      arrange(Specificity.index) %>%
      pull(gene)
    
    Branch.markers[[i]] <- c(pos.genes,neg.genes)
  }
  
  list.index.start <- min(unique(Datasrt@cluster.tree[[1]][["edge"]][,1])) -1
  list.index.end <- max(unique(Datasrt@cluster.tree[[1]][["edge"]][,1])) -list.index.start
  
  Branch.markers <- Branch.markers[c(rep(F,list.index.start), rep(T,list.index.end))] 
  
  return(Branch.markers)
}


##############################################################################
###### Return a list of ggplot object for each node markers ##################
##############################################################################

Node.markers.plots <- function(Datasrt,
                              Branch.markers.list) {
  p <- list()
  for (i in seq(Branch.markers.list)) {
    p[[i+1]] <- Cluster.dotplot(Datasrt,
                                Marker.genes = rev(Branch.markers.list[[i]]),
                                min.expression = 0.7 ,percent.mi= 0.2, maxdot.size = 5)  
  } 
  return(p)
}

##########################################################################
############### Plot nodes dot plots with all clusters ###################
##########################################################################

Plot.node.markers <- function(Datasrt,
                              node,
                              plot){
      tree <- Datasrt@cluster.tree[[1]]
      plotindx <- node - tree$Nnode
      return(plot[plotindx])
}

################################################################################
################ Plot nodes dotplots with only a subset of clusters ############
################################################################################

Plot.node.subset.markers <- function(Datasrt,
                                     node,
                                     Branch.markers.list){
  
  tree <- Datasrt@cluster.tree[[1]]
  
  clusters <- c()
  branches <- tree$edge[which(tree$edge[, 1] == node), 2]
  next.branche <- T
  
  while (next.branche == T) {
    next.branche  <- ifelse(!sum(branches > tree$Nnode+1) == 0, T, F)
    Next.nodes <- list()
    
    for (i in seq(length(branches))) {
      if (branches[i] <= tree$Nnode+1) {
        clusters <- c(clusters,tree$tip.label[branches[i]])
      }
      else{
        Next.nodes[[i]] <- tree$edge[which(tree$edge[, 1] == branches[i]), 2]
      } 
    } 
    branches <- unlist(Next.nodes)
    
  }
  
  tmp.Datasrt <- SubsetData(Datasrt, ident.use = clusters, subset.raw = T,  do.clean = F)
  
  
  plot.list <- Node.markers.plots(tmp.Datasrt,
                                 Branch.markers.list)
  
  Plot.node.markers(tmp.Datasrt,
                    node = node,
                    plot =  plot.list)
  
  
}

##### This function is for getting ratios at major cell type level and subtype level
calculate_ratios <- function(obj, order, cluster_name) {
    res <- data.frame(c = character(), cell = character(), ratio = double(), stringsAsFactors = FALSE)
    for (cluster in order) {
        tmp <- colnames(obj)[which(obj@meta.data$iter.clustering == cluster)]
        for (cell_type in unique(obj@meta.data[, cluster_name])) {
            tmp2 <- colnames(obj)[which(obj@meta.data[, cluster_name] == cell_type)]
            # Calculate the ratio
            ratio <- length(intersect(tmp, tmp2)) / length(tmp)
            res <- rbind(res, data.frame(c = cluster, num = length(tmp), cell = cell_type, ratio = ratio))
        }
    }
    # Set column types
    res$c <- factor(res$c, levels = order)
    res$ratio <- as.double(res$ratio)
    return(res)
}

