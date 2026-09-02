rm(list=ls())
library(igraph)
library(sbm)
require(ggplot2)
library(viridisLite)
library(stringr)
exportFig <- FALSE

where_fig = getwd()

# Unified theme for figures
theme_set(theme_bw(base_family = "LM Roman 10"))
graph_unit <- "px"
graph_width <- 1500
graph_height <- 1500


data("fungusTreeNetwork")

# Define symmetric adjacency matrix
tree_tree_weighted <- fungusTreeNetwork$tree_tree
row.names(tree_tree_weighted) <- colnames(tree_tree_weighted) <- str_sub(fungusTreeNetwork$tree_names, 1, 30)

##################################################################################
## ----tree_tree network plot data------------------------------------------------
##################################################################################
# as matrix 
g <- plotMyMatrix(tree_tree_weighted , dimLabels = list(row = 'Trees', col = 'Trees'), plotOptions = list(legend = TRUE,rowNames = TRUE,colNames = TRUE))
g <- g + theme(axis.text=element_text(size=8)) 
g
if(exportFig)
  ggsave(paste0(where_fig,"/FungusTree_Poisson_data.png"),
         plot = g ,
         scale = 1,
         width = 20,
         units = c("cm"))

#as network 

if(exportFig){png(paste0(where_fig,"/FungusTree_Poisson_data_network.png"), 
                  width = graph_width, 
                  height = graph_width, units=graph_unit, res = 150)}
nnodes <- nrow(tree_tree_weighted)
# Create graph
g <- graph_from_adjacency_matrix(tree_tree_weighted, mode = "undirected", weighted = TRUE, diag = FALSE)

# Clean labels (remove text after '(')
labels <- sub("\\s*\\(.*", "", V(g)$name)

# Get circular layout
layout_circle <- layout_in_circle(g) * 1.01 + matrix(rnorm(nnodes,0,1/100),nnodes,2,byrow = TRUE)
label_pos <- layout_circle 
# Combine layout and label positions to compute plot limits
all_x <- c(layout_circle[,1], label_pos[,1])
all_y <- c(layout_circle[,2], label_pos[,2])

# Set plot limits with some margin
xlim <- range(all_x) + c(-0.2, 0.2)
ylim <- range(all_y) + c(-0.2, 0.2)

# Compute angles for label orientation
angles <- atan2(layout_circle[, 2], layout_circle[, 1]) * 180 / pi

# Adjust angles for readability (flip upside down labels)
label_angles <- ifelse(angles < -90 | angles > 90, angles + 180, angles)

# Determine alignment (right/left) based on angle
label_adj <- ifelse(angles < -90 | angles > 90, 1, 0)

# Plot the graph without vertex labels
plot( g,
      layout = layout_circle,
      edge.width = E(g)$weight/5,
      vertex.label = NA,
      vertex.size = 4,
      vertex.color=viridis(n = 5, option = "C")[3],
      vertex.frame.color = viridis(n = 5, option = "C")[3],
      xlim = xlim,
      ylim = ylim,
      main = ""
)

# Add rotated labels manually
for (i in seq_along(labels)) {
  text(
    x = layout_circle[i, 1],
    y = layout_circle[i, 2],
    labels[i],
    srt = label_angles[i],
    adj = c(label_adj[i], 0.5),
    cex = 0.8,
    family = "LM Roman 10" 
  )
}
if(exportFig){dev.off()}




####################################################################################
###########################        INFERENCE    #######################################
######################################################################################"


mySimpleSBMPoisson <- tree_tree_weighted  %>%
  estimateSimpleSBM("poisson", directed = FALSE,
                    estimOptions = list(verbosity = 0 , plot = FALSE),
                    dimLabels = c('Trees'))


## ----simpleSBMfitPoisson plot1--------------------------------------------------
myColour =  viridis(n = 5, option = "C")[4]
results_model_select <- mySimpleSBMPoisson$storedModels
pPoisson <-  ggplot(results_model_select) + aes(x = nbBlocks, y = ICL)  + geom_line(color =myColour) + geom_point(shape = 24,size=4,color = myColour) +theme(axis.text=element_text(size=9))
pPoisson <- pPoisson + geom_vline(xintercept = mySimpleSBMPoisson$nbBlocks,linetype="dashed",color=myColour)
#pPoisson + geom_segment(aes(x = mySimpleSBMPoisson$nbBlocks, y = min(mySimpleSBMPoisson$storedModels$ICL) , xend = mySimpleSBMPoisson$nbBlocks, yend =  max(mySimpleSBMPoisson$ICL)))
pPoisson <- pPoisson + scale_x_discrete(limits=c("1", "2","3","4","5","6","7","8","9")) + xlab("Number of blocks")
pPoisson

if(exportFig)
  ggsave(paste0(where_fig,"/FungusTree_Poisson_ICL.png"),
         plot = pPoisson,
         width = graph_width,
         height=graph_height,
         units = graph_unit)

## ----simpleSBMfit plot1---------------------------------------------------------
g <- plot(mySimpleSBMPoisson, type = "expected", dimLabels = list(row = 'Trees', col = 'Trees'), plotOptions = list(legend = TRUE,rowNames = TRUE,colNames = TRUE))
g + theme(axis.text=element_text(size=8))

if(exportFig)
  ggsave(paste0(where_fig,"/FungusTree_Poisson_expected.png"),
         plot = g,
         scale = 1,
         width = 20,
         units = c("cm"))

g <- plot(mySimpleSBMPoisson, type = "data", dimLabels = list(row = 'Trees', col = 'Trees'), plotOptions = list(legend = TRUE,rowNames = TRUE,colNames = TRUE))
g <- g + theme(axis.text=element_text(size=8))

if(exportFig)
  ggsave(paste0(where_fig,"/FungusTree_Poisson_dataordered.png"),
         plot = g,
         scale = 1,
         width = 20,
         units = c("cm"))


