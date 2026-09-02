rm(list=ls())
library(igraph)
library(sbm)

require(ggplot2)
library(viridisLite)
library(stringr)
library(bipartite)

where_fig = getwd()

# Unified theme for figures
theme_set(theme_bw(base_family = "LM Roman 10"))
graph_unit <- "px"
graph_width <- 1500
graph_height <- 1500
exportFig <- FALSE



data("fungusTreeNetwork")
fungus_tree <- t(fungusTreeNetwork$fungus_tree)
row.names(fungus_tree)  <- sub("\\s*\\(.*", "",str_sub(fungusTreeNetwork$tree_names, 1, 30))
colnames(fungus_tree)  <- sub("\\s*\\(.*", "",str_sub(fungusTreeNetwork$fungus_names, 1, 30))


################### bipartite network



ntree= nrow(fungus_tree)
nfung= ncol(fungus_tree)

# Plot the bipartite network with vertical labels for the lower level (species)

par(mar = c(1, 1, 1, 1))

plotweb(
  fungus_tree,
  sorting = "normal",
  
  # Labels
  higher_labels = names(colSums(fungus_tree)),
  lower_labels = rownames(rowSums(fungus_tree)),
  text_size = 0.35,
  
  # Géométrie
  spacing = 0.50,
  box_size = 0.08,
  lab_distance = 0.03,
  # Divers
  horizontal = FALSE,
  plot_axes = FALSE,
  
)


top_trees  <- order(rowSums(fungus_tree), decreasing = TRUE)[1:10]
top_fungi  <- order(colSums(fungus_tree), decreasing = TRUE)[1:10]

fungus_tree_10 <- fungus_tree[top_trees, top_fungi]
par(mar = c(1, 1, 1, 1))
plotweb(
  fungus_tree_10,
  higher_labels = names(colSums(fungus_tree_10)),
  lower_labels = rownames(rowSums(fungus_tree_10)),
  sorting = "normal",
)




##################### incidence matrix plot
pBip <- plotMyMatrix(fungus_tree, dimLabels=list(col = 'Fungus',col = 'Tree'), plotOptions = list(legend = FALSE,rowNames = TRUE,colNames = TRUE))
pBip <- pBip + theme(axis.text.x = element_text( size = 5),
                     axis.text.y = element_text( size = 5))  
pBip

if(exportFig)
  ggsave(paste0(where_fig,"/FungusTree_Biparite_data.png"),
         plot = pBip ,
         scale = 1,
         width = 30,
         units = c("cm"))


## ----tree_fungi_bipartite network, eval=TRUE, echo = TRUE-----------------------
myBipartiteSBM <- estimateBipartiteSBM(
  netMat =fungus_tree,
  model = 'bernoulli',
  dimLabels=c(col = 'Fungus',row = 'Tree'),
  estimOptions = list(verbosity = 0, plot=FALSE))




## ----plot bipartite estim-------------------------------------------------------
plot(myBipartiteSBM, type = "data")
g <- plot(myBipartiteSBM, type = "expected", dimLabels = list(col = 'Fungus',row = 'Tree'), plotOptions = list(legend = FALSE,rowNames = TRUE,colNames = TRUE))
g <- g + theme(axis.text=element_text(size=5))
g
if(exportFig)
  ggsave(paste0(where_fig,"/FungusTree_Bipartite_expected.png"),
         plot = g,
         scale = 1,
         width = 30,
         units = c("cm"))

## ----plot bipartite expect------------------------------------------------------
g <- plot(myBipartiteSBM, type = "data", dimLabels = list(col = 'Fungus',row = 'Tree'), plotOptions = list(legend = FALSE,rowNames = TRUE,colNames = TRUE))
g <- g + theme(axis.text=element_text(size=5))
g
if(exportFig)
  ggsave(paste0(where_fig,"/FungusTree_Bipartite_dataordered.png"),
         plot = g,
         scale = 1,
         width = 30,
         units = c("cm"))

## ----plot BM network tree fungi,  echo=TRUE, eval = TRUE------------------------
plot(myBipartiteSBM, type = "meso", plotOptions = list(vertex.size=c(0.5,1) , edge.width  = 1))
