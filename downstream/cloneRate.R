# Test association with lineage

# CloneRate vignette

# Code
library(tidyverse)
library(parallel)
library(Matrix)
library(ape)
library(cloneRate)
library(ggplot2)


##


## 
path_nodes <- '/Users/ieo6983/Desktop/phylo_inference/data/nodes.csv'
path_edges <- '/Users/ieo6983/Desktop/phylo_inference/data/edges.csv'
# path_meta <- '/Users/IEO5505/Desktop/meta.csv'
# lineage_column <- 'aggregated_ct'
##


##


# Load data

# Reconstruct tree
nodes <- read.csv(path_nodes, row.names=1)[1:6] %>% distinct   # Avoid counting >1 nodes with >1 assigned mut
edges <- read.csv(path_edges, row.names=1)
tree <- list(
  edge = edges %>% dplyr::select(u, v) %>% as.matrix,
  Nnode = sum(!is.na(nodes$support)),
  tip.label = nodes %>% arrange(node) %>% drop_na(cell) %>% .$cell,
  edge.length = edges$branch_lenghts
)
class(tree) <- "phylo"

str(tree)
head(nodes)


##


# Convert the tree to an ultrametric tree
ultrametric_tree <- chronos(tree)
str(ultrametric_tree)
#plot(ultrametric_tree)


##


# Split clades
node_numbers <- c(unique(nodes$clonal_ancestor)[!is.na(unique(nodes$clonal_ancestor))])
clades <- list()

for (node in node_numbers) {
  clades[[as.character(node)]] <- extract.clade(ultrametric_tree, node)
}

# You can now access individual clades using their node numbers
plot.phylo(clades$`302`,
           direction = "downwards",
           show.tip.label = T, edge.width = 1.5
)
axisPhylo(side = 2, backward = FALSE, las = 1)
title(main = "Tree clone 302")


##


# Calculate growth rate using cloneRate main function 

# Full tree
resultsLengths <- internalLengths(ultrametric_tree)
# Clades 
resultsLengths <- internalLengths(clades)

#' @Warning: internalLengths() is not suitable for these trees. 
#' Using birthDeathMCMC() function is suggested 
#' 
#' For more @info about this: 
#' Why birthDeathMCMC: https://cran.r-project.org/web/packages/cloneRate/vignettes/cloneRate-simulate.html 
#' About birth-death models: https://bio.libretexts.org/Bookshelves/Evolutionary_Developmental_Biology/Phylogenetic_Comparative_Methods_(Harmon)/10%3A_Introduction_to_Birth-Death_Models 


##


# Calculate growth rate with birthDeathMCMC()
results <- cloneRate::birthDeathMCMC(clades, 
                                     maxGrowthRate = 4, # def
                                     alpha = 0.05) # def
View(results)


##


