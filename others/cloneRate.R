# Test association with lineage

# CloneRate vignette

# Code
library(tidyverse)
library(parallel)
library(Matrix)
library(ape)
library(PATH)


##


## 
# path_nodes <- '/Users/IEO5505/Desktop/nodes.csv'
# path_edges <- '/Users/IEO5505/Desktop/edges.csv'
# path_meta <- '/Users/IEO5505/Desktop/meta.csv'
# lineage_column <- 'aggregated_ct'
##


##


# Load data

# Reconstruct tree
nodes <- read.csv(path_nodes, row.names=1)[1:3] %>% distinct   # Avoid counting >1 nodes with >1 assigned mut
edges <- read.csv(path_edges, row.names=1)
tree <- list(
  edge = edges %>% select(u, v) %>% as.matrix,
  Nnode = sum(!is.na(nodes$support)),
  tip.label = nodes %>% arrange(node) %>% drop_na(cell) %>% .$cell,
  edge.length = edges$branch_lenghts
)
class(tree) <- "phylo"