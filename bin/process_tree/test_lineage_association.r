# Test association with lineage

# PATH (Landau pre-print) analysis

# Code
library(argparse)
library(tidyverse)
library(parallel)
library(Matrix)
library(ape)
library(PATH)


##


# Args
parser <- ArgumentParser(description="test_association with lineage.")
parser$add_argument("--nodes", required=TRUE, help="Nodes df. Default: nodes.csv")
parser$add_argument("--edges", required=TRUE, help="Edges df. Default: edges.csv")
parser$add_argument("--meta", required=TRUE, help="Meta df. Default: meta.csv")
parser$add_argument("--lineage_column", required=TRUE, help="Lineage column. Default: GBC")

# Parse 
args <- parser$parse_args()
path_nodes <- args$nodes
path_edges <- args$edges
path_meta <- args$meta
lineage_column <- args$lineage_column

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

# Load meta
meta <- read.csv(path_meta, row.names=1)
cov <- meta[tree$tip.label,][[lineage_column]]
cats <- cov %>% unique
cov <- factor(cov, levels=cats)
X <- state2mat.sparse(cov %>% as.numeric)

##

# PATH analysis
Winv <- inv.tree.dist(tree, node=TRUE, norm=FALSE) 
modxcor <- xcor(X, Winv)
Idf <- reshape2::melt(modxcor$Morans.I, value.name="I")
Zdf <- reshape2::melt(modxcor$Z.score, value.name="Z")
pdf <- reshape2::melt(modxcor$one.sided.pvalue, value.name="p")
dfs <- list(Idf, Zdf, pdf)
results <- Reduce(function(x, y) full_join(x, y, by=c("Var1", "Var2")), dfs)
mapping <- setNames(1:length(cats), cats)
results$Var1 <- sapply(results$Var1, function(x){ names(which(mapping==x)) })
results$Var2 <- sapply(results$Var2, function(x){ names(which(mapping==x)) })

# Save
write.csv(results, 'lineage_association.csv')


##
