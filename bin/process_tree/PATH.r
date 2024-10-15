
#!/usr/bin/env Rscript

library(argparse)
library(dplyr)
library(reshape2)
library(parallel)
library(Matrix)
library(ape)
library(PATH)


##


# cov <- 'GBC'
# sample <- 'MDA_PT'
# path_tree <- '/Users/IEO5505/Desktop/MI_TO/phylo_inference/scratch/final_tree.newick'
# path_meta <- '/Users/IEO5505/Desktop/MI_TO/phylo_inference/scratch/cells_meta.csv'

parser <- ArgumentParser(description = 'PATH wrapper')

parser$add_argument("path_tree", type = "character", help = "Path to .newick tree")
parser$add_argument("path_meta", type = "character", help = "Path to cells meta data")
parser$add_argument("lineage_column", type = "character", help = "Lineage column to evaluate phylocorr")
args <- parser$parse_args()


##


# Processing
tree <- ape::read.tree(args$path_tree)
cov_df <- read.csv(args$path_meta, row.names = 1)[tree$tip.label,]
x_factor <- as.factor(cov_df[[args$lineage_column]])
X <- state2mat.sparse(x_factor)

Winv <- inv.tree.dist(tree, node=TRUE, norm=FALSE)
modxcor <- xcor(X, Winv)
Idf <- reshape2::melt(modxcor$Morans.I, value.name = "I")
Zdf <- reshape2::melt(modxcor$Z.score, value.name = "Z")
pdf <- reshape2::melt(modxcor$one.sided.pvalue, value.name = "p")
dfs <- list(Idf, Zdf, pdf)

# Save
results <- Reduce(function(x, y) full_join(x, y, by = c("Var1", "Var2")), dfs)
results$Var1 <- x_factor[results$Var1]
results$Var2 <- x_factor[results$Var2]

write.csv(results, 'phylocorr.csv')


##

