# PATH (Landau pre-print) analysis

library(tidyverse)
library(parallel)
library(Matrix)
library(ape)
library(PATH)


##


# Set paths
cov <- 'GBC'
sample <- 'MDA_PT'
path_sample <- paste0('/Users/IEO5505/Desktop/mito_bench/results/phylo_inference/top_trees_MQuad/', sample)


##


# Load data
# Tree
tree <- ape::read.tree(paste0(path_sample, '/tree.newick'))
# Cat variable, to numeric
cov_df <- read.csv(paste0(path_sample, '/cov_df_', cov, '.csv'), row.names = 1)[tree$tip.label,]
table(cov_df)
X <- state2mat.sparse(cov_df)


##


# Phylo-correlation cat 

# Phylogenetic weight matrix
Winv <- inv.tree.dist(tree, node=TRUE, norm=FALSE)
modxcor <- xcor(X, Winv)
diag(modxcor$Z.score)
diag(modxcor$one.sided.pvalue)

# Format result
Idf <- reshape2::melt(modxcor$Morans.I, value.name = "I")
Zdf <- reshape2::melt(modxcor$Z.score, value.name = "Z")
pdf <- reshape2::melt(modxcor$one.sided.pvalue, value.name = "p")
dfs <- list(Idf, Zdf, pdf)
results <- Reduce(
  function(x, y) full_join(x, y, by = c("Var1", "Var2")), 
  dfs
  ) %>%
  mutate(sample=sample) %>% 
  rename(GBC1=Var1, GBC2=Var2)

# Save
write.csv(results, paste0(path_sample, '/phylocorr_', cov, '.csv'))


##
