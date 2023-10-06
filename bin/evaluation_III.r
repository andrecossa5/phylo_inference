#!/usr/bin/env Rscript

# R phylo script

# Code
library(tidyverse)
library(parallel)
library(reticulate)
library(ape)
library(picante)
library(phytools)

# use_condaenv("MI_TO")                        # Assume we are in our docker image
scipy_sparse <- import("scipy.sparse")

# Parse args
args <- commandArgs(trailingOnly = TRUE)
path_tree <- args[1]
path_input <- args[2]
ncores <- args[3]

# Paths
path_meta <- paste0(path_input, '/meta.csv')
path_variants <- paste0(path_input, '/variants.csv')
path_AD <- paste0(path_input, '/AD.npz')
path_DP <- paste0(path_input, '/DP.npz')

# Read data
tree <- ape::read.tree(path_tree)
meta <- read_csv(path_meta) %>% column_to_rownames(var = "index")
variants <- read_csv(path_variants, col_names=FALSE)$X1
AD <- scipy_sparse$load_npz(path_AD) %>% as.matrix()
DP <- scipy_sparse$load_npz(path_DP) %>% as.matrix()
afm <- as.data.frame(AD/DP)
colnames(afm) <- variants
row.names(afm) <- row.names(meta)

# Calculate phylogenetic signal of individual muts
L <- mclapply(
  colnames(afm), 
  function(mut) {
    x <- setNames(afm[tree$tip.label, mut], tree$tip.label)
    lambda <- phytools::phylosig(tree, x, method='lambda', test=TRUE)
    k <- phytools::phylosig(tree, x, method='K', test=TRUE)
    s <- c(lambda$lambda, k$K, lambda$P, k$P)
  }, 
  mc.cores=ncores
)
muts_df <- data.frame(setNames(L, colnames(afm)), check.names=FALSE) %>% t()
colnames(muts_df) <- c('lambda', 'K', 'lambda_p', 'k_p')

# Calculate phylogenetic clustering of clonal labels
clones_df <- picante::ses.mpd(table(meta$GBC, row.names(meta)) %>% as.matrix(), cophenetic(tree))

# Save all
write.csv(muts_df, 'muts.csv')
write.csv(clones_df, 'clones.csv')


##