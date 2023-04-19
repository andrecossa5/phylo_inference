# Phylo script for first trees with ape and ggtree

# Code
library(Rfast)
library(ape)
library(phangorn)
library(tidyverse)
library(ggtree)
library(TreeTools)
library(treeio)
library(viridis)
library(phytools)


##


################################################################ Utils
annotate_tree <- function(tree, df, cov='GBC') {

    if (is.character(df[[cov]])) {
        labels <- df[[cov]] %>% unique()
        L <- lapply( labels, function(x) { x <- df %>% filter(.data[[cov]] == x) %>% row.names() } )
        names(L) <- labels
        tree <- groupOTU(tree, L)

    } else if (is.numeric(df[[cov]])) {
        tree <- full_join(
            tree %>% as_tibble(), 
            df %>% select(cov) %>% rownames_to_column('label'), 
            by = 'label'
        ) %>% as.treedata()
    }
    return(tree)
}


##


prep_data <- function(path_phylo_data, path_main, sample, method) {
    df <- read.csv(paste0(path_phylo_data, sample, '_', method, '_50_50_cell_x_vars.csv'), row.names=1)
    X_char <- df %>% select(-GBC)
    cells <- row.names(X_char)
    colors <- read.csv(paste0(path_main, 'data/', sample, '_clones_colors.csv'), row.names=1)
    colors <- colors[row.names(colors) %in% unique(df$GBC), ]
    return(list(df, X_char, colors))
}


##


plot_tree_cov <- function(tree, cov='group', colors=NULL) {
    p <- ggtree(treeNJ, aes(color=.data[[cov]]), layout='circular') +
    theme(legend.position="none") + 
    scale_color_manual(values=colors, name=NULL)
    return(p)
}


##


bootstrapping_tree <- function(tree, X_char, FUN=function(xx) NJ(dist(xx)), B=100, ncores=8) {
    
    # Bootstrap a tree
    set.seed(1234)
    bp <- boot.phylo(tree, X_char, FUN, mc.cores=ncores, B=B, trees=TRUE)
    bp$BP[is.na(bp$BP)] <- 0

    # Annotate it
    df_boot <- tibble(node=1:Nnode(tree)+Ntip(tree), bootstrap=bp$BP/B)
    df_boot$support_cat <- cut(df_boot$bootstrap, c(0, 0.5, 0.75, 1), include.lowest=TRUE)
    tree <- full_join(tree, df_boot, by="node")

    return(tree)
}


##


consensus_tree <- function(tree, X_char, FUN=function(xx) NJ(dist(xx)), 
    B=100, ncores=8, p=0.5, method='average') {
    
    # Bootstrap a tree
    set.seed(1234)
    bp <- boot.phylo(tree, X_char, FUN, mc.cores=ncores, B=B, trees=TRUE)
    bp$BP[is.na(bp$BP)] <- 0

    # Create a list of trees from the bootstrap, and the consensus tree
    tree.list <- lapply(1:B, function(i) {bp$trees[[i]]})

    if (method == 'consensus') {
    con.tree <- Consensus(tree.list, p=0.5)
    } else if (method == 'average') {
        con.tree <- averageTree(tree.list, method="symmetric.difference")
    }

    return(con.tree)
}


##


plot_boot_support <- function(tree, colors=c("white", "#ce8312", "firebrick")) {

    mean.supp <- mean(tree@data$bootstrap[!is.na(tree@data$bootstrap)])
    p <- ggtree(tree, layout='circular') + 
    geom_text(aes(label=bootstrap), size=2) +
    geom_point2(aes(subset=!isTip, fill=support_cat), shape=21, size=5) + 
    theme_tree(legend.position=c(0.85, 0.9)) +
    scale_fill_manual(
        values=colors, 
        guide='legend', 
        name=paste0('Bootstrap support (mean: ', round(mean.supp, 2), ')'),
        breaks=c('[0,0.5]', '(0.5,0.75]', '(0.75,1]')
    ) 
    return(p)
}
################################################################


##


# Set paths 
path_main <- '/Users/IEO5505/Desktop/MI_TO/'
path_results <- paste0(path_main, '/results_and_plots/phylo/')
path_phylo_data <- paste0(path_main, '/data/phylo_input/')

?paste0

# Read charachter and euclidean distance matrices for each sample
for (sample in c('AML', 'MDA', 'PDX')) {

    if (sample == 'PDX') {method <- 'pegasus'} else {method <- 'MQuad'}
    L <- prep_data(path_phylo_data, path_main, sample, method)
    df <- L[[1]]; X_char <- L[[2]]; colors <- L[[3]]

    # Compute simple trees
    treeNJ  <- NJ(dist(X_char))
    treeUPGMA  <- upgma(dist(X_char))

    # Annotate with clones
    treeNJ <- annotate_tree(treeNJ, df)
    treeUPGMA <- annotate_tree(treeUPGMA, df)

    # Plot and save
    pdf(paste0(path_results, 'NJ_', sample, '_first.pdf'), width=15, height=15)
    p <- plot_tree_cov(treeNJ, colors=colors)
    print(p)
    dev.off()

    pdf(paste0(path_results, 'UPGMA_', sample, '_first.pdf'), width=15, height=15)
    p <- plot_tree_cov(treeUPGMA, colors=colors)
    print(p)
    dev.off()
}


##


# Top clone and its muts example: AML, TTACCCTGTGCGACCCAG clone
sample <- 'AML'; method <- 'MQuad'
L <- prep_data(path_phylo_data, path_main, sample, method)
df <- L[[1]]; X_char <- L[[2]]; colors <- L[[3]]

# Compute NJ tree
treeNJ  <- NJ(dist(X_char))

# Annotate with clone type and top important muts for clone ...
clone <- 'TTACCCTGTGCGACCCAG'
muts <- c('X2946_C', 'X15172_A', 'X8464_T')
df <- df %>% mutate(
    clone = case_when(
        GBC == clone ~ clone, 
        GBC != clone ~ 'other'
    )
)

# Plot 
pdf(paste0(path_results, 'NJ_', clone, '_example.pdf'), width=15, height=15)

for (x in c('clone', muts)) { 
    treeNJ <- annotate_tree(treeNJ, df, x) 
    if (x == 'clone') {
        p <- ggtree(treeNJ, aes(color=group), layout='circular') +
        theme(legend.position="none") +
        scale_color_manual(values=c('black', 'red'), name='Clone type')
    } else {
        p <- ggtree(treeNJ, aes(color=.data[[x]]), layout='circular') +
        theme(legend.position="none") +
        scale_color_viridis(name=x) 
    }
    print(p)
}

dev.off()


##


## Bootstrapping trees, boot.phylo version
for (sample in c('AML', 'MDA', 'PDX')) {

    if (sample == 'PDX') {method <- 'pegasus'} else {method <- 'MQuad'}

    L <- prep_data(path_phylo_data, path_main, sample, method)
    df <- L[[1]]; X_char <- L[[2]]; colors <- L[[3]]
    
    # Compute simple trees
    treeNJ  <- NJ(dist(X_char))
    treeUPGMA  <- upgma(dist(X_char))

    # Bootstrap both trees
    boot.NJ <- bootstrapping_tree(treeNJ, X_char, FUN=function(xx) NJ(dist(xx)))
    boot.UPGMA <- bootstrapping_tree(treeUPGMA, X_char, FUN=function(xx) upgma(dist(xx)))

    # Plot and save
    pdf(paste0(path_results, 'NJ_bootstap_', sample, '_first.pdf'), width=15, height=15)
    p <- plot_boot_support(boot.NJ)
    print(p)
    dev.off()

    pdf(paste0(path_results, 'UPGMA_bootstrap_', sample, '_first.pdf'), width=15, height=15)
    p <- plot_boot_support(boot.UPGMA)
    print(p)
    dev.off()
}


##


## Consensus trees
av.tree <- consensus_tree(treeNJ, X_char, p=0.7, method='average')
ggtree(av.tree, layout='circular')