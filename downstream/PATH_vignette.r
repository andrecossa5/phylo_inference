# PATH (Landau pre-print) analysis

library(tidyverse)
library(parallel)
library(Matrix)
library(ape)
library(fgsea)
library(PATH)


##


# Bug fix
parM_I <- function(d,w, break.length=100) {
  
  
  y <- c(seq(1,ncol(d), break.length), (ncol(d)+1) )
  y2 <- y[-1] - 1
  y1 <- y[-length(y)]
  
  mfunc <- function(b, e, d1=d, w1=w) {
    M <- xcor(d1[,b:e ], w1)
    x1 <- diag(M$Z)
    x2 <- diag(M$Morans.I)
    out <- data.frame("Z"=x1, "I"=x2)
    return(out)
  }
  
  out <- do.call(rbind,
                 mcmapply(function(b1,e1) mfunc(b=b1,e=e1, d, w), b1=y1, e1=y2,
                          mc.cores = detectCores(), SIMPLIFY = F))
  return(out)
}


##


xcor_gsea_I <- function(gene.data, weight.matrix, pathways, maxS=Inf, minS=5, nperm=10000) {
  z <- parM_I(gene.data, weight.matrix)
  z$gene <- rownames(z)
  
  z <- z %>% as_tibble() %>% group_by(gene) %>% dplyr::select(Z) %>% arrange(desc(Z)) %$% set_names(Z, gene)
  
  out <- fgsea(pathways, z, maxSize=maxS, minSize=minS, scoreType="std", nPermSimple=nperm)
  return(out)
}


##


# See data
tree115               # Phylo with $edge.length and $edge full.
modules115 %>% dim()  # df cell x gene module    
genes115 %>% dim()    # df cell x gene

row.names(modules115) == tree115$tip.label
# class(tree115)
# plot(tree115)


##

# Phylogenetic correlation

# Phylogenetic weight matrix
Winv <- inv.tree.dist(tree115, node=TRUE, norm=FALSE)
class(Winv)

# Compute phylogenetic correlations between GBM modules.
class(modules115) # Matrix, NB
modxcor <- xcor(modules115 %>% as.matrix(), Winv)

# Process xcor output for plotting and visualization.
Idf <- reshape2::melt(modxcor$Morans.I, value.name = "I")
Zdf <- reshape2::melt(modxcor$Z.score, value.name = "Z")
df <- full_join(Idf, Zdf, by=c("Var1", "Var2"))
df <- df %>% mutate(Var1 = as.factor(Var1), Var2 = as.factor(Var2))


##


# Correlation and transition probabilities categorical
table(colnames(modules115)[max.col(modules115)])
states4 <- max.col(modules115)
states3 <- ifelse(states4 == 1, 2, states4) - 1

# Add categorical state data to phylogeny
tree115$states <- states3                        # Add wichever categorical column!

# Corr with state3
M <- states3 %>% as.matrix()
rownames(M) <- tree115$tip.label
colnames(M) <- 'AA'
modxcor <- xcor(M, Winv)

# Use PATH to infer transition probabilities between stem-like, AC, and MES cells
Pinf <- PATH.inference(tree115, impute_branches = T, sample_rate_est = 10^-6)

# Format PATH output
Ptinf.df <- data.frame(Pinf$Pt)
rownames(Ptinf.df) <- c("Stem", "AC", "MES")


##


# Discovery of new modules from HVGs expression and a set of genes
genesets = msigdbr::msigdbr(species = "Homo sapiens", subcategory = "CGP")
pathwaysH = split(x = genesets$gene_symbol, f = genesets$gs_name)
pathwaysH <- lapply(pathwaysH, unique)
pathwaysH <- c(pathwaysH, as.list(neftel.genes))   # A list of gene names

# Correlation with gene sets
class(genes115)
class(Winv)
class(pathwaysH)
xcgsea <- xcor_gsea_I(genes115, Winv, pathwaysH) %>% 
  arrange(desc(-padj)) %>% 
  mutate(rank=rank(padj)) 


##



