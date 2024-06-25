# Prune Tree

# Workflow readapted from Weng et al., Nature 2024

# Code
library(ape)
library(phangorn)
library(tidyverse)
library(foreach)
library(doParallel)
library(doMC)
library(pryr)

# Source utils
get_script_dir <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    return(dirname(rstudioapi::getActiveDocumentContext()$path))
  } else {
    cmd_args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("--file=", cmd_args, value = TRUE)
    if (length(file_arg) > 0) {
      return(dirname(normalizePath(sub("--file=", "", file_arg))))
    } 
  }
}
source(file.path(get_script_dir(), "utils.r"))


##


# Workflow

# Args. TODO: move to command line
af.t <- 0.05
n.cores <- 8
prob.cut <- 0.3
MinCell <- 5


##


# Read and format inputs
tree <- ape::read.tree('tree.newick')
mtr <- read.csv('input_folder/AD.csv'), row.names = 1) %>% as.matrix() %>% t
depth <- read.csv('input_folder/DP.csv'), row.names = 1) %>% as.matrix() %>% t
mtr.bi <- matrix((mtr/(depth+0.000001) >= af.t) %>% as.numeric(), nrow=nrow(mtr), ncol=ncol(mtr))
row.names(mtr.bi) <- row.names(mtr)
colnames(mtr.bi) <- colnames(mtr)
mtr.bi.t<-t(mtr.bi) # each row is a cell
mtr.bi<-mtr.bi[,tree$tip.label]


##


# Assign variants, cut tree, drop tips 
var_assignment <- assign_variants(tree, mtr, mtr.bi, depth, n.cores=n.cores, assign.to='Child')
res_cut_tree <- cut_tree(tree, var_assignment, prob.cut=prob.cut, MinCell=MinCell)
tree <- drop.tip(tree, setdiff(tree$tip.label, res_cut_tree$cell_assignment$Cell))

# Write tree
write.csv(var_assignment, 'var_assignment.csv')
write.csv(res_cut_tree$cell_assignment, 'cell_assignment.csv')
write.tree(tree, 'tree_pruned.newick')


##


# Viz
# library(RColorBrewer)
# library(ggtree)
# library(ggtreeExtra)
# 
# cell_assignment$Clone_merge <- as.factor(cell_assignment$Clone_merge)
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# col_vector<-sample(col_vector,length(levels(cell_assignment$Clone_merge)),replace=T)
# 
# p<-ggtree(tree, layout="circular", branch.length='none') +
#   geom_fruit( 
#   data=cell_assignment, 
#   mapping=aes(y=Cell,x=2,fill=Clone_merge), geom=geom_tile,
#   pwidth=0.001, 
#   width=3, 
#   offset=0.05
#   )+
#   scale_fill_manual(values=col_vector)
# 
# pdf(paste0(input_path, 'tree.pdf'), height=10)
# print(p)
# dev.off()




