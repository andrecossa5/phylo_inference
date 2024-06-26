# Prune Tree

# Workflow readapted from Weng et al., Nature 2024

# Code
library(argparse)

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


# Args
parser <- ArgumentParser(description="prune_final_tree: assign mutations to a MT-tree branches, cut the tree in discrete clones.")
parser$add_argument("--tree", required=TRUE, help="Final tree. Default: final_tree.newick")
parser$add_argument("--AD", required=TRUE, help="AD matrix. Default: AD.csv")
parser$add_argument("--DP", required=TRUE, help="DP matrix. Default: DP.csv")
parser$add_argument("--af_t", required=TRUE, help="AF treshold to call a cell positive. Default: 0.05")
parser$add_argument("--ncores", required=TRUE, help="ncores. Default: 8")
parser$add_argument("--prob_cut", required=TRUE, help="prob_cut. Default: min probability to retain a branch-assigned mutation when tree cutting.")
parser$add_argument("--min_cell_clone", required=TRUE, help="Min n of cells for a clone. Default: 5")

# Parse 
args <- parser$parse_args()
path_tree <- args$tree
path_AD <- args$AD
path_DP <- args$DP
af.t <- as.numeric(args$af_t)
n.cores <- as.numeric(args$ncores)
prob.cut <- as.numeric(args$prob_cut)
MinCell <- as.numeric(args$min_cell_clone)


##


# Read and format inputs
tree <- ape::read.tree(path_tree)
mtr <- read.csv(path_AD, row.names=1) %>% as.matrix() %>% t
depth <- read.csv(path_DP, row.names=1) %>% as.matrix() %>% t
mtr.bi <- matrix((mtr/(depth+0.000001) >= af.t) %>% as.numeric(), nrow=nrow(mtr), ncol=ncol(mtr))
row.names(mtr.bi) <- row.names(mtr)
colnames(mtr.bi) <- colnames(mtr)
mtr.bi.t <- t(mtr.bi) # each row is a cell
mtr.bi <- mtr.bi[,tree$tip.label]

# Assign variants, cut tree, drop tips 
var_assignment <- assign_variants(tree, mtr, mtr.bi, depth, n.cores=n.cores, assign.to='Child')
res_cut_tree <- cut_tree(tree, var_assignment, prob.cut=prob.cut, MinCell=MinCell)
tree <- drop.tip(tree, setdiff(tree$tip.label, res_cut_tree$cell_assignment$Cell))

# Write tree
write.csv(var_assignment, 'var_assignment.csv')
write.csv(res_cut_tree$cell_assignment, 'cell_assignment.csv')
ape::write.tree(tree, 'final_tree.newick')


##

