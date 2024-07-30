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
parser <- ArgumentParser(description="prune_raw_tree: assign mutations to a MT-tree branches, cut the tree in discrete clones.")
parser$add_argument("--tree", required=TRUE, help="Raw tree. Default: raw_tree.newick")
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

main <- '/Users/IEO5505/Desktop/AML_clonal_reconstruction/results/phylo/barcodes/AML2/set1/'
path_tree <- paste0(main, 'tree.newick')
path_AD <- paste0(main, 'AD.csv')
path_DP <- paste0(main, 'DP.csv')
af.t <- 0.01

af.t <- as.numeric(args$af_t)
n.cores <- as.numeric(args$ncores)
prob.cut <- as.numeric(args$prob_cut)
MinCell <- as.numeric(args$min_cell_clone)


##


# Read and format inputs
tree <- ape::read.tree(path_tree)
mtr <- read.csv(path_AD, row.names=1) %>% as.matrix() %>% t %>% .[,tree$tip.label]
depth <- read.csv(path_DP, row.names=1) %>% as.matrix() %>% t %>% .[,tree$tip.label]
mtr.bi <- matrix((mtr/(depth+0.000001) >= af.t), nrow=nrow(mtr), ncol=ncol(mtr))
row.names(mtr.bi) <- row.names(mtr)
colnames(mtr.bi) <- colnames(mtr)
mtr.bi.t <- t(mtr.bi) # each row is a cell
muts <- colnames(mtr.bi.t)[which(colSums(mtr.bi.t)>0)]
mtr.bi <- mtr.bi[muts,]
mtr.bi.t <- mtr.bi.t[, muts]
depth <- depth[muts,]
mtr <- mtr[muts,]


# Assign variants, cut tree, drop tips 
# var_assignment <- assign_variants(tree, mtr, mtr.bi, depth, n.cores=n.cores, assign.to='Child')

# Get cell profiles
df <- reconstruct_genotype_summary(tree)
df_profile_mtx.t <- df$profile_int %>% do.call(rbind,.) %>% t

## Set expected global heteroplasmy level 
af.m <- mtr/(depth+0.000001)
Global.heteroplasmy <- af.m[!is.na(af.m)] %>% .[.!=0] %>% median  # this is a number ~ 0.03-0.05, an expected p for modeling dropout
if (Global.heteroplasmy>=0.1){Global.heteroplasmy <- 0.1}                         # Cut for numeric reasons...

## Generate a binomial-based probability distribution for each cell to model the dropout
Zero.p <- dbinom(0,depth,Global.heteroplasmy) # PROBABILITY OF OBSERVING AF=0, with a certain global AF
Zero.p[Zero.p==1] <- 0.99
Zero.logV <- log(Zero.p)  # A matrix of log(probability of zero given w/ mutation in the cell)  Variant vs cell
NonZero.logV <- log(1-Zero.p) # A matrix of log(probability of 1 or >=1 given w/ mutation in the cell)  Variant vs cell
Zero.logP <- log(0.95) 
NonZero.logP <- log(1-0.95) 

# Calculate LogLihelihoods for each variant


registerDoMC(n.cores)
Loglik <- foreach (i=1:nrow(mtr.bi), .combine='rbind') %dopar%  # To loop through all variants
{
    
    ## Make indicator matrices for 1-1 (Inside the clade, and is mutation) and 0-0 (outside of clade, and is not mutation) 
    x <- df_profile_mtx.t + mtr.bi[i,]
    
    x_11 <- matrix(nrow=nrow(x),ncol=ncol(x)) # cell x branch
    x_00 <- matrix(nrow=nrow(x),ncol=ncol(x)) # cell x branch
    x_11[which(x!=2)] <- 0                   
    x_11[which(x==2)] <- 1                   # Cell inside the branch with mutation 
    x_00[which(x!=0)] <- -1  
    x_00[which(x==0)] <- 1                   # Cell outside the branch and no mutation 
    x_00[which(x_00==-1)] <-0
    
    ## Make indicator matrices for 1-0 and 0-1 scenarios
    x_10 <- df_profile_mtx.t - mtr.bi[i,] 
    x_10[which(x_10==-1)] <- 0               # Cell inside the branch and no mutation 
    x_01 <- mtr.bi[i,] - df_profile_mtx.t
    x_01[which(x_01==-1)] <- 0                # Cell outside the branch with mutation 
    
    ## Compute the Loglikihood of mutation i across branches
    Loglik_11 <- x_11*NonZero.logV[i,]
    Loglik_10 <- x_10*Zero.logV[i,]
    Loglik_01 <- x_01*NonZero.logP
    Loglik_00 <- x_00*Zero.logP
    Loglik.i <- colSums(Loglik_11+Loglik_10+Loglik_01+Loglik_00)
    
    return(Loglik.i)
    
}

# Get final probabilities of assignment
row.names(Loglik) <- row.names(mtr.bi)
n <- dim(mtr.bi)[1]
edge_ml <- apply(Loglik,1,which.max) %>% as.numeric

p <- exp(Loglik)# -Loglik[(as.numeric(edge_ml)-1)*n+1:n])
p <- p/rowSums(p)
p[is.na(p)] = min(p, na.rm = TRUE)

assign.to <- 'Parent'
if (assign.to=="Parent"){edge.node=1} else {edge.node=2}
var_assignment <- data.frame(
  temp=tree$edge[,edge.node][apply(p,1,which.max)], # NBB: final assignment of each variant to the CHILD NODE of each branch
  p=apply(p,1,max)                                   
)
var_assignment <- setNames(var_assignment, c(assign.to, 'p')) %>% arrange(desc(p))

var_assignment <- cbind(
  var_assignment,
  apply(
    var_assignment %>% rownames_to_column(var='variant'), 1, function(x) { 
      node <- as.numeric(x[assign.to])
      mut <- x['variant']
      cells <- tree$tip.label[Descendants(tree, node, 'tips')[[1]]]
      other <- setdiff(colnames(mtr.bi), cells)
      af.1 <- af.m[mut, cells] %>% mean
      p.1 <- sum(mtr.bi[mut, cells]) / length(cells)
      af.0 <- af.m[mut, other] %>% mean
      p.0 <- sum(mtr.bi[mut, other]) / length(other)
      return(c(af1=af.1, af0=af.0, p1=p.1, p0=p.0))
    }
  ) %>% 
    t %>% as.data.frame %>% mutate_all(~ifelse(is.na(.), 0, .))
)

var_assignment %>% filter(p>=0.3)


##
