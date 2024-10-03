#!/usr/bin/python

# Build cassiopeia

########################################################################
 
# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='build_cassiopeia',
    description="Build cassiopeia trees"
)

# Add arguments
my_parser.add_argument(
    '--AD', 
    type=str,
    default='..',
    help='Path to AD.npz file, a table of AD counts. Default: .. .'
)

my_parser.add_argument(
    '--DP', 
    type=str,
    default=None,
    help='Path to DP.npz file, a table of DP counts. Default: .. .'
)

# Add arguments
my_parser.add_argument(
    '--cell_meta', 
    type=str,
    default='..',
    help='Path to cell_meta.csv file. Default: .. .'
)

my_parser.add_argument(
    '--char_meta', 
    type=str,
    default=None,
    help='Path to char_meta.csv file. Default: .. .'
)

my_parser.add_argument(
    '--dists', 
    type=str,
    default=None,
    help='Path to dist.npz file, storing pairwise cell-cell distances. Default: .. .'
)

my_parser.add_argument(
    '--path_bin', 
    type=str,
    default=None,
    help='Path to bin_ops.json file. Default: None.'
)

my_parser.add_argument(
    '--path_tree', 
    type=str,
    default=None,
    help='Path to tree_ops.json file. Default: None.'
)

my_parser.add_argument(
    '--bin_key', 
    type=str,
    default='default',
    help='Filtering option in bin_ops.json. Default: default.'
)

my_parser.add_argument(
    '--tree_key', 
    type=str,
    default='default',
    help='Filtering option in tree_ops.json. Default: default.'
)

my_parser.add_argument(
    '--n_cores', 
    type=int,
    default='8',
    help='n cores to use. Default: 8.'
)

my_parser.add_argument(
    '--boot_replicate', 
    type=str,
    default='original',
    help='Name of the boot replicate. Default: original.'
)


##


# Parse arguments
args = my_parser.parse_args()

path_AD = args.AD
path_DP = args.DP
path_cell_meta = args.cell_meta
path_char_meta = args.char_meta
path_dists = args.dists
path_bin = args.path_bin
path_tree = args.path_tree
bin_key = args.bin_key
tree_key = args.tree_key
n_cores = args.n_cores
boot_replicate = args.boot_replicate


##


########################################################################

# Preparing run: import code, prepare directories, set logger

# Code
from scipy.sparse import load_npz
from mito_utils.utils import *
from mito_utils.phylo import *


##


########################################################################

# Main
def main():

    # Prep input
    AD = load_npz(path_AD).A.T.astype(np.int16)
    DP = load_npz(path_DP).A.T.astype(np.int16)
    D = load_npz(path_dists).A
    cell_meta = pd.read_csv(path_cell_meta, index_col=0)
    char_meta = pd.read_csv(path_char_meta, index_col=0)

    # Params
    bin_method, bin_kwargs = process_bin_kwargs(path_bin, bin_key)
    tree_kwargs = process_kwargs(path_tree, tree_key)

    # Build tree
    tree = build_tree(AD=AD, DP=DP, D=D, meta=cell_meta, var_names=char_meta.index, 
                      bin_method=bin_method, binarization_kwargs=bin_kwargs, **tree_kwargs)

    # Write tree
    write_newick(tree, path=f'{boot_replicate}.newick')
    

    ##


#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################