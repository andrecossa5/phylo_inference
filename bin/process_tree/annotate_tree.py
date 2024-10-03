#!/usr/bin/python

# Build cassiopeia

########################################################################
 
# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='annotate_tree',
    description="Process and annotate tree."
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
    '--X_bin', 
    type=str,
    default=None,
    help='Path to X_bin.npz file, storing the binarized character matrix. Default: .. .'
)

my_parser.add_argument(
    '--dists', 
    type=str,
    default=None,
    help='Path to dist.npz file, storing pairwise cell-cell distances. Default: .. .'
)

my_parser.add_argument(
    '--tree', 
    type=str,
    default=None,
    help='Path to tree in .newick format. Default: .. .'
)


##


# Parse arguments
args = my_parser.parse_args()
path_tree = args.tree
path_AD = args.AD
path_DP = args.DP
path_X_bin = args.X_bin
path_cell_meta = args.cell_meta
path_char_meta = args.char_meta
path_dists = args.dists


##


########################################################################

# Preparing run: import code, prepare directories, set logger

# Code
import pickle
from scipy.sparse import load_npz
from mito_utils.utils import *
from mito_utils.phylo import *


##


########################################################################

# Main
def main():

    # Prep annot
    AD = load_npz(path_AD).A.T.astype(np.int16)
    DP = load_npz(path_DP).A.T.astype(np.int16)
    X_bin = load_npz(path_X_bin).A.T.astype(np.int16)
    D = load_npz(path_dists).A
    cell_meta = pd.read_csv(path_cell_meta, index_col=0)
    char_meta = pd.read_csv(path_char_meta, index_col=0)
    X_raw = pd.DataFrame(np.divide(AD, DP+.0000001), index=cell_meta.index, columns=char_meta.index)
    X_bin = pd.DataFrame(X_bin, index=cell_meta.index, columns=char_meta.index)
    D = pd.DataFrame(D, index=cell_meta.index, columns=cell_meta.index)

    # Load tree
    tree = read_newick(path_tree, X_raw=X_raw, X_bin=X_bin, D=D, meta=cell_meta)

    # Cut and annotate tree
    tree, _, _ = cut_and_annotate_tree(tree)

    # Write as pickle
    with open('annotated_tree.pickle', 'wb') as f:
        pickle.dump(tree, f)
    

    ##


#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################