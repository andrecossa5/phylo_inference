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
    '--afm', 
    type=str,
    default='..',
    help='Path to afm.h5ad file. Default: .. .'
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

path_afm = args.afm
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
from mito_utils.utils import *
from mito_utils.phylo import *


##


########################################################################

# Main
def main():

    # Read afm
    afm = sc.read(path_afm)

    # Tree
    tree_kwargs = process_kwargs(path_tree, tree_key)
    bin_method, bin_kwargs = process_bin_kwargs(path_bin, bin_key)

    tree = build_tree(afm, precomputed=True, ncores=n_cores, 
                      bin_method=bin_method, binarization_kwargs=bin_kwargs, **tree_kwargs)
    
    write_newick(tree, path=f'{boot_replicate}.newick')
    

    ##


#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################