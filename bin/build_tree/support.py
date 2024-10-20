#!/usr/bin/python

# TBE script

########################################################################
 
# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='support',
    description=
    """
    Evaluate internal nodes support.
    """
)

# obs_tree
my_parser.add_argument(
    '--trees', 
    type=str,
    default=None,
    help='Path to observed and bootstrapped trees. Default: None'
)

# obs_tree
my_parser.add_argument(
    '--n_cores', 
    type=int,
    default=8,
    help='n cores for the process. Default: 8'
)

my_parser.add_argument(
    '--support_method', 
    type=str,
    default='TBE',
    help='Support method. Default: TBE'
)


##


# Parse arguments
args = my_parser.parse_args()
path_list = args.trees.strip('[|]').split(', ')
support_method = args.support_method
n_cores = args.n_cores


##


########################################################################

# Preparing run: import code, prepare directories

# Code
from mito_utils.utils import *
from mito_utils.phylo import *

##

########################################################################

# Main
def main():

    # Read trees
    boot_trees = []
    for path in path_list:
        k = path.split('/')[-1].split('.')[0]
        if k == 'observed':
            obs_tree = read_newick(path)
        else:
            boot_trees.append(read_newick(path))

    # Add internal node support
    tree = calculate_supports(obs_tree, boot_trees, method=support_method, n_jobs=n_cores)

    # Write final annotated tree
    write_newick(tree, path='final_tree.newick')


    ##


#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################