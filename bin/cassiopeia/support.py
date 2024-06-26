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

# Add arguments

# Filtering
my_parser.add_argument(
    '--filtering_key', 
    type=str,
    default=None,
    help='Filtering_key in filtering_oprions.json. Default: weng2024.'
)

# Solver
my_parser.add_argument(
    '--solver', 
    type=str,
    default='UPMGA',
    help='''
    Solver. Default: UPMGA. Other options: 
    NJ, spectral, max_cut, greedy. See Cassiopeia (MW Jones et al., 2020) 
    for docs.
    '''
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


##


# Parse arguments
args = my_parser.parse_args()

filtering_key = args.filtering_key
solver = args.solver
trees = args.trees
n_cores = args.n_cores

########################################################################

# Preparing run: import code, prepare directories

# Code
from mito_utils.preprocessing import *
from mito_utils.utils import *
from mito_utils.phylo import *

########################################################################

# Main
def main():

    # Read trees
    boot_trees = []
    path_list = [ x.replace('[', '').replace(']', '') for x in trees.split(', ') ]
    for path in path_list:
        k = path.split('/')[-1].split('.')[0]
        if k == 'observed':
            obs_tree = read_newick(path)
        else:
            boot_trees.append(read_newick(path))

    # Compute branch supports: Transfer Bootstrap Expectations (TBE)
    tree = calculate_supports(obs_tree, tree_list=boot_trees, method='TBE', n_jobs=n_cores)

    # Write final annotated tree
    write_newick(tree, path='final_tree.newick')


    ##


#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################


