#!/usr/bin/python

# Evaluate I

########################################################################

# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='evaluate_II',
    description=
    """
    Evaluate trees consistency
    """
)

# Add arguments

# Sample name
my_parser.add_argument(
    '--sample_name', 
    type=str,
    default=None,
    help='Sample name. Default: None.'
)

# Filtering
my_parser.add_argument(
    '--filtering', 
    type=str,
    default=None,
    help='Filtering method. Default: miller2022.'
)

# Metric
my_parser.add_argument(
    '--metric', 
    type=str,
    default='cosine',
    help='Metric used for cell-cell distance matrix computation. Default: cosine.'
)

# Solver
my_parser.add_argument(
    '--solver', 
    type=str,
    default='UPMGA',
    help='''
    Tree building algorithm. Default: UPMGA. Other options: 
    NJ, spectral, max_cut, greedy. See Cassiopeia (MW Jones et al., 2020) 
    for docs.
    '''
)

# obs_tree
my_parser.add_argument(
    '--obs_tree', 
    type=str,
    default=None,
    help='Path to observed tree.'
)

# Parse arguments
args = my_parser.parse_args()

sample_name = args.sample_name
filtering = args.filtering
metric_l = args.metric
solver_l = args.solver
obs_tree_l = args.obs_tree

# sample_name = 'MDA_clones'
# filtering = 'miller2022'
# solver = '[NJ, max_cut, NJ]'
# metric = '[cosine, None, hamming]'
# obs_tree = '[/Users/IEO5505/Desktop/MI_TO/phylo_inference/work/87/25c31b923ee110ab35b3ab696ead5f/tree.newick, /Users/IEO5505/Desktop/MI_TO/phylo_inference/work/97/c8e7017edff8dbfcbdef1be037cbee/tree.newick, /Users/IEO5505/Desktop/MI_TO/phylo_inference/work/5c/2643a7e24213fb261bc004bcfff75c/tree.newick]]'

##

########################################################################

# Preparing run: import code, prepare directories

# Code
from itertools import combinations, chain
from mito_utils.preprocessing import *
from mito_utils.utils import *
from mito_utils.phylo import *

########################################################################

# Main
def main():

    # Read trees
    TREE_D = {}

    # Boot trees
    path_list = [ x.replace('[', '').replace(']', '') for x in obs_tree_l.split(', ') ]
    solvers = [ x.replace('[', '').replace(']', '') for x in solver_l.split(', ') ]
    metrics = [ x.replace('[', '').replace(']', '') for x in metric_l.split(', ') ]

    for path, solver, metric in zip(path_list, solvers, metrics):
        with open(path, 'r') as f:
            tree = f.readlines()[0]
            TREE_D[f'{solver}_{metric}'] = CassiopeiaTree(tree=tree)

    # Compute median RF distance
    L = []
    for x, y in combinations([ k for k in TREE_D ], 2):
        d, max_d = cs.critique.robinson_foulds(TREE_D[y], TREE_D[x])
        L.append(d/max_d)
    
    # Save
    s = ','.join([sample_name, filtering, str(np.median(L))])
    with open('results.txt', 'w') as f:
        f.write(s)

    ##


#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################