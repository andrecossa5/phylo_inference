#!/usr/bin/python

# Evaluate I

########################################################################

# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='evaluate_I',
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
    Bootstrap resampling strategy. Default: UPMGA. Other options: 
    NJ, spectral, max_cut, greedy. See Cassiopeia (MW Jones et al., 2020) 
    for docs.
    '''
)

# boot_trees
my_parser.add_argument(
    '--boot_trees', 
    type=str,
    default=None,
    help='String of paths to bootstrap trees to be collected and parsed.'
)

# boot_trees
my_parser.add_argument(
    '--boot_method', 
    type=str,
    default=None,
    help='Bootstrapping method used.'
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
metric = args.metric
solver = args.solver
boot_method = args.boot_method
boot_trees_s = args.boot_trees
obs_tree = args.obs_tree

# sample_name = 'MDA_clones'
# filtering = 'miller2022'
# solver = 'max_cut'
# metric = None 
# boot_method = 'counts_resampling'
# boot_trees_s = '[/Users/IEO5505/Desktop/MI_TO/phylo_inference/work/b2/bd1444eb504cf0166dd4f2a0770c6e/rep8.newick, /Users/IEO5505/Desktop/MI_TO/phylo_inference/work/7f/13c20ba54383545085a0df04ccf906/rep6.newick, /Users/IEO5505/Desktop/MI_TO/phylo_inference/work/26/a1aa3c18ee5c790a485d83b18164a8/rep3.newick, /Users/IEO5505/Desktop/MI_TO/phylo_inference/work/cc/8e3181458e3cfc925f5d5a1604d493/rep4.newick, /Users/IEO5505/Desktop/MI_TO/phylo_inference/work/47/df26e5767f218730edb69d132d27d0/rep2.newick, /Users/IEO5505/Desktop/MI_TO/phylo_inference/work/ef/a65d1de1b8ddbb5ea61d0014276d0f/rep7.newick, /Users/IEO5505/Desktop/MI_TO/phylo_inference/work/90/63e99791575d2551b9150bf6f20bbf/rep5.newick, /Users/IEO5505/Desktop/MI_TO/phylo_inference/work/65/3e2c684fb2c5f813d1ca291f43b1a3/rep1.newick, /Users/IEO5505/Desktop/MI_TO/phylo_inference/work/6a/cea48666c85eec729f3434e2ce544f/rep9.newick, /Users/IEO5505/Desktop/MI_TO/phylo_inference/work/3e/039826b57f5d2b35af065e5e4311df/rep10.newick]'
# obs_tree = '/Users/IEO5505/Desktop/MI_TO/phylo_inference/work/7c/8648de1abd8604a943361486b74a5b/tree.newick'
# ncores = 4

########################################################################

# Preparing run: import code, prepare directories

# Code
import pickle
from mito_utils.preprocessing import *
from mito_utils.utils import *
from mito_utils.phylo import *

########################################################################

# Main
def main():

    # Read trees
    TREE_D = {}

    # Boot trees
    path_list = [ x.replace('[', '').replace(']', '') for x in boot_trees_s.split(', ') ]
    for x in path_list:
        k = x.split('/')[-1].split('.')[0]
        with open(x, 'r') as f:
            tree = f.readlines()[0]
            TREE_D[k] = CassiopeiaTree(tree=tree)
    
    # Observed tree
    with open(obs_tree, 'r') as f:
        tree = f.readlines()[0]
    TREE_D['observed'] = CassiopeiaTree(tree=tree)

    # Save trees
    with open('trees.pickle', 'wb') as f:
        pickle.dump(TREE_D, f)
    
    # Compute supports
    df, _ = calculate_support(
        TREE_D['observed'], 
        tree_list=[ TREE_D[k] for k in TREE_D if k != 'observed'], 
        n=len(TREE_D)-1
    )

    # Compute median RF distance
    L = []
    for x in [ k for k in TREE_D if k != 'observed']:
        d, max_d = cs.critique.robinson_foulds(TREE_D['observed'], TREE_D[x])
        L.append(d/max_d)
    
    # Add specifics of the run
    df['median_RF'] = np.median(L)
    df['sample'] = sample_name
    df['filtering'] = filtering
    df['solver'] = solver
    df['metric'] = metric
    df['bootstrap'] = boot_method

    # Save
    cats = ['sample', 'filtering', 'bootstrap', 'solver', 'metric']
    conts = ['support', 'time', 'n_cells', 'median_RF']
    df[cats+conts].to_csv('extended_supports.csv')


    ##

#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################