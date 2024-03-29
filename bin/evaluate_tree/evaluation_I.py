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

# boot_trees
my_parser.add_argument(
    '--boot_support', 
    type=str,
    default=None,
    help='Support type.'
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
boot_support =  args.boot_support
boot_method = args.boot_method
boot_trees_s = args.boot_trees
obs_tree = args.obs_tree

# sample_name = 'MDA_clones'
# filtering = 'GT'
# solver = 'max_cut'
# metric = None 
# boot_method = 'counts_resampling'
# boot_trees_s = "[/Users/IEO5505/Desktop/MI_TO/phylo_inference/work/6f/3a59869760daffc096049e4e132c7d/rep5.newick, /Users/IEO5505/Desktop/MI_TO/phylo_inference/work/a3/c38059dd1d4bb74313b4914e0a3ecf/rep1.newick, /Users/IEO5505/Desktop/MI_TO/phylo_inference/work/f3/60e7de824875a075223800cd61faed/rep3.newick, /Users/IEO5505/Desktop/MI_TO/phylo_inference/work/f7/4934b2f9a1f26a3a97856972373318/rep2.newick, /Users/IEO5505/Desktop/MI_TO/phylo_inference/work/94/f2c630bc8660de02deb6a063f9e397/rep4.newick]"
# obs_tree = '/Users/IEO5505/Desktop/MI_TO/phylo_inference/work/85/f5fe9702f0864664ab3f066d0d4a56/tree.newick'
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
        
    # TBE
    tbe = calculate_supports(
        TREE_D['observed'], 
        tree_list=[ TREE_D[k] for k in TREE_D if k != 'observed'], 
        method='tbe_II', 
        n_jobs=8
    )
    # FBP
    fbp = calculate_supports(
        TREE_D['observed'], 
        tree_list=[ TREE_D[k] for k in TREE_D if k != 'observed'], 
        method='fbp', 
        n_jobs=8
    )
    # Merge
    df = (
        tbe.rename(columns={'support':'TBE'})
        .join(fbp[['support']].rename(columns={'support':'FBP'}))
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
    conts = ['time', 'n_cells', 'median_RF', 'TBE', 'FBP']
    df[cats+conts].to_csv('extended_supports.csv')


    ##

#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################