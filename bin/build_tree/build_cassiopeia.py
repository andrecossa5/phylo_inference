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
    '--path_pickles', 
    type=str,
    default=None,
    help='Path to pickles main folder. Default: None.'
)

my_parser.add_argument(
    '--sample', 
    type=str,
    default=None,
    help='Sample name. Default: None.'
)

my_parser.add_argument(
    '--job_id', 
    type=str,
    default=None,
    help='Job id. Default: None.'
)

my_parser.add_argument(
    '--solver', 
    type=str,
    default='UPMGA',
    help='Cassiopeia solver. Default: UPMGA.'
)

my_parser.add_argument(
    '--metric', 
    type=str,
    default='jaccard',
    help='Distance metric. Default: jaccard.'
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


########################################################################

# Code
import pickle
from mito_utils.utils import *
from mito_utils.phylo import *

########################################################################

# Main
def main():

    # Parse arguments
    args = my_parser.parse_args()

    path_afm = args.afm
    path_pickles = args.path_pickles
    sample = args.sample
    job_id = args.job_id
    solver = args.solver
    metric = args.metric
    n_cores = args.n_cores
    boot_replicate = args.boot_replicate

    # Modify if job_id pickle from tuning, if present
    if path_pickles is not None and job_id is not None:

        path_pickle = os.path.join(path_pickles, sample, f'{job_id}_stats.pickle')

        if os.path.exists(path_pickle):
            with open(path_pickle, 'rb') as f:
                d = pickle.load(f)
            tree_kwargs = d['options']['tree_kwargs']
        
        else:
            raise ValueError(f'{path_pickle} does not exists!')
    
    else:
        tree_kwargs = {'solver' : solver, 'metric': metric}

    # Build tree
    afm = sc.read(path_afm)
    tree = build_tree(afm, precomputed=True, ncores=n_cores, **tree_kwargs)
    
    # Write out
    write_newick(tree, path=f'rep_{boot_replicate}.newick')
    

    ##


#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################