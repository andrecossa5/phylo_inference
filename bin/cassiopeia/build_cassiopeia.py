#!/usr/bin/python

# Build cassiopeia

########################################################################
 
# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='build_cassiopeia',
    description=
    """
    Build cassiopeia trees.
    """
)

# Add arguments

# Path_main
my_parser.add_argument(
    '-p', 
    '--path_data', 
    type=str,
    default='..',
    help='Path to data folder. Default: .. .'
)

# Sample
my_parser.add_argument(
    '--sample', 
    type=str,
    default='sAML1',
    help='Sample name. Default: sAML1.'
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

# Solver
my_parser.add_argument(
    '--treshold_calling', 
    type=float,
    default=0.025,
    help='Treshold for variant allele calling. Default: .025.'
)


# Path meta
my_parser.add_argument(
    '--name', 
    type=str,
    default='tree',
    help='Tree name. Default: tree.'
)


# ncores
my_parser.add_argument(
    '--ncores', 
    type=int,
    default=4,
    help='N cores for pairwise distance calculation. Default: 4.'
)

# Path meta
my_parser.add_argument(
    '--path_meta', 
    type=str,
    default=None,
    help='Path to meta_cells .csv file. Default: None .'
)

# Path meta
my_parser.add_argument(
    '--path_priors', 
    type=str,
    default='None',
    help='Path to priors.csv file. Default: None .'
)

# Path meta
my_parser.add_argument(
    '--path_filtering', 
    type=str,
    default=None,
    help='Path to filtering_options.csv file. Default: None .'
)

# Sample
my_parser.add_argument(
    '--filtering_key', 
    type=str,
    default='weng2024',
    help='Filtering option in filtering_options.yml. Default: weng2024.'
)

# Sample
my_parser.add_argument(
    '--lineage_column', 
    type=str,
    default=None,
    help='Sample to use. Default: None.'
)

# Parse arguments
args = my_parser.parse_args()

sample = args.sample
path_data = args.path_data
metric = args.metric
solver = args.solver
ncores = args.ncores
name = args.name
path_priors = args.path_priors
path_filtering = args.path_filtering
path_meta = args.path_meta
filtering_key = args.filtering_key
lineage_column = args.lineage_column

########################################################################

# Preparing run: import code, prepare directories, set logger

# Code
import json
import numpy as np
from scipy.sparse import load_npz
from anndata import AnnData
from mito_utils.preprocessing import *
from mito_utils.utils import *
from mito_utils.phylo import *


##


########################################################################

# Main
def main():

    # Prep input

    # Reconstruct AFM
    AD = load_npz(os.path.join(path_data, 'AD.npz')).astype(np.int16).A
    DP = load_npz(os.path.join(path_data, 'DP.npz')).astype(np.int16).A
    if 'boot' in os.path.basename(path_data):
        cell_df = pd.read_csv(os.path.join(path_data, 'cells.csv'), index_col=0, header=None)
        var_df = pd.read_csv(os.path.join(path_data, 'variants.csv'), index_col=0, header=None)
    else:
        cell_df = pd.read_csv(os.path.join(path_data, 'meta.csv'), index_col=0)
        var_df = pd.read_csv(os.path.join(path_data, 'var_meta.csv'), index_col=0)
    a = AnnData(X=np.divide(AD, DP), obs=cell_df, var=var_df, dtype=np.float16)
    
    # Read af_confident_detection for binarization
    _, filtering_kwargs, _ = process_json(path_filtering, filtering_key)
    t = filtering_kwargs['af_confident_detection'] if 'af_confident_detection' in filtering_kwargs else .05
    
    # Compute variants weights, if necessary
    if os.path.exists(path_priors):
        priors = pd.read_csv(path_priors, index_col=0)
        priors = priors.loc[a.var_names, priors.columns!=sample].mean(axis=1)
        weights = 1-priors.values
    else:
        weights = None

    # Build tree
    tree = build_tree(a, solver=solver, metric=metric, weights=weights, t=t, ncores=ncores)

    # Write tree
    write_newick(tree, path=f'{name}.newick')
    

    ##


#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################