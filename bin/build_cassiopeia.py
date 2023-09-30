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

# Metric
my_parser.add_argument(
    '--metric', 
    type=str,
    default=None,
    help='Metric used for cell-cell distance matrix computation. Default: None.'
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
    '--name', 
    type=str,
    default='tree',
    help='Name for the saved newick object. Default: tree.'
)

# Solver
my_parser.add_argument(
    '--treshold_calling', 
    type=float,
    default=0.025,
    help='Treshold for variant allele calling. Default: .025.'
)

# ncores
my_parser.add_argument(
    '--ncores', 
    type=int,
    default=4,
    help='N cores for pairwise distance calculation. Default: 4.'
)

# Parse arguments
args = my_parser.parse_args()

path_data = args.path_data
metric = args.metric
solver = args.solver
name = args.name
t = args.treshold_calling
ncores = args.ncores

########################################################################

# Preparing run: import code, prepare directories, set logger

# Code
from scipy.sparse import load_npz
from anndata import AnnData
from mito_utils.preprocessing import *
from mito_utils.utils import *
from mito_utils.phylo import *

########################################################################

# Main
def main():

    # Read input 
    cells = pd.read_csv(os.path.join(path_data, 'cells.csv'), index_col=0, header=None)
    variants = pd.read_csv(os.path.join(path_data, 'variants.csv'), index_col=0, header=None)
    AD = load_npz(os.path.join(path_data, 'AD_boot.npz')).A
    DP = load_npz(os.path.join(path_data, 'DP_boot.npz')).A
    afm = AnnData(np.divide(AD, DP), obs=cells, var=variants, dtype=np.float16)
    afm = nans_as_zeros(afm)

    # Build tree
    tree = build_tree(afm, metric=metric, solver=solver, t=t, ncores=ncores)

    # Save tree in .netwick format
    with open(f'{name}.newick', 'w') as f:
        f.write(f'{tree.get_newick()}')
    

    ##


#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################




