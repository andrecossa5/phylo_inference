#!/usr/bin/python

# Calculate distances

########################################################################

# Libraries
import os 
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='calculate_distances',
    description=
    """
    Calc distances among cells AFs.
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
    default='cosine',
    help='Metric used for cell-cell distance matrix computation. Default: cosine.'
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
metric = args.metric if args.metric in ['UPMGA', 'NJ', 'spectral'] else 'cosine'
t = args.treshold_calling
ncores = args.ncores

########################################################################

# Preparing run: import code, prepare directories, set logger

# Code
import re
from scipy.sparse import load_npz, save_npz, csr_matrix
from anndata import AnnData
from mito_utils.preprocessing import *
from mito_utils.utils import *
from mito_utils.phylo import *

########################################################################

# Main
def main():

    # Read input 
    is_boot = any([ True for x in os.listdir(path_data) if bool(re.search('boot', x)) ])
    if is_boot:
        AD = load_npz(os.path.join(path_data, 'AD_boot.npz')).astype(np.int16).A
        DP = load_npz(os.path.join(path_data, 'DP_boot.npz')).astype(np.int16).A
        cells = pd.read_csv(
            os.path.join(path_data, 'cells.csv'), index_col=0, header=None
        )
        variants = pd.read_csv(
            os.path.join(path_data, 'variants.csv'), index_col=0, header=None
        )
    else:
        AD = load_npz(os.path.join(path_data, 'AD.npz')).astype(np.int16).A
        DP = load_npz(os.path.join(path_data, 'DP.npz')).astype(np.int16).A
        cells = pd.read_csv(
            os.path.join(path_data, 'meta.csv'), index_col=0
        )
        variants = pd.read_csv(
            os.path.join(path_data, 'variants.csv'), index_col=0, header=None
        )

    # Calc distance matrix
    afm = AnnData(np.divide(AD, DP), obs=cells, var=variants, dtype=np.float16)
    X = nans_as_zeros(afm).X.astype(np.float16)
    m = 'cosine' if metric is None else metric
    D = pair_d(X, ncores=ncores, metric=metric)
    D[np.isnan(D)] = 0

    # Save
    save_npz('dist.npz', csr_matrix(D))


    ##


#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################




