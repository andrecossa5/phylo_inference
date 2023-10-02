#!/usr/bin/python

# Evaluate III

########################################################################

# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='evaluate_III',
    description=
    """
    Evaluate trees external associations.
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

# input folder
my_parser.add_argument(
    '--input_folder', 
    type=str,
    default=None,
    help='Path to input_folder.'
)

# Parse arguments
args = my_parser.parse_args()

sample_name = args.sample_name
filtering = args.filtering
solver = args.solver
metric = args.metric if solver in ['UPMGA', 'NJ', 'spectral'] else 'cosine'
obs_tree = args.obs_tree
input_folder = args.input_folder

# sample_name = 'MDA_clones'
# filtering = 'miller2022'
# solver = 'NJ'
# metric = 'cosine'
# obs_tree = '/Users/IEO5505/Desktop/MI_TO/phylo_inference/work/ee/a57ede09b663d5a5f49364edae9db0/tree.newick'
# input_folder = '/Users/IEO5505/Desktop/MI_TO/phylo_inference/work/26/1a2298fb3097e777de326f82478e82/input_folder'


##

########################################################################

# Preparing run: import code, prepare directories

# Code
from scipy.sparse import load_npz
from mito_utils.preprocessing import *
from mito_utils.utils import *
from mito_utils.phylo import *

########################################################################

# Main
def main():

    # Read tree newick and AFM
    with open(obs_tree, 'r') as f:
        tree = f.readlines()[0]

    # Load meta and muts
    meta = pd.read_csv(os.path.join(input_folder, 'meta.csv'), index_col=0)
    variants_df = pd.read_csv(os.path.join(input_folder, 'variants.csv'), header=None, index_col=0)

    # Get AFM
    AD = load_npz(os.path.join(input_folder, 'AD.npz'))
    DP = load_npz(os.path.join(input_folder, 'DP.npz'))
    afm = AnnData(np.divide(AD.A, DP.A), obs=meta, var=variants_df, dtype=np.float16)
    afm = nans_as_zeros(afm)
    meta = meta.join(pd.DataFrame(afm.X, columns=variants_df.index, index=meta.index))

    if metric is None:
        m = 'cosine'
    else:
        m = metric
        
    D = pair_d(a=afm.X, metric=m)
    D[np.ix_(range(D.shape[0]), range(D.shape[0]))] = 0
    D = pd.DataFrame(D, index=meta.index, columns=meta.index)

    # Add to tree
    tree = CassiopeiaTree(tree=tree, cell_meta=meta, dissimilarity_map=D)

    # Calc moran'I for each variable
    morans_muts = cs.tl.compute_morans_i(tree, variants_df.index)
    couplings = cs.tl.compute_evolutionary_coupling(tree, meta_variable='GBC')

    # Save
    morans_muts.to_csv('morans_muts.csv')
    couplings.to_csv('couplings.csv')


    ##


#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################