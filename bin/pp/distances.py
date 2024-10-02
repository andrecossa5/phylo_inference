#!/usr/bin/python

# prep_MAESTER script

########################################################################

# Libraries 
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='distances',
    description=
    """
    Prepare input for tree building: character/distance matrices and sequences (.fasta) file.
    """
)

# Add arguments
my_parser.add_argument(
    '--AD', 
    type=str,
    default='..',
    help='Path to AD.npz file, a table of AD counts. Default: .. .'
)

my_parser.add_argument(
    '--DP', 
    type=str,
    default=None,
    help='Path to DP.npz file, a table of DP counts. Default: .. .'
)

my_parser.add_argument(
    '--path_bin', 
    type=str,
    default=None,
    help='Path to bin_ops.json file. Default: None.'
)

my_parser.add_argument(
    '--path_tree', 
    type=str,
    default=None,
    help='Path to tree_ops.json file. Default: None.'
)

my_parser.add_argument(
    '--bin_key', 
    type=str,
    default='default',
    help='Filtering option in bin_ops.json. Default: default.'
)

my_parser.add_argument(
    '--tree_key', 
    type=str,
    default='default',
    help='Filtering option in tree_ops.json. Default: default.'
)

my_parser.add_argument(
    '--n_cores', 
    type=int,
    default='8',
    help='n cores to use. Default: 8.'
)

my_parser.add_argument(
    '--boot_strategy', 
    type=str,
    default='feature_resampling',
    help='Bootstrap method. Default: feature_resampling.'
)

my_parser.add_argument(
    '--boot_replicate', 
    type=str,
    default='original',
    help='Name of the boot replicate. Default: original.'
)

my_parser.add_argument(
    '--frac_char_resampling', 
    type=float,
    default=.8,
    help='Fraction of characters to resample at each bootstrapping iteration. Default: .8.'
)


##


# Parse arguments
args = my_parser.parse_args()

path_AD = args.AD
path_DP = args.DP
path_bin = args.path_bin
path_tree = args.path_tree
bin_key = args.bin_key
tree_key = args.tree_key
n_cores = args.n_cores
boot_strategy = args.boot_strategy
boot_replicate = args.boot_replicate
frac_char_resampling = args.frac_char_resampling


##


########################################################################

# Preparing run: import code, prepare directories

# Code
from scipy.sparse import save_npz, load_npz, csr_matrix
from mito_utils.utils import *
from mito_utils.phylo import *

########################################################################

# Main
def main():

    # Parse kwargs options
    bin_method, bin_kwargs = process_bin_kwargs(path_bin, bin_key)
    tree_kwargs = process_kwargs(path_tree, tree_key)

    # Bootstrap
    AD_original = load_npz(path_AD).A.T.astype(np.int16)
    DP_original = load_npz(path_DP).A.T.astype(np.int16)

    if boot_replicate != 'observed':
        if boot_strategy == 'jacknife':
            AD, DP, _ = jackknife_allele_tables(AD_original, DP_original)
        elif boot_strategy == 'counts_resampling':
            AD, DP, _ = bootstrap_allele_counts(AD_original, DP_original)
        elif boot_strategy == 'feature_resampling':
            AD, DP, _ = bootstrap_allele_tables(AD_original, DP_original, frac_resampled=frac_char_resampling)
        else:
            raise ValueError(f'{boot_strategy} boot_strategy is not supported...')
    else:
        AD = AD_original
        DP = DP_original
    
    # Calculate distances
    D = compute_distances(AD=AD, DP=DP, metric=tree_kwargs['metric'], bin_method=bin_method, binarization_kwargs=bin_kwargs)

    # Save dist matrix 
    save_npz('dist.npz', coo_matrix(D))


    ##


#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################