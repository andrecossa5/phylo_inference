#!/usr/bin/python

# Bootstrap script

########################################################################

# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='bootstrap',
    description=
    """
    Prepare bootstrapped input for tree building: character and/or distance matrix + sequences .fasta.
    """
)

# Add arguments

# Path_main
my_parser.add_argument(
    '-p', 
    '--path_data', 
    type=str,
    default='..',
    help='Path to AD, DP folder. Default: .. .'
)

# Method
my_parser.add_argument(
    '--method', 
    type=str,
    default='jacknife',
    help='Bootstarpping strategy. Default: jacknife.'
)

# Method
my_parser.add_argument(
    '--feature_resampling_perc', 
    type=float,
    default=.9,
    help='Percentage of resampled features if method: feature_resampling. Default: feature_resampling.'
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
    '--boot_replicate', 
    type=str,
    default='observed',
    help='Name of the bootstrapped input. Default: observed.'
)

# Parse arguments
args = my_parser.parse_args()

path_i = args.path_data
method = args.method
feature_resampling_perc = args.feature_resampling_perc
path_filtering = args.path_filtering
filtering_key = args.filtering_key
boot_replicate = args.boot_replicate


##


########################################################################

# Preparing run: import code, prepare directories

# Code
import json
from scipy.sparse import load_npz, save_npz, csr_matrix
from anndata import AnnData
from mito_utils.phylo import *
from mito_utils.utils import *

########################################################################

# Main
def main():

    # AD and DP
    AD_original = load_npz(os.path.join(path_i, 'AD.npz')).astype(np.int16)
    DP_original = load_npz(os.path.join(path_i, 'DP.npz')).astype(np.int16)

    # Perturb AD and DP, if necessary
    if 'boot_replicate' != 'observed':
        if method == 'jacknife':
            AD_boot, DP_boot, sel_idx = jackknife_allele_tables(AD_original.A, DP_original.A)
        elif method == 'counts_resampling':
            AD_boot, DP_boot, sel_idx = bootstrap_allele_counts(AD_original.A, DP_original.A)
        elif method == 'feature_resampling':
            AD_boot, DP_boot, sel_idx = bootstrap_allele_tables(AD_original.A, DP_original.A, frac_resampled=feature_resampling_perc)
        else:
            raise ValueError(f'{method} is not supported...')
    else:
        AD_boot = AD_original.A
        DP_boot = DP_original.A
        sel_idx = [ x for x in np.arange(AD_boot.shape[1]) ]

    # Cells and variants
    cells = pd.read_csv(os.path.join(path_i, 'meta.csv'), index_col=0).index
    variants = pd.read_csv(os.path.join(path_i, 'variants.csv'), index_col=0, header=None).index
    variants = variants[sel_idx]

    # Save bootstrapped AD and DP
    os.mkdir('bootstrapped_input')
    os.chdir('bootstrapped_input')  
    cells.to_series().to_csv('cells.csv', index=False, header=None)
    variants.to_series().to_csv('variants.csv', index=False, header=None)
    save_npz('AD.npz', csr_matrix(AD_boot.astype(np.float16)))
    save_npz('DP.npz', csr_matrix(DP_boot.astype(np.float16)))

    # Prep fasta
    afm = AnnData(
        X=np.divide(AD_boot, DP_boot), 
        var=pd.DataFrame(index=variants), 
        dtype=np.float16
    )

    # Read t, if available
    _, filtering_kwargs, _ = process_json(path_filtering, filtering_key)
    t = filtering_kwargs['af_confident_detection'] if 'af_confident_detection' in filtering_kwargs else .05
    
    # Covert to .fasta and write
    seqs = AFM_to_seqs(afm, t=t)
    with open('sequences.fasta', 'w') as f:
        for k in seqs:
            f.write(f'>{k}\n')
            f.write(f'{seqs[k]}\n')

    ##
    
#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################