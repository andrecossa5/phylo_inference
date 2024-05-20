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

# Parse arguments
args = my_parser.parse_args()

path_i = args.path_data
method = args.method
path_filtering = args.path_filtering
filtering_key = args.filtering_key


##


########################################################################

# Preparing run: import code, prepare directories

# Code
import json
from scipy.sparse import load_npz, save_npz, csr_matrix
from anndata import AnnData
from mito_utils.phylo import *

########################################################################

# Main
def main():

    # AD and DP
    AD_original = load_npz(os.path.join(path_i, 'AD.npz')).astype(np.int16)
    DP_original = load_npz(os.path.join(path_i, 'DP.npz')).astype(np.int16)

    if method == 'jacknife':
        AD_boot, DP_boot, sel_idx = jackknife_allele_tables(AD_original.A, DP_original.A)
    elif method == 'counts_resampling':
        AD_boot, DP_boot, sel_idx = bootstrap_allele_counts(AD_original.A, DP_original.A)
    elif method == 'feature_resampling':
        AD_boot, DP_boot, sel_idx = bootstrap_allele_tables(AD_original.A, DP_original.A)

    # Cells and variants
    cells = pd.read_csv(os.path.join(path_i, 'meta.csv'), index_col=0).index
    variants = pd.read_csv(os.path.join(path_i, 'variants.csv'), index_col=0, header=None).index
    variants = variants[sel_idx]

    # Save bootstrapped AD and DP
    os.mkdir('bootstrapped_input')
    os.chdir('bootstrapped_input')  
    cells.to_series().to_csv('cells.csv', index=False, header=None)
    variants.to_series().to_csv('variants.csv', index=False, header=None)
    save_npz('AD_boot.npz', csr_matrix(AD_boot.astype(np.float16)))
    save_npz('DP_boot.npz', csr_matrix(DP_boot.astype(np.float16)))

    # Prep fasta
    afm = AnnData(
        X=np.divide(AD_boot, DP_boot), 
        var=pd.DataFrame(index=variants), 
        dtype=np.float16
    )

    # Read t, if available
    with open(path_filtering, 'r') as file:
        FILTERING_OPTIONS = json.load(file)

    if filtering_key in FILTERING_OPTIONS:
        d = FILTERING_OPTIONS[filtering_key]
        t = d['t'] if 't' in d else .05
    else:
        raise KeyError(f'{filtering_key} not in {path_filtering}!')
    
    # Covert to .fasta and write
    seqs = AFM_to_seqs(afm, t=t)
    with open('sequences.fasta', 'w') as f:
        for k in seqs:
            f.write(f'>seq_{k}\n')
            f.write(f'{seqs[k]}\n')

    ##
    
#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################