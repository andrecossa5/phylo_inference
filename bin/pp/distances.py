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
    '--afm', 
    type=str,
    default='..',
    help='Path to preprocessed afm.h5ad file. Default: .. .'
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

path_afm = args.afm
n_cores = args.n_cores
boot_strategy = args.boot_strategy
boot_replicate = args.boot_replicate
frac_char_resampling = args.frac_char_resampling


##


########################################################################

# Preparing run: import code, prepare directories

# Code
from scipy.sparse import save_npz, csr_matrix
from anndata import AnnData
from mito_utils.utils import *
from mito_utils.phylo import *

########################################################################

# Main
def main():

    # Load afm and extract AD and DP and options
    afm = sc.read(path_afm)
    metric = afm.uns['distance_calculations']['distances']['metric']
    bin_method = afm.uns['genotyping']['bin_method']
    bin_kwargs = { k:v for k,v in afm.uns['genotyping'].items() if k!='bin_method' }

    # Bootstrap
    AD_original = afm.layers['AD'].A
    DP_original = afm.layers['DP'].A

    if boot_replicate != 'observed':
        if boot_strategy == 'jacknife':
            AD, DP, idx = jackknife_allele_tables(AD_original, DP_original)
        elif boot_strategy == 'counts_resampling':
            AD, DP, idx = bootstrap_allele_counts(AD_original, DP_original)
        elif boot_strategy == 'feature_resampling':
            AD, DP, idx = bootstrap_allele_tables(AD_original, DP_original, frac_resampled=frac_char_resampling)
        else:
            raise ValueError(f'{boot_strategy} boot_strategy is not supported...')
    else:
        AD = AD_original
        DP = DP_original
    
    # Calculate distances
    AF = csr_matrix(np.divide(AD, (DP+.0000001)))
    AD = csr_matrix(AD)
    DP = csr_matrix(DP)
    afm_new = AnnData(X=AF, obs=afm.obs, var=afm.var.iloc[idx,:], uns=afm.uns, layers={'AD':AD, 'DP':DP})
    compute_distances(
        afm_new, 
        metric=metric, 
        bin_method=bin_method, 
        binarization_kwargs=bin_kwargs, 
        ncores=n_cores
    )

    # Save dist matrix 
    afm_new.write('afm_new.h5ad')


    ##


#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################