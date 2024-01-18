#!/usr/bin/python

# Prep input script

########################################################################

# Libraries 
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='prep_input',
    description=
    """
    Prepare input for tree building: character and/or distance matrix, sequences .fasta.
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
    default='MDA_clones',
    help='Sample to use. Default: MDA_clones.'
)

# Filter
my_parser.add_argument(
    '--filtering', 
    type=str,
    default='ludwig2019',
    help='Method to filter MT-SNVs. Default: ludwig2019.'
)

# ncores
my_parser.add_argument(
    '--ncores', 
    type=int,
    default=8,
    help='ncores filtering. Default: 8.'
)

# Treshold_calling
my_parser.add_argument(
    '--treshold_calling', 
    type=float,
    default=0.025,
    help='Treshold for variant allele calling. Default: .025.'
)

# GT_reference
my_parser.add_argument(
    '--GT_reference', 
    type=str,
    default='no_reference',
    help='Read GT reference genomic barcodes or not. Default: GBC.'
)

# Parse arguments
args = my_parser.parse_args()

path_data = args.path_data
sample = args.sample
filtering = args.filtering
ncores = args.ncores
t = args.treshold_calling
GT_reference = args.GT_reference


##


########################################################################

# Preparing run: import code, prepare directories

# Code
from scipy.sparse import save_npz
from mito_utils.preprocessing import *

# Paths and dirs
os.mkdir('input_folder')
os.chdir('input_folder')

########################################################################

# Main
def main():

    # Read AFM and filter vars
    with_GBC = True if GT_reference == 'GBC' else False
    afm = read_one_sample(path_data, sample, with_GBC=with_GBC)

    if filtering == 'GT':
        if with_GBC:
            _, a = filter_afm_with_gt(afm, min_cells_clone=5)
        else:
            raise ValueError('GT only for dataset with lentiviral barcoding Ground Truth!')
    else:
        _, a = filter_cells_and_vars(
            afm, filtering=filtering, path_=os.getcwd(), nproc=ncores
        )

    # Remove zeros and get AD, DP
    a = nans_as_zeros(a)
    AD, DP, _ = get_AD_DP(a)

    # Write meta, variants, AD and DP
    a.obs.to_csv('meta.csv')
    a.var_names.to_series().to_csv('variants.csv', index=False, header=None)
    save_npz('AD.npz', AD.T.astype(np.float16))
    save_npz('DP.npz', DP.T.astype(np.float16))


    ##
    
#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################









