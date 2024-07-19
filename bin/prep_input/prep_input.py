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
    default=None,
    help='Path to priors.csv file. Default: None .'
)

# Path meta
my_parser.add_argument(
    '--path_filtering', 
    type=str,
    default=None,
    help='Path to filtering_options.csv file. Default: None .'
)

# nmads MT-coverage
my_parser.add_argument(
    '--nmads', 
    type=int,
    default=7,
    help='nmads MT-coverage: Default: 7.'
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
    '--sample', 
    type=str,
    default='MDA_clones',
    help='Sample to use. Default: MDA_clones.'
)

# Sample
my_parser.add_argument(
    '--n_cores', 
    type=int,
    default='8',
    help='n cores to use. Default: 8.'
)

# Sample
my_parser.add_argument(
    '--lineage_column', 
    type=str,
    default=None,
    help='Sample to use. Default: None.'
)

# Sample
my_parser.add_argument(
    '--cell_file', 
    type=str,
    default="barcodes.txt",
    help='Sample to use. Default: barcodes.txt.'
)

# Parse arguments
args = my_parser.parse_args()

path_data = args.path_data
path_meta = args.path_meta
path_priors = args.path_priors
path_filtering = args.path_filtering
sample = args.sample
nmads = args.nmads
filtering_key = args.filtering_key
n_cores = args.n_cores
lineage_column = args.lineage_column
cell_file = args.cell_file


##


########################################################################

# Preparing run: import code, prepare directories

# Code
from scipy.sparse import save_npz
from mito_utils.utils import *
from mito_utils.preprocessing import *

# Paths and dirs
os.mkdir('input_folder')
os.chdir('input_folder')

########################################################################

# Main
def main():

    # Read AFM and add metadata
    meta = pd.read_csv(path_meta, index_col=0)
    with_GBC = False
    if lineage_column == 'GBC' and lineage_column in meta.columns:
        with_GBC = True
    afm = read_one_sample(path_data, sample=sample, with_GBC=with_GBC, cell_file=cell_file, nmads=nmads)
    sample_col = meta.columns[meta.columns.str.contains('sample')][0]
    meta = meta.loc[lambda x: meta[sample_col]==sample]
    afm.obs = afm.obs.join(meta[[ x for x in meta.columns if x not in afm.obs.columns ]])

    # Prep all filtering kwargs
    filtering, filtering_kwargs, kwargs = process_json(path_filtering, filtering_key)

    # Filter variants
    dataset_df, a = filter_cells_and_vars(
        afm, 
        sample_name=sample,
        filtering=filtering, 
        nproc=n_cores,
        filtering_kwargs=filtering_kwargs,
        lineage_column=lineage_column,
        path_priors=path_priors if os.path.exists(path_priors) else None,
        **kwargs
    )

    # Get AD, DP
    AD, DP, _ = get_AD_DP(a)

    # Write meta, variants, AD and DP
    a.obs.to_csv('meta.csv')
    a.var.to_csv('var_meta.csv')
    dataset_df.to_csv('dataset_df.csv')
    a.var_names.to_series().to_csv('variants.csv', index=False, header=None)
    save_npz('AD.npz', AD.T.astype(np.float16))
    save_npz('DP.npz', DP.T.astype(np.float16))


    ##
    
#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################

