#!/usr/bin/python

# prep_MAESTER script

########################################################################

# Libraries 
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='prep_MAESTER',
    description=
    """
    Prepare input for tree building: character/distance matrices and sequences (.fasta) file.
    """
)

# Add arguments

my_parser.add_argument(
    '--path_afm', 
    type=str,
    default='..',
    help='Path to <name>.h5ad file from mito_preprocessing pipelin. Default: .. .'
)

my_parser.add_argument(
    '--path_meta', 
    type=str,
    default=None,
    help='Path to meta_cells.csv file. Default: None .'
)

my_parser.add_argument(
    '--path_char_filtering', 
    type=str,
    default=None,
    help='Path to char_filtering_ops.json file. Default: None .'
)

my_parser.add_argument(
    '--path_cell_filtering', 
    type=str,
    default=None,
    help='Path to cell_filtering_ops.json file. Default: None.'
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
    '--char_filtering_key', 
    type=str,
    default='default',
    help='Filtering option in char_filtering_ops.json. Default: default.'
)

my_parser.add_argument(
    '--cell_filtering_key', 
    type=str,
    default='filter1',
    help='Filtering option in cell_filtering_ops.json. Default: default.'
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
    '--sample', 
    type=str,
    default='MDA_clones',
    help='Sample to use. Default: MDA_clones.'
)

my_parser.add_argument(
    '--job_id', 
    type=str,
    default='job_id',
    help='Job id. Default: None.'
)

my_parser.add_argument(
    '--lineage_column', 
    type=str,
    default=None,
    help='Lineage column for benchmarking, if necessary. Default: None.'
)

my_parser.add_argument(
    '--n_cores', 
    type=int,
    default='8',
    help='n cores to use. Default: 8.'
)

my_parser.add_argument(
    '--cell_file', 
    type=str,
    default="None",
    help='Path to subset of cells to utilize for the analysis. Default: None.'
)

my_parser.add_argument(
    '--path_dbSNP', 
    type=str,
    default=None,
    help='Path to dbSNP database. Default: None.'
)

my_parser.add_argument(
    '--path_REDIdb', 
    type=str,
    default=None,
    help='Path to REDIdb database. Default: None.'
)


##


# Parse arguments
args = my_parser.parse_args()

path_afm = args.path_afm
path_meta = args.path_meta
path_char_filtering = args.path_char_filtering
path_cell_filtering = args.path_cell_filtering
path_bin = args.path_bin
path_tree = args.path_tree
char_filtering_key = args.char_filtering_key
cell_filtering_key = args.cell_filtering_key
bin_key = args.bin_key
tree_key = args.tree_key
sample = args.sample
job_id = args.job_id
lineage_column = args.lineage_column
n_cores = args.n_cores
cell_file = args.cell_file if args.cell_file != "None" else None
path_dbSNP = args.path_dbSNP
path_REDIdb = args.path_REDIdb


##


########################################################################

# Preparing run: import code, prepare directories

# Code
from scipy.sparse import save_npz, load_npz
from mito_utils.utils import *
from mito_utils.preprocessing import *

########################################################################

# Main
def main():

    os.mkdir(f'{job_id}_pp')
    os.chdir(f'{job_id}_pp')

    # Parse kwargs options
    filtering, char_filtering_kwargs, kwargs = process_char_filtering_kwargs(path_char_filtering, char_filtering_key)
    cell_filtering_kwargs = process_kwargs(path_cell_filtering, cell_filtering_key)
    bin_method, bin_kwargs = process_bin_kwargs(path_bin, bin_key)
    tree_kwargs = process_kwargs(path_tree, tree_key)

    # Read AFM and add metadata
    afm = make_AFM(path_afm, cell_file=cell_file, path_meta=path_meta, sample=sample, **cell_filtering_kwargs)
    
    # Filter MT-SNVs and calculate metrics
    a, dataset_df, X_bin = filter_AFM(
        afm, 
        filtering=filtering, 
        filtering_kwargs=char_filtering_kwargs,
        bin_method=bin_method,
        binarization_kwargs=bin_kwargs,
        tree_kwargs=tree_kwargs,
        path_dbSNP=path_dbSNP,
        path_REDIdb=path_REDIdb,
        lineage_column=lineage_column,
        return_X_bin=True,
        **kwargs
    )

    # Get AD, DP
    AD, DP, _ = get_AD_DP(a)

    # Write AD and DP, cell and var meta, and dataset metric 
    save_npz('AD.npz', AD)
    save_npz('DP.npz', DP)
    save_npz('X_bin.npz', X_bin)
    a.obs.to_csv('cell_meta.csv')
    a.var.to_csv('char_meta.csv')
    dataset_df.to_csv('dataset_df.csv')
    seqs = AFM_to_seqs(a, X_bin=X_bin.A.T)
    with open('genotypes.fa', 'w') as f:
        for cell in seqs:
            f.write(f'>{cell}\n')
            f.write(f'{seqs[cell]}\n')
    
    # Write pp options
    kwargs['filtering'] = filtering
    to_frame_kwargs = lambda x: pd.Series(x).to_frame('value').reset_index(names='option').assign(job_id=job_id)
    to_frame_kwargs({**char_filtering_kwargs,**kwargs}).to_csv('char_filtering_ops.csv')
    to_frame_kwargs(cell_filtering_kwargs).to_csv('cell_filtering_ops.csv')
    to_frame_kwargs(bin_kwargs).to_csv('bin_ops.csv')
    to_frame_kwargs(tree_kwargs).to_csv('tree_ops.csv')


    ##


#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################