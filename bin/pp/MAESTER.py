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
    '--path_pickles', 
    type=str,
    default=None,
    help='Path to pickles main folder. Default: None.'
)

my_parser.add_argument(
    '--sample', 
    type=str,
    default=None,
    help='Sample name. Default: None.'
)

my_parser.add_argument(
    '--job_id', 
    type=str,
    default=None,
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
    default=1,
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
path_pickles = args.path_pickles
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
import pickle
from mito_utils.utils import *
from mito_utils.preprocessing import *

########################################################################

# Main
def main():

    # Load job_id_stats
    with open(os.path.join(path_pickles, sample, f'{job_id}_stats.pickle'), 'rb') as f:
        d = pickle.load(f)

    # Read AFM and add metadata
    afm = sc.read(path_afm)
    
    # Filter MT-SNVs and calculate metrics
    
    # TO DO: cell_subset = pd.read_csv(cell_file)[0].to_list() if cell_file is not None else None
    
    afm = sc.read(path_afm)
    afm = filter_cells(afm, **d['options']['cell_filter'])
    afm = filter_afm(
        afm,
        min_cell_number=d['options']['min_cell_number'],
        lineage_column=d['options']['lineage_column'],
        filtering=d['options']['filtering'],
        filtering_kwargs=d['options']['filtering_kwargs'],
        binarization_kwargs=d['options']['binarization_kwargs'],
        bin_method=d['options']['bin_method'],
        tree_kwargs=d['options']['tree_kwargs'],
        path_dbSNP=path_dbSNP, 
        path_REDIdb=path_REDIdb,
        spatial_metrics=True,
        compute_enrichment=True,
        max_AD_counts=2,
        ncores=n_cores,
        return_tree=False
    )

    # Add .uns for distance calculations
    afm.uns['distance_calculations'] = {}
    afm.uns['distance_calculations']['distances'] = {}
    afm.uns['distance_calculations']['distances']['metric'] = d['options']['tree_kwargs']['metric']
    afm.write('afm.h5ad')
            

    ##


#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################