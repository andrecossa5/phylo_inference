#!/usr/bin/python

# Final filter script

########################################################################

# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='final_filter',
    description=
    """
    Final filtering step to remove cells and variants before final tree reconstruction.
    """
)

# Add arguments

# Path_main
my_parser.add_argument(
    '-p', 
    '--path_input', 
    type=str,
    default='..',
    help='Path to AD, DP folder. Default: .. .'
)

# cell_assignment
my_parser.add_argument(
    '--filtering_key', 
    type=str,
    default=None,
    help='filtering_key in path_filtering .json file. Default: None .'
)

# cell_assignment
my_parser.add_argument(
    '--path_filtering', 
    type=str,
    default=None,
    help='path_filtering .json file to specify filtering options. Default: None .'
)

# cell_assignment
my_parser.add_argument(
    '--path_cell_assignment', 
    type=str,
    default=None,
    help='Path to cell_assignment.csv df. Default: None .'
)

# var_assignment
my_parser.add_argument(
    '--path_var_assignment', 
    type=str,
    default=None,
    help='Path to cell_assignment.csv df. Default: None .'
)

# var_assignment
my_parser.add_argument(
    '--min_assignment_prob', 
    type=float,
    default=.3,
    help='Minimum assignment probability for a variant. Default: .3.'
)

# Parse arguments
args = my_parser.parse_args()

path_i = args.path_input
filtering_key = args.filtering_key
path_filtering = args.path_filtering
path_cell_assignment = args.path_cell_assignment
path_var_assignment = args.path_var_assignment
min_assignment_prob = args.min_assignment_prob


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


    os.listdir(path_i)

    # Reconstruct AFM
    AD = load_npz(os.path.join(path_i, 'AD.npz')).astype(np.int16).A         # Cell x var
    DP = load_npz(os.path.join(path_i, 'DP.npz')).astype(np.int16).A
    meta = pd.read_csv(os.path.join(path_i, 'meta.csv'), index_col=0)
    var_meta = pd.read_csv(os.path.join(path_i, 'var_meta.csv'), index_col=0)
    afm = AnnData(X=np.divide(AD, DP), obs=meta, var=var_meta, dtype=np.float16)

    # Read cell and var assignments files
    filtered_cells = pd.read_csv(path_cell_assignment)['Cell'].to_list()
    vars_df = pd.read_csv(path_var_assignment, index_col=0)
    vars_df.index = vars_df.index.map(lambda x: x[1:].replace('.', '>'))    # R problems...
    filtered_vars = vars_df.query('p>=@min_assignment_prob').index.to_list()

    # Create tests
    test_cells = afm.obs_names.isin(filtered_cells)
    test_vars = afm.var_names.isin(filtered_vars)

    # Write filtered AD, DP, cells and variants.
    os.mkdir('filtered_input')
    os.chdir('filtered_input')  
    meta.loc[test_cells].to_csv('meta.csv')
    var_meta.loc[test_vars].to_csv('var_meta.csv')
    var_meta.loc[test_vars].index.to_series().to_csv('variants.csv', index=False, header=None)
    save_npz('AD.npz', csr_matrix(AD[np.ix_(test_cells, test_vars)].astype(np.float16)))
    save_npz('DP.npz', csr_matrix(DP[np.ix_(test_cells, test_vars)].astype(np.float16)))

    # Prep fasta

    # Read t, if available
    with open(path_filtering, 'r') as file:
        FILTERING_OPTIONS = json.load(file)

    if filtering_key in FILTERING_OPTIONS:
        d = FILTERING_OPTIONS[filtering_key]
        t = d['t'] if 't' in d else .05
    else:
        raise KeyError(f'{filtering_key} not in {path_filtering}!')
    
    # Covert to .fasta and write
    seqs = AFM_to_seqs(afm[test_cells, test_vars], t=t)
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