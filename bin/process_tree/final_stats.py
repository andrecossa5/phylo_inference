#!/usr/bin/python

# Final tree stats script

########################################################################

# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='final_stats',
    description=
    """
    Final tree last checks.
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
    '--path_filtering', 
    type=str,
    default=None,
    help='path_filtering .json file to specify filtering options. Default: None .'
)


# path_nodes
my_parser.add_argument(
    '--path_nodes', 
    type=str,
    default=None,
    help='Path to nodes.csv. Default: None.'
)

# path_edges
my_parser.add_argument(
    '--path_edges', 
    type=str,
    default=None,
    help='Path to edges.csv. Default: None.'
)

# ncores
my_parser.add_argument(
    '--ncores', 
    type=int,
    default=8,
    help='n cores execution. Default: 8.'
)


# Parse arguments
args = my_parser.parse_args()

path_i = args.path_input
path_nodes = args.path_nodes
path_edges = args.path_edges
ncores = args.ncores

# path_i = ''
# path_nodes = '...'
# path_edges = '...'
# ncores = 8

########################################################################

# Preparing run: import code, prepare directories

# Code
import json
from scipy.sparse import load_npz
from anndata import AnnData
from mito_utils.phylo import *

########################################################################

# Main
def main():

    # Reconstruct AFM: do we need it??
    AD = load_npz(os.path.join(path_i, 'AD.npz')).astype(np.int16).A                        # Cell x var
    DP = load_npz(os.path.join(path_i, 'DP.npz')).astype(np.int16).A
    meta = pd.read_csv(os.path.join(path_i, 'meta.csv'), index_col=0)
    var_meta = pd.read_csv(os.path.join(path_i, 'var_meta.csv'), index_col=0)
    afm = AnnData(X=np.divide(AD, DP), obs=meta, var=var_meta, dtype=np.float16)

    # Read annotated tree elements
    edges = pd.read_csv(path_edges, index_col=0)
    nodes = pd.read_csv(path_nodes, index_col=0)
    nodes['assigned_var'] = (
        nodes['assigned_var']
        .map(lambda x: x.replace('X', '').replace('.', '>') if pd.notna(x) else np.NaN)     # R problems...
    )

    # Build CassiopeiaTree
    tree = create_from_annot(nodes, edges)
    
    # General tree stats

    # Tree vs ch matrix distance correlation 

    # Ch matrix distance robustness to bootstrapping
 
    # Write stats
    os.system('touch final_tree_stats.csv')



    ##
    
#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################