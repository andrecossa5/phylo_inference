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

# metric
my_parser.add_argument(
    '--metric', 
    type=str,
    default='jaccard',
    help='Path to edges.csv. Default: jaccard.'
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
metric = args.metric
ncores = args.ncores

# path_i = '/Users/IEO5505/Desktop/MI_TO/phylo_inference/work/b1/0374586a971219073042563a8c7c7b/filtered_input'
# path_nodes = '/Users/IEO5505/Desktop/MI_TO/phylo_inference/work/b1/0374586a971219073042563a8c7c7b/nodes.csv'
# path_edges = '/Users/IEO5505/Desktop/MI_TO/phylo_inference/work/b1/0374586a971219073042563a8c7c7b/edges.csv'
# ncores = 8


########################################################################

# Preparing run: import code, prepare directories

# Code
import json
import pickle
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

    # Create CassiopeiaTree
    tree = create_from_annot(nodes, edges, afm)
    
    # General tree stats
    stats = {}
    stats['n_cells'] = tree.n_cell
    stats['n_vars'] = nodes['assigned_var'].dropna().unique().size
    stats['n_clones'] = nodes['clone'].dropna().unique().size
    stats['ratio_supported_unsupported_branches'] = stats['n_clones'] / len(tree.edges) 
    stats['median_support'] = nodes['support'].dropna().median()
    stats['std_support'] = nodes['support'].dropna().std()
    stats['median_time'] = np.mean(list(tree.get_times().values()))
    stats['mean_depth'] = tree.get_mean_depth_of_tree()

    # Spatial stats
    X_bin = np.where(tree.character_matrix>=.05, 1, 0)
    D = pairwise_distances(X_bin, metric=lambda x, y: np.sum(np.logical_and(x, y)))
    n_shared_muts = np.ma.masked_equal(D, np.diag(D)).mean(axis=1).data
    cell_conn = (D>0).sum(axis=1)-1
    stats['median_n_shared_muts'] = np.median(n_shared_muts)
    stats['median_connectedness'] = np.median(cell_conn)

    # Tree vs ch matrix distance correlation 
    stats['char_tree_dist_corr'] = calculate_corr_distances(tree)

    # Ch matrix distance robustness to bootstrapping
    n_samples = 100
    n_vars_sampling = round((tree.character_matrix.shape[1] / 100) * 80)
    L = []
    for _ in range(n_samples):
        vars_ = np.random.choice(tree.character_matrix.columns, size=n_vars_sampling, replace=False)
        D = pair_d(tree.character_matrix.loc[:,vars_].values, metric=metric, t=.05, ncores=ncores)
        L.append(D.flatten())
    stats['char_dist_boot'] = np.mean(np.corrcoef(np.array(L)))

    # Write final annotated tree as pickle, and its stats
    with open('annotated_tree.pickle', 'wb') as f:
        pickle.dump(tree, f)
    pd.Series(stats).to_frame('metric').to_csv('final_tree_stats.csv')


    ##
    
#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################