#!/usr/bin/python

# Build cassiopeia

########################################################################
 
# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='annotate_tree',
    description="Process and annotate tree."
)

my_parser.add_argument(
    '--tree', 
    type=str,
    default=None,
    help='Path to annotated_tree in .pickle format. Default: .. .'
)

my_parser.add_argument(
    '--lineage_column', 
    type=str,
    default=None,
    help='lineage_column. Default: None .'
)

my_parser.add_argument(
    '--job_id', 
    type=str,
    default=None,
    help='job_id. Default: None.'
)


##


# Parse arguments
args = my_parser.parse_args()
path_tree = args.tree
lineage_column = args.lineage_column
job_id = args.job_id


##


########################################################################

# Preparing run: import code, prepare directories, set logger

# Code
import pickle
from itertools import chain
from mito_utils.utils import *
from mito_utils.metrics import *
from mito_utils.phylo import *


##


########################################################################

# Main
def main():

    # Metrics
    metrics = {}

    # Load tree
    with open(path_tree, 'rb') as f:
        tree = pickle.load(f)
    
    # General 
    n_cells, n_vars = tree.layers['raw'].shape
    metrics['n_cells'] = n_cells
    metrics['n_vars'] = n_vars
    metrics['median_var_per_cell'] = np.median(np.sum(tree.layers['transformed']==1, axis=1))
    metrics['density'] = np.sum(tree.layers['transformed'].values==1) / np.product(tree.layers['transformed'].shape)

    # Bootstrap support
    df_support = get_internal_node_stats(tree)
    metrics['median_support'] = df_support['support'].median()
    metrics['median_support_mut_clades'] = df_support.loc[df_support['mut'],'support'].median()
    max_clade = np.percentile(df_support['clade_size'], 90)
    metrics['median_support_biggest_clades'] = df_support.query(f'clade_size>={max_clade}')['support'].median()

    # Time
    metrics['median_time'] = df_support['time'].median()

    # Others
    corr_dists, p_corr_dists = calculate_corr_distances(tree)
    metrics['corr_distances'] = corr_dists
    metrics['corr_distances_pvalue'] = p_corr_dists
    metrics['n_clones'] = tree.cell_meta['MT_clone'].nunique()
    metrics['median_n_cells_clone'] = tree.cell_meta['MT_clone'].value_counts().median()
    metrics['min_n_cells_clone'] = tree.cell_meta['MT_clone'].value_counts().min()
    metrics['max_n_cells_clone'] = tree.cell_meta['MT_clone'].value_counts().max()
    metrics['median_n_assigned_char_per_clone'] = np.median([ len(x.split(';')) for x in tree.cell_meta['clone_muts'].loc[lambda x: ~x.isna()].unique() ])
    metrics['total_assigned_char'] = len(list(chain.from_iterable([ x.split(';') for x in tree.cell_meta['clone_muts'].loc[lambda x: ~x.isna()].unique() ])))
    metrics['median_CI'] = np.median(CI(tree))
    metrics['median_RI'] = np.median(RI(tree))

    # Benchmark specifics
    if lineage_column is not None:
        test = ~tree.cell_meta['MT_clone'].isna()
        metrics['ARI'] = custom_ARI(tree.cell_meta.loc[test, lineage_column], tree.cell_meta.loc[test, 'MT_clone'])
        metrics['NMI'] = normalized_mutual_info_score(tree.cell_meta.loc[test, lineage_column], tree.cell_meta.loc[test, 'MT_clone'])

    # Collapse tree
    n = len(tree.internal_nodes)
    tree.collapse_mutationless_edges(True)
    n_ = len(tree.internal_nodes)
    metrics['frac_remaining_nodes_after_mutationless_edges_collapse'] = n_/n

    # Save
    to_frame = lambda x: pd.Series(x).to_frame('value').reset_index(names='metric').assign(job_id=job_id)
    to_frame(metrics).to_csv('tree_metrics.csv')


##


#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################