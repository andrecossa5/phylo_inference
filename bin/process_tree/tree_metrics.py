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
from sklearn.metrics import normalized_mutual_info_score
from mito_utils.utils import *
from mito_utils.metrics import *
from mito_utils.phylo import *
from mito_utils.clustering import *


##


########################################################################

# Main
def main():

    # Metrics
    metrics = {}

    # Load tree
    with open(path_tree, 'rb') as f:
        tree = pickle.load(f)
    
    # Bootstrap support
    supports = []
    times = []
    nodes = []
    for node in tree.internal_nodes:
        try:
            supports.append(tree.get_attribute(node, 'support'))
            times.append(tree.get_time(node))
            nodes.append(node)
        except:
            pass

    df_support = pd.DataFrame({'support':supports, 'node':nodes, 'time':times})
    metrics['median_support'] = df_support['support'].median()
    metrics['median_time'] = df_support['time'].median()
    metrics['median_support_upmost_nodes'] = df_support.query('time<=10')['support'].median()

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
    metrics['ratio_char_supported_nodes_from_ASR'] = n_/n

    # Save
    to_frame = lambda x: pd.Series(x).to_frame('value').reset_index(names='metric').assign(job_id=job_id)
    to_frame(metrics).to_csv('tree_metrics.csv')


##


#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################