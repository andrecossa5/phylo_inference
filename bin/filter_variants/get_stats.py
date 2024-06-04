#!/usr/bin/python

# Get MT-SNVs stats

########################################################################
 
# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='get_stats',
    description=
    """
    Filter an AFM an return variants- and dataset-level stats.
    """
)

# Add arguments

# Path_main
my_parser.add_argument(
    '--path_data', 
    type=str,
    default=None,
    help='Path to samples input folder. Default: None.'
)

# Sample name
my_parser.add_argument(
    '--sample_name', 
    type=str,
    default=None,
    help='Sample name. Default: None.'
)

# Path_filtering
my_parser.add_argument(
    '--path_filtering', 
    type=str,
    default=None,
    help='Path to variant_filtering .json file. Default: None.'
)

# Filterin_key
my_parser.add_argument(
    '--filtering_key', 
    type=str,
    default=None,
    help='Key to a variant_filtering combination in the variant_filtering .json file. Default: None.'
)

# lineage_column
my_parser.add_argument(
    '--lineage_column', 
    type=str,
    default=None,
    help='AFM.obs column to evaluate MT-SNVs enrichment: Default: None.'
)

# solver
my_parser.add_argument(
    '--solver', 
    type=str,
    default='UPMGA',
    help='Tree solver: Default: UPMGA.'
)

# solver
my_parser.add_argument(
    '--metric', 
    type=str,
    default='jaccard',
    help='Metric for cell-cell dissimilarity computation: Default: jaccard.'
)

# solver
my_parser.add_argument(
    '--ncores', 
    type=int,
    default=8,
    help='n cores for distances computation: Default: 8.'
)

# solver
my_parser.add_argument(
    '--path_meta', 
    type=str,
    default=None,
    help='Path to cells metadata: Default: None.'
)

# path_priors
my_parser.add_argument(
    '--path_priors', 
    type=str,
    default=None,
    help='Path to MT-SNVs priors (i.e., median prevalence across datasets): Default: None.'
)

# path_priors
my_parser.add_argument(
    '--spatial_metrics', 
    type=str,
    default="False",
    help='Do or do not compute spatial metrics. Default: False.'
)

# cell_file
my_parser.add_argument(
    '--cell_file', 
    type=str,
    default="barcodes.txt",
    help='Sample to use. Default: barcodes.txt.'
)



## 


# Parse arguments
args = my_parser.parse_args()
path_data = args.path_data
sample_name = args.sample_name
path_filtering = args.path_filtering
filtering_key = args.filtering_key
lineage_column = args.lineage_column
solver = args.solver
metric = args.metric
spatial_metrics = args.spatial_metrics
ncores = args.ncores
path_meta = args.path_meta
path_priors = args.path_priors if os.path.exists(args.path_priors) else None
cell_file = args.cell_file


##


########################################################################

# Preparing run: import code

# Code
import json
from mito_utils.preprocessing import *
from mito_utils.clustering import *
from mito_utils.phylo import *
from mito_utils.utils import *
from mito_utils.plotting_base import *
from mito_utils.phylo_plots import *
import warnings
warnings.simplefilter('ignore')


##


def main():

    # Read AFM and add metadata
    meta = pd.read_csv(path_meta, index_col=0)
    with_GBC = False
    if lineage_column == 'GBC' and lineage_column in meta.columns:
        with_GBC = True
    afm = read_one_sample(path_data, sample=sample_name, with_GBC=with_GBC, cell_file=cell_file)
    sample_col = meta.columns[meta.columns.str.contains('sample')][0]
    meta = meta.loc[lambda x: meta[sample_col]==sample_name]
    afm.obs = afm.obs.join(meta[[ x for x in meta.columns if x not in afm.obs.columns ]])

    # Prep all filtering kwargs
    filtering, filtering_kwargs, kwargs = process_json(path_filtering, filtering_key)

    # Filter variants
    dataset_df, a = filter_cells_and_vars(
        afm, 
        sample_name=sample_name,
        filtering=filtering, 
        filtering_kwargs=filtering_kwargs,
        spatial_metrics=spatial_metrics,
        lineage_column=lineage_column,
        path_priors=path_priors,
        **kwargs
    )

    # Write down dfs
    a.var.assign(filtering_key=filtering_key, sample=sample_name).to_csv(f'{sample_name}_{filtering_key}_vars_df.csv')
    dataset_df.assign(filtering_key=filtering_key, sample=sample_name).to_csv(f'{sample_name}_{filtering_key}_dataset_df.csv')
    
    # Dists
    priors = 1-a.var['prior'].values if 'prior' in a.var.columns else None
    t = kwargs['af_confident_detection'] if 'af_confident_detection' in kwargs else .05
    D = pair_d(a, metric=metric, weights=priors, t=t)
    order = leaves_list(linkage(D))

    fig, ax = plt.subplots(figsize=(5,5))
    ax.imshow(D[np.ix_(order, order)], cmap='viridis_r')
    format_ax(ax, title=f'Metric: {metric}, ncells: {a.shape[0]}, nvars: {a.shape[1]}', 
              xticks=[], yticks=[], xlabel='Cells', ylabel='Cells')
    fig.tight_layout()
    fig.savefig(f'{sample_name}_{filtering_key}_distances.png', dpi=500)

    # Phylo (raw, tree, full tree with no pruned edges)
    tree = build_tree(a, solver=solver, metric=metric, t=t, weights=priors)

    fig, ax = plt.subplots(figsize=(7,7))
    plot_tree(tree, ax=ax, internal_node_kwargs={'markersize':0})
    fig.tight_layout()
    fig.savefig(f'{sample_name}_{filtering_key}_tree_raw.png', dpi=1000)

    # Get only branch-supporting muts, rebuild and visualize tree
    muts = get_supporting_muts(tree, a, t=t)
    dataset_df, a = filter_cells_and_vars(a, variants=muts, path_priors=path_priors, sample_name=sample_name)
    priors = 1-a.var['prior'].values if 'prior' in a.var.columns else None
    tree = build_tree(a, solver=solver, metric=metric, t=t, weights=priors)
    tree.cell_meta = pd.DataFrame(a.X, index=a.obs_names, columns=muts)

    fig, ax = plt.subplots(figsize=(7,7))
    muts = get_supporting_muts(tree, a, t=t)
    plot_tree(
        tree, ax=ax, meta=muts, orient='down', 
        colorstrip_spacing=0.0001, colorstrip_width=1.5,
        internal_node_kwargs={'markersize':0}
    )
    fig.tight_layout()
    fig.savefig(f'{sample_name}_{filtering_key}_tree_muts.png', dpi=1000)


##


########################################################################

# Run
if __name__ == '__main__': 
    main()





