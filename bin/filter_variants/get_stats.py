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
path_priors = args.path_priors


##


########################################################################

# Preparing run: import code

# Code
import json
from mito_utils.preprocessing import *
from mito_utils.utils import *
import warnings
warnings.simplefilter('ignore')


##


def main():
    
    # Read AFM and add cells meta
    afm = read_one_sample(path_data, sample_name)
    meta = pd.read_csv(path_meta, index_col=0)
    sample_col = meta.columns[meta.columns.str.contains('sample')][0]
    meta = meta.loc[lambda x: meta[sample_col]==sample_name]
    afm.obs = afm.obs.join(meta[[ x for x in meta.columns if x not in afm.obs.columns ]])

    # Prep filtering kwargs
    with open(path_filtering, 'r') as file:
        FILTERING_OPTIONS = json.load(file)

    if filtering_key in FILTERING_OPTIONS:
        d = FILTERING_OPTIONS[filtering_key]
        filtering = d['filtering']
        filtering_kwargs = d['filtering_kwargs'] if 'filtering_kwargs' in d else {}
    else:
        raise KeyError(f'{filtering_key} not in {path_filtering}!')

    # Prep tree_kwargs
    af_confident_detection = filtering_kwargs['af_confident_detection'] if 'af_confident_detection' in filtering_kwargs else .05
    tree_kwargs = { 'solver' : solver, 'metric' : metric, 'ncores' : ncores, 't' : af_confident_detection }

    # Run
    t = Timer()
    t.start()

    dataset_df, a = filter_cells_and_vars(
        afm, 
        sample_name=sample_name,
        filtering=filtering, 
        filtering_kwargs=filtering_kwargs,
        tree_kwargs=tree_kwargs,
        lineage_column=lineage_column,
        spatial_metrics=True if spatial_metrics == "True" else False, 
        fit_mixtures=True, 
        path_priors=path_priors if os.path.exists(path_priors) else None 
    )
    a.var.assign(filtering_key=filtering_key, sample=sample_name).to_csv(f'{sample_name}_{filtering_key}_vars_df.csv')
    dataset_df.assign(filtering_key=filtering_key, sample=sample_name).to_csv(f'{sample_name}_{filtering_key}_dataset_df.csv')
    pd.Series(filtering_kwargs).to_frame().T.assign(filtering_key=filtering_key, sample=sample_name).to_csv(f'{sample_name}_{filtering_key}_job_df.csv')

    print(f'Job finished: {t.stop()}\n')


##


########################################################################

# Run
if __name__ == '__main__': 
    main()





