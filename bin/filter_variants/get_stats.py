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

# Filtering
my_parser.add_argument(
    '--min_site_cov', 
    type=int,
    default=5,
    help='Min mean site coverage. Default: 5.'
)

# Metric
my_parser.add_argument(
    '--min_var_quality', 
    type=int,
    default=30,
    help='Min mean base quality of MT-SNV (consensus) UMIs base qualities. Default: 30.'
)

# Solver
my_parser.add_argument(
    '--min_frac_negative', 
    type=float,
    default=.9,
    help='Min fraction of negative cells. Default: .9'
)

# boot_trees
my_parser.add_argument(
    '--min_n_positive', 
    type=int,
    default=2,
    help='Min n of positive cells. Default: 2.'
)

# boot_trees
my_parser.add_argument(
    '--low_confidence_af', 
    type=float,
    default=.01,
    help='AF value considered as "low confidence" detection. Default: .01'
)

# boot_trees
my_parser.add_argument(
    '--high_confidence_af', 
    type=float,
    default=.1,
    help='AF value considered as "high confidence" detection. Default: .1'
)

# obs_tree
my_parser.add_argument(
    '--min_prevalence_low_confidence_af', 
    type=float,
    default=.01,
    help='Min prevalence of the muts, detected at AF >= low_confidence_af. Default: .01'
)

# min_cells_high_confidence_af
my_parser.add_argument(
    '--min_cells_high_confidence_af', 
    type=int,
    default=2,
    help='Min n cells for which the variant has been detected at AF >= high_confidence_af. Default: 2'
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
    '--job_id', 
    type=int,
    default=None,
    help='Job id (numeric): Default: None.'
)


## 


# Parse arguments
args = my_parser.parse_args()
path_data = args.path_data
sample_name = args.sample_name
min_site_cov = args.min_site_cov
min_var_quality = args.min_var_quality
min_frac_negative = args.min_frac_negative
min_n_positive = args.min_n_positive
low_confidence_af = args.low_confidence_af
high_confidence_af = args.high_confidence_af
min_prevalence_low_confidence_af = args.min_prevalence_low_confidence_af
min_cells_high_confidence_af = args.min_cells_high_confidence_af
lineage_column = args.lineage_column
solver = args.solver
metric = args.metric
ncores = args.ncores
path_meta = args.path_meta
path_priors = args.path_priors
job_id = args.job_id


##


########################################################################

# Preparing run: import code

# Code
from mito_utils.preprocessing import *
from mito_utils.utils import *
import warnings
warnings.simplefilter('ignore')


##


def main():
    
    # Read AFM and add cells meta
    afm = read_one_sample(path_data, sample_name)
    meta = pd.read_csv(path_meta, index_col=0)
    meta = meta.query('sample_id==@sample_name')
    afm.obs = afm.obs.join(meta)

    # Prep kwargs
    filtering_kwargs = {
        'min_site_cov' : min_site_cov,
        'min_var_quality' : min_var_quality, 
        'min_frac_negative' : min_frac_negative,
        'min_n_positive' : min_n_positive, 
        'low_confidence_af' : low_confidence_af,
        'high_confidence_af' : high_confidence_af, 
        'min_prevalence_low_confidence_af' : min_prevalence_low_confidence_af,
        'min_cells_high_confidence_af' : min_cells_high_confidence_af
    }
    tree_kwargs = {
        'solver' : solver,
        'metric' : metric,
        'ncores' : ncores 
    }

    # Run
    t = Timer()
    t.start()

    vars_df, dataset_df, _ = filter_cells_and_vars(
        afm, 
        filtering='weng2024', 
        filtering_kwargs=filtering_kwargs,
        tree_kwargs=tree_kwargs,
        lineage_column=lineage_column,
        spatial_metrics=True, 
        fit_mixtures=False, 
        path_priors=path_priors
    )
    vars_df.assign(job_id=f'job_{job_id}').to_csv('vars_df.csv')
    dataset_df.assign(job_id=f'job_{job_id}').to_csv('dataset_df.csv')
    job_params = ','.join([f'job_{job_id}'] + list(filtering_kwargs.values()))
    with open('job_df.csv', 'w') as f:
        f.write(job_params)

    print(f'Job finished: {t.stop()}\n')


##


########################################################################

# Run
if __name__ == '__main__': 
    main()





