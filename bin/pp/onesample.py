#!/usr/bin/python

# onesample script

########################################################################

# Libraries 
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='onesample',
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
    '--job_id', 
    type=str,
    default='job_0',
    help='job id. Default: job_0.'
)

my_parser.add_argument(
    '--cell_filter', 
    type=str,
    default='filter2',
    help='Cell filtering method. Default: filter2.'
)

my_parser.add_argument(
    '--filtering', 
    type=str,
    default='MI_TO',
    help='Variant filtering method. Default: MI_TO.'
)

my_parser.add_argument(
    '--min_cell_number', 
    type=int,
    default=0,
    help='Min number of cell in <lineage_column> categories to retain them. Default: 10.'
)

my_parser.add_argument(
    '--min_var_quality', 
    type=int,
    default=30,
    help='Min Q30 phred score of a MT-SNV ADs. Default: 30.'
)

my_parser.add_argument(
    '--min_cov', 
    type=int,
    default=10,
    help='Minimum site coverage. Default: 10.'
)

my_parser.add_argument(
    '--min_frac_negative', 
    type=float,
    default=.2,
    help='Minimum fraction of -cells to consider a MT-SNV. Default: .2.'
)

my_parser.add_argument(
    '--min_n_positive', 
    type=int,
    default=5,
    help='Minimum number of +cells to consider a MT-SNV. Default: 2.'
)

my_parser.add_argument(
    '--af_confident_detection', 
    type=float,
    default=.01,
    help='Allelic Frequency of confident detection. Default: .01.'
)

my_parser.add_argument(
    '--min_n_confidently_detected', 
    type=int,
    default=2,
    help='Minimum number of confidently detected +cells to consider a MT-SNV. Default: 2.'
)

my_parser.add_argument(
    '--min_mean_AD_in_positives', 
    type=float,
    default=1.5,
    help='Minimum number of mean AD in +cells to consider a MT-SNV. Default: 1.5.'
)

my_parser.add_argument(
    '--min_mean_DP_in_positives', 
    type=float,
    default=20,
    help='Minimum number of mean DP in +cells to consider a MT-SNV. Default: 20.'
)

my_parser.add_argument(
    '--t_prob', 
    type=float,
    default=.7,
    help='Probability threshold for assigning cells to 0/1 mixture binomial components if bin_method=MI_TO. Default: .7.'
)

my_parser.add_argument(
    '--t_vanilla', 
    type=float,
    default=0,
    help='AF threshold to assigning cells to 0/1 genotypes if bin_method=MI_TO or vanilla. Default: 0.'
)

my_parser.add_argument(
    '--min_AD', 
    type=int,
    default=1,
    help='Min number of AD to assign a 0/1 genotype. Default: 1.'
)

my_parser.add_argument(
    '--bin_method', 
    type=str,
    default='MI_TO',
    help='Binarization method. Default: MI_TO.'
)

my_parser.add_argument(
    '--metric', 
    type=str,
    default='jaccard',
    help='Distance metric for cell-cell distance computation. Default: jaccard.'
)

my_parser.add_argument(
    '--solver', 
    type=str,
    default='NJ',
    help='Tree-building algorithm for fast lineage inference. Default: NJ.'
)

my_parser.add_argument(
    '--min_cell_prevalence', 
    type=float,
    default=1,
    help='Min number of AD to assign a 0/1 genotype. Default: 1.'
)

my_parser.add_argument(
    '--lineage_column', 
    type=str,
    default=None,
    help='Lineage column (i.e., GBC for benchmarks). Default: None.'
)

my_parser.add_argument(
    '--ncores', 
    type=int,
    default=8,
    help='n cores for computing distances. Default: 8.'
)

my_parser.add_argument(
    '--cell_file', 
    type=str,
    default=None,
    help='Path selected subset of cell barcodes to analyze (N.B. in <CB>_<sample_name> format). Default: None.'
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
job_id = args.job_id
cell_filter = args.cell_filter if args.cell_filter != "None" else None
filtering = args.filtering if args.filtering != "None" else None
min_cell_number = args.min_cell_number
min_cov = args.min_cov
min_var_quality = args.min_var_quality
min_frac_negative = args.min_frac_negative
min_n_positive = args.min_n_positive
af_confident_detection = args.af_confident_detection
min_n_confidently_detected = args.min_n_confidently_detected
min_mean_AD_in_positives = args.min_mean_AD_in_positives
min_mean_DP_in_positives = args.min_mean_DP_in_positives
t_prob = args.t_prob
t_vanilla = args.t_vanilla
min_AD = args.min_AD
min_cell_prevalence = args.min_cell_prevalence
bin_method = args.bin_method
metric = args.metric
metric = 'custom_MI_TO_jaccard' if bin_method == 'MI_TO' else metric
solver = args.solver
lineage_column = args.lineage_column
ncores = args.ncores
cell_file = args.cell_file if args.cell_file != "None" else None
path_dbSNP = args.path_dbSNP
path_REDIdb = args.path_REDIdb

##

# Parameter space
# grid = list(
#     product(
#         ['mito_preprocessing', 'maegatk'],  # pp_method
#         np.linspace(.005, .05, 5),          # af_confident_detection
#         np.arange(1, 3+1, 1),               # min_n_confidently_detected
#         np.linspace(1, 3, 5),               # min_mean_AD_in_positives
#         [1,2],                              # min_AD 
#         ['vanilla', 'MI_TO']                # bin_method 
#     )
# )
# ntot = len(grid)


##


########################################################################

# Preparing run: import code, prepare directories

# Code
import os
import pickle
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.phylo import *
from mito_utils.clustering import custom_ARI
from sklearn.metrics import normalized_mutual_info_score

########################################################################


def main():

    afm = sc.read(path_afm)
    afm = filter_cells(afm, cell_filter=cell_filter)

    if filtering == "MI_TO":
        filtering_kwargs = {
            'min_cov' : min_cov,
            'min_var_quality' : min_var_quality,
            'min_frac_negative' : min_frac_negative,
            'min_n_positive' : min_n_positive,
            'af_confident_detection' : af_confident_detection,
            'min_n_confidently_detected' : min_n_confidently_detected,
            'min_mean_AD_in_positives' : min_mean_AD_in_positives,       # 1.25,
            'min_mean_DP_in_positives' : min_mean_DP_in_positives
        }
    else:
        filtering_kwargs = {}

    binarization_kwargs = {
        't_prob':t_prob, 't_vanilla':t_vanilla, 'min_AD':min_AD, 'min_cell_prevalence':min_cell_prevalence
    }
    tree_kwargs = {'metric':metric, 'solver':solver}

    afm, tree = filter_afm(
        afm,
        min_cell_number=min_cell_number,
        lineage_column=lineage_column,
        filtering=filtering,
        filtering_kwargs=filtering_kwargs,
        binarization_kwargs=binarization_kwargs,
        bin_method=bin_method,
        tree_kwargs=tree_kwargs,
        path_dbSNP=path_dbSNP, 
        path_REDIdb=path_REDIdb,
        spatial_metrics=True,
        compute_enrichment=True,
        max_AD_counts=2,
        ncores=ncores,
        return_tree=True
    )
    tree, _, _ = cut_and_annotate_tree(tree)
    
    # Prep stats dictionary
    stats = {}
    
    # Options
    options = {}
    options['pp_method'] = afm.uns['pp_method']
    options['min_cell_number'] = min_cell_number
    options['lineage_column'] = lineage_column
    options['cell_filter'] = afm.uns['cell_filter']
    options['filtering'] = filtering
    options['filtering_kwargs'] = filtering_kwargs
    options['bin_method'] = bin_method
    options['binarization_kwargs'] = binarization_kwargs
    options['tree_kwargs'] = tree_kwargs
    stats['options'] = options

    # Metrics
    metrics = {}
    metrics['n_MT_clone'] = tree.cell_meta['MT_clone'].nunique()
    metrics['unassigned'] = tree.cell_meta['MT_clone'].isna().sum()
    metrics['corr'] = calculate_corr_distances(tree)[0]

    if lineage_column is not None:
        metrics[f'n_{lineage_column}_groups'] = tree.cell_meta[lineage_column].nunique()
        metrics['AUPRC'] = distance_AUPRC(afm.obsp['distances'].A, afm.obs[lineage_column])
        metrics['ARI'] = custom_ARI(tree.cell_meta[lineage_column], tree.cell_meta['MT_clone'])
        metrics['NMI'] = normalized_mutual_info_score(tree.cell_meta.dropna()[lineage_column], tree.cell_meta.dropna()['MT_clone'])

    if 'raw_basecalls_metrics' in afm.uns:
        metrics.update(afm.uns['raw_basecalls_metrics'])
        
    metrics.update(afm.uns['dataset_metrics'])
    metrics['n_dbSNP'] = afm.uns['char_filter']['n_dbSNP'] 
    metrics['n_REDIdb'] = afm.uns['char_filter']['n_REDIdb'] 
    stats['metrics'] = metrics

    # Cells and vars
    stats['cells'] = afm.obs_names
    stats['vars'] = afm.var_names

    # Save
    with open(f'tuning{job_id}_stats.pickle', 'wb') as f:
        pickle.dump(stats, f)


##


#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################
