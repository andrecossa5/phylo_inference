#!/usr/bin/python

# MAESTER script

########################################################################

# Libraries 
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='MiTo',
    description=
    """
    Prepare input for tree building, from MAESTER or RedeeM Allele Frequency Matrix.
    """
)

# Add arguments

my_parser.add_argument(
    '--path_afm', 
    type=str,
    default='.',
    help='Path to afm.h5ad file. Default: . .'
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
    '--cell_filter', 
    type=str,
    default='filter2',
    help='Cell filtering method. Default: filter2.'
)

my_parser.add_argument(
    '--filtering', 
    type=str,
    default=None,
    help='Variant filtering method. Default: None (i.e., all MT-SNVs variants will be retained).'
)

my_parser.add_argument(
    '--min_cell_number', 
    type=int,
    default=0,
    help='Min number of cell in <lineage_column> categories to retain them. Default: 0.'
)

my_parser.add_argument(
    '--min_cov', 
    type=int,
    default=10,
    help='Minimum mean coverage of the candidate variant site. Default: 10.'
)

my_parser.add_argument(
    '--min_var_quality', 
    type=int,
    default=30,
    help='Min phred score of a MT-SNV ADs. Default: 30.'
)

my_parser.add_argument(
    '--min_frac_negative', 
    type=float,
    default=.2,
    help='Minimum fraction of negative (i.e., AF==0) cells to consider a MT-SNV. Default: .2.'
)

my_parser.add_argument(
    '--min_n_positive', 
    type=int,
    default=2,
    help='Minimum number of positive (i.e., AF>0) cells to consider a MT-SNV. Default: 2.'
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
    help='Minimum number of confidently detected positive cells to consider a MT-SNV. Default: 2.'
)

my_parser.add_argument(
    '--min_mean_AD_in_positives', 
    type=float,
    default=1.5,
    help='Minimum number of mean AD in positive cells to consider a MT-SNV. Default: 1.5.'
)

my_parser.add_argument(
    '--min_mean_DP_in_positives', 
    type=float,
    default=20,
    help='Minimum number of mean DP in positive cells to consider a MT-SNV. Default: 20.'
)

my_parser.add_argument(
    '--t_prob', 
    type=float,
    default=.7,
    help='Probability threshold for assigning cells to 0/1 mixture binomial components if bin_method=MiTo. Default: .7.'
)

my_parser.add_argument(
    '--t_vanilla', 
    type=float,
    default=0,
    help='AF threshold to assigning cells to 0/1 genotypes if bin_method=MiTo or vanilla. Default: 0.'
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
    default='MiTo',
    help='Binarization method. Default: MiTo.'
)

my_parser.add_argument(
    '--min_cell_prevalence', 
    type=float,
    default=.1,
    help='Min cell prevalence to assign 0/1 genotype with the MiTo method. Default: .1.'
)

my_parser.add_argument(
    '--lineage_column', 
    type=str,
    default=None,
    help='Lineage column (i.e., GBC for benchmarks). Default: None.'
)

my_parser.add_argument(
    '--solver', 
    type=str,
    default='UPMGA',
    help='Cassiopeia solver. Default: UPMGA.'
)

my_parser.add_argument(
    '--metric', 
    type=str,
    default='jaccard',
    help='Distance metric. Default: jaccard.'
)

my_parser.add_argument(
    '--n_cores', 
    type=int,
    default=1,
    help='n cores to use. Default: 1.'
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


########################################################################

# Code
import pickle
from mito_utils.utils import *
from mito_utils.preprocessing import *

########################################################################

# Main
def main():

    # Parse arguments
    args = my_parser.parse_args()

    path_afm = args.path_afm
    path_pickles = args.path_pickles
    sample = args.sample
    job_id = args.job_id
    cell_filter = args.cell_filter 
    filtering = args.filtering
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
    solver = args.solver
    metric = args.metric
    lineage_column = args.lineage_column
    n_cores = args.n_cores
    path_dbSNP = args.path_dbSNP
    path_REDIdb = args.path_REDIdb


    ##


    # Handle params
    if path_pickles is not None and job_id is not None:

        path_pickle = os.path.join(path_pickles, sample, f'{job_id}_stats.pickle')

        if os.path.exists(path_pickle):
            with open(path_pickle, 'rb') as f:
                d = pickle.load(f)

            params = d['options']
            cell_filter = params['cell_filter']['cell_filter']
            min_cell_number = params['min_cell_number']
            lineage_column = params['lineage_column']
            filtering = params['filtering']
            filtering_kwargs = params['filtering_kwargs']
            tree_kwargs = params['tree_kwargs']
            binarization_kwargs = params['binarization_kwargs']
            bin_method = params['bin_method']
        
        else:
            raise ValueError(f'{path_pickle} does not exists!')
    
    else:

        filtering_kwargs = {
            'min_cov' : min_cov,
            'min_var_quality': min_var_quality,
            'min_frac_negative' : min_frac_negative,
            'min_n_positive' : min_n_positive,
            'af_confident_detection' : af_confident_detection,
            'min_n_confidently_detected' : min_n_confidently_detected,
            'min_mean_AD_in_positives' : min_mean_AD_in_positives,
            'min_mean_DP_in_positives' : min_mean_DP_in_positives 
        }
        binarization_kwargs = {
            't_prob' : t_prob, 
            't_vanilla' : t_vanilla,
            'min_AD' : min_AD,
            'min_cell_prevalence' : min_cell_prevalence,
        }
        tree_kwargs = {'solver':solver, 'metric':metric}

    # Filter matrix and calculate metrics
    afm = sc.read(path_afm)
    afm = filter_cells(afm, cell_filter=cell_filter)
    afm = filter_afm(
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
        ncores=n_cores,
        return_tree=False,
    )

    # Write out filtered matrix
    afm.write('afm.h5ad')
            

    ##


#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################