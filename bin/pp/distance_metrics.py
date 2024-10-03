#!/usr/bin/python

# prep_MAESTER script

########################################################################

# Libraries 
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='distance_metrics',
    description=
    """
    Prepare input for tree building: character/distance matrices and sequences (.fasta) file.
    """
)

# Add arguments
my_parser.add_argument(
    '--path_dists', 
    type=str,
    default=None,
    help='List of paths to all distances .npz files. Default: None'
)

my_parser.add_argument(
    '--replicates', 
    type=str,
    default=None,
    help='List of replicates identifiers mapping to all distances .npz files. Default: None'
)

my_parser.add_argument(
    '--path_meta', 
    type=str,
    default=None,
    help='Path to cells_meta.csv file. Default: None .'
)

my_parser.add_argument(
    '--job_id', 
    type=str,
    default=None,
    help='job_id identifier. Default: None .'
)

my_parser.add_argument(
    '--n_cores', 
    type=int,
    default=8,
    help='n cores to use. Default: 8.'
)

my_parser.add_argument(
    '--K', 
    type=int,
    default=10,
    help='k NN for NN metrics. Default: 10.'
)

my_parser.add_argument(
    '--lineage_column', 
    type=str,
    default=None,
    help='Lineage column for benchmark. Default: None.'
)



##


# Parse arguments
args = my_parser.parse_args()
path_dists = args.path_dists.strip('[|]').split(', ')
replicates = args.replicates.strip('[|]').split(', ')
path_meta = args.path_meta
K = args.K
job_id = args.job_id
lineage_column = args.lineage_column


##


########################################################################

# Preparing run: import code, prepare directories

# Code
from scipy.sparse import load_npz
from mito_utils.utils import *
from mito_utils.metrics import *
from mito_utils.phylo import *

########################################################################

# Main
def main():

    # Metrics d
    metrics = {}

    # Load distance matrices
    D = {}
    for k,v in zip(replicates, path_dists):
        D[k] = load_npz(v).A

    # Load cell meta
    meta = pd.read_csv(path_meta, index_col=0)

    if lineage_column is not None:
        
        labels = meta[lineage_column].astype(str)
    
        # kNN metrics
        idx = kNN_graph(D=D['observed'], k=K, from_distances=True)[0]
        _, _, acc_rate = kbet(idx, labels, only_score=False)
        median_entropy = NN_entropy(idx, labels)
        median_purity = NN_purity(idx, labels)
        metrics['kBET_rejection_rate'] = 1-acc_rate
        metrics['median_NN_entropy'] = median_entropy
        metrics['median_NN_purity'] = median_purity
    
        # AUPRC
        metrics['AUPRC'] = distance_AUPRC(D['observed'], labels)

    # Corr
    L = []
    for k in D:
        L.append(D[k].flatten())
    del D 
    metrics['corr'] = np.mean(np.corrcoef(np.array(L)))

    # Save
    to_frame = lambda x: pd.Series(x).to_frame('value').reset_index(names='metric').assign(job_id=job_id)
    to_frame(metrics).to_csv('distance_metrics.csv')


#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################


