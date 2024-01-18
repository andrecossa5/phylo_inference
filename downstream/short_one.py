#!/usr/bin/python

# Build cassiopeia

########################################################################
 
# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='short and at the point',
    description=
    """
    Build cassiopeia trees.
    """
)

# Add arguments

# Path_main
my_parser.add_argument(
    '-p', 
    '--path_data', 
    type=str,
    default='..',
    help='Path to data folder. Default: .. .'
)

# Sample
my_parser.add_argument(
    '--sample', 
    type=str,
    default='MDA_clones',
    help='Sample to use. Default: MDA_clones.'
)

# Filter
my_parser.add_argument(
    '--filtering', 
    type=str,
    default='ludwig2019',
    help='Method to filter MT-SNVs. Default: ludwig2019.'
)

# Solver
my_parser.add_argument(
    '--solver', 
    type=str,
    default='UPMGA',
    help='''
    Bootstrap resampling strategy. Default: UPMGA. Other options: 
    NJ, spectral, max_cut, greedy. See Cassiopeia (MW Jones et al., 2020) 
    for docs.
    '''
)

# ncores
my_parser.add_argument(
    '--frac_subsampling', 
    type=float,
    default=None,
    help='Fraction of cells to subsample. Default: None.'
)

# Metric
my_parser.add_argument(
    '--metric', 
    type=str,
    default='cosine',
    help='''
    Distance metric. Default: cosine
    '''
)

# boot_strategy
my_parser.add_argument(
    '--boot_strategy', 
    type=str,
    default='jacknife',
    help='Bootstarpping strategy. Default: jacknife.'
)

# ncores
my_parser.add_argument(
    '--nboot', 
    type=int,
    default=100,
    help='n_bootstrap. Default: 100.'
)

# ncores
my_parser.add_argument(
    '--t', 
    type=float,
    default=.025,
    help='t. Default: .025.'
)

# ncores
my_parser.add_argument(
    '--ncores', 
    type=int,
    default=4,
    help='N cores for pairwise distance calculation. Default: 4.'
)

# Parse arguments
args = my_parser.parse_args()
path_data = args.path_data
sample = args.sample
frac_subsampling = args.frac_subsampling
filtering = args.filtering
solver = args.solver
metric = args.metric
ncores = args.ncores
n_boot = args.nboot
t = args.t
boot_strategy = args.boot_strategy

# boot_strategy = 'jacknife'
# solver = 'UPMGA'
# metric = 'hamming'
# ncores = 8
# n_boot = 100
# t = .1

########################################################################

# Preparing run: import code, prepare directories, set logger
from mito_utils.preprocessing import *
from mito_utils.phylo import *

########################################################################

# Main
def main():

    # Read AFM and filter vars
    path_data = '/Users/IEO5505/Desktop/mito_bench/data'
    sample = 'MDA_clones'
    afm = read_one_sample(path_data, sample, with_GBC=True)

    if frac_subsampling is not None: 
        afm = filter_cells_coverage(afm)
        cells = (
            afm.obs.groupby('GBC')
            .apply(lambda x: x.sample(frac=frac_subsampling))
            .droplevel(0)
            .index
        )
        afm = remove_excluded_sites(afm[cells,:].copy())
    if filtering == 'GT':
        _, a = filter_afm_with_gt(afm, min_cells_clone=5)
    else:
        _, a = filter_cells_and_vars(
            afm, filtering=filtering, path_=os.getcwd(), nproc=ncores
        )

    # Remove zeros and get AD, DP
    a = nans_as_zeros(a)
    AD, DP, _ = get_AD_DP(a)

    TREE_D = {}
    kwargs = {'solver':solver, 'metric':metric, 'boot_strategy':boot_strategy}
    for k in [ f'rep_{i}' for i in range(n_boot) ]:
        tree = bootstrap_iteration(a, AD, DP, kwargs=kwargs)
        TREE_D[k] = tree

    # Support
    obs_tree = build_tree(a, t=t)
    df, _ = calculate_support(obs_tree, TREE_D.values(), n=n_boot)
    # RF
    rf_l = []
    for k in TREE_D:
        d, max_d = cs.critique.robinson_foulds(obs_tree, TREE_D[k])
        rf_l.append(d/max_d)

    # Describe
    print('\n')
    print(a)
    print('\n')
    print('Support \n')
    print(df['support'].describe())
    print('\n')
    print('Support by time\n')
    print(df.groupby('time')['support'].describe())
    print('\n')
    print('RF \n')
    print(pd.Series(rf_l).describe())


##


# Run
if __name__ == '__main__':
    main()
    