#!/usr/bin/python

##  Clones classification script

########################################################################

# Libraries
import sys 
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='phylo_cassiopeia',
    description=
        """
        Systematically testing the ability of phylogenetic reconstruction algorithms 
        to derive good MT-phylogenies.
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

# Metric
my_parser.add_argument(
    '--metric', 
    type=str,
    default='cosine',
    help='Metric used for cell-cell distance matrix computation. Default: cosine.'
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

# n boot
my_parser.add_argument(
    '--n_boot', 
    type=int,
    default=100,
    help='n bootstrap samples for support calculation. Default: 100.'
)

# ncombos
my_parser.add_argument(
    '--boot_strategy', 
    type=str,
    default='jacknife',
    help=
        '''
        Bootstrap resampling strategy. Default: jacknife. Other options: 
        features_resampling and counts_resampling
        '''
)

# ncombos
my_parser.add_argument(
    '--subsample', 
    type=int,
    default=None,
    help='N of cells to resample, if necessary. Default: None.'
)

# ncores
my_parser.add_argument(
    '--ncores', 
    type=int,
    default=8,
    help='ncores to use for model training. Default: 8.'
)

# skip
my_parser.add_argument(
    '--skip', 
    action='store_true',
    help='Skip analysis. Default: False.'
)

# Parse arguments
args = my_parser.parse_args()

path_data = args.path_data
sample = args.sample
filtering = args.filtering
metric = args.metric
solver = args.solver
n_boot = args.n_boot
boot_strategy = args.boot_strategy
subsample = args.subsample 
ncores = args.ncores 

# path_main = '/Users/IEO5505/Desktop/mito_bench'
# sample = 'MDA_clones'
# filtering = 'GT'
# metric = 'cosine'
# solver = 'NJ'
# n_boot = 50
# boot_strategy = 'jacknife'
# subsample = 0
# ncores = 8
# matplotlib.use('macOSX')


########################################################################

# Preparing run: import code, prepare directories, set logger
if not args.skip:

    # Code
    import cassiopeia as cs
    from mito_utils.phylo_plots import *
    from mito_utils.diagnostic_plots import sturges
    from mito_utils.preprocessing import *
    from mito_utils.utils import *
    from mito_utils.distances import *
    from mito_utils.it_diagnostics import *
    from mito_utils.phylo import *
    from mito_utils.plotting_base import *

    # Job key
    job_key = f'{sample}_{filtering}_{metric}_{solver}_{boot_strategy}_{n_boot}'
    
    # Set paths and folders
    # path_data = os.path.join(path_data, 'data') 
    # path_output = os.path.join(path_main, 'results', 'phylo', 'output')
    # path_viz = os.path.join(path_main, 'results', 'phylo', 'visualization')
    # for x in [path_output, path_viz]:
    #     make_folder(x, sample, overwrite=False)
    #     make_folder(os.path.join(x, sample), job_key, overwrite=True)

    os.mkdir(job_key)
    os.chdir(job_key)
    os.mkdir('output')
    os.mkdir('viz')
    
    # Update
    path_output = os.path.join(os.getcwd(), 'output')
    path_viz = os.path.join(os.getcwd(), 'viz')

    # Set logger
    logger = set_logger(path_output, f'log.txt')

########################################################################

# Main
def main():

    T = Timer()
    T.start()

    # Load data
    t = Timer()
    t.start()

    logger.info(
        f""" 
        Execute phylogeny inference: \n
        --sample {sample} 
        --filtering {filtering} 
        --metric {metric} 
        --solver {solver}
        --boot_strategy {boot_strategy} 
        --subsample {subsample}
        --n_boot {n_boot}
        """
    )

    # Read AFM and filter 
    afm = read_one_sample(path_data, sample, with_GBC=True)

    if filtering == 'GT':
        _, a = filter_afm_with_gt(afm, min_cells_clone=5)
    else:
        _, a = filter_cells_and_vars(afm, filtering=filtering, path_=path_output)

    # Subsample if necessary
    if subsample is not None and subsample>0:
        subsample_cells = a.obs.sample(subsample).index
        _, a = filter_cells_and_vars(a, cells=subsample_cells)

    # AD,DP
    a = nans_as_zeros(a)
    AD, DP, ad_vars = get_AD_DP(a)

    # Plot selected variants matrix aggregated by clone
    fig, ax = plt.subplots(figsize=(5,5))
    X = (
        a.obs
        .join(pd.DataFrame(a.X, index=a.obs_names, columns=a.var_names))
        .groupby('GBC').agg('median')
    )
    plot_heatmap(
        X, ax=ax, title=f'{filtering}, {sample}',
        label='AF', rank_diagonal=True, vmin=.01, vmax=.1
    )
    format_ax(
        ax, title=f'{sample}, {filtering} filtering', 
        xlabel='Variants', ylabel='Clones', xticks=[], yticks=[]
    )
    fig.tight_layout()
    fig.savefig(os.path.join(path_viz, 'aggregated_AF.png'), dpi=300)
    
    logger.info(f'Cells and variants filtered: {t.stop()}')


    ##


    # Bootstrap trees
    logger.info('Starting bootstrap...')

    kwargs = {
        'metric' : metric, 
        'solver' : solver,
        'solver_kwargs' : {}
    }
    boot_kwargs = {**kwargs, **{'boot_strategy' : boot_strategy}}

    tree_list = []
    for i in range(n_boot):
        t.start()
        tree = bootstrap_iteration(a, AD, DP, kwargs=boot_kwargs.copy())
        print(f'Bootstrap sample {i}')
        tree_list.append(tree)
        logger.info(f'Bootstrap sample {i}: {t.stop()}')


    ##


    # Calculate observed branches supports
    t.start()

    obs_tree = build_tree(a, X=None, **kwargs)
    supports_df, bootstrap_record = calculate_support(obs_tree, tree_list, n=len(tree_list))
    supports_df.describe()
    bootstrap_record.describe()

    # Most expanded clones and more abundant clones support
    exp_clones = get_expanded_clones(obs_tree, t=.05, min_depth=3)
    top_n = supports_df.sort_values('n_cells', ascending=False).index[:25]
    supports_df['expanded_clones'] = supports_df.index.map(
        lambda x: 'expanded' if x in exp_clones else 'other'
    )
    supports_df.loc[exp_clones].sort_values('support', ascending=False)
    supports_df.loc[top_n].sort_values('support', ascending=False)

    # Add leaves metadata and save dfs
    obs_tree.cell_meta = pd.DataFrame(
        a.obs.loc[obs_tree.leaves, "GBC"].astype(str)
    )
    obs_tree.cell_meta[a.var_names] = a[obs_tree.leaves, :].X.toarray()

    # Save outputs
    supports_df.to_csv(os.path.join(path_output, 'supports_df.csv'))
    bootstrap_record.to_csv(os.path.join(path_output, 'bootstrap_record.csv'))
    with open(os.path.join(path_output, 'tree.pickle'), 'wb') as f:
        pickle.dump(obs_tree, f)

    logger.info(f'Finished support calculation and saving data: {t.stop()}')


    ##

    # Plot some boot-trees per variant
    plot_boot_trees(a, tree_list, a.var_names[:5], path_viz)

    ##

    # Boot support and clade type diagnostics
    t.start()

    fig = plot_bootstrap_relationships(supports_df, figsize=(10,5)) 
    fig.savefig(os.path.join(path_viz, 'bootstrap_relationship.png'), dpi=400)
    
    fig = plot_bootstrap_record(bootstrap_record, figsize=(6,5))
    fig.savefig(os.path.join(path_viz, f'bootstrap_record.png'), dpi=400)

    ##

    # Main phylo plots
    fig = plot_main(a, obs_tree, supports_df, a.var_names, bootstrap_record, sample)
    fig.savefig(os.path.join(path_viz, f'main.png'), dpi=400)

    fig = plot_main_support(obs_tree, supports_df)
    fig.savefig(os.path.join(path_viz, f'main_supports.png'), dpi=400)

    logger.info(f'Finished plotting: {t.stop()}')

    ##


    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    if not args.skip:
        main()

#######################################################################









