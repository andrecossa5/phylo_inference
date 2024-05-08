"""
Visualiza phylo output.
"""

import os
import pickle
from scipy.sparse import load_npz
from anndata import AnnData
from mito_utils.preprocessing import *
from mito_utils.phylo_plots import *
from mito_utils.heatmaps_plots import *
from mito_utils.plotting_base import *
matplotlib.use('macOSX')


##


def prep_sample(sample):

    # Filter top
    d = (
        # df
        df.query('filtering!="GT"')
        #.query('sample==@sample and time<=10 and time>=5')
        .groupby(['sample', 'filtering', 'bootstrap', 'solver', 'metric'])
        ['TBE'].median().sort_values(ascending=False)
        .reset_index()
        .iloc[0,:].to_dict()
    )
    print(d)

    # SET PARAMs
    filtering = d['filtering']
    solver = d['solver']
    metric = d['metric']             # default, cosine
    bootstrap = d['bootstrap']
    path_data = os.path.join(path_results, 'output', sample, filtering, 'input_folder')
    path_trees = os.path.join(
        path_results, 'output', sample, filtering, 
        solver, metric, bootstrap, 'trees.pickle'
    )

    # Read tree
    with open(path_trees, 'rb') as f:
        trees = pickle.load(f)

    # Get observed tree
    obs_tree = trees['observed']

    # Read cells, variants and build the afm
    variants = pd.read_csv(os.path.join(path_data, 'variants.csv'), index_col=0, header=None)
    meta = pd.read_csv(os.path.join(path_data, 'meta.csv'), index_col=0)
    AD = load_npz(os.path.join(path_data, 'AD.npz')).astype(np.int16).A
    DP = load_npz(os.path.join(path_data, 'DP.npz')).astype(np.int16).A
    afm = AnnData(AD/DP, obs=meta, var=variants)

    # Create annot dfs
    obs_tree.cell_meta = pd.DataFrame(
        afm.obs.loc[obs_tree.leaves, "GBC"].astype('str')
        ).join(
            pd.DataFrame(
            afm[obs_tree.leaves,:].X, index=obs_tree.leaves, columns=afm.var_names)
    )
    # obs_tree.cell_meta['GBC'].value_counts(normalize=True)
    df_supports = (
        df.query('sample==@sample and filtering==@filtering and solver==@solver and bootstrap==@bootstrap and metric==@metric')
        .loc[obs_tree.internal_nodes, ['TBE', 'FBP', 'median_RF']]
    )

    return afm, obs_tree, df_supports, variants


##


# Paths
path_main = '/Users/IEO5505/Desktop/mito_bench/'
path_results = os.path.join(path_main, 'results/phylo_inference')
path_top_trees = os.path.join(path_results, 'top_trees_MQuad')
samples = ['AML_clones', 'MDA_clones', 'MDA_PT', 'MDA_lung']

# Plot one good tree example per sample
df = pd.read_csv(os.path.join(path_results, 'output/supports_df.csv'), index_col=0)

# Save individual sample top tree for PATH analysis
# cov = 'GBC'
# for sample in samples:
#     make_folder(path_top_trees, sample, overwrite=False)
#     path_sample = os.path.join(path_top_trees, sample)
#     afm, obs_tree, df_supports, variants = prep_sample(sample)
#     counts = afm.obs['GBC'].value_counts().to_frame('n_cells')
#     counts['to_numeric'] = range(1,counts.shape[0]+1)
#     GT = afm.obs['GBC'].map(counts['to_numeric'].to_dict()).to_frame('GBC')
#     first_to_discard = counts.query('n_cells<10')['to_numeric'].min()
#     counts[counts['to_numeric']>=first_to_discard] = first_to_discard
#     counts.to_csv(os.path.join(path_sample, f'mapping_{cov}.csv'))
#     GT.loc[GT['GBC']>=first_to_discard] = first_to_discard
#     GT.to_csv(os.path.join(path_sample, f'cov_df_{cov}.csv'))
#     with open(os.path.join(path_sample, 'tree.newick'), 'w') as f:
#         f.write(obs_tree.get_newick(record_branch_lengths=True))


##

# Colors 
with open(os.path.join(path_main, 'data', 'clones_colors_sc.pickle'), 'rb') as f:
    colors = pickle.load(f)

# Viz 4 of them
vmin_annot = .001
vmax_annot = .1

fig, axs = plt.subplots(1,4,figsize=(14,4.5))

for ax, sample in zip(axs, samples):
    afm, obs_tree, df_supports, variants = prep_sample(sample)
    sample_colors = { k:colors[k] for k in colors if k in afm.obs['GBC'].values }
    plot_tree(
       obs_tree, ax=ax, meta=['GBC'],
        colorstrip_width=4,
        colorstrip_spacing=0.001,
        extend_branches=True,
        orient=90,
        categorical_cmap_annot=sample_colors,
        meta_internal_nodes=df_supports,
        cov_internal_nodes='TBE',
        internal_node_kwargs={'markersize':3, 'markeredgecolor':None},
        internal_node_vmin=.5, internal_node_vmax=.9, 
    )
    format_ax(
        ax=ax, 
        title=f'''
        {sample} \n TBE {df_supports['TBE'].median():.2f} (+-{df_supports['TBE'].std():.2f}), FBP {df_supports['FBP'].median():.2f} (+-{df_supports['FBP'].std():.2f})
        ''',
        title_size=11
    )

# Format and save
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'top_trees_MQuad.pdf'), dpi=500)


##


# One at the time with variants
sample = 'MDA_clones'
vmin_annot = .001
vmax_annot = .1

fig, ax = plt.subplots(figsize=(7,6.5))
afm, obs_tree, df_supports, variants = prep_sample(sample)
sample_colors = { k:colors[k] for k in colors if k in afm.obs['GBC'].values }
plot_tree(
   obs_tree, ax=ax, meta=['GBC', *variants.index.to_list()],
    colorstrip_width=1,
    colorstrip_spacing=0.001,
    orient='down',
    categorical_cmap_annot=sample_colors,
    extend_branches=True,
    meta_internal_nodes=df_supports,
    cov_internal_nodes='TBE',
    internal_node_kwargs={'markersize':4, 'markeredgecolor':None},
    internal_node_vmin=.5, internal_node_vmax=.9, 
)
format_ax(ax=ax, title=sample, title_size=11)
add_cbar(
    df_supports['TBE'], palette='Spectral_r', ax=ax, label='TBE',
    ticks_size=8, label_size=9, vmin=.5, vmax=.9, 
    layout=( (1.03,.25,.03,.5), 'right', 'vertical' )
)
add_cbar(
    afm.X.flatten(), palette='mako', ax=ax, label='AF',
    ticks_size=8, label_size=9, layout=( (1.2,.25,.03,.5), 'right', 'vertical' ),
    vmax=vmax_annot, vmin=vmin_annot
)

# Format and save
# fig.subplots_adjust(left=.05, right=.7, top=.8, bottom=.1)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'{sample}.pdf'), dpi=500)


##





