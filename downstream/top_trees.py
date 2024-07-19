"""
Visualize top trees.
"""

import os
import pickle
from scipy.sparse import load_npz
from anndata import AnnData
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.phylo_plots import *
from mito_utils.heatmaps_plots import *
from mito_utils.plotting_base import *
matplotlib.use('macOSX')


##


# def prep_sample(sample):
# 
#     # Filter top
#     d = (
#         # df
#         df.query('filtering!="GT"')
#         #.query('sample==@sample and time<=10 and time>=5')
#         .groupby(['sample', 'filtering', 'bootstrap', 'solver', 'metric'])
#         ['TBE'].median().sort_values(ascending=False)
#         .reset_index()
#         .iloc[0,:].to_dict()
#     )
#     print(d)
# 
#     # SET PARAMs
#     filtering = d['filtering']
#     solver = d['solver']
#     metric = d['metric']             # default, cosine
#     bootstrap = d['bootstrap']
#     path_data = os.path.join(path_results, 'output', sample, filtering, 'input_folder')
#     path_trees = os.path.join(
#         path_results, 'output', sample, filtering, 
#         solver, metric, bootstrap, 'trees.pickle'
#     )
# 
#     # Read tree
#     with open(path_trees, 'rb') as f:
#         trees = pickle.load(f)
# 
#     # Get observed tree
#     obs_tree = trees['observed']
# 
#     # Read cells, variants and build the afm
#     variants = pd.read_csv(os.path.join(path_data, 'variants.csv'), index_col=0, header=None)
#     meta = pd.read_csv(os.path.join(path_data, 'meta.csv'), index_col=0)
#     AD = load_npz(os.path.join(path_data, 'AD.npz')).astype(np.int16).A
#     DP = load_npz(os.path.join(path_data, 'DP.npz')).astype(np.int16).A
#     afm = AnnData(AD/DP, obs=meta, var=variants)
# 
#     # Create annot dfs
#     obs_tree.cell_meta = pd.DataFrame(
#         afm.obs.loc[obs_tree.leaves, "GBC"].astype('str')
#         ).join(
#             pd.DataFrame(
#             afm[obs_tree.leaves,:].X, index=obs_tree.leaves, columns=afm.var_names)
#     )
#     # obs_tree.cell_meta['GBC'].value_counts(normalize=True)
#     df_supports = (
#         df.query('sample==@sample and filtering==@filtering and solver==@solver and bootstrap==@bootstrap and metric==@metric')
#         .loc[obs_tree.internal_nodes, ['TBE', 'FBP', 'median_RF']]
#     )
# 
#     return afm, obs_tree, df_supports, variants


##


# Paths
path_main = '/Users/IEO5505/Desktop/mito_bench/'
path_results = os.path.join(path_main, 'results/phylo')


samples = ['AML_clones', 'MDA_clones', 'MDA_PT', 'MDA_lung']

trees = traverse_and_extract_flat(path_results, file_name='annotated_tree.pickle')
nodes = traverse_and_extract_flat(path_results, file_name='nodes.csv')


##

# Colors 
with open(os.path.join(path_main, 'data', 'clones_colors_sc.pickle'), 'rb') as f:
    colors = pickle.load(f)

# Viz 4 of them
vmin_annot = .001
vmax_annot = .1



key = ('MDA_lung', 'MQuad')



# Further annotate the tree
df_nodes = (
    nodes[key].loc[lambda x: ~x['cell'].isna()]
    [['cell', 'clonal_ancestor', 'clone']]
    .drop_duplicates().set_index('cell')
)
df_nodes['clone'] = df_nodes['clone'].astype('str')
tree = trees[key]
tree.cell_meta = tree.cell_meta.join(df_nodes)


from sklearn.metrics import normalized_mutual_info_score, adjusted_rand_score

normalized_mutual_info_score(tree.cell_meta['clone'], tree.cell_meta['GBC'])
adjusted_rand_score(tree.cell_meta['clone'], tree.cell_meta['GBC'])

lenti_cmap = create_palette(tree.cell_meta, 'GBC', 'tab20')
mt_cmap = create_palette(tree.cell_meta, 'clone', sc.pl.palettes.godsnot_102)

fig, ax = plt.subplots(figsize=(6,6))
plot_tree(
   tree, ax=ax, meta=['clone', 'GBC'],
    colorstrip_width=3,
    colorstrip_spacing=0.001,
    extend_branches=True,
    orient=90,
    categorical_cmaps={'GBC':lenti_cmap, 'clone':mt_cmap},
    cov_internal_nodes='support',
    internal_node_kwargs={'markersize':3, 'markeredgecolor':None},
    internal_node_vmin=.5, internal_node_vmax=.9, 
)
plt.show()


tree.cell_meta.shape




































##


# # One at the time with variants
# sample = 'MDA_clones'
# vmin_annot = .001
# vmax_annot = .1
# 
# fig, ax = plt.subplots(figsize=(7,6.5))
# afm, obs_tree, df_supports, variants = prep_sample(sample)
# sample_colors = { k:colors[k] for k in colors if k in afm.obs['GBC'].values }
# plot_tree(
#    obs_tree, ax=ax, meta=['GBC', *variants.index.to_list()],
#     colorstrip_width=1,
#     colorstrip_spacing=0.001,
#     orient='down',
#     categorical_cmap_annot=sample_colors,
#     extend_branches=True,
#     meta_internal_nodes=df_supports,
#     cov_internal_nodes='TBE',
#     internal_node_kwargs={'markersize':4, 'markeredgecolor':None},
#     internal_node_vmin=.5, internal_node_vmax=.9, 
# )
# format_ax(ax=ax, title=sample, title_size=11)
# add_cbar(
#     df_supports['TBE'], palette='Spectral_r', ax=ax, label='TBE',
#     ticks_size=8, label_size=9, vmin=.5, vmax=.9, 
#     layout=( (1.03,.25,.03,.5), 'right', 'vertical' )
# )
# add_cbar(
#     afm.X.flatten(), palette='mako', ax=ax, label='AF',
#     ticks_size=8, label_size=9, layout=( (1.2,.25,.03,.5), 'right', 'vertical' ),
#     vmax=vmax_annot, vmin=vmin_annot
# )
# 
# # Format and save
# # fig.subplots_adjust(left=.05, right=.7, top=.8, bottom=.1)
# fig.tight_layout()
# fig.savefig(os.path.join(path_results, f'{sample}.pdf'), dpi=500)
# 
# 
# ##





