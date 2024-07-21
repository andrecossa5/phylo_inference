"""
Visualize top trees.
"""

import os
import pickle
from scipy.sparse import load_npz
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.phylo_plots import *
from mito_utils.heatmaps_plots import *
from mito_utils.plotting_base import *
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score
matplotlib.use('macOSX')


##



def _annotate_tree(trees, nodes, key=None): 

    # Further annotate the cell_meta
    df_nodes = (
        nodes[key].loc[lambda x: ~x['cell'].isna()]
        [['cell', 'clonal_ancestor']]
        .drop_duplicates().set_index('cell')
    )
    df_nodes['clonal_ancestor'] = df_nodes['clonal_ancestor'].astype('int').astype('str')
    tree = trees[key].copy()
    tree.cell_meta = tree.cell_meta.join(df_nodes)

    # Further annotate internal nodes
    df_nodes = nodes[key].loc[lambda x: ~x['assigned_var'].isna() & x['cell'].isna()]
    df_nodes['node'] = df_nodes['node'].astype('str')
    n_vars = df_nodes.groupby('node').size()
    for node in tree.depth_first_traverse_nodes():
        nvar = n_vars.loc[node] if node in n_vars.index else 0
        tree.set_attribute(node, attribute_name='nvar', value=nvar)
    
    return tree, n_vars


##

# Paths
path_main = '/Users/IEO5505/Desktop/mito_bench/'
path_results = os.path.join(path_main, 'results/phylo')

# Load trees and nodes annot
trees = traverse_and_extract_flat(path_results, file_name='annotated_tree.pickle')
nodes = traverse_and_extract_flat(path_results, file_name='nodes.csv')

# Colors 
with open(os.path.join(path_main, 'data', 'clones_colors_sc.pickle'), 'rb') as f:
    clone_colors = pickle.load(f)


##


# Fig 1
samples = ['AML_clones', 'MDA_clones', 'MDA_PT', 'MDA_lung']
variant_set = 'GT_enriched'

fig, axs = plt.subplots(1,4,figsize=(15,4))

for ax, sample in zip(axs.ravel(), samples):

    tree, _ = _annotate_tree(trees, nodes, key=(sample, variant_set))
    ARI = adjusted_mutual_info_score(tree.cell_meta['clonal_ancestor'], tree.cell_meta['GBC'])
    NMI = adjusted_rand_score(tree.cell_meta['clonal_ancestor'], tree.cell_meta['GBC'])
    mt_cmap = create_palette(tree.cell_meta, 'clonal_ancestor', sc.pl.palettes.godsnot_102)

    plot_tree(
       tree, ax=ax, meta=['clonal_ancestor'],
        colorstrip_width=1.5,
        colorstrip_spacing=0.001,
        extend_branches=True,
        orient=90,
        categorical_cmaps={'clonal_ancestor':mt_cmap}, # 'GBC':clone_colors, 
        cov_internal_nodes='support', # 'nvar',
        # internal_node_subset=n_vars.index,
        internal_node_kwargs={'markersize':5, 'markeredgecolor':'k'},
        internal_node_vmin=40, internal_node_vmax=70, 
    )
    ax.set(title=f'ncells: {len(tree.leaves)}, NMI: {NMI:.2f}, ARI: {ARI:.2f}')

fig.tight_layout()
fig.savefig(os.path.join(path_results, f'{variant_set}_trees_support.png'), dpi=1000)


## 


# MDA_PT only
fig, ax = plt.subplots(figsize=(6,6))

sample = 'MDA_PT'
variant_set = 'set1'

# Annotate with vars
tree, n_vars = _annotate_tree(trees, nodes, key=(sample, variant_set))
AD = load_npz(os.path.join(path_results, sample, variant_set, 'filtered_input', 'AD.npz')).A
DP = load_npz(os.path.join(path_results, sample, variant_set, 'filtered_input', 'DP.npz')).A
cells = pd.read_csv(os.path.join(path_results, sample, variant_set, 'filtered_input', 'meta.csv'), index_col=0).index
variants = pd.read_csv(os.path.join(path_results, sample, variant_set, 'filtered_input', 'variants.csv'), index_col=0, header=None).index
afm = pd.DataFrame(AD+0.0001/DP, columns=variants, index=cells)
vois = nodes[(sample, variant_set)].loc[lambda x: ~x['assigned_var'].isna() & x['cell'].isna()]['assigned_var'].map(lambda x: x.replace('X','').replace('.', '>')).values
tree.cell_meta = tree.cell_meta.join(afm.loc[:,vois])

ARI = adjusted_mutual_info_score(tree.cell_meta['clonal_ancestor'], tree.cell_meta['GBC'])
NMI = adjusted_rand_score(tree.cell_meta['clonal_ancestor'], tree.cell_meta['GBC'])
mt_cmap = create_palette(tree.cell_meta, 'clonal_ancestor', 'tab10')
lenti_cmap = create_palette(tree.cell_meta, 'GBC', 'dark')

plot_tree(
   tree, ax=ax, meta=['clonal_ancestor', 'GBC'] + vois.tolist(),
    colorstrip_width=1.5,
    colorstrip_spacing=0.001,
    extend_branches=True,
    orient='down',
    categorical_cmaps={'clonal_ancestor':mt_cmap, 'GBC':lenti_cmap}, 
    continuous_cmap='mako',
    vmin_annot=0, vmax_annot=.1,
    cov_internal_nodes='nvar',
    internal_node_subset=n_vars.index,
    internal_node_kwargs={'markersize':5, 'markeredgecolor':'k', 'c':'darkred'},
    internal_node_vmin=0, internal_node_vmax=1, 
)
ax.set(title=f'n cells: {len(tree.leaves)}, n assigned MT-SNVs: {vois.size} \n NMI: {NMI:.2f}, ARI: {ARI:.2f}')

fig.tight_layout()
fig.savefig(os.path.join(path_results, f'{sample}_{variant_set}_tree_vars.png'), dpi=1000)


##

##