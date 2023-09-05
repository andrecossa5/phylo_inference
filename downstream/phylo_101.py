"""
Tree building with cassiopea
"""

# Code
import os
import warnings
from sklearn.metrics import normalized_mutual_info_score
import matplotlib
from matplotlib.colors import ListedColormap
warnings.filterwarnings('ignore')

from mito_utils.preprocessing import *
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.dimred import *
from mito_utils.distances import *
from mito_utils.clustering import *
from mito_utils.it_diagnostics import *
from mito_utils.it_iterations import *
from mito_utils.plotting_base import *
from mito_utils.embeddings_plots import *
from mito_utils.heatmaps_plots import plot_heatmap
matplotlib.use('macOSX')


# Args
path_main = '/Users/IEO5505/Desktop/mito_bench'
sample = 'MDA_PT'

# Paths
path_data = os.path.join(path_main, 'data') 
path_output = os.path.join(path_main, 'results', 'unsupervised_clones', 'output')
path_viz = os.path.join(path_main, 'results', 'unsupervised_clones', 'visualization')
path_tmp = os.path.join(path_main, 'results', 'unsupervised_clones', 'downstream_files')


##

# Load
afm = read_one_sample(path_data, sample, with_GBC=True)
a_cells = filter_cells_coverage(afm)
a_cells = filter_baseline(a_cells)

t = .75
gt_l = [
    rank_clone_variants(
        a_cells, var='GBC', group=g, rank_by='custom_perc_tresholds',
        min_clone_perc=t, max_perc_rest=.25
    ).assign(clone=g)
    for g in a_cells.obs['GBC'].unique()
]
df_gt = pd.concat(gt_l).join(summary_stats_vars(a_cells))

# Filter 
vois_df = (
    df_gt
    .query('n_cells_clone>10')
    .sort_values('log2_perc_ratio', ascending=False)
    .loc[:, 
        [
            'median_AF_clone', 'median_AF_rest', 'perc_clone', 
            'perc_rest', 'log2_perc_ratio', 'n_cells_clone', 'clone'
        ]
    ]
)
vois = vois_df.index.unique()
cells = a_cells.obs['GBC'].loc[lambda x: x.isin(vois_df['clone'])].index


# UMAP and visualization
a_good = filter_cells_and_vars(a_cells, cells=cells, variants=vois)[0]
a_good = nans_as_zeros(a_good)
X_umap = reduce_dimensions(a_good, 'UMAP', metric='cosine', n_comps=2)

df_ = (
    pd.DataFrame(X_umap[0], columns=X_umap[1], index=cells)
    .join(a_cells.obs)
)
df_['GBC'] = df_['GBC'].astype('str')

fig, ax = plt.subplots(figsize=(4.8,5))
draw_embeddings(
    df_, cat='GBC', ax=ax, 
    axes_kwargs={'legend':False},
    title=f'MDA_PT, only recoverable cells and GT variants ({100 * cells.size / a_cells.shape[0]:.2f}%)',
    legend_kwargs={'loc':'center left', 'bbox_to_anchor':(1,.5)}
)
fig.tight_layout()
plt.show()


##


# UMAP and visualization
a_miller = filter_cells_and_vars(a_cells, filtering='miller2022')[0]
a_miller = nans_as_zeros(a_miller)
X_umap = reduce_dimensions(a_miller, 'UMAP', metric='cosine', n_comps=2)

df_ = (
    pd.DataFrame(X_umap[0], columns=X_umap[1], index=a_miller.obs_names)
    .join(a_cells.obs)
)
df_['GBC'] = df_['GBC'].astype('str')

fig, ax = plt.subplots(figsize=(4.8,5))
draw_embeddings(
    df_, cat='GBC', ax=ax, 
    axes_kwargs={'legend':False},
    title=f'MDA_PT, all cells, miller2022 variants',
    legend_kwargs={'loc':'center left', 'bbox_to_anchor':(1,.5)}
)
fig.tight_layout()
plt.show()
















import cassiopeia as cs


M = pd.DataFrame(
    np.where(a_good.X>.1, 1, 0),
    index=a_good.obs_names,
    columns=a_good.var_names
)

D = pd.DataFrame(
    pair_d(a_good, metric='cosine'),
    index=a_good.obs_names,
    columns=a_good.obs_names
)

tree = cs.data.CassiopeiaTree(
    character_matrix=M, 
    dissimilarity_map=D, 
    cell_meta=a_good.obs
)
solver = cs.solver.UPGMASolver()
# solver = cs.solver.NeighborJoiningSolver(add_root=True)
# solver = cs.solver.ILPSolver()
# solver = cs.solver.VanillaGreedySolver()
solver.solve(tree)



dir(tree)
tree

tree.cell_meta = pd.DataFrame(
    a_cells.obs.loc[tree.leaves, "GBC"].astype(str)
)
to_plot = vois_df.sort_values('n_cells_clone', ascending=False).index[:6]
tree.cell_meta[to_plot] = a_good[tree.leaves, to_plot].X.toarray()

fig, ax = plt.subplots(figsize=(8,8))
cs.pl.plot_matplotlib(
    tree, meta_data=[to_plot[0], "GBC"], 
    categorical_cmap=ListedColormap(sc.pl.palettes.godsnot_102), 
    continuous_cmap='mako',
    add_root=True,
    ax=ax
)
format_ax(ax, title=f'Jaccard distance, ILP')
# format_ax(ax, title=f'MDA_PT, all cells, miller2022 variants')
plt.show()


len(tree.leaves) / a_cells.shape[0]

















































from Bio.Phylo.Consensus import bootstrap









probability_threshold = 0.05

cs.tl.compute_expansion_pvalues(tree, min_clade_size=(0.05 * tree.n_cell), min_depth=6)
expanding_nodes = []
for node in tree.depth_first_traverse_nodes():
    if tree.get_attribute(node, "expansion_pvalue") < probability_threshold:
        expanding_nodes.append(node)


make_folder(path_viz, 'trees', overwrite=True)

for i in range(len(expanding_nodes)):
    fig, ax = plt.subplots(figsize=(9,9))
    cs.pl.plot_matplotlib(
        tree, meta_data=["GBC"], 
        categorical_cmap=ListedColormap(sc.pl.palettes.godsnot_102), 
        clade_colors={expanding_nodes[i]: "red"},
        add_root=True,
        ax=ax
    )
    fig.tight_layout()
    fig.savefig(os.path.join(path_viz, 'trees', f'clone_{i}.png'))




# Association


cs.tl.compute_morans_i(tree, meta_columns=[to_plot[1]])



from scipy.stats import f_oneway


# Perform the phylogenetic ANOVA (one-way ANOVA)
anova_result = f_oneway(*tree.cell_meta['GBC'].astype('category').cat.codes)

# Print the F-statistic and p-value
print("F-statistic:", anova_result.statistic)
print("p-value:", anova_result.pvalue)




