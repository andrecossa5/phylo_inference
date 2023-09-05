"""
Tree building and consistency checks.
"""

# Code
import cassiopeia as cs
from matplotlib.colors import ListedColormap
from mito_utils.preprocessing import *
from mito_utils.utils import *
from mito_utils.distances import *
from mito_utils.it_diagnostics import *
from mito_utils.phylo import *
from mito_utils.plotting_base import *
matplotlib.use('macOSX')


# Args
path_main = '/Users/IEO5505/Desktop/mito_bench'
sample = 'MDA_clones'


##


# Paths
path_data = os.path.join(path_main, 'data') 
path_output = os.path.join(path_main, 'results', 'phylo', 'output')
path_viz = os.path.join(path_main, 'results', 'phylo', 'visualization')
path_tmp = os.path.join(path_main, 'results', 'phylo', 'downstream_files')


##

# Load afm and filter with GT
afm = read_one_sample(path_data, sample, with_GBC=True)
a_cells, a = filter_afm_with_gt(afm, min_cells_clone=5)

# Create allele tables
AD, DP, ad_vars = get_AD_DP(a)

# Bootstrap trees
tree_list = []
for i in range(100):

    print(f'Bootstrap sample {i}')

    AD_, DP_ = bootstrap_allele_tables(AD.A.T, DP.A.T)
    X = np.divide(AD_, DP_)
    tree = build_tree(a, X=X)
    tree.cell_meta = pd.DataFrame(
        a.obs.loc[tree.leaves, "GBC"].astype(str)
    )
    tree.cell_meta[a.var_names] = a[tree.leaves, :].X.toarray()

    tree_list.append(tree)


##


# Plot some bootstrapped trees
# fig = plt.figure(figsize=(15,15))
# 
# for i in range(16):
# 
#     tree = tree_list[i]
#     ax = fig.add_subplot(4,4,i+1)
#     cs.pl.plot_matplotlib(
#         tree, meta_data=[a.var_names[0], "GBC"], 
#         categorical_cmap=ListedColormap(sc.pl.palettes.godsnot_102), 
#         continuous_cmap='mako',
#         add_root=True,
#         ax=ax
#     )
#     format_ax(ax, title=f'{sample}, bootstrap {i+1}')
# 
# # Save fig
# fig.tight_layout()
# fig.savefig(os.path.join(path_viz, 'bootstrap.png'), dpi=300)


##


# Calculate observed branches supports
obs_tree = build_tree(a, X=None)
obs_tree


## ....


##





















## 




# Tree analysis
# probability_threshold = 0.05
# cs.tl.compute_expansion_pvalues(tree, min_clade_size=(0.05 * tree.n_cell), min_depth=6)
# expanding_nodes = []
# for node in tree.depth_first_traverse_nodes():
#     if tree.get_attribute(node, "expansion_pvalue") < probability_threshold:
#         expanding_nodes.append(node)
# 
# # Association tests
# cs.tl.compute_morans_i(tree, meta_columns=[a.var_names[1]])
# 
# # Perform the phylogenetic ANOVA (one-way ANOVA)
# from scipy.stats import f_oneway
# anova_result = f_oneway(*tree.cell_meta['GBC'].astype('category').cat.codes)
# 
# # Print the F-statistic and p-value
# print("F-statistic:", anova_result.statistic)
# print("p-value:", anova_result.pvalue)
