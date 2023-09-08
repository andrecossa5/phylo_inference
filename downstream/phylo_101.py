"""
Tree building and consistency checks.
"""

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
a_cells, a = filter_cells_and_vars(afm, filtering='miller2022', path_=os.getcwd())

# a_cells, a = filter_afm_with_gt(afm, min_cells_clone=5)

# Create AD, DP tables
# subsample_cells = a.obs.groupby('GBC').apply(lambda x: x.sample(5)).droplevel(0).index
# subsample_cells = a.obs.sample(200).index
# _, a = filter_cells_and_vars(a, cells=subsample_cells)

a = nans_as_zeros(a)
AD, DP, ad_vars = get_AD_DP(a)

# AD.A.T
# fig, ax = plt.subplots()
# X = (
#     a.obs
#     .join(pd.DataFrame(a.X, index=a.obs_names, columns=a.var_names))
#     .groupby('GBC').agg('median')
#     #.columns
# )
# ax.imshow(X)
# add_cbar(ax=ax, x=X.values.flatten(), vmin=.01, vmax=.1, label='AF')
# fig.tight_layout()
# plt.show()


# Bootstrap trees
tree_list = []
for i in range(100):
    print(f'Bootstrap sample {i}')
    # AD_, DP_, sel_idx = jackknife_allele_tables(AD.A.T, DP.A.T)
    AD_, DP_ = bootstrap_allele_counts(AD.A.T, DP.A.T)
    X = np.divide(AD_, DP_)
    X[np.isnan(X)] = 0
    tree = build_tree(
        a, 
        X=X, 
        metric='cosine',
        # variants=a.var_names[sel_idx],
        solver=cs.solver, 
        # kwargs={'add_root':True}
    )
    tree_list.append(tree)


##


# # Plot some bootstrapped trees
# for var in a.var_names:
# 
#     fig = plt.figure(figsize=(15,15))
# 
#     for i in range(16):
# 
#         tree = tree_list[i*3]
#         tree.cell_meta = pd.DataFrame(
#             a.obs.loc[tree.leaves, "GBC"].astype(str)
#         )
#         tree.cell_meta[a.var_names] = a[tree.leaves, :].X.toarray()
# 
#         ax = fig.add_subplot(4,4,i+1)
#         cs.pl.plot_matplotlib(
#             tree, 
#             orient='down',
#             meta_data=[var, 'GBC'], 
#             categorical_cmap=ListedColormap(sc.pl.palettes.godsnot_102), 
#             continuous_cmap='mako',
#             add_root=True,
#             ax=ax
#         )
#         format_ax(ax, title=f'{sample}, bootstrap {i+1}')
# 
#     # Save fig
#     fig.tight_layout()
#     fig.savefig(os.path.join(path_viz, f'{var}_bootstrap.png'), dpi=400)


##


# Calculate observed branches supports
obs_tree = build_tree(
    a, X=None, 
    solver=cs.solver.NeighborJoiningSolver, 
    metric='cosine',
    kwargs={'add_root':True}
)
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


##


# Clades and bootstrap diagnostics
fig, axs = plt.subplots(1,2,figsize=(10,5))

colors = { 'expanded' : 'r', 'other' : 'k' }

scatter(supports_df, 'time', 'support', ax=axs[0], by='expanded_clones', c=colors, s=50)
sns.regplot(data=supports_df, x='time', y='support', scatter=False, ax=axs[0])
pho = np.corrcoef(supports_df['time'], supports_df['support'])[0,1]
format_ax(
    axs[0], xlabel='Molecular time', ylabel='Support', 
    title=f'Pearson\'s r: {pho:.2f}, median support: {supports_df["support"].median():.2f}',
    reduced_spines=True
)

scatter(supports_df, 'n_cells', 'support', ax=axs[1], by='expanded_clones', c=colors, s=50)
sns.regplot(data=supports_df, x='n_cells', y='support', scatter=False, ax=axs[1])
pho = np.corrcoef(supports_df['n_cells'], supports_df['support'])[0,1]
format_ax(
    axs[1], xlabel='n cells', ylabel='Support', 
    title=f'Pearson\'s r: {pho:.2f}, median support: {supports_df["support"].median():.2f}',
    reduced_spines=True
)
add_legend(
    ax=axs[1], label='Clade', colors=colors, 
    ticks_size=8, artists_size=7, label_size=9,
    bbox_to_anchor=(.95,.95), loc='upper right'
)

fig.tight_layout()
fig.savefig(os.path.join(path_viz, f'{sample}_boot_rel.png'), dpi=400)


# Bootstrap record
fig, ax = plt.subplots(figsize=(6,5))

n = sturges(bootstrap_record)
hist(bootstrap_record, 'perc_clades_found', ax=ax, c='k', n=n)
m = bootstrap_record['perc_clades_found'].median()
n = bootstrap_record.shape[0]
format_ax(
    ax, ylabel='n bootstrap samples', xlabel='% clades found',
    title=f'{n} bootstrap samples \n median % of clades found per replicate: {m:.2f}', 
    reduced_spines=True
)
fig.tight_layout()
fig.savefig(os.path.join(path_viz, f'{sample}_boot_record.png'), dpi=400)


##


# Annotate tree for plotting
obs_tree.cell_meta = pd.DataFrame(
    a.obs.loc[obs_tree.leaves, "GBC"].astype(str)
)
obs_tree.cell_meta[a.var_names] = a[obs_tree.leaves, :].X.toarray()


##


# Main
fig, ax = plt.subplots(figsize=(13,8))

mask = (supports_df['time']<=5) # & (supports_df['support']>.5) 
clones_colors = create_palette(obs_tree.cell_meta, 'GBC', 'tab10')

plot_tree(
    obs_tree,
    meta_data=[*a.var_names.to_list(),'GBC'], 
    orient='down',
    internal_node_kwargs={
        's': 8, 
        'by': 'support',
        'meta': supports_df,
        'plot' : True,
        'annot' : False,
        'continuous_cmap' : 'YlOrBr'
       # 'mask' : mask,
    },
    branch_kwargs={'c':'k', 'linewidth':.7},
    colorstrip_spacing=.1, 
    categorical_cmap=clones_colors,
    ax=ax,
)

format_ax(
    ax=ax, 
    title=f'''
        {sample}, {bootstrap_record.shape[0]} bootstrap samples \n
        {a.shape[0]} cells, {supports_df.shape[0]} internal nodes. Support: {supports_df['support'].median():.2f} (+-{supports_df['support'].std():.2f})
        '''
)
add_legend(
    'Clones', colors=clones_colors, ax=ax, bbox_to_anchor=(-.3,.5), loc='center left',
    label_size=10, artists_size=9, ticks_size=9
)
add_cbar(
    supports_df['support'], palette='YlOrBr', ax=ax, label='Support',
    ticks_size=7, label_size=9
)
add_cbar(
    a.X.flatten(), palette='mako', ax=ax, label='AF',
    ticks_size=7, label_size=9, layout=( (1.2,.25,.03,.5), 'right' ),
    vmax=.1, vmin=.01
)

fig.subplots_adjust(left=.22, right=.78, top=.8)
fig.savefig(os.path.join(path_viz, f'{sample}_boot_supports.png'), dpi=400)


##


# Annotated supports
fig, ax = plt.subplots(figsize=(13,13))

plot_tree(
    obs_tree,
    orient=90,
    internal_node_kwargs={
        's': 0, 
        'by': 'support',
        'meta': supports_df,
        'plot' : False,
        'annot' : True,
    },
    branch_kwargs={'c':'k', 'linewidth':.7},
    colorstrip_spacing=.1, 
    annot_size=7,
    ax=ax,
)

fig.tight_layout()
fig.savefig(os.path.join(path_viz, f'{sample}_boot_supports.png'), dpi=400)






##









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
