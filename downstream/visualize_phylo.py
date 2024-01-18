"""
Visualiza phylo output.
"""

import os
from scipy.sparse import load_npz
from anndata import AnnData
from mito_utils.preprocessing import *
from mito_utils.phylo_plots import *
from mito_utils.heatmaps_plots import *
from mito_utils.plotting_base import *
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/mito_bench/'
path_results = os.path.join(path_main, 'results/phylo/output_in_vitro')
path_viz = os.path.join(path_main, 'results/phylo/visualization')

# Read results
os.listdir(path_results)


##


# Tree viz
df = pd.read_csv(os.path.join(path_results, 'supports_df.csv'), index_col=0)

# 1. How robust are trees to bootstrapping?

# Individual nodes
fig, axs = plt.subplots(1,2, figsize=(10, 5))
sns.kdeplot(data=df.reset_index(), x="support", fill=True, hue='bootstrap', ax=axs[0])
axs[0].set(title=f'Internal nodes support distribution')
m = df.query('bootstrap == "jacknife"')['support'].median()
std = df.query('bootstrap == "jacknife"')['support'].std()
axs[0].text(.25, .7, f'jacknife support: {m:.2f} (+-{std:.2f})', transform=axs[0].transAxes)
m = df.query('bootstrap == "counts_resampling"')['support'].median()
std = df.query('bootstrap == "counts_resampling"')['support'].std()
axs[0].text(.25, .65, f'counts_resampling support: {m:.2f} (+-{std:.2f})', transform=axs[0].transAxes)
axs[0].axvline(x=.75, c='grey', linestyle='dashed')

sns.kdeplot(data=df.reset_index(), x="median_RF", fill=True, hue='bootstrap', ax=axs[1])
axs[1].set(title=f'median RF distances distribution')
m = df.query('bootstrap == "jacknife"')['median_RF'].median()
std = df.query('bootstrap == "jacknife"')['median_RF'].std()
axs[1].text(.05, .7, f'jacknife median_RF: {m:.2f} (+-{std:.2f})', transform=axs[1].transAxes)
m = df.query('bootstrap == "counts_resampling"')['median_RF'].median()
std = df.query('bootstrap == "counts_resampling"')['median_RF'].std()
axs[1].text(.05, .65, f'counts_resampling median_RF: {m:.2f} (+-{std:.2f})', transform=axs[1].transAxes)
fig.tight_layout()


# Take only jacknife results
df = df.query('bootstrap == "jacknife"')

# Nodes, molecular time relationship
fig, axs = plt.subplots(1,2, figsize=(10, 5))
sns.regplot(x='time', y='support', data=df, scatter=False, ax=axs[0])
axs[0].plot(df['time'], df['support'], 'ko', markersize=1)
format_ax(axs[0], reduced_spines=True, xlabel='Molecular time', ylabel='Support',
    title=f'Pearson\'s correlation: {np.corrcoef(df["time"], df["support"])[0,1]:.2f}'
)
sns.regplot(x='n_cells', y='support', data=df, scatter=False, ax=axs[1])
axs[1].plot(df['n_cells'], df['support'], 'ko', markersize=1)
format_ax(axs[1], reduced_spines=True, xlabel='n cells', ylabel='Support',
    title=f'Pearson\'s correlation: {np.corrcoef(df["n_cells"], df["support"])[0,1]:.2f}'
)
fig.tight_layout()
plt.show()


# Fig
fig, axs = plt.subplots(1,2, figsize=(12, 5))

colors = create_palette(df, 'filtering', 'tab10')
order = df.groupby('solver')['support'].median().sort_values(ascending=False).index
df['solver'] = pd.Categorical(df['solver'], categories=order)
box(df, x='solver', y='support', by='filtering', c=colors, ax=axs[0], order=order)

order = df.groupby('solver')['median_RF'].median().sort_values().index
df['solver'] = pd.Categorical(df['solver'], categories=order)
box(df, x='solver', y='median_RF', by='filtering', c=colors, ax=axs[1])

add_legend(
    label='MT-SNV filtering method', ax=axs[0], colors=colors, 
    label_size=10, artists_size=9, ticks_size=8,
    loc='lower left', bbox_to_anchor=(.6,1.05), 
    ncols=df['filtering'].unique().size
)
fig.subplots_adjust(top=.8)
plt.show()


##


# 2. How consistent are trees, across solvers?
df = pd.read_csv(os.path.join(path_results, 'filtering_df.csv'), index_col=0)
df['median_RF'] = np.round(df['median_RF'], 2)
df_ = df.loc[lambda x: x.index.str.contains('MDA')]

fig, ax = plt.subplots(figsize=(7,5))
bar(df_.sort_values('median_RF'),
    x='filtering', y='median_RF', s=.5, 
    c='#D2CBCB', edgecolor='k', ax=ax
)
ax.set_ylim((0.95, 1.005))
format_ax(
    ax=ax, title='RF distance across solvers', xticks=df_.sort_values('median_RF')['filtering'],
    xlabel='MT-SNVs filtering', ylabel='Median RF distance'
)
ax.spines[['left', 'right', 'top']].set_visible(False)
plt.show()


##


# 3. Are gt clones associated with tree topologies??
df = pd.read_csv(os.path.join(path_results, 'clones_df.csv'), index_col=0)
df_ = (
    df
    .groupby(['sample', 'filtering', 'solver', 'metric'])
    .apply(lambda x: (np.sum(x['mpd.obs.p']<=.05) / x['mpd.obs.p'].size))
    .reset_index(name='perc_significant')
)

##

fig, ax = plt.subplots(figsize=(6, 5))

colors = create_palette(df_, 'filtering', 'tab10')
order = df_.groupby('solver')['perc_significant'].median().sort_values(ascending=False).index
df_['solver'] = pd.Categorical(df_['solver'], categories=order)
box(df_, x='solver', y='perc_significant', by='filtering', c=colors, ax=ax, order=order)
format_ax(
    ax, title='Phylogenetic clustering of clonal labels', 
    ylabel='% of clones with significant phylogenetic clustering',
    xlabel='Solver'
)
add_legend(
    label='MT-SNV filtering method', ax=ax, colors=colors, 
    label_size=10, artists_size=9, ticks_size=8,
    loc='center left', bbox_to_anchor=(1,.5)
)

fig.tight_layout()
plt.show()

##

fig, ax = plt.subplots(figsize=(6, 5))

df_to_plot = (
    df_.groupby('filtering')['perc_significant'].median()
    .reset_index().sort_values('perc_significant', ascending=False)
)
df_to_plot['perc_significant'] = np.round(df_to_plot['perc_significant'], 2)
bar(df_to_plot, x='filtering', y='perc_significant', s=.5, c='#D2CBCB', edgecolor='k', ax=ax)
format_ax(
    ax=ax, title='Phylogenetic clustering of clonal labels', xticks=df_to_plot['filtering'],
    xlabel='MT-SNVs filtering', ylabel='% clones with significant phylogenetic clustering'
)
ax.spines[['left', 'right', 'top']].set_visible(False)
plt.show()


##


# 4. Are input mutations significantly associated with tree topologies??
df = pd.read_csv(os.path.join(path_results, 'muts_df.csv'), index_col=0)

df_ = (
    df
    .groupby(['sample', 'filtering', 'solver', 'metric'])
    .apply(lambda x: (np.sum(x['k_p']<=.05) / x['k_p'].size))
    .reset_index(name='perc_significant')
)

fig, ax = plt.subplots(figsize=(6, 5))

colors = create_palette(df_, 'filtering', 'tab10')
order = df_.groupby('solver')['perc_significant'].median().sort_values(ascending=False).index
df_['solver'] = pd.Categorical(df_['solver'], categories=order)
box(df_, x='solver', y='perc_significant', by='filtering', c=colors, ax=ax, order=order)
format_ax(
    ax, title='Phylogenetic signal of input MT-SNVs', 
    ylabel='% of MT-SNVS with significant phylogenetic signal',
    xlabel='Solver'
)
add_legend(
    label='MT-SNV filtering method', ax=ax, colors=colors, 
    label_size=10, artists_size=9, ticks_size=8,
    loc='center left', bbox_to_anchor=(1,.5)
)

fig.tight_layout()
plt.show()


##


# Plot one good tree example
df = pd.read_csv(os.path.join(path_results, 'supports_df.csv'), index_col=0)
sample = 'MDA_clones'

(
    df.query('sample==@sample and time<=10 and time>=5')
    .groupby(['sample', 'filtering', 'bootstrap', 'solver', 'metric'])
    ['support'].median().sort_values(ascending=False)
    #.describe()
)

filtering = 'MQuad'
solver = 'UPMGA'
metric = 'cosine'
bootstrap = 'jacknife'
path_tree = os.path.join(
    path_results, sample, filtering, 
    solver, metric, bootstrap, 'trees.pickle'
)
path_data = os.path.join(path_results, sample, filtering, 'input_folder')

##

# Read tree
with open(path_tree, 'rb') as f:
    tree = pickle.load(f)

# Get observed tree
obs_tree = tree['observed']

# Read cells, variants and build the afm
variants = pd.read_csv(os.path.join(path_data, 'variants.csv'), index_col=0, header=None)
meta = pd.read_csv(os.path.join(path_data, 'meta.csv'), index_col=0)
AD = load_npz(os.path.join(path_data, 'AD.npz')).astype(np.int16).A
DP = load_npz(os.path.join(path_data, 'DP.npz')).astype(np.int16).A
afm = AnnData(AD/DP, obs=meta, var=variants)

# Create annot dfs
obs_tree.cell_meta = pd.DataFrame(
    afm.obs.loc[obs_tree.leaves, "GBC"].astype(str)
    ).join(
        pd.DataFrame(
        afm[obs_tree.leaves,:].X, index=obs_tree.leaves, columns=afm.var_names)
)

df_supports = (
    df.query('sample==@sample and filtering==@filtering and solver==@solver and bootstrap==@bootstrap and metric==@metric')
    .loc[obs_tree.internal_nodes, ['support', 'median_RF']]
)

# Viz
fig, ax = plt.subplots(figsize=(7, 7))

vmin_annot = .001
vmax_annot = .1
plot_tree(
    obs_tree,
    meta_data=['GBC'] + afm.var_names.to_list(),
    orient='down',
    internal_node_kwargs={
        's': 20, 
        'by': 'support',
        'meta': df_supports,
        'plot' : True,
        'annot' : False,
        'continuous_cmap' : 'YlOrBr',
        'vmin' : .5,
        'vmax' : .7
        # 'mask' : mask,
    },
    colorstrip_spacing=.05,
    colorstrip_width=1,
    vmin_annot=vmin_annot, vmax_annot=vmax_annot,
    ax=ax,
)
add_cbar(df_supports['support'], palette='YlOrBr', ax=ax, label='Support', vmin=.5, vmax=.7)
add_cbar(
    afm.X.flatten(), palette='mako', ax=ax, label='AF', vmin=vmin_annot, vmax=vmax_annot,
    layout=( (-.08,.25,.03,.5), 'left' )
)
m_s = df_supports['support'].median()
std_s = df_supports['support'].std()
m_RF = df_supports['median_RF'].unique()[0]

t = f'sample: {sample}; filtering: {filtering},\n solver:{solver}; bootstrap: {bootstrap} \n Support: {m_s:.2f} (+-{std_s:.2f}); Median RF: {m_RF:.2f}'
format_ax(ax=ax, title=t, title_size=9)

fig.tight_layout()
fig.savefig(os.path.join(path_viz, 'tree_.png'), dpi=500)


##








