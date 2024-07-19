"""
Visualiza phylo output.
"""

import os
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.phylo_plots import *
from mito_utils.heatmaps_plots import *
from mito_utils.plotting_base import *
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/mito_bench/'
path_results = os.path.join(path_main, 'results/phylo')


##

base_dir = path_results
file_name = 'nodes.csv'



nodes_d = traverse_and_extract_flat(base_dir, file_name=file_name)
df = pd.concat([   
    df.assign(sample=sample, variant_set=variant_set) for \
    (sample, variant_set), df in nodes_d.items() 
])
    
(
    df[~df['support'].isna()]
    .groupby(['sample', 'variant_set']).apply(lambda x: x['support'].sort_values().tail(50).median()) #np.percentile(x['support'], 75))
    .reset_index()
)



sns.kdeplot()
plt.show()

df['support'].describe()
df.loc[lambda x: ~x['support'].isna() & ~x['assigned_var'].isna()]['support'].median()

df.






# 1. How robust are trees to bootstrapping?
df.describe()
df.groupby(['sample', 'filtering', 'solver', 'bootstrap']).median()['TBE'].sort_values(ascending=False)
df.groupby('sample').median()['TBE'].sort_values(ascending=False)

# Reformat for plotting
df = df.melt(id_vars=df.columns[~df.columns.isin(['TBE', 'FBP'])], var_name='support')
df['Support type'] =  df['support'] + ", " + df['bootstrap']

# Individual nodes
fig, ax = plt.subplots(figsize=(7, 3.5))
sns.kdeplot(data=df, x="value", fill=True, hue='Support type', ax=ax, legend=True)
ax.set(title=f'Internal nodes support distribution')
ax.set(xlabel='support')
ax.spines[['right', 'top']].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'support_distribution.pdf'), dpi=500)


##


# Take only jacknife results
df_ = df.loc[df['Support type'] == "TBE, jacknife"]

# Nodes, molecular time relationship
fig, axs = plt.subplots(1,2, figsize=(10, 5))
sns.regplot(x='time', y='value', data=df_, scatter=False, ax=axs[0])
axs[0].plot(df_['time'], df_['value'], 'ko', markersize=1)
format_ax(axs[0], reduced_spines=True, xlabel='Molecular time', ylabel='TBE',
    title=f'Pearson\'s correlation: {np.corrcoef(df_["value"], df_["value"])[0,1]:.2f}'
)
sns.regplot(x='n_cells', y='value', data=df_, scatter=False, ax=axs[1])
axs[1].plot(df_['n_cells'], df_['value'], 'ko', markersize=1)
format_ax(axs[1], reduced_spines=True, xlabel='n cells', ylabel='TBE',
    title=f'Pearson\'s correlation: {np.corrcoef(df_["n_cells"], df["value"])[0,1]:.2f}'
)
fig.tight_layout()
plt.show()


# Fig
fig, axs = plt.subplots(1,2, figsize=(4.5, 4.5))

colors = {
    'GT': '#1f77b4', 'MQuad': '#ff7f0e'
}

df_ = df.loc[df['Support type'] == "TBE, jacknife"]
order = df_.groupby('solver')['value'].median().sort_values(ascending=False).index
df_['solver'] = pd.Categorical(df_['solver'], categories=order)
box(df_, x='solver', y='value', by='filtering', c=colors, ax=axs[0], order=order)
format_ax(ax=axs[0], ylabel='TBE', reduced_spines=True)

order = df_.groupby('solver')['median_RF'].median().sort_values().index
df_['solver'] = pd.Categorical(df_['solver'], categories=order)
box(df_, x='solver', y='median_RF', by='filtering', c=colors, ax=axs[1])
add_legend(
    label='MT-SNV filtering method', ax=axs[0], colors=colors, 
    label_size=10, artists_size=9, ticks_size=8,
    loc='lower left', bbox_to_anchor=(.6,1.05), 
    ncols=df['filtering'].unique().size
)
format_ax(ax=axs[1], reduced_spines=True, ylabel='RF distance')
fig.subplots_adjust(top=.8, wspace=.55, left=.2, right=.9)
fig.savefig(os.path.join(path_results, 'RF_and_TBE_solvers_filtering.pdf'), dpi=500)


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
# df = pd.read_csv(os.path.join(path_results, 'clones_df.csv'), index_col=0)
# df_ = (
#     df
#     .groupby(['sample', 'filtering', 'solver', 'metric'])
#     .apply(lambda x: (np.sum(x['mpd.obs.p']<=.05) / x['mpd.obs.p'].size))
#     .reset_index(name='perc_significant')
# )
# 
# ##
# 
# fig, ax = plt.subplots(figsize=(6, 5))
# 
# colors = create_palette(df_, 'filtering', 'tab10')
# order = df_.groupby('solver')['perc_significant'].median().sort_values(ascending=False).index
# df_['solver'] = pd.Categorical(df_['solver'], categories=order)
# box(df_, x='solver', y='perc_significant', by='filtering', c=colors, ax=ax, order=order)
# format_ax(
#     ax, title='Phylogenetic clustering of clonal labels', 
#     ylabel='% of clones with significant phylogenetic clustering',
#     xlabel='Solver'
# )
# add_legend(
#     label='MT-SNV filtering method', ax=ax, colors=colors, 
#     label_size=10, artists_size=9, ticks_size=8,
#     loc='center left', bbox_to_anchor=(1,.5)
# )
# 
# fig.tight_layout()
# plt.show()

##

# fig, ax = plt.subplots(figsize=(6, 5))
# 
# df_to_plot = (
#     df_.groupby('filtering')['perc_significant'].median()
#     .reset_index().sort_values('perc_significant', ascending=False)
# )
# df_to_plot['perc_significant'] = np.round(df_to_plot['perc_significant'], 2)
# bar(df_to_plot, x='filtering', y='perc_significant', s=.5, c='#D2CBCB', edgecolor='k', ax=ax)
# format_ax(
#     ax=ax, title='Phylogenetic clustering of clonal labels', xticks=df_to_plot['filtering'],
#     xlabel='MT-SNVs filtering', ylabel='% clones with significant phylogenetic clustering'
# )
# ax.spines[['left', 'right', 'top']].set_visible(False)
# plt.show()


##


# 4. Are input mutations significantly associated with tree topologies??
# df = pd.read_csv(os.path.join(path_results, 'muts_df.csv'), index_col=0)
# 
# df_ = (
#     df
#     .groupby(['sample', 'filtering', 'solver', 'metric'])
#     .apply(lambda x: (np.sum(x['k_p']<=.05) / x['k_p'].size))
#     .reset_index(name='perc_significant')
# )
# 
# fig, ax = plt.subplots(figsize=(6, 5))
# 
# colors = create_palette(df_, 'filtering', 'tab10')
# order = df_.groupby('solver')['perc_significant'].median().sort_values(ascending=False).index
# df_['solver'] = pd.Categorical(df_['solver'], categories=order)
# box(df_, x='solver', y='perc_significant', by='filtering', c=colors, ax=ax, order=order)
# format_ax(
#     ax, title='Phylogenetic signal of input MT-SNVs', 
#     ylabel='% of MT-SNVS with significant phylogenetic signal',
#     xlabel='Solver'
# )
# add_legend(
#     label='MT-SNV filtering method', ax=ax, colors=colors, 
#     label_size=10, artists_size=9, ticks_size=8,
#     loc='center left', bbox_to_anchor=(1,.5)
# )
# 
# fig.tight_layout()
# plt.show()

##