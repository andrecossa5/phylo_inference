"""
Collapse pipeline output.
"""

# Code
import os
import re
import pandas as pd
import seaborn as sns
from mito_utils.plotting_base import *
from matplotlib.gridspec import GridSpec

from mito_utils.phylo_plots import plot_tree
from mito_utils.preprocessing import *

matplotlib.use('macOSX')


## Helper

def get_files_l(path, pattern=None, print_all=False, from_tokey=None):

    path = path_data
    L = []
    d = {}
    for root, _, files in os.walk(path):
        for file in files:
            if print_all:
                print(os.path.join(root, file))
            else:
                if bool(re.search(pattern, file)):
                    if from_tokey is not None:
                        key = root.split('/')[-from_tokey:]
                        if len(key)>1:
                            key = tuple(key)
                        else:
                            key = key[0]
                        d[key] = os.path.join(root, file)
                    else:
                        L.append(os.path.join(root, file))

    if from_tokey is not None:
        return d
    else:          
        return L


# Paths
path_data = '/Users/IEO5505/Desktop/mito_bench/results/phylo/out_in_vitro'


##


# Read tree-level outputs
get_files_l(path_data, print_all=True)
df = pd.concat([ 
    pd.read_csv(x, index_col=0) 
    for x in get_files_l(os.path.join(path_data), pattern='supp')
])

# df.iloc[:,8].sort_values()[:100]

# (
#     df.query('bootstrap == "jacknife" and filtering == "MQuad_optimized"')
#     .groupby(['sample', 'filtering', 'solver', 'metric'])
#     ['median_RF'].describe()
#     ['min'].sort_values(ascending=False)
#     .describe()
# )

# Plot tree diagnostic metrics

# Fig
fig, ax = plt.subplots(figsize=(5, 5))
m = df["median_RF"].median()
std = df["median_RF"].std()
sns.kdeplot(data=df.reset_index(), x="median_RF", fill=True, hue='bootstrap', ax=ax)
ax.set(title=f'median_RF, observed vs bootstrap trees: {m:.2f}+-{std:.2f}')
# ax.axvline(x=.75, c='grey', linestyle='dashed')
fig.tight_layout()

plt.show()


##


# Fig
fig, axs = plt.subplots(1,2, figsize=(12, 5))

colors = create_palette(df, 'filtering', 'tab10')
box(df, x='solver', y='support', by='filtering', c=colors, ax=axs[0])
box(df, x='solver', y='median_RF', by='filtering', c=colors, ax=axs[1])
add_legend(label='MT-SNV filtering method', ax=axs[0], colors=colors, label_size=10, artists_size=9, ticks_size=8,
           loc='lower left', bbox_to_anchor=(.6,1.05), ncols=df['filtering'].unique().size)
fig.subplots_adjust(top=.8)
plt.show()


# Fig
fig, ax = plt.subplots(figsize=(5, 4.5))
order = df.groupby('filtering')['support'].median().sort_values(ascending=False).index
box(df, x='filtering', y='support', c='white', ax=ax, order=order)
fig.tight_layout()
plt.show()


##


# Plot some good trees, and save for R session
sample = 'MDA_clones'

(
    df.query('sample==@sample')
    .groupby(['sample', 'filtering', 'bootstrap', 'solver', 'metric'])
    ['support'].median().sort_values(ascending=False)
    #.describe()
)

filtering = 'GT'
solver = 'UPMGA'
metric = 'hamming'
bootstrap = 'jacknife'
path_tree = os.path.join(path_data, sample, filtering, solver, metric, bootstrap, 'trees.pickle')
path_afm = '/Users/IEO5505/Desktop/mito_bench/data/'

# Read 
with open(path_tree, 'rb') as f:
    tree = pickle.load(f)
# Get observed
obs_tree = tree['observed']
tree_nw = obs_tree.get_newick(record_branch_lengths=True)

# Get cell meta

# deltaBIC but more...
afm = read_one_sample(path_afm, sample=sample, with_GBC=True)

a_cells = filter_cells_coverage(afm)
a_cells = filter_baseline(a_cells)
df = fit_MQuad_mixtures(a_cells, nproc=8, path_=os.getcwd())

variants = (
    df
    .join(summary_stats_vars(a_cells))
    .sort_values('deltaBIC', ascending=False)
    .head(100)
    .query('fr_positives>=.01 and num_cells_nonzero_AD>=1')
    .index
)

a = a_cells[:, variants].copy()
a = remove_excluded_sites(a)


plt.plot(np.sort(a[:,variants[-1]].X.flatten()), 'ko')
plt.show()








from mito_utils.distances import pair_d
from mito_utils.phylo import AFM_to_seqs


seqs = AFM_to_seqs(a, .1)

np.unique(list(seqs.values())).size / len(seqs)





red = pd.DataFrame((np.sum(D==0, axis=0)/D.shape[0]).flatten(), columns=['red'])


(D==0).all()

sns.kdeplot('red', data=red)
plt.show()





# Save data
path_ = '/Users/IEO5505/Desktop/mito_bench/results/phylo/to_R'
df_afm = pd.DataFrame(a.X, index=a.obs_names, columns=a.var_names)
df_afm = df_afm.loc[obs_tree.leaves]
df_afm.to_csv(os.path.join(path_, 'afm.csv'))
a.obs.loc[obs_tree.leaves].to_csv(os.path.join(path_, 'meta.csv'))
with open(os.path.join(path_, 'tree.newick'), 'w') as f:
    f.write(tree_nw)


##


# Create annot dfs
obs_tree.cell_meta = pd.DataFrame(
    a_cells.obs.loc[obs_tree.leaves, "GBC"].astype(str)
)
df_supports = (
    pd.DataFrame(
    index=obs_tree.internal_nodes
    ).join(
        df.query('sample==@sample and filtering==@filtering and solver==@solver and bootstrap==@bootstrap and metric==@metric')
        [['support', 'median_RF']]
    )# .drop_duplicates()
)

# Viz
fig, ax = plt.subplots(figsize=(5, 5))

vmin = .5
vmax = .9
plot_tree(
    obs_tree,
    meta_data=['GBC'],
    internal_node_kwargs={
        's': 35, 
        'by': 'support',
        'meta': df_supports,
        'plot' : True,
        'annot' : False,
        'continuous_cmap' : 'YlOrBr',
        # 'mask' : mask,
    },
    vmin=vmin, vmax=vmax,
    ax=ax,
)
add_cbar(df_supports['support'], palette='YlOrBr', ax=ax, label='Support', vmin=vmin, vmax=vmax)
m_s = df_supports['support'].median()
std_s = df_supports['support'].std()
m_RF = df_supports['median_RF'].unique()[0]

t = f'sample: {sample}; filtering: {filtering},\n solver:{solver}; bootstrap: {bootstrap} \n Support: {m_s:.2f} (+-{std_s:.2f}); Median RF: {m_RF:.2f}'
format_ax(ax=ax, title=t, title_size=9)

fig.tight_layout()

plt.show()










# matplotlib.use('macOSX')








