"""
Viz PATH output.
"""

import os
from mito_utils.plotting_base import *
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/mito_bench/'
path_results = os.path.join(path_main, 'results/phylo_inference')


##



# Read PATH output
samples = ['MDA_clones', 'AML_clones', 'MDA_lung', 'MDA_PT']

# Viz 4 heatmap phylo correlations
fig, axs = plt.subplots(1,4,figsize=(10,3))
for i, sample in enumerate(samples):
    df = pd.read_csv(os.path.join(path_results, 'top_trees_MQuad', sample, 'phylocorr_GBC.csv'), index_col=0)
    X = df[['GBC1', 'GBC2', 'Z']].pivot_table(index='GBC1', columns='GBC2').values
    axs[i].imshow(X, cmap='inferno', vmin=-10, vmax=30)
    if i==0:
        format_ax(axs[i], title=sample, xlabel='Lentiviral clones', ylabel='Lentiviral clones', xticks=[], yticks=[])
    else:
        format_ax(axs[i], title=sample, xlabel=None, ylabel=None, xticks=[], yticks=[])
    if i==3:
        add_cbar(X.flatten(), ax=axs[i], palette='inferno', vmin=-10, vmax=30, layout='outside', label='Phylo correlation')

fig.tight_layout()
fig.suptitle('MQuad')
# fig.savefig(os.path.join(path_results, 'PATH_GBCs.pdf'), dpi=500)
fig.savefig(os.path.join(path_results, 'PATH_MQuad_GBCs.pdf'), dpi=500)


##


# All
L = []
for i, sample in enumerate(samples):
    df = pd.read_csv(os.path.join(path_results, 'top_trees_MQuad', sample, 'phylocorr_GBC.csv'), index_col=0)
    L.append(df)
df = pd.concat(L)
df['corr_type'] = np.where(df['GBC1'] == df['GBC2'], 'auto', 'cross')

# Viz
fig, ax = plt.subplots(figsize=(3,3.5))
box(df, 'corr_type', 'Z', ax=ax, c='grey', with_stats=True, pairs=[['auto', 'cross']])
format_ax(ax, title='MQuad', ylabel='z-scored phylogenetic correlation', reduced_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'MQuad_phylocorr.pdf'), dpi=500)