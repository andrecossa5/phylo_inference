"""
Viz PATH output.
"""

import os
from itertools import product
from mito_utils.plotting_base import *
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/mito_bench/'
path_results = os.path.join(path_main, 'results/phylo')


##



# Read PATH output
samples = ['MDA_clones', 'AML_clones', 'MDA_lung', 'MDA_PT']
sets = os.listdir(os.path.join(path_results, samples[0]))

# Viz 4 heatmap phylo correlations
for variant_set in sets:

    fig, axs = plt.subplots(1,4,figsize=(10,3))
    for i, sample in enumerate(samples):
        df = pd.read_csv(os.path.join(path_results, sample, variant_set, 'lineage_association.csv'), index_col=0)
        X = df[['Var1', 'Var2', 'Z']].pivot_table(index='Var1', columns='Var2').values
        vmin = np.percentile(X.flatten(), 5)
        vmax = np.percentile(X.flatten(), 95)
        axs[i].imshow(X, cmap='inferno', vmin=vmin, vmax=vmax)
        if i==0:
            format_ax(axs[i], title=sample, xlabel='Lentiviral clones', ylabel='Lentiviral clones', xticks=[], yticks=[])
        else:
            format_ax(axs[i], title=sample, xlabel=None, ylabel=None, xticks=[], yticks=[])
        if i==3:
            add_cbar(X.flatten(), ax=axs[i], palette='inferno', vmin=vmin, vmax=vmax, layout='outside', label='Phylo correlation')

    fig.tight_layout()
    fig.suptitle(variant_set)
    fig.savefig(os.path.join(path_results, f'{variant_set}_GBCs.png'), dpi=1000)


##


# All sets, all samples
L = []
sets = os.listdir(os.path.join(path_results, sample))
for variant_set, sample in product(sets, samples):
    df = pd.read_csv(os.path.join(path_results, sample, variant_set, 'lineage_association.csv'), index_col=0)
    L.append(df.assign(sample=sample, variant_set=variant_set))
df = pd.concat(L)
df['corr_type'] = np.where(df['Var1'] == df['Var2'], 'auto', 'cross')


df_ = (
    df
    .groupby(['sample', 'variant_set', 'corr_type'])
    ['Z'].median().reset_index()
    .pivot(index=['sample', 'variant_set'], values='Z', columns=['corr_type'])
    .reset_index()
    .assign(delta=lambda x: np.abs(x['auto']-x['cross']))
)

# Viz
fig, ax = plt.subplots(figsize=(3,3.5))
sample = 'AML_clones'
set_ = 'set2'
box(df.query('sample==@sample and variant_set==@set_'), 'corr_type', 'Z', ax=ax, c='grey', with_stats=True, pairs=[['auto', 'cross']])
format_ax(ax, title=variant_set, ylabel='z-scored phylogenetic correlation', reduced_spines=True)
fig.tight_layout()
plt.show()


# fig.savefig(os.path.join(path_results, f'{variant_set}_phylocorr.png'), dpi=1000)