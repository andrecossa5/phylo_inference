"""
Visualize tree topology/ GT clonal labels association.
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

# For each variant set and sample, visualize phylo correlations
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
    df['corr_type'] = np.where(df['Var1'] == df['Var2'], 'auto', 'cross')
    df = (
        df.query('corr_type=="auto"').drop(['Var2', 'corr_type'], axis=1).rename(columns={'Var1':'GBC'}).set_index('GBC')
        .join(
            pd.read_csv(os.path.join(path_results, sample, variant_set, 'filtered_input', 'meta.csv'), index_col=0)
            ['GBC'].value_counts(normalize=True).to_frame('prevalence')
        )
    )
    L.append(df.assign(sample=sample, variant_set=variant_set))

# Concat
df = pd.concat(L)

# Set variant set order
df_perc_significant = (
    df.groupby(['variant_set', 'sample'])
    .apply( lambda x: np.sum((x['Z']>10) & (x['p'] <= 0.05) & (x['prevalence'] >= .01)) / (x['prevalence'] >= .01).sum() )
    .to_frame('perc_significant').reset_index()
    .groupby('variant_set')
    ['perc_significant'].median().sort_values()
)
df['variant_set'] = pd.Categorical(df['variant_set'], categories=df_perc_significant.index)


##


# Viz
df['sample'] = pd.Categorical(df['sample'], categories=['AML_clones', 'MDA_clones', 'MDA_PT', 'MDA_lung'])
sample_colors = create_palette(df, 'sample', 'tab10')

fig, ax = plt.subplots(figsize=(5.5,5.5))

for i, variant_set in enumerate(df['variant_set'].unique()):
    subset = df.query('variant_set==@variant_set')
    jitter = 0.1 * np.random.randn(len(subset))
    ax.scatter(
        x=subset['Z'],
        y=np.ones(len(subset)) * i + jitter,
        s=subset['prevalence'] * 100,
        c=[ sample_colors[x] for x in subset['sample'] ] , 
        alpha=0.7
    )

format_ax(ax, title='', xlabel='z-scored phylogenetic correlation', ylabel='Variant set', 
        yticks=df_perc_significant.index, reduced_spines=True)
add_legend('Sample', ax=ax, colors=sample_colors, loc='upper right', bbox_to_anchor=(1,1), label_size=9, artists_size=9, ticks_size=8)

fig.tight_layout()
fig.savefig(os.path.join(path_results, f'summary_GBC_phylocorr.png'), dpi=1000)


##


# Viz
fig, ax = plt.subplots(figsize=(5.5,5))
df_ = (df_perc_significant*100).round()[::-1].to_frame('perc').reset_index()
bar(df_, x='variant_set', y='perc', s=.75, ax=ax, edgecolor='k', c='grey')
format_ax(
    ax, title='', xlabel='Variant set', 
    ylabel='(Median) % of clones (prevalence>1%) \n with positive and significative phylocorr', 
    xticks=df_['variant_set'], reduced_spines=True, rotx=90
)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'perc_sign_GBCs.png'), dpi=1000)


##



