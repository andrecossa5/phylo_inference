"""
Visualiza phylo output.
"""

import os
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.phylo_plots import *
from mito_utils.heatmaps_plots import *
from mito_utils.plotting_base import *
from mito_utils.metrics import *
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/mito_bench/'
path_results = os.path.join(path_main, 'results/phylo')


##


# # Collect metrics
# 
# # 1. precomputed
# d = traverse_and_extract_flat(path_results, file_name='final_tree_stats.csv')
# df_stats = pd.concat([   
#     df.assign(sample=sample, variant_set=variant_set)
#     .reset_index()
#     .rename(columns={'index':'metric', 'metric':"value"}) \
#     for (sample, variant_set), df in d.items() 
# ])
#     
# ##
# 
# # 2. Lineage association, continuos
# d = traverse_and_extract_flat(path_results, file_name='lineage_association.csv')
# 
# L = []
# for (sample, variant_set), df in d.items():
#     df['corr_type'] = np.where(df['Var1'] == df['Var2'], 'auto', 'cross')
#     df = (
#         df.query('corr_type=="auto"').drop(['Var2', 'corr_type'], axis=1).rename(columns={'Var1':'GBC'}).set_index('GBC')
#         .join(
#             pd.read_csv(os.path.join(path_results, sample, variant_set, 'filtered_input', 'meta.csv'), index_col=0)
#             ['GBC'].value_counts(normalize=True).to_frame('prevalence')
#         )
#     )
#     L.append(df.assign(sample=sample, variant_set=variant_set))
# 
# df_perc_significant = (
#     pd.concat(L)
#     .groupby(['variant_set', 'sample'])
#     .apply( lambda x: np.sum((x['Z']>10) & (x['p'] <= 0.05) & (x['prevalence'] >= .01)) / (x['prevalence'] >= .01).sum() )
#     .to_frame('value').reset_index()
#     .assign(metric='perc_significant')
# )
# 
# ##
# 
# # 3. n cells analyzed, nGT clones, NMI, ARI
# cells_meta = pd.read_csv(os.path.join(path_main, 'data', 'cells_meta.csv'), index_col=0)
# d = traverse_and_extract_flat(path_results, file_name='nodes.csv')
# 
# L = []
# for (sample, variant_set), df_nodes in d.items():
#     meta_sample = cells_meta.query('sample == @sample')
#     cell_sample = df_nodes.loc[lambda x: ~x['cell'].isna()][['cell', 'clonal_ancestor']].drop_duplicates().set_index('cell')
#     perc_cells_analyzed = cell_sample.shape[0] / meta_sample.shape[0]
#     nGT_clones = meta_sample['GBC'].nunique()
#     cells = list(set(meta_sample.index) & set(cell_sample.index))
#     x = meta_sample.loc[cells, 'GBC'].values
#     y = cell_sample.loc[cells, 'clonal_ancestor']
#     assert x.size == y.size
#     ARI = adjusted_rand_score(x,y)
#     NMI = adjusted_mutual_info_score(x,y)
#     L.append([perc_cells_analyzed, nGT_clones, NMI, ARI, sample, variant_set])
# 
# df_last_metrics = (
#     pd.DataFrame(L, columns=['perc_cells_analyzed', 'nGT_clones', 'NMI', 'ARI', 'sample', 'variant_set'])
#     .melt(id_vars=['sample', 'variant_set'], var_name='metric', value_name='value')
# )
# 
# 
# ##
# 
# 
# # Concat
# phylo_df = pd.concat([ 
#     df_stats.reset_index(drop=True), 
#     df_perc_significant.reset_index(drop=True), 
#     df_last_metrics.reset_index(drop=True) 
# ])
# phylo_df.to_csv(os.path.join(path_results, 'phylo_df.csv'))


##


# Read df metrics
df = pd.read_csv(os.path.join(path_results, 'phylo_df.csv'), index_col=0)
interesting_metrics = [ 'char_tree_dist_corr', 'char_dist_boot', 'median_support', 'perc_cells_analyzed', 'perc_significant', 'NMI' ]
d_rename = {
    'char_tree_dist_corr' : 'Character- vs tree-based \n distance correlation',
    'char_dist_boot' : 'Character distance correlation \n upon bootstrapping',
    'median_support' : 'Median internal nodes support \n upon bootstrapping',
    'perc_cells_analyzed' : 'Fraction of cells in the final tree',
    'perc_significant' : 'Fraction of clones (>1% prevalence) \n with significant and positive phylocorr',
    'NMI' : 'Normalized mutual information'
}


df = df.query('metric in @interesting_metrics')

# Colors 
df['sample'] = pd.Categorical(df['sample'], categories=['AML_clones', 'MDA_clones', 'MDA_PT', 'MDA_lung'])
sample_colors = create_palette(df, 'sample', 'tab10')

# Here we go
fig, axs = plt.subplots(2,3,figsize=(10,8))

for ax, metric in zip(axs.ravel(), interesting_metrics):
    df_ = df.query('metric==@metric')
    print(df_.shape)
    order_ = df_.groupby('variant_set')['value'].median().sort_values(ascending=False).index
    df_['variant_set'] = pd.Categorical(df_['variant_set'], categories=order_)
    box(df_, x='variant_set', y='value', c='white', ax=ax)
    strip(df_, x='variant_set', y='value', by='sample', c=sample_colors, ax=ax)
    format_ax(ax, ylabel='value', xlabel='', rotx=90, title=d_rename[metric], reduced_spines=True, title_size=11)

add_legend(label='Sample', ax=axs[0,2], colors=sample_colors, loc='upper right', bbox_to_anchor=(1,1),
           artists_size=9, label_size=9, ticks_size=7)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'phylo_performance.png'), dpi=1000)


##