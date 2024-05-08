"""
Longitudinal analysis top clones.
"""

import os
import pickle
import anndata
from mito_utils.preprocessing import *
from mito_utils.phylo_plots import *
from mito_utils.phylo import *
from mito_utils.heatmaps_plots import *
from mito_utils.plotting_base import *
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/mito_bench'
path_data = os.path.join(path_main, 'data')
path_variants = os.path.join(path_main, 'results/supervised_clones/variants.pickle')
path_results = os.path.join(path_main, 'results/phylo_inference')
path_top_trees = os.path.join(path_results, 'top_trees')

# Colors
with open(os.path.join(path_data, 'clones_colors_sc.pickle'), 'rb') as f:
    colors = pickle.load(f)


##
    

# Read variants
with open(path_variants, 'rb') as f:
    variants = pickle.load(f)
variants = list(set(variants[('MDA_PT', 'MQuad')]) | set(variants[('MDA_lung', 'MQuad')]))

# Read AFMs
PT = read_one_sample(path_data, 'MDA_PT')
_, PT = filter_cells_and_vars(PT, variants=variants)
lung = read_one_sample(path_data, 'MDA_lung')
_, lung = filter_cells_and_vars(lung, variants=variants)

# Subsample and merge
weights = PT.obs['GBC'].value_counts(normalize=True)
np.random.seed(1)
PT = PT[PT.obs.sample(n=1000, replace=False, weights=PT.obs['GBC'].map(weights)).index,:]
weights = lung.obs['GBC'].value_counts(normalize=True)
np.random.seed(1)
lung = lung[lung.obs.sample(n=1000, replace=False, weights=lung.obs['GBC'].map(weights)).index,:]
afm = anndata.concat([PT, lung])
afm.obs['origin'] = np.where(afm.obs['sample']=='MDA_PT', 'PT', 'lung')

# Build UPMGA tree, all clones
afm = nans_as_zeros(afm)
obs_tree = build_tree(afm)

# Get covariate dfs
cov = 'origin'
path_sample = os.path.join(path_top_trees, 'couple')
assert all([ x in afm.obs_names for x in obs_tree.leaves ])
afm = afm[obs_tree.leaves,:].copy()
obs_tree.cell_meta[cov] = afm.obs['sample'].map(lambda x: x.split('_')[-1])
counts = afm.obs[cov].value_counts().to_frame('n_cells')
counts['to_numeric'] = range(1,counts.shape[0]+1)
cov_df = afm.obs[cov].map(counts['to_numeric'].to_dict()).to_frame('cov')
cov_df.to_csv(os.path.join(path_sample, f'cov_df.csv'))
counts.to_csv(os.path.join(path_sample, f'mapping.csv'))
with open(os.path.join(path_sample, f'tree.newick'), 'w') as f:
    f.write(obs_tree.get_newick(record_branch_lengths=True))


##


# All couple, subsamples
origin_colors = { 
    "PT" : sns.color_palette("Spectral_r")[0], 
    "lung" : sns.color_palette("Spectral_r")[-1] 
}

# One clone at the time
afm = anndata.concat([PT, lung])
afm.obs['origin'] = np.where(afm.obs['sample']=='MDA_PT', 'PT', 'lung')
clones = (
    afm.obs.groupby(['GBC', 'origin'])
    .size().to_frame('n').reset_index()
    .pivot_table(index='GBC', values='n', columns='origin')
    .dropna()
    .assign(mean=lambda x: x.sum(axis=1)/2)
    .sort_values('mean', ascending=False)
    .query('PT>=10 and lung>=10')
    .index
)

# Single clone

# Viz
fig = plt.figure(figsize=(9,3.5))

for i, clone in enumerate(clones):

    ax = fig.add_subplot(1,3,i+1)
    afm = anndata.concat([PT, lung])
    afm.obs['origin'] = np.where(afm.obs['sample']=='MDA_PT', 'PT', 'lung')
    cells = afm.obs.query('GBC==@clone').index
    afm = nans_as_zeros(afm[cells])
    obs_tree = build_tree(afm[cells])

    plot_tree(
        obs_tree, ax=ax,
        meta=['origin'], 
        categorical_cmap_annot=origin_colors, 
        colorstrip_width=4,
        colorstrip_spacing=0.001,
        extend_branches=True,
        internal_node_kwargs={'markersize':0},
        orient=90,
    )
    ax.set(title=clone)

fig.tight_layout()
fig.savefig(os.path.join(path_top_trees, 'couple', 'top_longitudinal_clones_met.pdf'), dpi=500)


##