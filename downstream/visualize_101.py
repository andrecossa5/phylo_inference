"""
Viz cassiopeia results. 
"""

# Code
from mito_utils.plotting_base import *
matplotlib.use('macOSX')


# Args
path_main = '/Users/IEO5505/Desktop/mito_bench'

##

# Paths
path_output = os.path.join(path_main, 'results', 'phylo')
path_viz = os.path.join(path_main, 'results', 'phylo', 'visualization')


##


# Read out_supports
df = pd.read_csv(os.path.join(path_output, 'output_df.csv'), index_col=0)
df.groupby('run').agg('median').sort_values('support')[::-1]

# Format
df[['filtering', 'metric']] = df['run'].str.split('_', expand=True).iloc[:, 2:4]
df['n_boot'] = df['run'].map(lambda x: x.split('_')[-1])
df['n_boot'] = df['run'].map(lambda x: x.split('_')[-1])
df['sample'] = df['run'].map(lambda x: '_'.join(x.split('_')[:2]))
patterns = ['greedy', 'UPMGA', 'NJ', 'max_cut', 'spectral']
df['solver'] = np.select([ df['run'].str.contains(x) for x in patterns ], patterns)
patterns = ['jacknife', 'features_resampling', 'counts_resampling']
df['boot_strategy'] = np.select([ df['run'].str.contains(x) for x in patterns ], patterns)

# Viz and quantify
fig, ax = plt.subplots(figsize=(6,5))
box(df, 'sample', 'support', ax=ax, c='white')
# strip(df, 'solver', 'support', ax=ax, c='black')
plt.show()

# Take best ones
(
    df.groupby('run')
    # .apply(lambda x: x.query('expanded_clones!="other"')['support'].median())
    .apply(lambda x: x.query('time<5')['support'].min())
    .sort_values(ascending=False)
    .head(20)
)

##