"""
Mutations redundancy script.
"""

# Code
import matplotlib
from mito_utils.preprocessing import *
from mito_utils.phylo import *
from mito_utils.plotting_base import *
from mito_utils.diagnostic_plots import *
from matplotlib.gridspec import GridSpec
matplotlib.use('macOSX')


##


# Read AFM and filter vars
path_main = '/Users/IEO5505/Desktop/mito_bench'
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results', 'phylo', 'visualization')

sample = 'MDA_clones'
filtering = 'GT'

# Filter AFM
afm = read_one_sample(path_data, sample, with_GBC=True)

if filtering == 'GT':
    _, a = filter_afm_with_gt(afm, min_cells_clone=5)
elif filtering == 'base_filter':
    a_cells = filter_cells_coverage(afm)
    a = filter_baseline(a_cells)
    a = remove_excluded_sites(a) 
else:
    _, a = filter_cells_and_vars(
        afm, filtering=filtering, path_=os.getcwd(), nproc=8
    )

a = nans_as_zeros(a)
AD, DP, _ = get_AD_DP(a)


##


# Redundancy across sites
L = []
for t in np.linspace(.001, .1, 10):
    L.append(
        a.obs
        .assign(haplotype=AFM_to_seqs(a, t=t))
        .groupby(['haplotype', 'GBC'])
        .size().loc[lambda x:x>0].reset_index()
        .groupby('haplotype').size()
        .to_frame('n clones')
        .assign(t=np.round(t, 3))
    )
df_ = pd.concat(L).reset_index()


##


# MT-genetic diversity plot
fig = plt.figure(figsize=(12, 4.5), constrained_layout=True)
gs = GridSpec(1, 3, figure=fig, width_ratios=[2,2,3.5])

# n unique haplotypes
ax = fig.add_subplot(gs[0])
scatter(
    df_.groupby('t').size().to_frame('n_unique_haplotypes').reset_index(),
    x='t', y='n_unique_haplotypes', c='k', s=100, ax=ax
)
format_ax(
    ax=ax, xticks_size=8, xlabel='AF treshold', ylabel='n distinct haplotypes', 
    title='n MT-haplotypes', reduced_spines=True
)
ax.grid(axis='y')

# Redundancy within cells
ax = fig.add_subplot(gs[1])
scatter(
    df_.groupby('t').size().map(lambda x:x/a.shape[0]).to_frame('ratio').reset_index(),
    x='t', y='ratio', c='k', s=100, ax=ax
)
format_ax(
    ax=ax, xticks_size=8, xlabel='AF treshold', ylabel='n distinct haplotypes / n cells', 
    title='MT-haplotype redundancy', reduced_spines=True
)
ax.grid(axis='y')

# Sharedness among clones
ax = fig.add_subplot(gs[2])
box(df_, x='t', y='n clones', c='white', ax=ax)
strip(df_, x='t', y='n clones', ax=ax, c='k') #, by='haplotype')
format_ax(
    ax=ax, xticks_size=8, xlabel='AF treshold', ylabel='n clones with the {x} haplotype', 
    title='MT-haplotype sharedness', reduced_spines=True
)

fig.suptitle(f'MT-genetic diversity: {filtering} (n={a.shape[1]} MT-SNVs)')
fig.savefig(os.path.join(path_results, f'{filtering}_muts_sharedness.png'), dpi=300)


##


sample = 'MDA_clones'


# Filter AFM
afm = read_one_sample(path_data, sample, with_GBC=True)

filters = ['base_filter', 'MQuad', 'GT']
colors = ['grey', 'green', 'r']

# Get AFMs
AFMs = {}
for filtering in filters:
    if filtering == 'GT':
        _, a = filter_afm_with_gt(afm, min_cells_clone=5)
    elif filtering == 'base_filter':
        a_cells = filter_cells_coverage(afm)
        a = filter_baseline(a_cells)
        a = remove_excluded_sites(a) 
    else:
        _, a = filter_cells_and_vars(
            afm, filtering=filtering, path_=os.getcwd(), nproc=8
        )
    AFMs[filtering] = a

# Figure
fig, axs = plt.subplots(1,3, figsize=(15,5))

j = 0
for filtering, color, in zip(filters, colors):

    a = nans_as_zeros(AFMs[filtering])
    to_plot = a.X.copy()
    to_plot[np.isnan(to_plot)] = 0

    for i in range(to_plot.shape[1]):
        x = to_plot[:, i]
        x = np.sort(x)
        axs[j].plot(x, '-', color=color, linewidth=0.5)
    axs[j].set_ylim((-0.05,1.05))
    format_ax(
        ax=axs[j], 
        title=f'Ranked AF distributions: {filtering} (n={a.shape[1]})', 
        xlabel='Cell rank', ylabel='AF'
    )
    j += 1

fig.tight_layout()
fig.savefig(os.path.join(path_results, f'muts_AF_spectrum.png'), dpi=300)