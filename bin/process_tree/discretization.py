"""
Genotype imputation for MT-SNVs.
"""

import os
from scipy.stats import binom, betabinom, poisson, nbinom
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.plotting_base import *


##


# Paths 
sample = 'MDA_clones_100'
path_main = '/Users/IEO5505/Desktop/mito_bench'
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results', 'discretization')
make_folder(path_results, sample, overwrite=False)


##


# AFM
afm = read_one_sample(path_data, sample=sample, with_GBC=True, cell_file='barcodes.txt', nmads=10)
meta = pd.read_csv(os.path.join(path_data, 'cells_meta.csv'), index_col=0)
afm.obs = afm.obs.join(meta[[ x for x in meta.columns if x not in afm.obs.columns ]])

# Filter variants: MI_TO
dataset_df, a = filter_cells_and_vars(
    afm, 
    sample_name='MDA_clones_100',
    filtering='MI_TO', 
    lineage_column='GBC',
    path_priors=None
)

# Discretization

# 1. Which muts are we looking at??
fig, axs = plt.subplots(1,2,figsize=(8,4.5))

ax = axs[0]
for x in a.var_names:
    sns.kdeplot(a[:,x].X.flatten(), ax=ax, fill=True)
ax.set(title=f'Allelic frequency distribution')

ax = axs[1]
AD, DP, var_names = get_AD_DP(a)
max_ADs = AD.A.max(axis=1)
n_positives = (a.X>.05).sum(axis=0)
ax.plot(n_positives, max_ADs, 'ko')
format_ax(ax=ax, xlabel='# positive cells', ylabel='Max alternative allele \n UMI counts (across cells)', title=f'MT-SNVs metrics diagnostics')

fig.suptitle(f'{sample}: n={a.shape[1]} muts')
fig.tight_layout()
plt.show()


##


# 2. Which distribution do these mutation fit?? 

# NB, Bin, mixtures
fig, axs = plt.subplots(1,7,figsize=(15,2.5))

for i,x in enumerate(var_names):
    ax = axs.ravel()[i]
    ax.plot(DP.A.T[:,i], AD.A.T[:,i], 'ko')
    ax.set(title=x)

fig.tight_layout()
plt.show()


# Given vector of total UMIs
total_umis_vector = DP.A.T[:,1].astype('int')

# Parameters
mutation_rate = 0.05  # Binomial model
alpha, beta = 2, 8  # Beta distribution parameters for Beta-Binomial
lambda_poisson = 50  # Poisson rate parameter
mean_nb, dispersion_nb = 50, 5  # Negative Binomial parameters

# Simulating mutated UMIs for each value in total_umis_vector
np.random.seed(0)  # For reproducibility

# Binomial Model
fitted_binomial = binom.rvs(n=total_umis_vector, p=mutation_rate)

# Beta-Binomial Model
fitted_betabinom = betabinom.rvs(n=total_umis_vector, a=alpha, b=beta)

# Poisson Model (assuming lambda is the same across all)
fitted_poisson = poisson.rvs(mu=lambda_poisson, size=len(total_umis_vector))

# Negative Binomial Model
if dispersion_nb > mean_nb:
    r_nb = mean_nb**2 / (dispersion_nb - mean_nb)
    p_nb = mean_nb / dispersion_nb
else:
    r_nb = dispersion_nb**2 / (mean_nb - dispersion_nb)
    p_nb = dispersion_nb / mean_nb
fitted_nbinom = nbinom.rvs(n=r_nb, p=1-p_nb, size=len(total_umis_vector))

# Compile the results into a DataFrame
fitted_values = {
    'Total UMIs': total_umis_vector,
    'Fitted Binomial': fitted_binomial,
    'Fitted Beta-Binomial': fitted_betabinom,
    'Fitted Poisson': fitted_poisson,
    'Fitted Negative Binomial': fitted_nbinom
}

fitted_df = pd.DataFrame(fitted_values)
fitted_df['observed'] = total_umis_vector
fitted_df.iloc[:,1:].corr()['observed']






