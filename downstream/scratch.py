"""
Input to fractal, script.
"""

# Code
from mito_utils.preprocessing import *
from mito_utils.utils import *
from mito_utils.phylo import *


# Args
path_main = '/Users/IEO5505/Desktop/mito_bench'
sample = 'MDA_clones'


##


# Paths
path_data = os.path.join(path_main, 'data') 
path_output = os.path.join(path_main, 'results', 'phylo', 'output')
path_viz = os.path.join(path_main, 'results', 'phylo', 'visualization')
path_tmp = os.path.join(path_main, 'results', 'phylo', 'downstream_files')


##

# Load afm and filter with GT
afm = read_one_sample(path_data, sample, with_GBC=True)


pd.Series((~np.isnan(afm.X)).sum(axis=1)/afm.X.shape[1]).describe()


X = filter_baseline(afm).X
sns.kdeplot((X>0).sum(axis=0))
plt.show()




sns.histplot(, bins=1000)
plt.show()

import matplotlib
matplotlib.use('macOSX')

a_cells, a = filter_cells_and_vars(afm, filtering='MQuad', path_=os.getcwd())


a = nans_as_zeros(a)
AD, DP, ad_vars = get_AD_DP(a)


# From AD, DP, ad_vars to fasta
d = AFM_to_seqs(a, t=.1, method='simple_treshold')


# Phylinsic method







