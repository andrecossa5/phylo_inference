"""
First phylo trees, prep input
"""

import os
import pickle
from Cellula.preprocessing._metrics import *
from Cellula.plotting._plotting_base import *
from Cellula.plotting._colors import *
from MI_TO.preprocessing import *


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/'
path_data = path_main + '/data/'
path_distances = path_main + '/results_and_plots/distances/'
path_class_performance = path_main + '/results_and_plots/classification_performance/'
path_results = path_main + '/results_and_plots/phylo/'
path_phylo_data = path_main + '/data/phylo_input/'


##


# Prep files for each sample

## MDA
sample = 'MDA'
path_MDA_input = path_class_performance + 'top_3/MDA/MQuad_50_50/'

# Cell x var matrix
with open(path_MDA_input + 'cell_x_var_hclust.pickle', 'rb') as f:
    cell_x_var_d = pickle.load(f)

cells = cell_x_var_d['cells']
variants = cell_x_var_d['vars']
afm = read_one_sample(path_main, sample='MDA')
X_cells_vars = afm[cells, variants].X
df = pd.DataFrame(X_cells_vars, columns=variants, index=cells)
df['GBC'] = afm.obs.loc[cells, 'GBC']
df.to_csv(path_phylo_data + 'MDA_MQuad_50_50_cell_x_vars.csv')


##


## AML
sample = 'AML'
path_AML_input = path_class_performance + 'top_3/AML/MQuad_50_50/'



pd.read_excel(path_clones + 'clones_AML_MQuad_50_50_xgboost_f1.xlsx', index_col=0).query(
    'comparison == "TTACCCTGTGCGACCCAG_vs_rest"'
)






# Cell x var matrix
with open(path_AML_input + 'cell_x_var_hclust.pickle', 'rb') as f:
    cell_x_var_d = pickle.load(f)

cells = cell_x_var_d['cells']
variants = cell_x_var_d['vars']
afm = read_one_sample(path_main, sample='AML')
X_cells_vars = afm[cells, variants].X
df = pd.DataFrame(X_cells_vars, columns=variants, index=cells)
df['GBC'] = afm.obs.loc[cells, 'GBC']
df.to_csv(path_phylo_data + 'AML_MQuad_50_50_cell_x_vars.csv')


##


## PDX
sample = 'PDX'
path_PDX_input = path_class_performance + 'top_3/PDX/pegasus_50_50/'

# Cell x var matrix
with open(path_PDX_input + 'cell_x_var_hclust.pickle', 'rb') as f:
    cell_x_var_d = pickle.load(f)

cells = cell_x_var_d['cells']
variants = cell_x_var_d['vars']
afm = read_one_sample(path_main, sample='PDX')
X_cells_vars = afm[cells, variants].X
df = pd.DataFrame(X_cells_vars, columns=variants, index=cells)
df['GBC'] = afm.obs.loc[cells, 'GBC']
df.to_csv(path_phylo_data + 'PDX_pegasus_50_50_cell_x_vars.csv')


##


# Change colors and make into .csv files
#os.listdir(path_data)
#
def change_cols(t):
    r, g, b = [ int(x*155) for x in t]
    return rgb2hex(r, g, b)

for x in ['AML', 'MDA', 'PDX']:
    with open(path_data + f'{x}_clones_colors.pkl', 'rb') as f:
        colors = pickle.load(f)
        colors = { k : change_cols(v) for k, v in colors.items() }
        pd.Series(colors).to_frame().to_csv(path_data + f'{x}_clones_colors.csv')


from colormap import rgb2hex

colors


colors