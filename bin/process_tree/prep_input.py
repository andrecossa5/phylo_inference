"""
Prep input for tree processing.
"""

import os
import sys
from mito_utils.phylo import *


##


# Args
path_input = sys.argv[1]
path_tree = sys.argv[2]
path_input = '/Users/IEO5505/Desktop/MI_TO/phylo_inference/input_folder'
path_tree = '/Users/IEO5505/Desktop/MI_TO/phylo_inference/tree.newick'


##


# Read tree and set node labels
tree = read_newick(path_tree)




write_newick

