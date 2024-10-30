#!/usr/bin/python

# Build cassiopeia

########################################################################
 
# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='annotate_tree',
    description="Process and annotate tree."
)

# Add arguments
my_parser.add_argument(
    '--afm', 
    type=str,
    default=None,
    help='Path to afm in afm.h5ad format. Default: .. .'
)

my_parser.add_argument(
    '--tree', 
    type=str,
    default=None,
    help='Path to tree in .newick format. Default: .. .'
)


##


# Parse arguments
args = my_parser.parse_args()
path_tree = args.tree
path_afm = args.afm

##


########################################################################

# Preparing run: import code, prepare directories, set logger

# Code
import pickle
from mito_utils.utils import *
from mito_utils.phylo import *


##


########################################################################

# Main
def main():

    # Prep annot
    afm = sc.read(path_afm)
    cell_meta = afm.obs
    char_meta = afm.var
    X_raw = pd.DataFrame(afm.X.A, index=cell_meta.index, columns=char_meta.index)
    X_bin = pd.DataFrame(afm.layers['bin'].A, index=cell_meta.index, columns=char_meta.index)
    D = pd.DataFrame(afm.obsp['distances'].A, index=cell_meta.index, columns=cell_meta.index)

    # Load tree
    tree = read_newick(path_tree, X_raw=X_raw, X_bin=X_bin, D=D, meta=cell_meta)

    # Cut and annotate tree
    tree, _, _ = MiToTreeAnnotator(tree)

    # Write as pickle
    with open('annotated_tree.pickle', 'wb') as f:
        pickle.dump(tree, f)
    

    ##


#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################