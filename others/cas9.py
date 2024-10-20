"""
Cas9 comparision.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
from cassiopeia.preprocess import lineage_utils
from mito_utils.preprocessing import *
from mito_utils.phylo import *


# Paths
path_data = '/Users/IEO5505/Desktop/MI_TO/cas9/KPTracer-Data'

# Read data
allele_table = pd.read_csv(
    os.path.join(path_data, "KPTracer.alleleTable.FINAL.txt"), sep="\t", index_col=0
)
all_tumors = allele_table["Tumor"].unique()
primary_nt_tumors = [
    tumor
    for tumor in all_tumors
    if "NT" in tumor and tumor.split("_")[2].startswith("T")
]
primary_nt_allele_table = allele_table[allele_table["Tumor"].isin(primary_nt_tumors)]

# Priors
indel_priors = cs.pp.compute_empirical_indel_priors(
    allele_table, grouping_variables=["intBC", "MetFamily"]
)
indel_priors.sort_values(by="count", ascending=False).head()

# One tumor
tumor = "3726_NT_T1"
tumor_allele_table = primary_nt_allele_table[primary_nt_allele_table["Tumor"] == tumor]
n_cells = tumor_allele_table["cellBC"].nunique()
n_intbc = tumor_allele_table["intBC"].nunique()

# To character matrix
(
    character_matrix,
    priors,
    state_to_indel,
) = cs.pp.convert_alleletable_to_character_matrix(
    tumor_allele_table, allele_rep_thresh=0.9, mutation_priors=indel_priors
)

# Tree and plot
tree = cs.data.CassiopeiaTree(character_matrix=character_matrix, priors=priors)
greedy_solver = cs.solver.UPGMASolver()
greedy_solver.solve(tree)

cs.pl.plot_matplotlib(tree, orient="down", allele_table=tumor_allele_table)
plt.show()



