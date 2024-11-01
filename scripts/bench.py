#!/usr/bin/python

import sys
import pickle
from sklearn.metrics import silhouette_score, normalized_mutual_info_score
from kneed import KneeLocator
from vireoSNP import BinomMixtureVB
from CCLONE.cluster.NMF import get_wNMF_matrices, NMF_weighted, orth_score
from mito_utils.utils import *
from mito_utils.kNN import kNN_graph
from mito_utils.clustering import leiden_clustering
from mito_utils.metrics import custom_ARI
from mito_utils.phylo import build_tree, MiToTreeAnnotator


##


# Paths
# path_afm = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/results/MI_TO_bench/phylo_inference/distance/MDA_PT/job2/afm.h5ad'
# maxK = 50


##


def main():

    # Read afm
    # afm = sc.read(path_afm)
    afm = sc.read(sys.argv[1])
    maxK = int(sys.argv[2])

    # Get 
    D = {}

    # Leiden
    D['leiden'] = {}
    scores = []
    resolutions = np.linspace(0.5,2.5,50)
    for res in resolutions:
        _, _, conn = kNN_graph(D=afm.obsp['distances'].A, k=15, from_distances=True)
        labels = leiden_clustering(conn, res=res)
        silhouette_avg = silhouette_score(afm.obsp['distances'].A, labels, metric='precomputed')
        scores.append(silhouette_avg)

    labels = leiden_clustering(conn, res=resolutions[np.argmax(scores)])
    test = np.isnan(labels)
    D['leiden']['labels'] = labels
    D['leiden']['NA'] = test.sum() / labels.size
    D['leiden']['ARI'] = custom_ARI(afm.obs['GBC'][~test], labels[~test])
    D['leiden']['NMI'] = normalized_mutual_info_score(afm.obs['GBC'][~test], labels[~test])


    ## 


    # Vireo
    D['vireoSNP'] = {}
    _ELBO_mat = []
    for k in range(2,maxK+1):
        print(f'Clone n: {k}')
        _model = BinomMixtureVB(n_var=afm.shape[1], n_cell=afm.shape[0], n_donor=k)
        _model.fit(afm.layers['AD'].T, afm.layers['DP'].T, min_iter=30, max_iter=500, max_iter_pre=250, n_init=50, random_seed=1234)
        _ELBO_mat.append(_model.ELBO_inits)

    _ELBO_mat = np.row_stack(_ELBO_mat)
    x = range(2,maxK+1)
    y = np.median(_ELBO_mat, axis=1)
    knee = KneeLocator(x, y).find_knee()[0]
    n_clones = knee

    _model = BinomMixtureVB(n_var=afm.shape[1], n_cell=afm.shape[0], n_donor=n_clones)
    _model.fit(afm.layers['AD'].T, afm.layers['DP'].T, min_iter=30, n_init=50, max_iter=500, max_iter_pre=250, random_seed=1234)

    clonal_assignment = _model.ID_prob
    idx = clonal_assignment.argmax(axis=1)
    labels = np.zeros(idx.shape)
    for i,clone_idx in enumerate(idx):
        labels[i] = clone_idx if clonal_assignment[i,clone_idx]>=.7 else np.nan

    test = np.isnan(labels)
    D['vireoSNP']['labels'] = labels
    D['vireoSNP']['NA'] = test.sum() / labels.size
    D['vireoSNP']['ARI'] = custom_ARI(afm.obs['GBC'][~test], labels[~test])
    D['vireoSNP']['NMI'] = normalized_mutual_info_score(afm.obs['GBC'][~test], labels[~test])       


    ##


    # MI_TO
    D['MI_TO'] = {}
    tree = build_tree(afm, precomputed=True, solver='UPMGA')
    tree, _,_ = MiToTreeAnnotator(tree)
    labels = tree.cell_meta['MT_clone'].values
    test = tree.cell_meta['MT_clone'].isna().values
    D['MI_TO']['labels'] = labels
    D['MI_TO']['NA'] = test.sum() / labels.size
    D['MI_TO']['ARI'] = custom_ARI(afm.obs['GBC'][~test], labels[~test])
    D['MI_TO']['NMI'] = normalized_mutual_info_score(afm.obs['GBC'][~test], labels[~test])


    ##


    # CClone
    afm.layers['ALT'] = afm.layers['AD'].A
    afm.layers['REF'] = afm.layers['DP'].A - afm.layers['AD'].A
    afm.X = afm.X.A

    D['CClone'] = {}

    scores = []
    for k in range(2,maxK+1):

        print(f'k={k}')
        afm = get_wNMF_matrices(afm, mode='10X')
        C, _ = NMF_weighted(afm, k=k, max_cycles=100, parallel=True, n_jobs=-1)
        labels = np.argmax(C, axis=1)
        scores.append(orth_score(C))

    k = list(range(2,maxK+1))[np.argmin(scores)]
    C, _ = NMF_weighted(afm, k=k, max_cycles=100, parallel=True, n_jobs=-1)
    labels = np.argmax(C, axis=1)
    test = np.isnan(labels)
    D['CClone']['labels'] = labels
    D['CClone']['NA'] = test.sum() / labels.size
    D['CClone']['ARI'] = custom_ARI(afm.obs['GBC'][~test], labels[~test])
    D['CClone']['NMI'] = normalized_mutual_info_score(afm.obs['GBC'][~test], labels[~test])  
    
    # Save
    with open(os.path.join(os.path.dirname(sys.argv[1]), 'bench_clonal_recontruction.pickle'), 'wb') as f:
        pickle.dump(D, f)


##


# Run
if __name__ == '__main__':
    main()


##