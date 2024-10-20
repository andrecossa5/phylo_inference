#!/usr/bin/python

import sys
import pickle
from sklearn.metrics import silhouette_score, normalized_mutual_info_score
from kneed import KneeLocator
from vireoSNP import BinomMixtureVB
from CCLONE.cluster.NMF import bootstrap_wNMF
from mito_utils.utils import *
from mito_utils.kNN import kNN_graph
from mito_utils.clustering import leiden_clustering
from mito_utils.clustering import custom_ARI
from mito_utils.phylo import build_tree, cut_and_annotate_tree


##


# Paths
# path_afm = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/results/MI_TO_bench/phylo_inference/distance/MDA_PT/job2/afm.h5ad'
# maxK = 10


##


def main():

    # Read afm
    afm = sc.read(sys.argv[1])
    maxK = sys.argv[2]

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

    D['leiden']['labels'] = leiden_clustering(conn, res=resolutions[np.argmax(scores)])
    D['leiden']['ARI'] = custom_ARI(afm.obs['GBC'], labels)
    D['leiden']['NMI'] = normalized_mutual_info_score(afm.obs['GBC'], labels)


    ## 


    # Vireo
    D['vireoSNP'] = {}
    _ELBO_mat = []
    for k in range(2,maxK):
        print(f'Clone n: {k}')
        _model = BinomMixtureVB(n_var=afm.shape[1], n_cell=afm.shape[0], n_donor=k)
        _model.fit(afm.layers['AD'].T, afm.layers['DP'].T, min_iter=30, max_iter=500, max_iter_pre=250, n_init=50, random_seed=1234)
        _ELBO_mat.append(_model.ELBO_inits)

    x = range(2,maxK)
    y = np.median(_ELBO_mat, axis=1)
    knee = KneeLocator(x, y).find_knee()[0]
    n_clones = knee

    _model = BinomMixtureVB(n_var=afm.shape[1], n_cell=afm.shape[0], n_donor=n_clones)
    _model.fit(afm.layers['AD'].T, afm.layers['DP'].T, min_iter=30, n_init=50, max_iter=500, max_iter_pre=250, random_seed=1234)

    clonal_assignment = _model.ID_prob
    df_ass = pd.DataFrame(
        clonal_assignment, 
        index=afm.obs_names, 
        columns=range(clonal_assignment.shape[1])
    )

    labels = []
    for i in range(df_ass.shape[0]):
        cell_ass = df_ass.iloc[i,:]
        try:
            labels.append(np.where(cell_ass>.7)[0][0])
        except:
            labels.append('unassigned')

    D['vireoSNP']['labels'] = labels
    D['vireoSNP']['ARI'] = custom_ARI(afm.obs['GBC'], labels)
    D['vireoSNP']['NMI'] = normalized_mutual_info_score(afm.obs['GBC'], labels)         


    ##


    # MI_TO
    tree = build_tree(afm, bin_method='vanilla', metric='jaccard', precomputed=True, solver='UPMGA')
    tree, _,_ = cut_and_annotate_tree(tree)
    D['MI_TO']['labels'] = tree.cell_meta['MT_clone'].to_list()
    D['MI_TO']['ARI'] = custom_ARI(tree.cell_meta['GBC'], tree.cell_meta['MT_clone'])
    D['MI_TO']['NMI'] = normalized_mutual_info_score(tree.cell_meta.dropna()['GBC'], tree.cell_meta.dropna()['MT_clone'])      


    ##


    # CClone
    afm.layers['ALT'] = afm.layers['AD'].A
    afm.layers['REF'] = afm.layers['DP'].A - afm.layers['AD'].A
    afm.X = afm.X.A

    scores = []
    for k in range(2,maxK):
        print(f'k={k}')
        afm = bootstrap_wNMF(afm, k, n_bootstrap=1, max_cycles=100, parallel=True, n_jobs=-1, mode='10X')
        labels = np.argmax(afm.obsm['C'], axis=1)
        silhouette_avg = silhouette_score(afm.obsp['distances'].A, labels, metric='precomputed')
        scores.append(silhouette_avg)

    k = list(range(2,maxK))[np.argmax(scores)]
    afm = bootstrap_wNMF(afm, k, n_bootstrap=1, max_cycles=100, parallel=True, n_jobs=-1, mode='10X')
    labels = np.argmax(afm.obsm['C'], axis=1).tolist()

    D['CClone']['labels'] = labels
    D['CClone']['ARI'] = custom_ARI(afm.obs['GBC'], labels)
    D['CClone']['NMI'] = normalized_mutual_info_score(afm.obs['GBC'], labels)   

    
    # Save
    with open(os.path.join(os.path.dirname(sys.argv[1]), 'bench_clonal_recontruction.pickle'), 'wb') as f:
        pickle.dump(D, f)


##


# Run
if __name__ == '__main__':
    main()


##