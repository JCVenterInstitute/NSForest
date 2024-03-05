
### Libraries ###
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import fbeta_score
from sklearn.metrics import precision_score
from sklearn.metrics import confusion_matrix
import time
from tqdm import tqdm
import os
import mydecisiontreeevaluation

def DecisionTree(adata, cluster_header, medians_header, binary_scores_header, 
                 marker_genes_dict = {}, beta = 0.5, exact_genes_eval = False, 
                 output_folder = "", outputfilename = ""):
    
    """\
    Calculating sklearn.metrics's fbeta_score, sklearn.metrics's prevision_score, sklearn.metrics's confusion_matrix for each `genes_eval` combination. 
    Returning set of genes and scores with highest score sum. 

    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    cluster_header
        Column in `adata`'s `.obs` representing cell annotation.
    medians_header
        Column in `adata`'s `.varm` storing median expression matrix. 
    binary_scores_header
        Column in `adata`'s `.varm` storing binary score matrix.
    marker_genes_dict
        Dictionary containing marker genes for cell annotations (clusterName: list of markers)
    beta
        Beta value in sklearn.metrics's fbeta_score. 
    exact_genes_eval
        Whether to use myDecisionTreeEvaluation on various combinations of `genes_eval`. 
    output_folder
        Output folder. 
    outputfilename
        Prefix for all output files. 
    """
    ##-----
    ## prepare adata
    ##-----
    print("Preparing data...")
    start_time = time.time()
    ## densify X from sparse matrix format
    adata.X = adata.to_df()
    ## categorial cluster labels
    adata.obs[cluster_header] = adata.obs[cluster_header].astype('category')
    ## dummy/indicator for one vs. all Random Forest model
    df_dummies = pd.get_dummies(adata.obs[cluster_header]) #cell-by-cluster
    print("--- %s seconds ---" % (time.time() - start_time))

    # Getting pre-calculated cluster_median matrix
    cluster_medians = adata.varm[medians_header].T
    # Getting pre-calculated binary_score matrix
    binary_scores = adata.varm[binary_scores_header].T
    
    ############################## START iterations ######################################
    cluster_list = list(marker_genes_dict.keys())
    n_clusters = len(cluster_list)
    
    print("Number of clusters to evaluate: " + str(n_clusters))
    df_supp = df_markers = df_results = pd.DataFrame()
    start_time = time.time()
    
    for cl in cluster_list[:]:
        ct = list(cluster_list).index(cl) + 1
        print(f"{ct} out of {n_clusters}:")

        ## cluster in iteration
        print(f"\t{cl}")
        print(f"\tmarker genes to be evaluated: {marker_genes_dict[cl]}")
        
        ##=== reset parameters for this iteration!!! (for taking care of special cases) ===##
        marker_genes = []
        # marker_genes = marker_genes_dict[cl]
        for marker in marker_genes_dict[cl]: 
            if marker in list(cluster_medians.columns): 
                marker_genes.append(marker)
            else: 
                print(f"cannot find {marker} in cluster_medians, check if marker exists in h5ad file")

        top_gene_medians = cluster_medians.loc[cl,marker_genes]
        top_rf_genes_positive = top_gene_medians[top_gene_medians>0]
        n_positive_genes = sum(top_gene_medians>0)
        
        ##=== special cases: ===##
        if n_positive_genes == 0:
            print("\t" + "No positive genes for evaluation. Skipped. ")
            continue

        ##===##
        
        ## Binary scoring step: rank genes by binary score
        top_binary_genes = pd.Series(binary_scores.loc[cl, marker_genes], index=top_rf_genes_positive.index).sort_values(ascending=False)

        ## Evaluation step: calculate F-beta score for gene combinations
        markers, scores, score_max = mydecisiontreeevaluation.myDecisionTreeEvaluation(adata, df_dummies, cl, marker_genes, beta, exact_genes_eval = exact_genes_eval)
        print("\t" + str(markers))
        print("\t" + str(score_max))

        ## return supplementary table as csv
        binary_genes_list = top_binary_genes.index.to_list()
        df_supp_cl = pd.DataFrame({'clusterName': cl,
                                   'binary_genes': binary_genes_list,
                                   'cluster_median': top_gene_medians[binary_genes_list],
                                   'binary_score': top_binary_genes[binary_genes_list]}).sort_values('binary_score', ascending=False)
        df_supp = pd.concat([df_supp,df_supp_cl]).reset_index(drop=True)
        df_supp.to_csv(output_folder + outputfilename + "_supplementary.csv", index=False)

        ## return markers table as csv
        df_markers_cl = pd.DataFrame({'clusterName': cl, 'markerGene': markers, 'score': scores[0]})
        df_markers = pd.concat([df_markers, df_markers_cl]).reset_index(drop=True)
        df_markers.to_csv(output_folder + outputfilename + "_markers.csv", index=False)

        ## return final results as dataframe
        dict_results_cl = {'clusterName': cl,
                           'clusterSize': int(scores[4]+scores[5]),
                           'f_score': scores[0],
                           'PPV': scores[1],
                           'TN': int(scores[2]),
                           'FP': int(scores[3]),
                           'FN': int(scores[4]),
                           'TP': int(scores[5]),
                           'marker_count': len(markers),
                           'NSForest_markers': [markers],
                           'binary_genes': [df_supp_cl['binary_genes'].to_list()] #for this, order is the same as the supp order
                           }
        df_results_cl = pd.DataFrame(dict_results_cl)
        df_results = pd.concat([df_results,df_results_cl]).reset_index(drop=True)
        df_results.to_csv(output_folder + outputfilename + "_results.csv", index=False)

    print("--- %s seconds ---" % (time.time() - start_time))
    ### END iterations ###
    
    return(df_results)
