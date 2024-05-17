
import os
import time
import pandas as pd
from nsforest.nsforesting import mydecisiontreeevaluation
from nsforest.nsforesting import calculate_fraction

def DecisionTree(adata, cluster_header, markers_dict, medians_header = "medians_", 
                 beta = 0.5, combinations = False, use_mean = False,
                 save_supplementary = False, output_folder = "", outputfilename_prefix = ""): 
    """\
    Calculating sklearn.metrics's fbeta_score, precision_score, recall_score, and confusion_matrix for `genes_eval`. 

    Parameters
    ----------
        adata: AnnData
            Annotated data matrix.
        cluster_header: str
            Column in `adata.obs` storing cell annotation.
        markers_dict: dict
            Dictionary containing genes for each `cluster_header` (clusterName: list of markers)
        medians_header: str (default: "medians_{cluster_header}")
            Key in `adata.varm` storing median expression matrix. 
        beta: float (default: 0.5)
            `beta` parameter in sklearn.metrics's fbeta_score. 
        combinations: bool (default: True)
            Whether to find the combination of `genes_eval` with the highest fbeta_score. 
        use_mean: bool (default: False)
            Whether to use the mean (vs median) for minimum gene expression threshold. 
        save_supplementary: bool (default: False)
            Whether to save additional supplementary csvs. 
        output_folder: str (default: "")
            Output folder. Created if doesn't exist. 
        outputfilename_prefix: str (default: "")
            Prefix for all output files. 
    
    Returns
    -------
    df_results: pd.DataFrame 
        NS-Forest results. Includes classification metrics (f_score, PPV, recall, onTarget). 
    """
    # default medians_header
    if medians_header == "medians_": medians_header = "medians_" + cluster_header
    # Creating directory if does not exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print(f"Creating new directory...\n{output_folder}")

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
    
    ############################## START iterations ######################################
    cluster_list = list(markers_dict.keys())
    n_clusters = len(cluster_list)
    
    print("Number of clusters to evaluate: " + str(n_clusters))
    df_results = pd.DataFrame()
    start_time = time.time()
    
    for cl in cluster_list[:]:
        ct = list(cluster_list).index(cl) + 1
        print(f"{ct} out of {n_clusters}:")
        print(f"\t{cl}")
        print(f"\tmarker genes to be evaluated: {markers_dict[cl]}")
        
        ##=== reset parameters for this iteration!!! (for taking care of special cases) ===##
        markers = []
        for marker in markers_dict[cl]: 
            if marker in list(adata.var_names): 
                markers.append(marker)
            else: 
                print(f"cannot find {marker} in adata.var_names, excluding from DecisionTree.")

        if len(markers) == 0: continue
        
        ## Evaluation step: calculate F-beta score for gene combinations
        markers, scores = mydecisiontreeevaluation.myDecisionTreeEvaluation(adata, df_dummies, cl, markers, beta, combinations = combinations)
        print("\t" + str(markers))
        print("\t" + "fbeta: " + str(scores[0]))
        print("\t" + "PPV: " + str(scores[1]))
        print("\t" + "recall: " + str(scores[2]))

        ## return final results as dataframe
        dict_results_cl = {'clusterName': cl,
                           'clusterSize': int(scores[5]+scores[6]),
                           'f_score': scores[0],
                           'PPV': scores[1],
                           'recall': scores[2],
                           'TN': int(scores[3]),
                           'FP': int(scores[4]),
                           'FN': int(scores[5]),
                           'TP': int(scores[6]),
                           'marker_count': len(markers),
                           'markers': [markers] 
                           }
        df_results_cl = pd.DataFrame(dict_results_cl)
        df_results = pd.concat([df_results,df_results_cl]).reset_index(drop=True)
        df_results.to_csv(output_folder + outputfilename_prefix + "_results.csv", index=False)

    markers_dict = dict(zip(df_results["clusterName"], df_results["markers"]))
    on_target_ratio = calculate_fraction.markers_onTarget(adata, cluster_header, markers_dict, use_mean = use_mean, save_supplementary = save_supplementary, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)
    df_results = df_results.merge(on_target_ratio, on = "clusterName", how = "left")
    df_results.to_csv(f"{output_folder}{outputfilename_prefix}_results.csv", index=False)
    print(f"Saving final results table as...\n{output_folder}{outputfilename_prefix}_results.csv")
    print("--- %s seconds ---" % (time.time() - start_time))
    ### END iterations ###
    
    return df_results

def add_fraction(adata, df_results, cluster_header, medians_header = "medians_", use_mean = False, save_supplementary = False, output_folder = "", outputfilename_prefix = ""): 
    """\
    Calculating sklearn.metrics's fbeta_score, sklearn.metrics's prevision_score, sklearn.metrics's confusion_matrix for each `genes_eval` combination. 
    Returning set of genes and scores with highest score sum. 

    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    df_results: pd.DataFrame
        NS-Forest results. Contains classification metrics (f_score, PPV, recall, onTarget). 
    cluster_header
        Column in `adata`'s `.obs` representing cell annotation.
    medians_header: str
        Key in `adata`'s `.varm` storing median expression matrix. 
    use_mean
        Whether to use the mean or median for minimum gene expression threshold.
    output_folder
        Output folder. 
    outputfilename_prefix
        Prefix for all output files. 
    
    Returns
    -------
    df_results: pd.DataFrame of the NS-Forest results. Contains classification metrics (f_score, PPV, recall, onTarget). 
    """

    # default medians_header
    if medians_header == "medians_": medians_header = "medians_" + cluster_header

    markers_dict = dict(zip(df_results["clusterName"], df_results["markers"]))
    on_target_ratio = calculate_fraction.markers_onTarget(adata, markers_dict, cluster_header, medians_header, use_mean, save_supplementary, output_folder, outputfilename_prefix)
    if "fraction" in list(df_results.columns): del df_results["fraction"]
    if "onTarget" in list(df_results.columns): del df_results["onTarget"]
    df_results = df_results.merge(on_target_ratio, on = "clusterName", how = "left")
    df_results.to_csv(f"{output_folder}{outputfilename_prefix}_results.csv", index=False)
    print(f"Saving final results table as...\n{output_folder}{outputfilename_prefix}_results.csv")
    return df_results
