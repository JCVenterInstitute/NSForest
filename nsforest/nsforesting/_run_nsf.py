
import os
import time
import pandas as pd
from nsforest.nsforesting import myrandomforest
from nsforest.nsforesting import mydecisiontreeevaluation
from nsforest.nsforesting import calculate_fraction

def NSForest(adata, cluster_header, medians_header = "medians_", binary_scores_header = "binary_scores_", 
             cluster_list = [], gene_selection = "BinaryFirst_high",
             n_trees = 1000, n_jobs = -1, beta = 0.5, n_top_genes = 15, n_binary_genes = 10, n_genes_eval = 6,
             save_supplementary = False, output_folder = "", outputfilename_prefix = ""):
    """\
    Performs the main NS-Forest algorithm to find a list of NS-Forest markers for each `cluster_header`. 

    Parameters
    ----------
        adata: AnnData
            Annotated data matrix.
        cluster_header: str
            Column in `adata.obs` storing cell annotation.
        medians_header: str (default: "medians_{cluster_header}")
            Key in `adata.varm` storing median expression matrix. 
        binary_scores_header: str (default: "binary_scores_{cluster_header}")
            Key in `adata.varm` storing binary score matrix.
        cluster_list: list (default: all clusters)
            For subsetting by specified cell annotations. Used for parallelizing NSForest. 
        gene_selection: str (default: "BinaryFirst_high")
            Level of filtering genes by binary score. Options: [None, "BinaryFirst_high", "BinaryFirst_moderate", "BinaryFirst_low"]. None includes all genes. BinaryFirst_high includes genes with binary scores > 2 std. BinaryFirst_moderate includes genes with binary scores > 1 std. BinaryFirst_low includes genes with binary scores > median. 
        n_trees: int (default: 1000)
            `n_estimators` parameter in sklearn.ensemble's RandomForestClassifier. 
        n_jobs: int (default: -1)
            `n_jobs` parameter in sklearn.ensemble's RandomForestClassifier. 
        beta: float (default: 0.5)
            `beta` parameter in sklearn.metrics's fbeta_score. 
        n_top_genes: int (default: 15)
            Taking the top `n_top_genes` genes ranked by sklearn.ensemble's RandomForestClassifier as input for sklearn.tree's DecisionTreeClassifier. 
        n_binary_genes: int (default: 10)
            Taking the top `n_binary_genes` genes ranked by binary score for supplementary table output. 
        n_genes_eval: int (default: 6)
            Taking the top `n_genes_eval` genes ranked by binary score as input for sklearn.tree's DecisionTreeClassifier. 
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

    # default medians_header and binary_scores_header
    if medians_header == "medians_": medians_header = "medians_" + cluster_header
    if binary_scores_header == "binary_scores_": binary_scores_header = "binary_scores_" + cluster_header
    # Creating directory if does not exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print(f"Creating new directory...\n{output_folder}")

    ##-----
    ## prepare adata
    ##-----
    print("Preparing adata...")
    start_time = time.time()
    ## densify X from sparse matrix format
    adata.X = adata.to_df()
    ## categorial cluster labels
    adata.obs[cluster_header] = adata.obs[cluster_header].astype('category')
    ## dummy/indicator for one vs. all Random Forest model
    df_dummies = pd.get_dummies(adata.obs[cluster_header]) #cell-by-cluster
    ## get number of clusters
    n_total_clusters = len(df_dummies.columns)
    print("--- %s seconds ---" % (time.time() - start_time))

    # Getting cluster_median matrix
    cluster_medians = adata.varm[medians_header].T
    # Getting binary_score matrix
    binary_scores = adata.varm[binary_scores_header].T

    ####################### SET THRESHOLD/CUTOFF VALUE FOR BINARY FIRST TO USE #############################
    print("Pre-selecting genes based on binary scores...")
    if gene_selection == "BinaryFirst_mild":
        # if this line below is too slow, could also convert binary_scores df to np array and use np statistics
        median = binary_scores.stack().median()
        threshold = median
        print("\t", "Threshold (median):", threshold)
        
    elif gene_selection == "BinaryFirst_moderate":
        mean = binary_scores.stack().mean()
        stddev = binary_scores.stack().std()
        threshold = mean + stddev
        print("\t", "Threshold (mean + 1 * std):", threshold)
    
    elif gene_selection == "BinaryFirst_high":
        mean = binary_scores.stack().mean()
        stddev = binary_scores.stack().std()
        threshold = mean + 2 * stddev
        print("\t", "Threshold (mean + 2 * std):", threshold)
    
    else: 
        threshold = 0
        print("\t", "Threshold (all genes (default)):", threshold)
    
    ## create dummy matrix based on if binary score is greater than threshold
    binary_dummies = binary_scores >= threshold
    binary_dummies = binary_dummies*1
    binary_dummies = pd.DataFrame(binary_dummies, index=binary_scores.index, columns=binary_scores.columns) #cluster-by-gene
    ## average number of genes selected
    binary_dummies_sum = binary_dummies.stack().sum()
    print("\tAverage number of genes after gene_selection in each cluster:", binary_dummies_sum/n_total_clusters)
    ## write number of genes selected for each cluster
    if save_supplementary: 
        print(f"Saving number of genes selected per cluster as...\n{output_folder}{outputfilename_prefix}_gene_selection.csv")
        pd.DataFrame(binary_dummies.sum(axis=1)).to_csv(output_folder + outputfilename_prefix + "_gene_selection.csv") # , header=["n_genes_selected"]
    
    ############################## START iterations ######################################
    if cluster_list == []:
        cluster_list = df_dummies.columns
    n_clusters = len(cluster_list)
    
    print ("Number of clusters to evaluate: " + str(n_clusters))
    df_supp = df_markers = df_results = pd.DataFrame()
    start_time = time.time()
    
    for cl in cluster_list: 
        ct = list(cluster_list).index(cl) + 1
        print(f"{ct} out of {n_clusters}:")
        print(f"\t{cl}")
        
        ##=== reset parameters for this iteration!!! (for taking care of special cases) ===##
        n_binary_genes_cl = n_binary_genes
        n_genes_eval_cl = n_genes_eval

        top_rf_genes = myrandomforest.myRandomForest(adata, df_dummies, cl, n_trees, n_jobs, n_top_genes, binary_dummies)      
        ## filter out negative genes by thresholding median>0 ==> to prevent dividing by 0 in binary score calculation
        top_gene_medians = cluster_medians.loc[cl,top_rf_genes.index]
        top_rf_genes_positive = top_gene_medians[top_gene_medians>0]
        n_positive_genes = sum(top_gene_medians>0)
        
        ##=== special cases: ===##
        if n_positive_genes == 0:
            print("\t" + "No positive genes for evaluation. Skipped. Optionally, consider increasing n_top_genes.")
            continue

        if n_positive_genes < n_binary_genes:
            print("\t" + f"Only {n_positive_genes} out of {n_top_genes} top Random Forest features with median > 0 will be further evaluated.")
            n_binary_genes_cl = n_positive_genes
            n_genes_eval_cl = min(n_positive_genes, n_genes_eval)
        ##===##
        
        ## Binary scoring step: rank genes by binary score
        top_binary_genes = pd.Series(binary_scores.loc[cl], index=top_rf_genes_positive.index).sort_values(ascending=False)

        ## Evaluation step: calculate F-beta score for gene combinations
        genes_eval = top_binary_genes.index[:n_genes_eval_cl].to_list()
        markers, scores = mydecisiontreeevaluation.myDecisionTreeEvaluation(adata, df_dummies, cl, genes_eval, beta)
        print("\t" + str(markers))
        print("\t" + "fbeta: " + str(scores[0]))
        print("\t" + "PPV: " + str(scores[1]))
        print("\t" + "recall: " + str(scores[2]))

        ## return supplementary table as csv
        binary_genes_list = top_binary_genes.index[:n_binary_genes_cl].to_list()
        df_supp_cl = pd.DataFrame({'clusterName': cl,
                                'binary_genes': binary_genes_list,
                                'rf_feature_importance': top_rf_genes[binary_genes_list],
                                'cluster_median': top_gene_medians[binary_genes_list],
                                'binary_score': top_binary_genes[binary_genes_list]}).sort_values('binary_score', ascending=False)
        df_supp = pd.concat([df_supp,df_supp_cl]).reset_index(drop=True)
        if save_supplementary: 
            df_supp.to_csv(output_folder + outputfilename_prefix + "_supplementary.csv", index=False)

        ## return markers table as csv
        df_markers_cl = pd.DataFrame({'clusterName': cl, 'markerGene': markers, 'score': scores[0]})
        df_markers = pd.concat([df_markers, df_markers_cl]).reset_index(drop=True)
        if save_supplementary: 
            df_markers.to_csv(output_folder + outputfilename_prefix + "_markers.csv", index=False)

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
                           'NSForest_markers': [markers],
                           'binary_genes': [df_supp_cl['binary_genes'].to_list()] #for this, order is the same as the supp order
                           }
        df_results_cl = pd.DataFrame(dict_results_cl)
        df_results = pd.concat([df_results,df_results_cl]).reset_index(drop=True)
        df_results.to_csv(f"{output_folder}{outputfilename_prefix}_results.csv", index=False)

    ### END iterations ###

    if save_supplementary: 
        print(f"Saving supplementary table as...\n{output_folder}{outputfilename_prefix}_supplementary.csv")
        print(f"Saving markers table as...\n{output_folder}{outputfilename_prefix}_markers.csv")

    print(f"Saving results table as...\n{output_folder}{outputfilename_prefix}_results.csv")
    markers_dict = dict(zip(df_results["clusterName"], df_results["NSForest_markers"]))
    on_target_ratio = calculate_fraction.markers_onTarget(adata, cluster_header, markers_dict, medians_header, save_supplementary = save_supplementary, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)
    df_results = df_results.merge(on_target_ratio, on = "clusterName", how = "left")
    df_results.to_csv(f"{output_folder}{outputfilename_prefix}_results.csv", index=False)
    print(f"Saving final results table as...\n{output_folder}{outputfilename_prefix}_results.csv")

    print("--- %s seconds ---" % (time.time() - start_time))
    
    return df_results
