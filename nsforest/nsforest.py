
### Libraries ###
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import fbeta_score
from sklearn.metrics import precision_score
from sklearn.metrics import confusion_matrix
import itertools
import time
from tqdm import tqdm
import os
import myrandomforest
import mydecisiontreeevaluation
import utils

# v4.0 includes new parameter, "gene_selection," which determines whether BinaryFirst is used or not and its cutoff value
def NSForest(adata, cluster_header, cluster_list=None,
             medians_header=None, binary_scores_header=None, gene_selection=None,
             n_trees=1000, n_jobs=-1, beta=0.5, n_top_genes=15, n_binary_genes=10, n_genes_eval=6,
             output_folder="NSForest_outputs/", outputfilename="default"):
    
    ## set up output folder
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    ##-----
    ## prepare adata
    ##-----
    print("Preparing data...")
    first_start_time = time.time()
    ## densify X from sparse matrix format
    adata.X = adata.to_df()
    ## categorial cluster labels
    adata.obs[cluster_header] = adata.obs[cluster_header].astype('category')
    ## dummy/indicator for one vs. all Random Forest model
    df_dummies = pd.get_dummies(adata.obs[cluster_header]) #cell-by-cluster
    ## get number of clusters
    n_total_clusters = len(df_dummies.columns)
    print("--- %s seconds ---" % (time.time() - first_start_time))
    
    ##-----
    ## calculate medians ==> cluster_medians: rows = clusters, cols = genes
    ##-----
    if medians_header == None:
        print("Calculating medians...")
        start_time = time.time()
        # TODO: saved cluster_medians to save time
        if os.path.exists(output_folder + outputfilename + "_cluster_medians.csv"): 
            cluster_medians = pd.read_csv(output_folder + outputfilename + "_cluster_medians.csv", index_col = 0) 
        else: 
            ## get dataframes for X and cluster in a column
            df_X = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names) #cell-by-gene
            clusters = adata.obs[cluster_header]
            df_X_clusters = pd.concat([df_X, clusters], axis=1)
            ## get cluster medians
            cluster_medians = df_X_clusters.groupby([cluster_header], observed = True).median() #cluster-by-gene
            ## delete to free up memories
            del df_X, clusters, df_X_clusters
            cluster_medians.to_csv(output_folder + outputfilename + "_cluster_medians.csv")
            print("saving file as " + output_folder + outputfilename + "_cluster_medians.csv")
        print("--- %s seconds ---" % (time.time() - start_time))
    else:
        print("Getting pre-calculated medians...")
        start_time = time.time()
        cluster_medians = adata.varm[medians_header].transpose() #cluster-by-gene
        print("--- %s seconds ---" % (time.time() - start_time))
    
    ##-----
    ## calculate binary scores based on cluster_medians ==> binary_scores: rows = clusters, cols = genes
    ##-----
    if binary_scores_header == None:
        print("Calculating binary scores...")
        start_time = time.time()
        # TODO: saved binary_scores to save time
        if os.path.exists(output_folder + outputfilename + "_binary_scores.csv"): 
            binary_scores = pd.read_csv(output_folder + outputfilename + "_binary_scores.csv", index_col = 0) # index_col = "cluster"
        else: 
            ## calculate binary scores based on cluster_medians for each cluster
            binary_scores = []
            for cl in cluster_medians.index:
                ## get binary scores for all genes in each row (cluster) in df
                binary_scores_cl = [sum(np.maximum(0,1-cluster_medians[i]/cluster_medians.loc[cl,i]))/(n_total_clusters-1) for i in cluster_medians.columns]
                binary_scores.append(binary_scores_cl)
            ## binary scores matrix and handle nan
            binary_scores = pd.DataFrame(binary_scores, index=cluster_medians.index, columns=cluster_medians.columns).fillna(0) #cluster-by-gene
            binary_scores.to_csv(output_folder + outputfilename + "_binary_scores.csv")
            print("saving file as " + output_folder + outputfilename + "_binary_scores.csv")
        print("--- %s seconds ---" % (time.time() - start_time))
    else:
        print("Getting pre-calculated binary scores...")
        start_time = time.time()
        binary_scores = adata.varm[binary_scores_header].transpose() #cluster-by-gene
        print("--- %s seconds ---" % (time.time() - start_time))
    
    ####################### SET THRESHOLD/CUTOFF VALUE FOR BINARY FIRST TO USE #############################
    ####### "None" = do not use BinaryFirst to filter out any genes, "_mild" = median, "_moderate" = mean + 1 std, "_high" = mean + 2 std
        
    if gene_selection == "BinaryFirst_mild":
        print("Pre-select genes based on binary scores:")
        # if this line below is too slow, could also convert binary_scores df to np array and use np statistics
        median = binary_scores.stack().median()
        threshold = median
        print("\t", "Threshold (median):", threshold)
        
    elif gene_selection == "BinaryFirst_moderate":
        print("Pre-select genes based on binary scores:")
        # if this line below is too slow, could also convert binary_scores df to np array and use np statistics
        mean = binary_scores.stack().mean()
        stddev = binary_scores.stack().std()
        threshold = mean + stddev
        print("\t", "Threshold (mean + 1 * std):", threshold)
    
    elif gene_selection == "BinaryFirst_high":
        print("Pre-select genes based on binary scores:")
        # if this line below is too slow, could also convert binary_scores df to np array and use np statistics
        mean = binary_scores.stack().mean()
        stddev = binary_scores.stack().std()
        threshold = mean + 2 * stddev
        print("\t", "Threshold (mean + 2 * std):", threshold)
    
    else: 
    # if gene_selection == None:
        print("Pre-select genes based on binary scores: all genes (default)")
        threshold = 0
        print("\t", "Threshold:", threshold)
    
    ## create dummy matrix based on if binary score is greater than threshold
    binary_dummies = binary_scores >= threshold
    binary_dummies = binary_dummies*1
    binary_dummies = pd.DataFrame(binary_dummies, index=binary_scores.index, columns=binary_scores.columns) #cluster-by-gene
    ## average number of genes selected
    binary_dummies_sum = binary_dummies.stack().sum()
    print("\t Average number of genes after gene_selection in each cluster:", binary_dummies_sum/n_total_clusters)
    ## write number of genes selected for each cluster
    pd.DataFrame(binary_dummies.sum(axis=1)).to_csv(output_folder + outputfilename + "_gene_selection.csv") # , header=["n_genes_selected"]
    print("saving file as " + output_folder + outputfilename + "_gene_selection.csv")
    
    ############################## START iterations ######################################
    if cluster_list == None:
        cluster_list = df_dummies.columns
    n_clusters = len(cluster_list)
    
    print ("Number of clusters to evaluate: " + str(n_clusters))
    ct = 0
    df_supp = df_markers = df_results = pd.DataFrame()
    start_time = time.time()
    
    for cl in cluster_list:
        ct+=1
        print(f"{ct} out of {n_clusters}:")

        ## cluster in iteration
        print(f"\t{cl}")
        
        ##=== reset parameters for this iteration!!! (for taking care of special cases) ===##
        n_binary_genes_cl = n_binary_genes
        n_genes_eval_cl = n_genes_eval

        # Run RF
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
        markers, scores, score_max = mydecisiontreeevaluation.myDecisionTreeEvaluation(adata, df_dummies, cl, genes_eval, beta)
        print("\t" + str(markers))
        print("\t" + str(score_max))

        ## return supplementary table as csv
        binary_genes_list = top_binary_genes.index[:n_binary_genes_cl].to_list()
        df_supp_cl = pd.DataFrame({'clusterName': cl,
                                   'binary_genes': binary_genes_list,
                                   'rf_feature_importance': top_rf_genes[binary_genes_list],
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

    print("--- %s seconds ---" % (time.time() - first_start_time))
    ### END iterations ###
    
    return(df_results)
