
import numpy as np
import pandas as pd
import itertools
import statistics
import scanpy as sc
import nsforest as ns
import nsforest.preprocessing as pp

# Calculate ratio of diagonal/total expression
def markers_onTarget(adata, markers_dict, cluster_header, use_mean = False, output_folder = "", outputfilename_prefix = ""):
    all_markers = list(set(itertools.chain.from_iterable(markers_dict.values())))
    adata_eval = adata[:,all_markers]
    # cluster_medians = adata.varm[medians_header].T
    if use_mean: print("using mean")
    else: print("using median")
    cluster_medians = ns.pp.get_medians(adata_eval, cluster_header, use_mean = use_mean) #gene-by-cluster
    cluster_medians_genesum = cluster_medians.sum(axis=1) #i.e. rowsum 

    df_ontarget_supp = pd.DataFrame()
    for cl in markers_dict.keys():
        ontarget_per_gene = []
        markers = markers_dict[cl]
        for gg in markers:
            if cluster_medians_genesum[gg] != 0:
                ontarget_gg = cluster_medians.loc[gg,cl] / cluster_medians.loc[gg,].sum()
            else: 
                ontarget_gg = 0
            ontarget_per_gene.append(ontarget_gg)
        if use_mean:
            ontarget = statistics.mean(ontarget_per_gene)
        else: 
            ontarget = statistics.median(ontarget_per_gene) #slighly in favors good markers

        ## return ontarget table as csv
        df_ontarget_cl = pd.DataFrame({'clusterName': cl, 'markerGene': markers, 
                                       'onTarget_per_gene': ontarget_per_gene, 'onTarget': ontarget})
        df_ontarget_supp = pd.concat([df_ontarget_supp, df_ontarget_cl]).reset_index(drop=True) 
        df_ontarget_supp.to_csv(output_folder + outputfilename_prefix + "_markers_onTarget_supp.csv", index=False)
    
    df_ontarget = df_ontarget_supp[['clusterName', 'onTarget']].drop_duplicates().reset_index(drop=True)
    df_ontarget.to_csv(output_folder + outputfilename_prefix + "_markers_onTarget.csv", index=False)
    return df_ontarget

# Calculate ratio of diagonal/total expression
def on_target_fraction_angela(adata, markers_dict, cluster_header, medians_header, output_folder, outputfilename_prefix): 
    """\
    Calculating the on-target fraction. 

    Parameters
    ----------
    adata
        Annotated data matrix.
    nsf_results_df
        Output dataframe of NSForest. 
    cluster_header
        Column in `adata`'s `.obs` representing cell annotation.
    medians_header
        Column in `adata`'s `.varm` storing median expression matrix. 
    output_folder
        Output folder. 
    outputfilename_prefix
        Prefix for all output files. 
    """
    cluster_medians = adata.varm[medians_header].transpose()
    
    target_clusters, markers, marker_target_exp, marker_total_exp, marker_fraction_values = [], [], [], [], []
    for key in markers_dict.keys():
        for value in markers_dict[key]:
            # append target cluster
            target_clusters.append(key)
            markers.append(value)
            # get marker expression in target cluster
            target_exp = cluster_medians.loc[key, value]
            marker_target_exp.append(target_exp)
            # get total expression for maker in all clusters
            total_exp = cluster_medians.loc[:, value].sum()
            marker_total_exp.append(total_exp)
            # get on-target fraction for this marker
            on_target_fraction = target_exp / total_exp
            marker_fraction_values.append(on_target_fraction)

    # make new df with target cluster, marker, target exp, total exp, and fraction
    marker_fraction_df = pd.DataFrame({'target_cluster': target_clusters, 'markerGene': markers, 'target_exp': marker_target_exp, 'total_exp': marker_total_exp, 'fraction': marker_fraction_values})
    marker_fraction_df.to_csv(output_folder + outputfilename_prefix + "_marker_fractions.csv", index=False)

    median_fractions = []
    for cluster in np.unique(adata.obs[cluster_header]):
        # get rows where target_cluster is cluster
        cluster_markers = marker_fraction_df[marker_fraction_df['target_cluster'] == cluster]
        # get median fraction value for this cluster
        median_fraction = cluster_markers['fraction'].median()
        median_fractions.append(median_fraction)

    fractions_df = pd.DataFrame({'clusterName': np.unique(adata.obs[cluster_header]), 'fraction': median_fractions})
    return fractions_df
