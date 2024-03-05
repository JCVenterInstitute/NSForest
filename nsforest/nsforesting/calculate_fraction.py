
import numpy as np
import pandas as pd
import scanpy as sc

# Calculate ratio of diagonal/total expression
def on_target_fraction(adata, markers_dict, cluster_header, medians_header, output_folder, outputfilename): 
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
    outputfilename
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
    marker_fraction_df.to_csv(output_folder + outputfilename + "_marker_fractions.csv", index=False)

    median_fractions = []
    for cluster in np.unique(adata.obs[cluster_header]):
        # get rows where target_cluster is cluster
        cluster_markers = marker_fraction_df[marker_fraction_df['target_cluster'] == cluster]
        # get median fraction value for this cluster
        median_fraction = cluster_markers['fraction'].median()
        median_fractions.append(median_fraction)

    fractions_df = pd.DataFrame({'clusterName': np.unique(adata.obs[cluster_header]), 'fraction': median_fractions})
    return fractions_df
