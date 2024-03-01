
import pandas as pd
import scanpy as sc
import utils

# Calculate ratio of diagonal/total expression
def on_target_ratio(adata, nsf_results_df, cluster_header, medians_header, output_folder = "", outputfilename = ""): 

    """\
    Calculating the on-target ratios. 

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

    # cluster_header = "ann_finest_level"
    # medians_header = "medians_" + cluster_header
    # adata = sc.read_h5ad("demo_data/hlca_core_ann_finest_level_precalculated.h5ad")
    # nsf_results_df = pd.read_csv("results_hlca_core.csv")

    markers_dict = dict(zip(nsf_results_df["clusterName"], nsf_results_df["NSForest_markers"]))
    markers = []
    for key in markers_dict.keys(): 
        markers.extend(markers_dict[key])

    dend_header = "dendrogram_" + cluster_header
    dend_order = adata.uns[dend_header].get('categories_ordered')

    cluster_medians = adata.varm[medians_header].transpose()

    target_clusters, marker_target_exp, marker_total_exp, marker_fraction_values = [], [], [], []

    # loop over clusters in dend_order
    for key in markers_dict.keys():
        for value in markers_dict[key]:
            # append target cluster
            target_clusters.append(key)

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
    for cluster in dend_order:
        # get rows where target_cluster is cluster
        cluster_markers = marker_fraction_df[marker_fraction_df['target_cluster'] == cluster]
        # get median fraction value for this cluster
        median_fraction = cluster_markers['fraction'].median()
        median_fractions.append(median_fraction)

    marker_fraction_med_df = pd.DataFrame({'cluster': dend_order, 'median_fraction': median_fractions})
    filename = output_folder + outputfilename + "_median_fractions.csv"
    print("Saving median fractions as...", filename)
    marker_fraction_med_df.to_csv(filename, index = False)

    return

def calculate_diagonals(nsf_markers_df, cluster_medians, output_folder = "", outputfilename = ""): 
    
    """\
    Calculating the diagonals. 

    Parameters
    ----------
    nsf_results_df
        Output dataframe of NSForest. 
    medians_header
        Column in `adata`'s `.varm` storing median expression matrix. 
    output_folder
        Output folder. 
    outputfilename
        Prefix for all output files. 
    
    """

    marker_list = nsf_markers_df['markerGene']
    clusters = pd.unique(nsf_markers_df['clusterName'])

    df = pd.DataFrame()
    total_target_exp = 0
    total_subclade_exp = 0
    for m in range(len(marker_list)):
        # get median expression of marker in its target cluster
        marker = marker_list[m]
        target_cluster = nsf_markers_df.loc[m, 'clusterName']
        target_exp = cluster_medians.loc[target_cluster, marker]
        total_target_exp += target_exp

        # get median expression values of this marker m in all other clusters
        for cluster in clusters:
            if cluster == target_cluster: # don't include expression in target cluster! (so this ratio could be  > 1!)
                continue
            else:
                total_subclade_exp += cluster_medians.loc[cluster, marker]

        df_marker = pd.DataFrame({'marker': [marker],
                                'total_target_exp': [total_target_exp],
                                'total_subclade_exp': [total_subclade_exp],
                                'difference': [total_subclade_exp - total_target_exp],
                                'ratio': [total_target_exp/total_subclade_exp]})
        df = pd.concat([df, df_marker]).reset_index(drop=True)

    filename = output_folder + outputfilename + "_diagonal.csv"
    print("Saving diagonals as...", filename)
    df.to_csv(filename, index = False)

    return df
