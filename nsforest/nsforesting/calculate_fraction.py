
import numpy as np
import pandas as pd
import itertools
import statistics
import scanpy as sc
import nsforest as ns
import nsforest.preprocessing as pp

def markers_onTarget_percluster(adata, cluster_header, markers_dict, use_mean = False, save_supplementary = False, output_folder = "", outputfilename_prefix = ""):
    """\
    Calculating the on-target fraction of each gene list. 

    Parameters
    ----------
        adata: AnnData
            Annotated data matrix.
        cluster_header: str
            Column in `adata.obs` storing cell annotation.
        markers_dict: dict
            Dictionary containing genes for each `cluster_header` (clusterName: list of markers)
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
    df_ontarget: pd.DataFrame
        AnnData with on-target fraction per `cluster_header`. 
    """
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
            ontarget = statistics.median(ontarget_per_gene) #slightly favors good markers

        ## return ontarget table as csv
        df_ontarget_cl = pd.DataFrame({'clusterName': cl, 'markerGene': markers, 
                                       'onTarget_per_gene': ontarget_per_gene, 'onTarget': ontarget})
        df_ontarget_supp = pd.concat([df_ontarget_supp, df_ontarget_cl]).reset_index(drop=True) 
        if save_supplementary: 
            df_ontarget_supp.to_csv(output_folder + outputfilename_prefix + "_markers_onTarget_supp.csv", index=False)
    
    df_ontarget = df_ontarget_supp[['clusterName', 'onTarget']].drop_duplicates().reset_index(drop=True)
    if save_supplementary: 
        df_ontarget.to_csv(output_folder + outputfilename_prefix + "_markers_onTarget.csv", index=False)

    return df_ontarget

def markers_onTarget(adata, cluster_header, markers_dict, use_mean = False, save_supplementary = False, output_folder = "", outputfilename_prefix = ""):
    """\
    Calculating the on-target fraction of each gene list. 

    Parameters
    ----------
        adata: AnnData
            Annotated data matrix.
        cluster_header: str
            Column in `adata.obs` storing cell annotation.
        markers_dict: dict
            Dictionary containing genes for each `cluster_header` (clusterName: list of markers)
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
    df_ontarget: pd.DataFrame
        AnnData with on-target fraction per `cluster_header`. 
    """
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
            ontarget = statistics.median(ontarget_per_gene) #slightly favors good markers

        ## return ontarget table as csv
        df_ontarget_cl = pd.DataFrame({'clusterName': cl, 'markerGene': markers, 
                                       'onTarget_per_gene': ontarget_per_gene, 'onTarget': ontarget})
        df_ontarget_supp = pd.concat([df_ontarget_supp, df_ontarget_cl]).reset_index(drop=True) 
        if save_supplementary: 
            df_ontarget_supp.to_csv(output_folder + outputfilename_prefix + "_markers_onTarget_supp.csv", index=False)
    
    df_ontarget = df_ontarget_supp[['clusterName', 'onTarget']].drop_duplicates().reset_index(drop=True)
    if save_supplementary: 
        print(f"Saving supplementary table as...\n{output_folder}{outputfilename_prefix}_markers_onTarget_supp.csv")
        df_ontarget_supp.to_csv(output_folder + outputfilename_prefix + "_markers_onTarget_supp.csv", index=False)
        print(f"Saving supplementary table as...\n{output_folder}{outputfilename_prefix}_markers_onTarget.csv")
        df_ontarget.to_csv(output_folder + outputfilename_prefix + "_markers_onTarget.csv", index=False)

    return df_ontarget
