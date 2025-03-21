
import time
import numpy as np
import pandas as pd
from tqdm import tqdm # may have to play with "import tqdm" vs "from tqdm import tqdm"
import matplotlib.pyplot as plt
import scanpy as sc
import logging

def dendrogram(adata, cluster_header, *, plot = False, save = False, figsize = (12, 2), 
               output_folder = "", outputfilename_suffix = "", **kwargs): 
    """\
    Generating a dendrogram from the AnnData object. 

    Parameters
    ----------
        adata: AnnData
            Annotated data matrix.
        cluster_header: str
            Column in `adata.obs` storing cell annotation. Passed into scanpy's dendrogram as `groupby`.
        plot: bool (default: False)
            Whether to use sc.pl.dendrogram instead of sc.tl.dendrogram. 
        save: bool | str (default: False)
            Whether to save plot in `output_folder`. If string, choose the type of file to save as ('png'(default), 'svg', 'pdf).
        figsize: tuple (default: (12, 2))
            figure.figsize for plt.rc_context. 
        output_folder: str (default: "")
            Output folder. Created if doesn't exist. 
        outputfilename_suffix: str (default: "")
            Suffix for all output files. 
        kwargs: dictionary (default: None)
            Additional parameters to pass to sc.tl.dendrogram or sc.pl.dendrogram.
    
    Returns
    -------
    does not return anything. Adds `adata.uns["dendrogram_{cluster_header}"]` to passed in adata. 
    """
    if not plot: # default no plot, no save
        sc.tl.dendrogram(adata, cluster_header, **kwargs)
    if save: 
        if save == True: 
            save = "png"
        if save not in ['pdf', 'png', 'svg']: 
            print("warning: `save` must be one of the following: 'pdf', 'png', 'svg'")
            print("saving as png")
            save = "png"
        sc.settings.figdir = output_folder
        save = f"_{outputfilename_suffix}.{save}"
    with plt.rc_context({"figure.figsize": figsize}): 
        sc.pl.dendrogram(adata, cluster_header, save = save, **kwargs)
    return

def get_medians(adata, cluster_header, use_mean = False): 
    """\
    Calculating the median (mean) expression per gene for each `cluster_header`. 

    Parameters
    ----------
        adata: AnnData
            Annotated data matrix.
        cluster_header: str
            Column in `adata.obs` storing cell annotation.
        use_mean: bool (default: False)
            Whether to use the mean (vs median) for minimum gene expression threshold. 
    
    Returns
    -------
    cluster_medians: pd.DataFrame
        Gene-by-cluster median (mean) expression dataframe. 
    """
    cluster_medians = pd.DataFrame()
    for cl in tqdm(sorted(set(adata.obs[cluster_header])), desc="Calculating medians (means) per cluster"):
        adata_cl = adata[adata.obs[cluster_header]==cl,]
        if use_mean: 
            medians_cl = adata_cl.to_df().mean()
        else: 
            medians_cl = adata_cl.to_df().median()
        cluster_medians = pd.concat([cluster_medians, pd.DataFrame({cl: medians_cl})], axis=1) #gene-by-cluster
    return cluster_medians

def prep_medians(adata, cluster_header, use_mean = False, positive_genes_only = True):
    """\
    Calculating the median expression matrix. Subsetting `adata` if `positive_genes_only` = True. 

    Parameters
    ----------
        adata: AnnData
            Annotated data matrix.
        cluster_header: str
            Column in `adata.obs` storing cell annotation.
        use_mean: bool (default: False)
            Whether to use the mean (vs median) for minimum gene expression threshold. 
        positive_genes_only: bool (default: True)
            Whether to subset AnnData to only have genes with median/mean expression greater than 0. 
    
    Returns
    -------
    adata: AnnData
        AnnData with median expression values stored in `adata.varm["medians_{cluster_header}"]`. 
    """
    print("Calculating medians...")
    start_time = time.time()
    ## get medians
    if use_mean: 
        print("use_mean is True. Using the mean expression of each gene per cluster. ")
    cluster_medians = get_medians(adata, cluster_header, use_mean) #gene-by-cluster
    ## attach calculated medians to adata
    print("Saving calculated medians as adata.varm.medians_" + cluster_header)
    adata.varm['medians_' + cluster_header] = cluster_medians #gene-by-cluster
    print("--- %s seconds ---" % (time.time() - start_time))
    print("median:", cluster_medians.stack().median())
    print("mean:", cluster_medians.stack().mean())
    print("std:", cluster_medians.stack().std())

    if positive_genes_only:
        ## select only genes with median > 0
        genes_selected = cluster_medians.index[cluster_medians.sum(axis=1)>0].to_list()
        print(f"Only positive genes selected. {len(genes_selected)} positive genes out of {adata.n_vars} total genes")
        ## subset data with only positive genes
        adata = adata[:,genes_selected].copy()
    return adata

def prep_binary_scores(adata, cluster_header, medians_header = "medians_"):
    """\
    Calculating the binary scores of each gene per `cluster_header`. 

    Parameters
    ----------
        adata: AnnData
            Annotated data matrix.
        cluster_header: str
            Column in `adata.obs` storing cell annotation.
        medians_header: str (default: "medians_{cluster_header}")
            Key in `adata.varm` storing median expression matrix. 
    
    Returns
    -------
    adata: AnnData
        AnnData with binary scores stored in `adata.varm["binary_scores_{cluster_header}"]`. 
    """
    # default medians_header
    if medians_header == "medians_": medians_header = "medians_" + cluster_header

    print("Calculating binary scores...")
    start_time = time.time()
    ## get medians
    cluster_medians = adata.varm[medians_header].transpose() #cluster-by-gene
    n_total_clusters = cluster_medians.shape[0]
    
    ## calculate binary scores based on cluster_medians for all genes per cluster
    binary_scores = []
    for cl in tqdm(cluster_medians.index, desc="Calculating binary scores per cluster"):
        ## get binary scores for all genes in each row (cluster) in df
        binary_scores_cl = [sum(np.maximum(0,1-cluster_medians[i]/cluster_medians.loc[cl,i]))/(n_total_clusters-1) for i in cluster_medians.columns]
        binary_scores.append(binary_scores_cl)

    ## binary scores matrix and handle nan
    binary_scores = pd.DataFrame(binary_scores, index=cluster_medians.index, columns=cluster_medians.columns).fillna(0) #cluster-by-gene
    ## attach pre-calculated binary scores to adata
    print("Saving calculated binary scores as adata.varm.binary_scores_" + cluster_header)
    adata.varm['binary_scores_' + cluster_header] = binary_scores.transpose() #gene-by-cluster
    print("--- %s seconds ---" % (time.time() - start_time))
    print("median:", binary_scores.stack().median())
    print("mean:", binary_scores.stack().mean())
    print("std:", binary_scores.stack().std())
    
    return adata
