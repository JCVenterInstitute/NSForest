
import time
import numpy as np
import pandas as pd
from tqdm import tqdm # may have to play with "import tqdm" vs "from tqdm import tqdm"
import matplotlib.pyplot as plt
import scanpy as sc

def dendrogram(adata, cluster_header, save = ".png"): 
    """\
    Generating a scanpy dendrogram from anndata. 

    Parameters:
    -----------
        adata: AnnData
            Annotated data matrix. 
        cluster_header: str
            Key in `adata.obs` storing cell type. 
        save: str (default: "dendrogram.png")
            Dendrogram file name suffix.  
    
    Returns:
    ========
        adata: AnnData
            AnnData with `cluster_header`_dendrogram in `adata.uns`. 
    """
    with plt.rc_context({"figure.figsize": (12, 1)}):
        sc.pl.dendrogram(adata, cluster_header, orientation = "top", save = save)
    return adata

def get_medians(adata, cluster_header, use_mean = False): 
    """\
    Calculating the median expression per gene for each cluster. 

    Parameters:
    -----------
        adata: AnnData
            Annotated data matrix. 
        cluster_header: str
            Key in `adata.obs` storing cell type. 
        use_mean: bool (default: False)
            Whether to use the mean or median for minimum gene expression threshold. 
    
    Returns:
    ========
        cluster_medians: pd.DataFrame
            gene-by-cluster median (mean) expression matrix. 
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
    Calling get_medians(). Filtering for positive genes if specified. 

    Parameters:
    -----------
        adata: AnnData
            Annotated data matrix. 
        cluster_header: str
            Key in `adata.obs` storing cell type. 
        use_mean: bool (default: False)
            Whether to use the mean or median for minimum gene expression threshold. 
        positive_genes_only: bool (default: True)
            Whether to subset `adata` to only have genes with median (mean) expression greater than 0. 
    
    Returns:
    ========
        adata: AnnData
            AnnData with `cluster_header`_medians in `adata.varm`. 
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

    if positive_genes_only:
        ## select only genes with median > 0
        genes_selected = cluster_medians.index[cluster_medians.sum(axis=1)>0].to_list()
        print(f"Only positive genes selected. {len(genes_selected)} positive genes out of {adata.n_vars} total genes")
        ## subset data with only positive genes
        adata = adata[:,genes_selected].copy()
    return adata

def prep_binary_scores(adata, cluster_header, medians_header):
    """\
    Calculating the binary scores for each gene per cluster. 

    Parameters:
    -----------
        adata: AnnData
            Annotated data matrix. 
        cluster_header: str
            Key in `adata.obs` storing cell type. 
        medians_header: str
            Key in `adata.varm` storing median expression matrix. 
    
    Returns:
    ========
        adata: AnnData
            AnnData with `cluster_header`_binary_scores in `adata.varm`. 
    """
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
