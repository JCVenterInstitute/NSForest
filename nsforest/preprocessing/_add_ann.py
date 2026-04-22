
import time
import numpy as np
import pandas as pd
from tqdm import tqdm # may have to play with "import tqdm" vs "from tqdm import tqdm"
import matplotlib.pyplot as plt
import scanpy as sc
import logging

def dendrogram(adata, cluster_header, *, tl_kwargs = {}, pl_kwargs = {}, save = False, figsize = None, 
               output_folder = "", outputfilename_suffix = ""): 
    """\
    Generating a dendrogram from the AnnData object. 

    Parameters
    ----------
        adata: AnnData
            Annotated data matrix.
        cluster_header: str
            Column in `adata.obs` storing cell annotation. Passed into scanpy's dendrogram as `groupby`.
        tl_kwargs: dict
            Additional parameters to pass to sc.tl.dendrogram.
        pl_kwargs: dict
            Additional parameters to pass to sc.pl.dendrogram.
        save: bool | str (default: False)
            Whether to save plot in `output_folder`. If string, choose the type of file to save as ('png'(default), 'svg', 'pdf).
        figsize: tuple (default: (12, 2))
            figure.figsize for plt.rc_context. 
        output_folder: str (default: "")
            Output folder. Created if doesn't exist. 
        outputfilename_suffix: str (default: "")
            Suffix for all output files. 
    
    Returns
    -------
    does not return anything. Adds `adata.uns["dendrogram_{cluster_header}"]` to passed in adata. 
    """
    if save: 
        if save == True: 
            save = "png"
        if save not in ['pdf', 'png', 'svg']: 
            print("warning: `save` must be one of the following: 'pdf', 'png', 'svg'")
            print("saving as png")
            save = "png"
        sc.settings.verbosity = 0
        sc.settings.figdir = output_folder
        save = f"_{outputfilename_suffix}.{save}"
        print(f"Saving dendrogram as...\n{output_folder}{save}")
    if not adata.obsm and "X_pca" not in adata.obsm: 
        sc.pp.pca(adata)
    sc.tl.dendrogram(adata, cluster_header, use_rep="X_pca", **tl_kwargs)
    if not figsize: 
        n_clusters = len(adata.obs[cluster_header].unique())
        fig_width = max(6.4, n_clusters/5)
        fig_height = max([len(z) for z in adata.obs[cluster_header].unique()]) / 30 + 1
        figsize = (fig_width, fig_height)
    with plt.rc_context({"figure.figsize": figsize}): 
        # sc.pl.dendrogram(adata, cluster_header, save = save, **pl_kwargs)
        if save: sc.pl.dendrogram(adata, cluster_header, **pl_kwargs).figure.savefig(save)
        else: sc.pl.dendrogram(adata, cluster_header, **pl_kwargs)
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
    if use_mean: value = "mean"
    else: value = "median"
    for cl in tqdm(sorted(set(adata.obs[cluster_header])), desc=f"Calculating {value}s per cluster"):
        adata_cl = adata[adata.obs[cluster_header]==cl,]
        if use_mean: 
            medians_cl = adata_cl.to_df().mean()
        else: 
            medians_cl = adata_cl.to_df().median()
        cluster_medians = pd.concat([cluster_medians, pd.DataFrame({cl: medians_cl})], axis=1) #gene-by-cluster
    return cluster_medians

def prep_medians(adata, cluster_header, use_mean = False, positive_genes_only = True, plot = False):
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
    start_time = time.time()

    ## get medians
    if use_mean: print("`use_mean` is True. Using the mean expression of each gene per cluster.")
    cluster_medians = get_medians(adata, cluster_header, use_mean) #gene-by-cluster
    adata.varm['medians_' + cluster_header] = cluster_medians
    print("Saving medians as adata.varm.medians_" + cluster_header)

    print("median:", round(cluster_medians.stack().median(), 3))
    print("mean:", round(cluster_medians.stack().mean(), 3))
    print("std:", round(cluster_medians.stack().std(), 3))

    if positive_genes_only:
        ## select only genes with median > 0
        genes_selected = cluster_medians.index[cluster_medians.sum(axis=1)>0].to_list()
        print(f"Only positive genes selected. {len(genes_selected)} positive genes out of {adata.n_vars} total genes")
        ## subset data with only positive genes
        adata = adata[:,genes_selected].copy()
    
    print("--- %s seconds ---" % (time.time() - start_time))

    

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
    start_time = time.time()
    
    ## get medians
    if medians_header == "medians_": medians_header = "medians_" + cluster_header
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
    adata.varm['binary_scores_' + cluster_header] = binary_scores.transpose() 
    print("Saving binary scores as adata.varm.binary_scores_" + cluster_header)
    print("median:", round(binary_scores.stack().median(), 3))
    print("mean:", round(binary_scores.stack().mean(), 3))
    print("std:", round(binary_scores.stack().std(), 3))
    print("--- %s seconds ---" % (time.time() - start_time))

    return adata

def plot_varm(adata, varm_key, nonzero = False, scale = None, figsize = (6, 4), show = True, save = False, output_folder = ""): 
    """\
    Plotting histogram of median expression per gene per cluster.  

    Parameters:
    -----------
        adata: AnnData
            Annotated data matrix.
        varm_key: str
            Key in `adata.varm` storing calculated medians or binary scores.
        nonzero: bool
            Whether to remove zeros from histogram.
        scale: str
            How to scale the y-axis. 
        figsize: tuple
            Width and height of plot. 
        show: bool
            Whether to show the plot. 
        save: bool | str (default: False)
            Whether to save plot. If string, choose the type of file to save as ("png", "svg", "pdf").
        output_folder: str (default: "")
            Output folder for output files. 
    
    Returns:
    ========
    fig: matplotlib.pyplot.figure
        Histogram of adata.varm[varm_key]
    """
    plt.figure(figsize = figsize)

    # If removing zero values
    if nonzero: 
        values = adata.varm[varm_key].unstack().replace(0, np.nan).dropna(how='all')
    else: 
        values = adata.varm[varm_key].unstack()

    plt.hist(values, bins = 100)

    # If y-axis is log scaled
    if scale == "log": plt.yscale("log")
    plt.xlabel(varm_key)

    # Adding title
    if "medians" in varm_key: title = "Median expression"
    elif "binary_score" in varm_key: title = "Binary score"
    else: title = f"adata.varm[{varm_key}]"
    if nonzero: title = title + " (non-zeros)"
    plt.title(title)
    
    if save: 
        if save == True: 
            save = "png"
        elif save not in ['pdf', 'png', 'svg']: 
            print("warning: `save` must be one of the following: 'pdf', 'png', 'svg'")
            print("saving as png")
            save = "png"
        filename = f'{output_folder}histogram_{varm_key}.{save}'
        print(f"Saving adata.varm[{varm_key}] as histogram as...\n{filename}")
        plt.savefig(filename, bbox_inches='tight')
    if show: plt.show()
    return

def spaceTx_genefilter(adata, lower_percentile = 0.1, upper_percentile = 0.99, min_txLength = 700, species = "human", species_dict = None, gencode_folder = "gencode_annotation"): 

    """\
    Filtering genes for spatial gene probe panel design. 

    Parameters
    ----------
        adata: AnnData
            Annotated data matrix.
        lower_percentile: float (default: 0.1)
            Lower quartile percentile to filter non-0 median gene expression. 
        upper_percentile: float (default: 0.99)
            Upper quartile percentile to filter non-0 median gene expression. 
        min_txLength: int (default: 700)
            Minimum transcript length. 
        species: ["human", "mouse", "other"] (default: "human")
            Species relating to gencode_annotation. 
    
    Returns
    -------
    adata: AnnData
        Subset AnnData based on `lower_percentile`, `upper_percentile`, `min_txLength`.
    """
    ### FILTER 1: EXPRESSION ###
    # Filtering non-zero median expression of all genes
    
    expr = adata.to_df()
    ## get non-zero median expression values
    expr_nonZeroMedian = expr.apply(lambda xx: xx.loc[xx!=0].median())
    
    ## check
    expr_ZeroMedian = expr.apply(lambda xx: xx.loc[xx==0].median())
    if expr_ZeroMedian.sum() != 0: 
        print("warning: expr_ZeroMedian.sum() != 0")#should be 0
    
    limits = expr_nonZeroMedian.quantile([lower_percentile, upper_percentile])
    print(f"Non-zero median expression percentile limits: \n{limits}")
    
    ## plot histogram
    ax = expr_nonZeroMedian.hist(bins=50, grid=False)
    # Add a vertical line at specified percentiles
    ax.axvline(limits.iloc[0], color='red', linestyle='dashed', label=f'{int(lower_percentile*100)}% percentile')
    ax.axvline(limits.iloc[1], color='red', linestyle='dotted', label=f'{int(upper_percentile*100)}% percentile')
    ax.legend()
    plt.title("Non-zero median expression value")
    plt.show()

    ## select genes that pass the expression filter
    ind_selected_expr = (expr_nonZeroMedian > limits.iloc[0]) & (expr_nonZeroMedian < limits.iloc[1])
    print(f'FILTER 1: {sum(ind_selected_expr)} out of {adata.n_vars} total genes passed the expression filter (lower_percentile = {lower_percentile}, upper_percentile = {upper_percentile}).\n')
    
    ### FILTER 2: TRANSCRIPT LENGTH ###
    # Filtering by transcript length

    # Chosing gencode csv
    if not species_dict: 
        species_dict = {"human": f"{gencode_folder}/gencode.v47.annotation_txLength.csv", 
                        "mouse": f"{gencode_folder}/gencode.vM36.annotation_txLength.csv"}
    if species=="other": 
        print('FILTER 2: Filter based on transcript length is omitted.')
        ind_selected_final = ind_selected_expr
    elif species not in species_dict: 
        print("ERROR: attempting to filter transcript length. Add species in species_dict and gencode_annotation folder")
    else: 
        annotation_txLength = pd.read_csv(species_dict[species], low_memory = False)
        # column check
        if "ENSEMBL_ID" not in list(annotation_txLength.columns) or "tx_length" not in list(annotation_txLength.columns): 
            print("ERROR: some required column names are missing: ENSEMBL_ID, tx_length")
            return
        tx_length = adata.var.copy().reset_index()
        if tx_length[tx_length.duplicated("index", keep = False)].shape[0] != 0: 
            print("ERROR: duplicate adata.var_names detected")
            return
        if "ENS" in list(tx_length["index"])[0]: col = "ENSEMBL_ID"
        else: col = "gene_name"
        tx_length = tx_length.rename(columns = {"index": col})
        tx_length = pd.merge(tx_length, annotation_txLength[[col, "tx_length"]], on = col, how = "left").sort_values("tx_length")
        tx_length = tx_length.drop_duplicates(col, keep = "first")
        tx_length.index = tx_length[col]
        ind_selected_txLength = tx_length["tx_length"] > min_txLength
        ## select genes that pass the txLength filter
        print(f'FILTER 2: {sum(ind_selected_txLength)} out of {adata.n_vars} total genes passed the transcript length filter (min_txLength = {min_txLength}).')
        ind_selected_final = ind_selected_expr & ind_selected_txLength

    ## subset adata for final selected genes
    print(f'\nFINAL SELECTION: {sum(ind_selected_final)} out of {adata.n_vars} total genes passed both filters.')
    adata_prep = adata[:,ind_selected_final].copy()
    
    return adata_prep
