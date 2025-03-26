
import pandas as pd
import scanpy as sc
import plotly.express as px

def boxplot(df, col, save = False, output_folder = "", outputfilename_prefix = ""): 
    """\
    Generating plotly boxplot. The `hover_name` is "clusterName". 

    Parameters:
    -----------
        df: pd.DataFrame
            NS-Forest results containing "clusterName" and `col` columns. 
        col: str
            Column in `df` to create the boxplot from. 
        save: bool | str (default: False)
            Whether to save plot. If string, choose the type of file to save as ("png", "jpeg", "webp", "svg", "pdf").
        output_folder: str (default: "")
            Output folder for output files. 
        outputfilename_prefix: str (default: "")
            Prefix for all output files. 
    
    Returns:
    ========
    fig: plotly.graph_objects.Figure
        Boxplot of `col` values. 
    """
    fig = px.box(df, y=col, points='all', range_y=[-.05,1.05],
                 title=f"{col} median = {round(df[col].median(),3)}",
                 width=400, height=500, hover_name='clusterName')
    if save: 
        if save in [True, "html"]: 
            filename = output_folder + outputfilename_prefix + f"_boxplot_{col}.html"
            fig.write_html(filename)
        elif save in ["png", "jpeg", "webp", "svg", "pdf"]: 
            filename = output_folder + outputfilename_prefix + f"_boxplot_{col}.{save}"
            fig.write_image(filename)
        else: 
            print(f"ERROR: invalid file extension: {save}")
        print("Saving...\n", filename)
        
    return fig

def scatter_w_clusterSize(df, col, save = False, output_folder = "", outputfilename_prefix = ""): 
    """\
    Generating plotly scatterplot with "clusterSize". The `hover_name` is "clusterName". 

    Parameters:
    -----------
        df: pd.DataFrame
            NS-Forest results containing "clusterName", "clusterSize", and `col` columns. 
        col: str
            Column in `df` to plot against "clusterSize". 
        save: bool | str (default: False)
            Whether to save plot. If string, choose the type of file to save as ("png", "jpeg", "webp", "svg", "pdf").
        output_folder: str (default: "")
            Output folder for output files. 
        outputfilename_prefix: str (default: "")
            Prefix for all output files. 
    
    Returns:
    ========
    fig: plotly.graph_objects.Figure
        Scatterplot of `col` values against `clusterSize`. 
    """
    fig = px.scatter(df, x='clusterSize', y=col, range_y=[-.05,1.05],
                     width=700, height=500, hover_name='clusterName')
    if save in [True, "html"]: 
        filename = output_folder + outputfilename_prefix + f"_scatter_{col}.html"
        print("Saving...\n", filename)
        fig.write_html(filename)
    elif save in ["png", "jpeg", "webp", "svg", "pdf"]: 
        filename = output_folder + outputfilename_prefix + f"_scatter_{col}.{save}"
        fig.write_image(filename)
    else: 
        print(f"ERROR: invalid file extension: {save}")
    return fig 

def dotplot(adata, markers, cluster_header, *, dendrogram = True, save = False, 
            output_folder = "", outputfilename_suffix = "", **kwargs): 
    """\
    Generating scanpy dotplot of `adata` with input marker list. 

    Parameters:
    -----------
        adata: AnnData
            Annotated data matrix.
        markers: list | dict
            List of markers to show in dotplot. Dictionary of markers per cell type to group by. 
        cluster_header: str
            Column in `adata.obs` storing cell annotation.
        dendrogram: bool | list (default: True)
            Whether to use dendrogram from `adata.uns["dendrogram_{cluster_header}"]`. If list, dendrogram order. 
        save: bool | str (default: False)
            Whether to save plot. If string, choose the type of file to save as ("png", "jpeg", "webp", "svg", "pdf").
        output_folder: str (default: ".")
            Output folder. Created if doesn't exist. 
        outputfilename_suffix: str (default: "")
            Prefix for all output files. 
        kwargs: dictionary (default: None)
            Additional parameters to pass to sc.pl.dotplot.
    """
    if save: 
        sc.settings.figdir = output_folder
        if save == True: 
            save = "png"
        if save in ["png", "jpeg", "webp", "svg", "pdf"]: 
            save = outputfilename_suffix + f".{save}"
        else: 
            print(f"ERROR: invalid file extension: {save}")
            save = False
    if isinstance(dendrogram, bool): # gene_symbols = gene_symbols, use_raw = False, standard_scale = "var", 
        sc.pl.dotplot(adata, markers, cluster_header, dendrogram = dendrogram, save = save, **kwargs)
    elif isinstance(dendrogram, list): 
        sc.pl.dotplot(adata, markers, cluster_header, categories_order = dendrogram, save = save, **kwargs)
    return 

def stackedviolin(adata, markers, cluster_header, *, dendrogram = True, save = False, 
                  output_folder = "", outputfilename_suffix = "", **kwargs): 
    """\
    Generating scanpy stacked_violin of `adata` with input marker list. 

    Parameters:
    -----------
        adata: AnnData
            Annotated data matrix.
        markers: list/dict
            List of markers to show in dotplot. Dictionary of markers per cell type to group by. 
        cluster_header: str
            Key in `adata.obs` storing cell type.
        dendrogram: bool/list (default: True)
            Whether to use dendrogram from `adata.uns["dendrogram_{cluster_header}"]`. Dendrogram order. 
        save: bool, str (default: False)
            Whether to save plot. Set as option for type of image to save ("html", "svg", "png", etc)
        output_folder: str (default: ".")
            Output folder. Created if doesn't exist. 
        outputfilename_suffix: str (default: "")
            Prefix for all output files. 
        kwargs: dictionary (default: None)
            Additional parameters to pass to sc.pl.stacked_violin.
    """
    if save: 
        sc.settings.figdir = output_folder
        if save in [True, "png"]: 
            save = outputfilename_suffix + ".png"
        elif save in ["png", "jpeg", "webp", "svg", "pdf"]: 
            save = outputfilename_suffix + f".{save}"
        else: 
            print(f"ERROR: invalid file extension: {save}")
            save = False
    if isinstance(dendrogram, bool): # gene_symbols = gene_symbols, use_raw = False, standard_scale = "var", 
        sc.pl.stacked_violin(adata, markers, cluster_header, dendrogram = dendrogram, save = save, **kwargs)
    elif isinstance(dendrogram, list): 
        sc.pl.stacked_violin(adata, markers, cluster_header, categories_order = dendrogram, save = save, **kwargs)
    return

def matrixplot(adata, markers, cluster_header, *, dendrogram = True, save = False, 
               output_folder = "", outputfilename_suffix = "", **kwargs): 
    """\
    Generating scanpy matrixplot of `adata` from `markers`. 

    Parameters:
    -----------
        adata: AnnData
            Annotated data matrix.
        markers: list/dict
            List of markers to show in dotplot. Dictionary of markers per cell type to group by. 
        cluster_header: str
            Key in `adata.obs` storing cell type.
        dendrogram: bool/list (default: True)
            Whether to use dendrogram from `adata.uns["dendrogram_{cluster_header}"]`. Dendrogram order. 
        save: bool, str (default: False)
            Whether to save plot. Set as option for type of image to save ("html", "svg", "png", etc)
        output_folder: str (default: ".")
            Output folder. Created if doesn't exist. 
        outputfilename_suffix: str (default: "")
            Prefix for all output files. 
        kwargs: dictionary (default: None)
            Additional parameters to pass to sc.pl.matrixplot.
    """
    if save: 
        sc.settings.figdir = output_folder
        if save in [True, "png"]: 
            save = outputfilename_suffix + ".png"
        elif save in ["png", "jpeg", "webp", "svg", "pdf"]: 
            save = outputfilename_suffix + f".{save}"
        else: 
            print(f"ERROR: invalid file extension: {save}")
            save = False
    if isinstance(dendrogram, bool): # heatmap # gene_symbols = gene_symbols, use_raw = False, standard_scale = "var", 
        sc.pl.matrixplot(adata, markers, cluster_header, dendrogram = dendrogram, save = save, **kwargs)
    elif isinstance(dendrogram, list): 
        sc.pl.matrixplot(adata, markers, cluster_header, categories_order = dendrogram, save = save, **kwargs)
    return
