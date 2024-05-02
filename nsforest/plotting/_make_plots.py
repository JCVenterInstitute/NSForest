
import pandas as pd
import scanpy as sc
import plotly.express as px

def boxplot(df, col, save = False, output_folder = "", outputfilename_prefix = ""): 
    """\
    Generating plotly boxplot. 

    Parameters:
    -----------
        df: pd.DataFrame
            NS-Forest results containing `clusterName` and `col`. 
        col: str
            Column in `df` to create the boxplot. 
        save: bool (default: False)
            Whether to save html file. 
        output_folder: str (default: "")
            Output folder for output files. 
        outputfilename_prefix: str (default: "")
            Prefix for all output files. 
    
    Returns:
    ========
        fig: plotly.graph_objects.Figure
            Boxplot of `col` values
    """
    fig = px.box(df, y=col, points='all', range_y=[-.05,1.05],
                 title=f"{col} median = {round(df['f_score'].median(),3)}",
                 width=400, height=500, hover_name='clusterName')
    if save: 
        filename = output_folder + outputfilename_prefix + f"_boxplot_{col}.html"
        print("Saving...\n", filename)
        fig.write_html(filename)
    return fig

def scatter_w_clusterSize(df, col, save = False, output_folder = "", outputfilename_prefix = ""): 
    """\
    Generating plotly scatterplot with `clusterSize`. 

    Parameters:
    -----------
        df: pd.DataFrame
            NS-Forest results containing `clusterName`, `clusterSize`, and `col`. 
        col: str
            Column in `df` to plot against `clusterSize`. 
        save: bool (default: False)
            Whether to save html file. 
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
    if save: 
        filename = output_folder + outputfilename_prefix + f"_scatter_{col}.html"
        print("Saving...\n", filename)
        fig.write_html(filename)
    return fig 

def dotplot(adata, markers, cluster_header, dendrogram = True, save = False, outputfilename_prefix = ""): 
    """\
    Generating scanpy dotplot of anndata with input marker list. 

    Parameters:
    -----------
        adata: AnnData
            Annotated data matrix.
        markers: list/dict
            List of markers to show in dotplot. Dictionary of markers per cell type to group by. 
        cluster_header: str
            Key in `adata.obs` storing cell type.
        dendrogram: bool/list (default: True)
            Whether to use user-defined dendrogram. Dendrogram order. 
        save: bool (default: False)
            Whether to save png file. 
        outputfilename_prefix: str (default: "")
            Prefix for all output files. 
    """
    if save: 
        print("Saving...\n", outputfilename_prefix + ".png")
        save = outputfilename_prefix + ".png"
    if isinstance(dendrogram, bool): 
        sc.pl.dotplot(adata, markers, cluster_header, use_raw = False, standard_scale = "var", dendrogram = dendrogram, save = save)
    elif isinstance(dendrogram, list): 
        sc.pl.dotplot(adata, markers, cluster_header, use_raw = False, standard_scale = "var", categories_order = dendrogram, save = save)
    return 

def stackedviolin(adata, markers, cluster_header, dendrogram = True, save = False, outputfilename_prefix = ""): 
    """\
    Generating scanpy stacked_violin of anndata with input marker list. 

    Parameters:
    -----------
        adata: AnnData
            Annotated data matrix.
        markers: list/dict
            List of markers to show in dotplot. Dictionary of markers per cell type to group by. 
        cluster_header: str
            Key in `adata.obs` storing cell type.
        dendrogram: bool/list (default: True)
            Whether to use user-defined dendrogram. Dendrogram order. 
        save: bool (default: False)
            Whether to save png file. 
        outputfilename_prefix: str (default: "")
            Prefix for all output files. 
    """
    if save: 
        print("Saving...\n", outputfilename_prefix + ".png")
        save = outputfilename_prefix + ".png"
    if isinstance(dendrogram, bool): 
        sc.pl.stacked_violin(adata, markers, cluster_header, use_raw = False, standard_scale = "var", dendrogram = dendrogram, save = save)
    elif isinstance(dendrogram, list): 
        sc.pl.stacked_violin(adata, markers, cluster_header, use_raw = False, standard_scale = "var", categories_order = dendrogram, save = save)
    return

def matrixplot(adata, markers, cluster_header, dendrogram = True, save = False, outputfilename_prefix = ""): 
    """\
    Generating scanpy matrixplot of anndata with input marker list. 

    Parameters:
    -----------
        adata: AnnData
            Annotated data matrix.
        markers: list/dict
            List of markers to show in dotplot. Dictionary of markers per cell type to group by. 
        cluster_header: str
            Key in `adata.obs` storing cell type.
        dendrogram: bool/list (default: True)
            Whether to use user-defined dendrogram. Dendrogram order. 
        save: bool (default: False)
            Whether to save png file. 
        outputfilename_prefix: str (default: "")
            Prefix for all output files. 
    """
    if save: 
        print("Saving...\n", outputfilename_prefix + ".png")
        save = outputfilename_prefix + ".png"
    if isinstance(dendrogram, bool): 
        # sc.pl.heatmap(adata, markers, cluster_header, use_raw = False, standard_scale = "var", dendrogram = dendrogram, save = save)
        sc.pl.matrixplot(adata, markers, cluster_header, use_raw=False, standard_scale='var', dendrogram = dendrogram, save = save)
    elif isinstance(dendrogram, list): 
        sc.pl.matrixplot(adata, markers, cluster_header, use_raw=False, standard_scale='var', categories_order = dendrogram, save = save)
    return
