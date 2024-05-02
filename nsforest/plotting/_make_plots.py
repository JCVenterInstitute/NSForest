
import pandas as pd
import scanpy as sc
import plotly.express as px

def boxplot(df, col, save = False, output_folder = "", outputfilename_prefix = ""): 
    """\
    Generating plotly boxplot of specified column in df. 

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
    Generating plotly scatterplot of specified column by `clusterSize. 

    Parameters:
    -----------
        df: pd.DataFrame
            NS-Forest results containing `clusterName`, `clusterSize`, and `col`. 
        col: str
            Column in `df` to create the scatterplot. 
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
    fig = px.scatter(df, x='clusterSize', y=col, range_y=[-.05,1.05],
                     width=700, height=500, hover_name='clusterName')
    if save: 
        filename = output_folder + outputfilename_prefix + f"_scatter_{col}.html"
        print("Saving...\n", filename)
        fig.write_html(filename)
    return fig 

def boxplot_fscore(df, save = False, output_folder = "", outputfilename_prefix = ""): 
    """\
    Generating plotly boxplot of f_score. 

    Parameters:
    -----------
        df: pd.DataFrame
            NS-Forest results containing `f_score` column. 
        save: bool (default: False)
            Whether to save html file. 
        output_folder
            Output folder. 
        outputfilename_prefix
            Prefix for all output files. 
    
    Returns:
    ========
    fig: plotly figure
    """
    fig = px.box(df, y='f_score', points='all', range_y=[-.05,1.05],
                 title=f"F-beta score median = {round(df['f_score'].median(),3)}",
                 width=400, height=500, hover_name='clusterName')
    if save: 
        filename = output_folder + outputfilename_prefix + "_boxplot_fscore.html"
        print("Saving...\n", filename)
        fig.write_html(filename)
    return fig

def boxplot_ppv(df, save = False, output_folder = "", outputfilename_prefix = ""): 
    """\
    Generating plotly boxplot of PPV. 

    Parameters:
    -----------
    df: pd.DataFrame
        NS-Forest results containing PPV column. 
    save: bool
        Whether to save html file. 
    output_folder
        Output folder. 
    outputfilename_prefix
        Prefix for all output files. 
    
    Returns:
    ========
    fig: plotly figure
    """
    fig = px.box(df, y='PPV', points='all', range_y=[-.05,1.05],
                 title=f"Positive predictive value median = {round(df['PPV'].median(),3)}",
                 width=400, height=500, hover_name='clusterName')
    if save: 
        filename = output_folder + outputfilename_prefix + "_boxplot_ppv.html"
        print("Saving...\n", filename)
        fig.write_html(filename)
    return fig

def boxplot_ontarget(df, save = False, output_folder = "", outputfilename_prefix = ""): 
    """\
    Generating plotly boxplot of onTarget fraction. 

    Parameters:
    -----------
    df: pd.DataFrame
        NS-Forest results containing onTarget column. 
    save: bool
        Whether to save html file. 
    output_folder
        Output folder. 
    outputfilename_prefix
        Prefix for all output files. 
    
    Returns:
    ========
    fig: plotly figure
    """
    fig = px.box(df, y='onTarget', points='all', range_y=[-.05,1.05],
                 title=f"On-target fraction = {round(df['onTarget'].median(),3)}",
                 width=400, height=500, hover_name='clusterName')
    if save: 
        filename = output_folder + outputfilename_prefix + "_boxplot_ontarget.html"
        print("Saving...\n", filename)
        fig.write_html(filename)
    return fig

def scatter_w_clusterSize_fscore(df, save = False, output_folder = "", outputfilename_prefix = ""): 
    """\
    Generating plotly scatterplot of f_score by clusterSize. 

    Parameters:
    -----------
    df: pd.DataFrame
        NS-Forest results containing f_score and clusterSize columns. 
    save: bool
        Whether to save html file. 
    output_folder
        Output folder. 
    outputfilename_prefix
        Prefix for all output files. 
    
    Returns:
    ========
    fig: plotly figure
    """
    fig = px.scatter(df, x='clusterSize', y='f_score', range_y=[-.05,1.05],
                     width=700, height=500, hover_name='clusterName')
    if save: 
        filename = output_folder + outputfilename_prefix + "_scatter_fscore.html"
        print("Saving...\n", filename)
        fig.write_html(filename)
    return fig 

def scatter_w_clusterSize_ppv(df, save = False, output_folder = "", outputfilename_prefix = ""): 
    """\
    Generating plotly scatterplot of PPV by clusterSize. 

    Parameters:
    -----------
    df: pd.DataFrame
        NS-Forest results containing PPV and clusterSize columns. 
    save: bool
        Whether to save html file. 
    output_folder
        Output folder. 
    outputfilename_prefix
        Prefix for all output files. 
    
    Returns:
    ========
    fig: plotly figure
    """
    fig = px.scatter(df, x='clusterSize', y='PPV', range_y=[-.05,1.05],
                     width=700, height=500, hover_name='clusterName')
    if save: 
        filename = output_folder + outputfilename_prefix + "_scatter_ppv.html"
        print("Saving...\n", filename)
        fig.write_html(filename)
    return fig 

def scatter_w_clusterSize_ontarget(df, save = False, output_folder = "", outputfilename_prefix = ""): 
    """\
    Generating plotly scatterplot of onTarget fraction by clusterSize. 

    Parameters:
    -----------
    df: pd.DataFrame
        NS-Forest results containing onTarget and clusterSize columns. 
    save: bool
        Whether to save html file. 
    output_folder
        Output folder. 
    outputfilename_prefix
        Prefix for all output files. 
    
    Returns:
    ========
    fig: plotly figure
    """
    fig = px.scatter(df, x='clusterSize', y='onTarget', range_y=[-.05,1.05], 
                     width=700, height=500, hover_name='clusterName')
    if save: 
        filename = output_folder + outputfilename_prefix + "_scatter_ontarget.html"
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
    markers: dict
        clusterName: list of markers
    dendrogram: bool/list
        Whether to use user-defined dendrogram. Dendrogram order. 
    save: bool
        Whether to save png file. 
    outputfilename_prefix
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
    markers: dict
        clusterName: list of markers
    dendrogram: bool/list
        Whether to use user-defined dendrogram. Dendrogram order. 
    save: bool
        Whether to save png file. 
    outputfilename_prefix
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
    markers: dict
        clusterName: list of markers
    dendrogram: bool/list
        Whether to use user-defined dendrogram. Dendrogram order. 
    save: bool
        Whether to save png file. 
    outputfilename_prefix
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
