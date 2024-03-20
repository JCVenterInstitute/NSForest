
import pandas as pd
import scanpy as sc
import plotly.express as px

def boxplot_fscore(df, save = False, output_folder = "", outputfilename = ""): 
    fig = px.box(df, y='f_score', points='all', range_y=[-.05,1.05],
                 title=f"F-beta score median = {round(df['f_score'].median(),3)}",
                 width=400, height=500, hover_name='clusterName')
    if save: 
        filename = output_folder + outputfilename + "_boxplot_fscore.html"
        print("Saving...\n", filename)
        fig.write_html(filename)
    return fig

def boxplot_ppv(df, save = False, output_folder = "", outputfilename = ""): 
    fig = px.box(df, y='PPV', points='all', range_y=[-.05,1.05],
                 title=f"Positive predictive value median = {round(df['PPV'].median(),3)}",
                 width=400, height=500, hover_name='clusterName')
    if save: 
        filename = output_folder + outputfilename + "_boxplot_ppv.html"
        print("Saving...\n", filename)
        fig.write_html(filename)
    return fig

def boxplot_fraction(df, save = False, output_folder = "", outputfilename = ""): 
    fig = px.box(df, y='fraction', points='all', range_y=[-.05,1.05],
                 title=f"On-target fraction = {round(df['fraction'].median(),3)}",
                 width=400, height=500, hover_name='clusterName')
    if save: 
        filename = output_folder + outputfilename + "_boxplot_fraction.html"
        print("Saving...\n", filename)
        fig.write_html(filename)
    return fig

def scatter_w_clusterSize_fscore(df, save = False, output_folder = "", outputfilename = ""): 
    fig = px.scatter(df, x='clusterSize', y='f_score', range_y=[-.05,1.05],
                     width=700, height=500, hover_name='clusterName')
    if save: 
        filename = output_folder + outputfilename + "_scatter_fscore.html"
        print("Saving...\n", filename)
        fig.write_html(filename)
    return fig 

def scatter_w_clusterSize_ppv(df, save = False, output_folder = "", outputfilename = ""): 
    fig = px.scatter(df, x='clusterSize', y='PPV', range_y=[-.05,1.05],
                     width=700, height=500, hover_name='clusterName')
    if save: 
        filename = output_folder + outputfilename + "_scatter_ppv.html"
        print("Saving...\n", filename)
        fig.write_html(filename)
    return fig 

def scatter_w_clusterSize_fraction(df, save = False, output_folder = "", outputfilename = ""): 
    fig = px.scatter(df, x='clusterSize', y='fraction', range_y=[-.05,1.05], 
                     width=700, height=500, hover_name='clusterName')
    if save: 
        filename = output_folder + outputfilename + "_scatter_fraction.html"
        print("Saving...\n", filename)
        fig.write_html(filename)
    return fig 

def dotplot(adata, markers, cluster_header, dendrogram = True, save = False, output_folder = "", outputfilename = ""): 
    # sc._settings.ScanpyConfig(figdir = './' + output_folder)
    if save: 
        print("Saving...\n", outputfilename + ".png")
        save = outputfilename + ".png"
    if isinstance(dendrogram, bool): 
        sc.pl.dotplot(adata, markers, cluster_header, use_raw = False, standard_scale = "var", dendrogram = dendrogram, save = save)
    elif isinstance(dendrogram, list): 
        sc.pl.dotplot(adata, markers, cluster_header, use_raw = False, standard_scale = "var", categories_order = dendrogram, save = save)
    return 

def stackedviolin(adata, markers, cluster_header, dendrogram = True, save = False, output_folder = "", outputfilename = ""): 
    if save: 
        print("Saving...\n", outputfilename + ".png")
        save = outputfilename + ".png"
    if isinstance(dendrogram, bool): 
        sc.pl.stacked_violin(adata, markers, cluster_header, use_raw = False, standard_scale = "var", dendrogram = dendrogram, save = save)
    elif isinstance(dendrogram, list): 
        sc.pl.stacked_violin(adata, markers, cluster_header, use_raw = False, standard_scale = "var", categories_order = dendrogram, save = save)
    return

def matrixplot(adata, markers, cluster_header, dendrogram = True, save = False, output_folder = "", outputfilename = ""): 
    if save: 
        print("Saving...\n", outputfilename + ".png")
        save = outputfilename + ".png"
    if isinstance(dendrogram, bool): 
        # sc.pl.heatmap(adata, markers, cluster_header, use_raw = False, standard_scale = "var", dendrogram = dendrogram, save = save)
        sc.pl.matrixplot(adata, markers, cluster_header, use_raw=False, standard_scale='var', dendrogram = dendrogram, save = save)
    elif isinstance(dendrogram, list): 
        sc.pl.matrixplot(adata, markers, cluster_header, use_raw=False, standard_scale='var', categories_order = dendrogram, save = save)
        # print("I think scanpy.pl.heatmap is missing categories_order functionality.")
        # print("https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.heatmap.html")
        # sc.pl.heatmap(adata, markers, cluster_header, use_raw = False, standard_scale = "var", categories_order = dendrogram, save = outputfilename + ".png")
    return
