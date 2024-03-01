
import pandas as pd
import plotly.express as px
import scanpy as sc
import utils

def plot_results(NSForest_results = pd.DataFrame(), output_folder = "", outputfilename = ""): 

    fig = px.box(NSForest_results, y='f_score', points='all', range_y=[-.05,1.05],
                title=f"F-beta score median = {round(NSForest_results['f_score'].median(),3)}",
                width=400, height=500, hover_name='clusterName')
    filename = output_folder + outputfilename + "_boxplot_fscore.html"
    print("Saving...", filename)
    fig.write_html(filename)

    fig = px.scatter(NSForest_results, x='clusterSize', y='f_score', range_y=[-.05,1.05],
                 width=700, height=500, hover_name='clusterName')
    filename = output_folder + outputfilename + "_scatter_fscore.html"
    print("Saving...", filename)
    fig.write_html(filename)

    fig = px.box(NSForest_results, y='PPV', points='all', range_y=[-.05,1.05],
             title=f"Positive predictive value median = {round(NSForest_results['PPV'].median(),3)}",
             width=400, height=500, hover_name='clusterName')
    filename = output_folder + outputfilename + "_boxplot_ppv.html"
    print("Saving...", filename)
    fig.write_html(filename)

    fig = px.scatter(NSForest_results, x='clusterSize', y='PPV', range_y=[-.05,1.05],
                 width=700, height=500, hover_name='clusterName')
    filename = output_folder + outputfilename + "_scatter_ppv.html"
    print("Saving...", filename)
    fig.write_html(filename)

    return

def plot_scanpy(adata, cluster_header, nsf_results_df = pd.DataFrame(), gene_symbols = "", output_folder = "", outputfilename = ""): 

    sc._settings.ScanpyConfig(figdir = './' + output_folder)
    
    markers_dict = {}
    markers = dict(zip(nsf_results_df["clusterName"], nsf_results_df["NSForest_markers"]))
    dend_header = "dendrogram_" + cluster_header
    dend_order = adata.uns[dend_header].get('categories_ordered')
    for cluster in dend_order: 
        if cluster in markers: 
            values = list(set(markers[cluster]) & set(adata.var[gene_symbols])) 
            markers_dict[cluster] = values

    print("markers_dict:", markers_dict)

    # TODO: Hard coding fix
    # adata.raw.var[gene_symbols] = adata.raw.var.index

    sc.pl.dotplot(adata, markers_dict, cluster_header, gene_symbols = gene_symbols, use_raw = False, standard_scale = "var", dendrogram = True, save = outputfilename + ".png")

    sc.pl.stacked_violin(adata, markers_dict, cluster_header, gene_symbols = gene_symbols, use_raw = False, standard_scale = "var", dendrogram = True, save = outputfilename + ".png")

    return
