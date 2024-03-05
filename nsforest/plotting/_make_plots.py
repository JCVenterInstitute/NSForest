
import pandas as pd
import scanpy as sc
import plotly.express as px
from nsforest.plotting import utils

# TODO: make better
def boxplot(NSForest_results = pd.DataFrame(), output_folder = "", outputfilename = ""): 

    fig = px.box(NSForest_results, y='f_score', points='all', range_y=[-.05,1.05],
                title=f"F-beta score median = {round(NSForest_results['f_score'].median(),3)}",
                width=400, height=500, hover_name='clusterName')
    filename = output_folder + outputfilename + "_boxplot_fscore.html"
    print("Saving...", filename)
    fig.write_html(filename)

    fig = px.box(NSForest_results, y='PPV', points='all', range_y=[-.05,1.05],
             title=f"Positive predictive value median = {round(NSForest_results['PPV'].median(),3)}",
             width=400, height=500, hover_name='clusterName')
    filename = output_folder + outputfilename + "_boxplot_ppv.html"
    print("Saving...", filename)
    fig.write_html(filename)

    fig = px.box(NSForest_results, y='fraction', points='all', range_y=[-.05,1.05],
             title=f"On-target fraction = {round(NSForest_results['fraction'].median(),3)}",
             width=400, height=500, hover_name='clusterName')
    filename = output_folder + outputfilename + "_boxplot_fraction.html"
    print("Saving...", filename)
    fig.write_html(filename)

    return

def scatter_w_clusterSize(NSForest_results = pd.DataFrame(), output_folder = "", outputfilename = ""): 

    fig = px.scatter(NSForest_results, x='clusterSize', y='f_score', range_y=[-.05,1.05],
                 width=700, height=500, hover_name='clusterName')
    filename = output_folder + outputfilename + "_scatter_fscore.html"
    print("Saving...", filename)
    fig.write_html(filename)

    fig = px.scatter(NSForest_results, x='clusterSize', y='PPV', range_y=[-.05,1.05],
                 width=700, height=500, hover_name='clusterName')
    filename = output_folder + outputfilename + "_scatter_ppv.html"
    print("Saving...", filename)
    fig.write_html(filename)

    fig = px.scatter(NSForest_results, x='clusterSize', y='fraction', range_y=[-.05,1.05],
                 width=700, height=500, hover_name='clusterName')
    filename = output_folder + outputfilename + "_scatter_fraction.html"
    print("Saving...", filename)
    fig.write_html(filename)

    return

def dotplot(adata, cluster_header, dendrogram, results_df = pd.DataFrame(), output_folder = "", outputfilename = ""): 

    # sc._settings.ScanpyConfig(figdir = './' + output_folder)
    if isinstance(dendrogram, list): # https://stackoverflow.com/questions/23482668/sorting-by-a-custom-list-in-pandas
        results_df["clusterName"] = results_df["clusterName"].astype("category")
        results_df["clusterName"] = results_df["clusterName"].cat.set_categories(dendrogram)
    markers_dict = utils.prepare_markers(results_df, "clusterName", "NSForest_markers")
    # print("markers_dict:", markers_dict)

    if isinstance(dendrogram, bool): 
        sc.pl.dotplot(adata, markers_dict, cluster_header, use_raw = False, standard_scale = "var", dendrogram = dendrogram, save = outputfilename + ".png")
    elif isinstance(dendrogram, list): 
        sc.pl.dotplot(adata, markers_dict, cluster_header, use_raw = False, standard_scale = "var", categories_order = dendrogram, save = outputfilename + ".png")

    return

def stackedviolin(adata, cluster_header, dendrogram, results_df = pd.DataFrame(), output_folder = "", outputfilename = ""): 

    if isinstance(dendrogram, list): # https://stackoverflow.com/questions/23482668/sorting-by-a-custom-list-in-pandas
        results_df["clusterName"] = results_df["clusterName"].astype("category")
        results_df["clusterName"] = results_df["clusterName"].cat.set_categories(dendrogram)
    markers_dict = utils.prepare_markers(results_df, "clusterName", "NSForest_markers")
    # print("markers_dict:", markers_dict)

    if isinstance(dendrogram, bool): 
        sc.pl.stacked_violin(adata, markers_dict, cluster_header, use_raw = False, standard_scale = "var", dendrogram = dendrogram, save = outputfilename + ".png")
    elif isinstance(dendrogram, list): 
        sc.pl.stacked_violin(adata, markers_dict, cluster_header, use_raw = False, standard_scale = "var", categories_order = dendrogram, save = outputfilename + ".png")

    return

def heatmap(adata, cluster_header, dendrogram, results_df = pd.DataFrame(), output_folder = "", outputfilename = ""): 

    if isinstance(dendrogram, list): # https://stackoverflow.com/questions/23482668/sorting-by-a-custom-list-in-pandas
        results_df["clusterName"] = results_df["clusterName"].astype("category")
        results_df["clusterName"] = results_df["clusterName"].cat.set_categories(dendrogram)
    markers_dict = utils.prepare_markers(results_df, "clusterName", "NSForest_markers")
    # print("markers_dict:", markers_dict)

    if isinstance(dendrogram, bool): 
        sc.pl.heatmap(adata, markers_dict, cluster_header, use_raw = False, standard_scale = "var", dendrogram = dendrogram, save = outputfilename + ".png")
    elif isinstance(dendrogram, list): 
        print("I think scanpy.pl.heatmap is missing categories_order functionality. ")
        print("https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.heatmap.html")
        # sc.pl.heatmap(adata, markers_dict, cluster_header, use_raw = False, standard_scale = "var", categories_order = dendrogram, save = outputfilename + ".png")

    return
