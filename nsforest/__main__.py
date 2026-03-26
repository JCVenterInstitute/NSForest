
import sys
import os
sys.path.insert(0, os.path.abspath("./"))
sys.path.insert(0, os.path.abspath("./nsforest/nsforesting"))
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import nsforest as ns
from nsforest import utils

"""\
If running locally: python3 nsforest
If HPC parallelizing: python3 nsforest -c <cluster>
"""

def main():

    # Preparing kwargs
    args = utils._parse_args(sys.argv[1:], description=__doc__)

    data_folder = "demo_data/"
    file = data_folder + "adata_layer1.h5ad"
    adata = sc.read_h5ad(file)

    cluster_header = "cluster"
    output_folder = "outputs_layer1/"

    #### Preprocessing ####
    ns.pp.dendrogram(adata, cluster_header, save = True, output_folder = output_folder, outputfilename_suffix = cluster_header)

    # Only run these 2 lines if running ns.nsforesting.NSForest
    adata = ns.pp.prep_medians(adata, cluster_header)
    adata = ns.pp.prep_binary_scores(adata, cluster_header)

    filename = file.replace(".h5ad", "_preprocessed.h5ad")
    print(f"Saving new anndata object as...\n{filename}")
    adata.write_h5ad(filename)

    #### Running NSForest ####
    outputfilename_prefix = cluster_header
    if args.c: 
        results = ns.nsforesting.NSForest(adata, cluster_header, cluster_list = [args.c], save_supplementary = True, save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)
    else: 
        results = ns.nsforesting.NSForest(adata, cluster_header, save_supplementary = True, save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)

    # #### Evaluating markers ####
    # outputfilename_prefix = "marker_eval"
    # markers = pd.read_csv("outputs_layer1/cluster_results.csv")
    # markers_dict = utils.prepare_markers(markers, "clusterName", "NSForest_markers")
    # results = ns.ev.DecisionTree(adata, cluster_header, markers_dict, combinations = False, use_mean = False, save_supplementary = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)

    # Loading results csv
    # if parallelizing, remember to combine the separate output csv's into a single results csv before plotting
    # results = pd.read_csv("../outputs_layer1/cluster_results.csv")
    # results["markers"] = utils.str_to_list(results["markers"])

    #### Plotting ####
    # If you want to specify dendrogram order
    dendrogram = list(adata.uns["dendrogram_" + cluster_header]["categories_ordered"])
    if dendrogram: # https://stackoverflow.com/questions/23482668/sorting-by-a-custom-list-in-pandas
        results["clusterName"] = results["clusterName"].astype("category")
        results["clusterName"] = results["clusterName"].cat.set_categories(dendrogram)
        results = results.sort_values("clusterName")
    
    markers_dict = utils.prepare_markers(results, "clusterName", "NSForest_markers")
    
    # # Scanpy plots
    # ns.pl.dotplot(adata, markers_dict, cluster_header, dendrogram = dendrogram, save = True, output_folder = output_folder, outputfilename_suffix = outputfilename_prefix)
    # ns.pl.stackedviolin(adata, markers_dict, cluster_header, dendrogram = dendrogram, save = True, output_folder = output_folder, outputfilename_suffix = outputfilename_prefix)
    # ns.pl.matrixplot(adata, markers_dict, cluster_header, dendrogram = dendrogram, save = True, output_folder = output_folder, outputfilename_suffix = outputfilename_prefix)

    # # Plotly plots
    # ns.pl.boxplot(results, ["f_score", "precision", "recall", "onTarget"], save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)
    # ns.pl.scatter_w_clusterSize(results, "f_score", save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)
    # ns.pl.scatter_w_clusterSize(results, "precision", save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)
    # ns.pl.scatter_w_clusterSize(results, "recall", save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)
    # ns.pl.scatter_w_clusterSize(results, "onTarget", save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)

main()
