
import sys
import os
sys.path.insert(0, os.path.abspath("./"))
sys.path.insert(0, os.path.abspath("./nsforest/nsforesting"))
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import nsforest as ns
from nsforest import utils
from nsforest import preprocessing as pp
from nsforest import nsforesting
from nsforest import evaluating
from nsforest import plotting as pl

"""\
If running locally: python3 nsforest
If HPC parallelizing: python3 nsforest -c <cluster>
"""

def main():

    # Preparing kwargs
    args = utils._parse_args(sys.argv[1:], description=__doc__)

    file = "../demo_data/adata_layer1.h5ad" 
    adata = sc.read_h5ad(file)
    adata.var_names = adata.var["feature_name"]

    cluster_header = "cluster"
    medians_header = "medians_" + cluster_header
    binary_scores_header = "binary_scores_" + cluster_header

    #### Preprocessing ####
    ns.pp.dendrogram(adata, cluster_header) 
    ns.pp.prep_medians(adata, cluster_header, use_mean = False, positive_genes_only = True) 
    ns.pp.prep_binary_scores(adata, cluster_header, medians_header) 

    filename = file.replace(".h5ad", "_preprocessed.h5ad")
    print(f"Saving new anndata object as...\n{filename}")
    adata.write_h5ad(filename)

    # Start here if you already preprocessed the anndata
    adata = sc.read_h5ad(file.replace(".h5ad", "_preprocessed.h5ad"))

    outputfilename = cluster_header
    output_folder = "../outputs_layer1/"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print(f"Creating new directory...\n{output_folder}")

    #### Running NSForest ####
    # Adding this line to actually subset to positive genes if running NSForest
    adata = ns.pp.prep_medians(adata, cluster_header, use_mean = False, positive_genes_only = True) 
    if args.c: 
        nsforesting.NSForest(adata, cluster_header, medians_header, binary_scores_header, cluster_list = [args.c], output_folder = output_folder, outputfilename = outputfilename)
    else: 
        results = nsforesting.NSForest(adata, cluster_header, medians_header, binary_scores_header, output_folder = output_folder, outputfilename = outputfilename)

    #### Evaluating markers ####
    outputfilename = "markers"
    markers = pd.read_csv("../outputs_layer1/cluster_results.csv")
    markers_dict = utils.prepare_markers(markers, "clusterName", "NSForest_markers")
    adata = sc.read_h5ad(file.replace(".h5ad", "_preprocessed.h5ad"))
    results = evaluating.DecisionTree(adata, cluster_header, medians_header, markers_dict, combinations = False, output_folder = output_folder, outputfilename = outputfilename)

    # Loading results csv
    # if parallelizing, remember to combine the separate output csv's into a single results csv before plotting
    # results = pd.read_csv("../outputs_layer1/cluster_results.csv")
    # results["markers"] = utils.str_to_list(results["markers"])

    #### Plotting ####
    ns.pl.boxplot_fscore(results, True, output_folder, outputfilename)
    ns.pl.boxplot_ppv(results, True, output_folder, outputfilename)
    ns.pl.boxplot_ontarget(results, True, output_folder, outputfilename)
    ns.pl.scatter_w_clusterSize_fscore(results, True, output_folder, outputfilename)
    ns.pl.scatter_w_clusterSize_ppv(results, True, output_folder, outputfilename)
    ns.pl.scatter_w_clusterSize_ontarget(results, True, output_folder, outputfilename)

    # If you want to specify dendrogram order
    dendrogram = list(adata.uns["dendrogram_" + cluster_header]["categories_ordered"])
    if dendrogram: # https://stackoverflow.com/questions/23482668/sorting-by-a-custom-list-in-pandas
        results["clusterName"] = results["clusterName"].astype("category")
        results["clusterName"] = results["clusterName"].cat.set_categories(dendrogram)
        results = results.sort_values("clusterName")
    
    markers = []
    for val in results["NSForest_markers"]: # markers NSForest_markers binary_markers
        markers.extend(val)
    
    ns.pl.dotplot(adata, markers, cluster_header, dendrogram, True, output_folder, outputfilename)
    ns.pl.stackedviolin(adata, markers, cluster_header, dendrogram, True, output_folder, outputfilename)
    ns.pl.matrixplot(adata, markers, cluster_header, True, True, output_folder, outputfilename)

main()
