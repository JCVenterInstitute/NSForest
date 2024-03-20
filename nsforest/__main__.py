
import sys
import os
sys.path.insert(0, os.path.abspath("./"))
sys.path.insert(0, os.path.abspath("./nsforest/nsforesting"))
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import nsforest as ns
from nsforest import preprocessing as pp
from nsforest.preprocessing import _add_ann
from nsforest import nsforesting
from nsforest import utils
from nsforest import evaluating
from nsforest.evaluating import _run_markers
from nsforest import plotting as pl
from nsforest.plotting import _make_plots

def main():

    # Loading h5ad
    file = "demo_data/hlca_core.h5ad"
    cluster_headers = ["ann_finest_level", "ann_level_1", "ann_level_2", "ann_level_3", "ann_level_4", "ann_level_5"]
    for cluster_header in cluster_headers: 
        adata = sc.read_h5ad(file)
        medians_header = "medians_" + cluster_header
        binary_scores_header = "binary_scores_" + cluster_header
        ns.pp.prep_medians(adata, cluster_header, positive_genes_only = False)
        ns.pp.prep_binary_scores(adata, cluster_header, medians_header)
        ns.pp.dendrogram(adata, cluster_header) # using all genes or positive genes?? rn using positive genes
        adata.write_h5ad(file.replace(".h5ad", "_" + cluster_header + ".h5ad"))
    # file = "demo_data/hlca_core_preprocessed.h5ad"
    # file = "demo_data/hlca_core_ann_finest_level_precalculated.h5ad"
    # adata = sc.read_h5ad(file)
    # adata.var_names = adata.var["feature_name"]
    # print(adata)

    # cluster_header = "ann_finest_level"
    # medians_header = "medians_" + cluster_header
    # binary_scores_header = "binary_scores_" + cluster_header
    # outputfilename = cluster_header
    # output_folder = "test/"

    # markers = pd.read_csv("inputs_hlca_core_hlca_markers/marker_genes_hlca_core_ann_finest_level.csv")
    # markers_dict = utils.prepare_markers(markers, "cell_annotation", "marker_gene")
    # results = evaluating.DecisionTree(adata, cluster_header, markers_dict, combinations = False, output_folder = output_folder, outputfilename = outputfilename)
    # results.to_csv(output_folder + outputfilename + "_results.csv", index = False)

main()
