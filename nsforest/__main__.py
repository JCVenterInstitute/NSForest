
_epilog = """
Examples

python3 nsforest -a arguments.csv -p
python3 nsforest -a arguments_layer1.csv

"""

import sys
import numpy as np
import pandas as pd
import scanpy as sc
import utils
import nsforest
import decisiontreewithmarkerlist
import plotting
import calculate_diagonals
import os

def main():

    # Preparing kwargs
    args = utils._parse_args(sys.argv[1:], epilog=_epilog, description=__doc__)
    kwargs = utils.get_kwargs(args.a)

    # If parallelizing NSForest
    if args.c: 
        kwargs["cluster_list"] = [args.c.replace("\r", "").replace('"', "").replace("'", "")]
        if not os.path.exists(kwargs["output_folder"]): 
            os.makedirs(kwargs["output_folder"])
        kwargs["output_folder"] = kwargs["output_folder"] + "nsforest_percluster/"
        kwargs["outputfilename"] = kwargs["cluster_list"][0].replace(" ", "_").replace("/", "_")
    
    # Printing inputs
    print("Input values...")
    print(kwargs)

    # Setting up folder directories
    if not os.path.exists(kwargs["output_folder"]):
        os.makedirs(kwargs["output_folder"])

    # Loading h5ad file
    adata = sc.read_h5ad(kwargs["h5ad"])
    print("adata...")
    print(adata)
    
    # Calculating median gene expression matrix if no header given
    to_save = False
    if "medians_header" not in kwargs: 
        adata = utils.preprocessing_medians(adata, kwargs["cluster_header"])
        kwargs["medians_header"] = 'medians_' + kwargs["cluster_header"]
        to_save = True
    # Calculating binary score gene expression matrix if no header given
    if "binary_scores_header" not in kwargs: 
        adata = utils.preprocessing_binary(adata, kwargs["cluster_header"], kwargs["medians_header"])
        kwargs["binary_scores_header"] = 'binary_scores_' + kwargs["cluster_header"]
        to_save = True
    # Saving new h5ad
    if to_save: 
        filename = kwargs["h5ad"].replace(".h5ad", "_" + kwargs["cluster_header"] + "_precalculated.h5ad")
        print("Saving new anndata object as...\n", filename)
        adata.write_h5ad(filename)
    else: 
        print(f'Using pre-calculated cluster_medians and binary_scores in {kwargs["h5ad"]} from {kwargs["medians_header"]} {kwargs["binary_scores_header"]}')
    print(f'\tanndata.varm[{kwargs["medians_header"]}] and anndata.varm[{kwargs["binary_scores_header"]}]')

    if "compute_decision_trees" in kwargs and kwargs["compute_decision_trees"]: 
        # Using nsforest to generate marker genes
        if "marker_genes_csv" not in kwargs: 
            nsf_kwargs = utils.subset_kwargs("nsf.NSF", kwargs)
            nsforest.NSForest(adata, **nsf_kwargs)
        # Inputting marker genes file to create decision trees
        else: 
            # Changing feature_name to ensembl_id
            mapping = pd.read_csv(kwargs["ensembl_mapping"], header = None)
            mapping = dict(zip(mapping[1], mapping[0]))
            kwargs["marker_genes_dict"] = utils.prepare_markers(kwargs["marker_genes_csv"], 
                                                                kwargs["marker_genes_csv_clustercol"], 
                                                                kwargs["marker_genes_csv_markercol"], 
                                                                kwargs["cluster_header"], 
                                                                mapping)
            dt_kwargs = utils.subset_kwargs("dtwml.DT", kwargs)
            decisiontreewithmarkerlist.DecisionTree(adata, **dt_kwargs)
    
    if "combine_results" in kwargs and kwargs["combine_results"]: 
        # input_folder, output_folder, outputfilename
        utils.combine_results(kwargs["output_folder"] + "nsforest_percluster/", 
                              kwargs["output_folder"], 
                              kwargs["outputfilename"])

    # Getting output marker list
    if "calculate_diagonals" in kwargs and kwargs["calculate_diagonals"]: 
        # 'clusterName', 'markerGene', 'score'
        nsf_markers_df = pd.read_csv(kwargs["output_folder"] + kwargs["outputfilename"] + "_markers.csv")
        if "subclade_list" in kwargs: 
            subclade_list = utils.str_to_list(kwargs["subclade_list"])
            nsf_markers_df = nsf_markers_df[nsf_markers_df['clusterName'].isin(subclade_list)].reset_index(drop=True)
        cluster_medians = adata.varm[kwargs["medians_header"]].transpose()
        # calculate_diagonals(nsf_markers_df, cluster_medians, output_folder = "", outputfilename = "")
        calculate_diagonals.calculate_diagonals(nsf_markers_df, cluster_medians, kwargs["output_folder"], kwargs["outputfilename"])
        
        # 'clusterName', 'clusterSize', 'f_score', 'PPV', 'TN', 'FP', 'FN', 'TP', 'marker_count', 'NSForest_markers', 'binary_genes'
        nsf_results_df = pd.read_csv(kwargs["output_folder"] + kwargs["outputfilename"] + "_results.csv")
        nsf_results_df["NSForest_markers"] = utils.str_to_list(nsf_results_df["NSForest_markers"])
        dend_header = "dendrogram_" + kwargs["cluster_header"]
        if dend_header not in adata.uns: 
            if kwargs["cluster_header"] in ["ann_finest_level", "Granular cell type"]: 
                sc.set_figure_params(scanpy=True, fontsize=4)
            sc.pl.dendrogram(adata, kwargs["cluster_header"], orientation = "top", save = "_" + kwargs["outputfilename"] + ".png")
            print("Saving new anndata object as...\n", kwargs["h5ad"])
            adata.write_h5ad(kwargs["h5ad"])
        # on_target_ratio(adata, nsf_results_df, cluster_header, medians_header, output_folder = "", outputfilename = "")
        calculate_diagonals.on_target_ratio(adata, nsf_results_df, kwargs["cluster_header"], kwargs["medians_header"], kwargs["output_folder"], kwargs["outputfilename"])
    
    # Plotting results
    if "plot_results" in kwargs and kwargs["plot_results"]: 
        print(kwargs["output_folder"] + kwargs["outputfilename"] + "_results.csv")
        nsf_results_df = pd.read_csv(kwargs["output_folder"] + kwargs["outputfilename"] + "_results.csv")
        nsf_results_df["NSForest_markers"] = utils.str_to_list(nsf_results_df["NSForest_markers"])
        # Changing ensembl_id to feature_name
        if "ENSG" in list(nsf_results_df["NSForest_markers"])[0][0]: 
            mapping = pd.read_csv(kwargs["ensembl_mapping"], header = None)
            mapping = dict(zip(mapping[0], mapping[1]))
            values = []
            for markers in nsf_results_df["NSForest_markers"]: 
                markers = [mapping[marker] for marker in markers]
                values.append(markers)
            nsf_results_df["NSForest_markers"] = values
        # plot_results(NSForest_results = pd.DataFrame(), output_folder = None, outputfilename = None)
        plotting.plot_results(nsf_results_df, kwargs["output_folder"], kwargs["outputfilename"])
        # plot_scanpy(adata, cluster_header, nsf_results_df = pd.DataFrame(), output_folder = "", outputfilename = "")
        plotting.plot_scanpy(adata, kwargs["cluster_header"], nsf_results_df, kwargs["gene_symbols"], kwargs["output_folder"], kwargs["outputfilename"])

main()
