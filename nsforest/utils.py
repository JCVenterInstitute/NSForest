
import pandas as pd
import argparse

# Simple argument parsing
def _parse_args(args, **kwargs): 
    parser = argparse.ArgumentParser(**kwargs, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-c", metavar = "cluster", type=str, 
                        help="parallelizing, one cluster at a time")
    return parser.parse_args()
  
def str_to_list(values): 
    """\
    Converting str representation of list to list

    Parameters:
    -----------
        values: str 
            String representation of list
    
    Returns:
    ========
        values: list
    """
    values = [val.replace("[", "").replace("]", "").replace(", ", ",").replace("'", "").replace('"', "").split(",") for val in values]
    return values

def prepare_markers(results, col_cluster, col_marker): 
    """\
    Converting cell type to marker dictionary. 

    Parameters:
    -----------
        results: pd.DataFrame
            DataFrame with `col_cluster` and `col_marker`. 
        col_cluster
            Column name in `results` storing cell type.
        col_marker
            Column name in `results` storing markers.
    
    Returns:
    ========
        marker_genes_dict: dict
            Dictionary of clusterName: marker_genes
    """
    # marker_genes_csv = pd.read_csv(filename) 
    not_found = []
    marker_genes_dict = {}
    # if markers are represented as list
    if "[" in list(results[col_marker])[0]: # "['NTNG1', 'EYA4']"
        marker_genes_dict = dict(zip(results[col_cluster], str_to_list(results[col_marker])))
    else: # 'NTNG1'\n'EYA4'
        for cluster, marker in zip(results[col_cluster], results[col_marker]): 
            if cluster not in marker_genes_dict: 
                marker_genes_dict[cluster] = [marker]
            else: 
                marker_genes_dict[cluster].append(marker)
    # print out too
    if len(not_found) > 1: 
        print("WARNING: input markers not found in anndata\nMarkers:", not_found)
        # df = pd.DataFrame()
        # df[col_marker] = not_found
        # df.to_csv(output_folder + outputfilename_prefix + "_not_found.csv", index = False, header = False)

    return marker_genes_dict 
