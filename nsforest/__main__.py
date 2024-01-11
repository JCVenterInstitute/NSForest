
_epilog = ""

import sys
import os
import scanpy as sc
from . import utils
from . import nsforest

def main():

    # Reading in arguments
    args = utils._parse_args(sys.argv[1:],
                            epilog=_epilog,
                            description=__doc__)
    
    args_1 = sys.argv[1] # filename
    args_2 = sys.argv[2] # cluster
    print(args_1, args_2)

    if os.path.isfile(args_1): 
        file = args_1
    else: 
        print("filepath not found")
    cluster_header = args_2

    adata = sc.read_h5ad(file)

    nsforest.NSForest(adata, cluster_header = cluster_header, gene_selection = "BinaryFirst_high", n_trees=10)

    # python3 nsforest.py demo_data/adata_layer1.h5ad

main()