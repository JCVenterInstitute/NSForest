#!/usr/bin/env python

import sys
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
CODE_PATH = "./" # location of NSForest folder
sys.path.insert(0, os.path.abspath(CODE_PATH))
from nsforest import ns, nsforesting, utils, NSFOREST_VERSION

data_folder = f"{CODE_PATH}/demo_data/"
file = data_folder + "adata_layer1.h5ad"
adata = sc.read_h5ad(file)

cluster_header = "cluster"

output_folder = f"{CODE_PATH}/outputs_layer1/"

ns.pp.dendrogram(adata, cluster_header, save = True, output_folder = output_folder, outputfilename_suffix = cluster_header)

adata_subset = adata.copy()

adata_subset = ns.pp.prep_medians(adata_subset, cluster_header)

adata_subset = ns.pp.prep_binary_scores(adata_subset, cluster_header)

outputfilename_prefix = cluster_header
results = nsforesting.NSForest(adata_subset, cluster_header, save = True, save_supplementary = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)
