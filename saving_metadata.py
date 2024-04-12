
import numpy as np
import pandas as pd
import scanpy as sc

kwargs = {}
# kwargs["h5ad"] = "demo_data/hlca_core.h5ad"
# kwargs["cluster_header"] = "ann_finest_level"
# kwargs["cluster_header"] = "ann_level_1"
# kwargs["cluster_header"] = "ann_level_2"
# kwargs["input_folder"] = "inputs_hlca_core_hlca_markers"

kwargs["h5ad"] = "demo_data/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad"
kwargs["cluster_header"] = "Granular cell type"
kwargs["input_folder"] = "inputs_gtex"

adata = sc.read_h5ad(kwargs["h5ad"])
print(adata)

# Saving list of clusters as csv
df = pd.DataFrame()
df["cluster"] = np.unique(adata.obs[kwargs["cluster_header"]])
kwargs["cluster_header"] = kwargs["cluster_header"].lower().replace(" ", "_")
filename = f'{kwargs["input_folder"]}/clusters_{kwargs["cluster_header"]}.csv'
print("Saving list of clusters as...\n", filename)
df.to_csv(filename, index = False, header = False)

# Saving ensembl_id and feature_name mapping as csv
df = adata.var.reset_index()
filename = f'{kwargs["input_folder"]}/ensembl_id_and_feature_name.csv'
print("Saving ensembl_id and feature_name mapping as...\n", filename)
# df[["ensembl_id", "feature_name"]].to_csv(filename, index = False)
df[["gene_ids", "gene_name"]].to_csv(filename, index = False)

# # Saving cell type ontology term mapping as csv
# df = adata.obs[["cell_type_ontology_term_id", "ann_finest_level", "ann_level_1", "ann_level_2", "ann_level_3", "ann_level_4", "ann_level_5"]]
# order = ["ann_level_1", "ann_level_2", "ann_level_3", "ann_level_4", "ann_level_5", "ann_finest_level"]
# df = df[~df.duplicated(keep = "first")].sort_values(order)
# filename = f'{kwargs["input_folder"]}/cell_type_ontology_terms.csv'
# print("Saving ensembl_id and feature_name mapping as...\n", filename)
# df.to_csv(filename, index = False)
