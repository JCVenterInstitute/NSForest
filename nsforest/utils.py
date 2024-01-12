
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm

def _parse_args(args, **kwargs): 
  parser = argparse.ArgumentParser(**kwargs,
                                    formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument("-a", metavar = "arguments", type=str, 
                      help="arguments csv, if you want to specify more nsforest arguments than just filename and cluster")
  
  parser.add_argument("-f", metavar = "filename", type=str, 
                      help="anndata object")
  
  parser.add_argument("-c", metavar = "cluster", type=str, 
                      help="specify the header name of clusters in adata.obs")
  
  return parser.parse_args()

def get_medians(adata, cluster_header):
    cluster_medians = pd.DataFrame()
    for cl in tqdm(sorted(set(adata.obs[cluster_header])), desc="Calculating medians per cluster"):
        adata_cl = adata[adata.obs[cluster_header]==cl,]
        medians_cl = adata_cl.to_df().median()
        cluster_medians = pd.concat([cluster_medians, pd.DataFrame({cl: medians_cl})], axis=1) #gene-by-cluster
    return cluster_medians

def preprocessing_medians(adata, cluster_header, positive_genes_only=True):
    ## get medians
    cluster_medians = get_medians(adata, cluster_header) #gene-by-cluster
    ## attach pre-calculated medians to adata
    print("Saved pre-calculated medians in adata.varm.")
    adata.varm['medians_' + cluster_header] = cluster_medians #gene-by-cluster
    
    if positive_genes_only:
        print("Only positive genes are selected.")
        ## select only genes with median > 0
        genes_selected = cluster_medians.index[cluster_medians.sum(axis=1)>0].to_list()
        ## subset data with only positive genes
        adata = adata[:,genes_selected].copy()
    return adata

def preprocessing_binary(adata, cluster_header, medians_header):
    ## get medians
    cluster_medians = adata.varm[medians_header].transpose() #cluster-by-gene
    n_total_clusters = cluster_medians.shape[0]
    
    ## calculate binary scores based on cluster_medians for all genes
    binary_scores = []
    for cl in tqdm(cluster_medians.index, desc="Calculating binary scores per cluster"):
        ## get binary scores for all genes in each row (cluster) in df
        binary_scores_cl = [sum(np.maximum(0,1-cluster_medians[i]/cluster_medians.loc[cl,i]))/(n_total_clusters-1) for i in cluster_medians.columns]
        binary_scores.append(binary_scores_cl)

    ## binary scores matrix and handle nan
    binary_scores = pd.DataFrame(binary_scores, index=cluster_medians.index, columns=cluster_medians.columns).fillna(0) #cluster-by-gene
    
    ## attach pre-calculated binary scores to adata
    print("Saved pre-calculated binary scores in adata.varm.")
    adata.varm['binary_scores_' + cluster_header] = binary_scores.transpose() #gene-by-cluster
    ##
    print("median:", binary_scores.stack().median())
    print("mean:", binary_scores.stack().mean())
    print("std. dev.:", binary_scores.stack().std())
    
    return adata
