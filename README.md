# Necessary and Sufficient Forest (NS-Forest) for Cell Type Marker Determination from cell type clusters

## Getting Started

Install python 3.6 or above
Download NSForest_v3.py file


### Prerequisites

* The  is a python script written in python 3.6. Required libraries: Numpy, Pandas, Sklearn, graphviz, numexpr, scanpy
* scanpy object (adata) with at least one column containing the cluster assignments. Default slot set to adata.obs["louvain"]; however parameter is tunable in function call.


### Description

Necessary and Sufficient Forest is a method that takes cluster results from single cell/nuclei RNAseq experiments 
and generates lists of minimal markers needed to define each “cell type cluster”. 
 
The method begins by re-encoding the cluster labels into binary classifications, and Random Forest models are generated comparing each 
cluster versus all. The top fifteen genes are then reranked using a score measuring how binary they are, e.g., a gene with expression in
the target cluster but no expression in the other clusters would have a high binary score. Expression cutoffs for the top six genes ranked
by binary score are then determined by generating individual decision trees and extracting the decision path information. Then all combinations 
of the top six most binary genes are evaluated using f-beta score as an objective function (the beta value default set at 0.5, which weights the 
f-measure score more toward precision as opposed to recall). 


See code for detailed comments. 


## Versioning

This is version 3.0 The earlier releases were described in the below publications.  

Version 2

Aevermann BD, Zhang Y, Novotny M, Keshk M, Bakken TE, Miller JA, Hodge RD, Lelieveldt B, Lein ES, Scheuermann RH. A machine learning method for the discovery of minimum marker gene combinations for cell-type identification from single-cell RNA sequencing. Genome Res. 2021 Jun 4:gr.275569.121. doi: 10.1101/gr.275569.121. Epub ahead of print. PMID: 34088715.

version 1.3/1.0:

Aevermann BD, Novotny M, Bakken T, Miller JA, Diehl AD, Osumi-Sutherland D, Lasken RS, Lein ES, Scheuermann RH.
Cell type discovery using single-cell transcriptomics: implications for ontological representation. 
Hum Mol Genet. 2018 May 1;27(R1):R40-R47. doi: 10.1093/hmg/ddy100.


## Authors

* Brian Aevermann baeverma@jcvi.org and Richard Scheuermann RScheuermann@jcvi.org


## License

This project is licensed under the MIT License - see the https://opensource.org/licenses/MIT for details

## Acknowledgments

* BICCN
* Allen Institute of Brain Science
* Chan Zuckerberg Initiative 
* California Institute for Regenerative Medicine 

