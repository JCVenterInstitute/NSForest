# Necessary and Sufficient Forest (NS-Forest) for Cell Type Marker Determination from cell type clusters

## Getting Started

Install Jupyter notebook and python 2.7

### Prerequisites

* The script is a Jupyter notebook in python 2.7. Required libraries: Numpy, Pandas, Sklearn, graphviz, numexpr
* The input data is a tab delimited expression Cell x Gene matrix with one column containing the cluster assignments 
* The cluster-label column must be named "Clusters" and the labels must be non-numeric (if currently numbers, please add "Cl" or any text would work). 
* The gene identifiers used must avoid special characters such as ./-/@ or beginning with numbers (I prefix identifiers beginning with numbers and substitute all special characters with "_")


### Description

Necessary and Sufficient Forest is a method that takes cluster results from single cell/nuclei RNAseq experiments 
and generates lists of minimal markers needed to define each “cell type cluster”. 
 
The method begins by re-encoding the cluster labels into binary classifications, and Random Forest models are generated comparing each 
cluster versus all. The top fifteen genes are then reranked using a score measuring how binary they are, e.g., a gene with expression in
the target cluster but no expression in the other clusters would have a high binary score. Expression cutoffs for the top six genes ranked
by binary score are then determined by generating individual decision trees and extracting the decision path information. Then all permutations 
of the top six most binary genes are evaluated using f-beta score as an objective function (the beta value default set at 0.5, which weights the 
f-measure score more toward precision as opposed to recall). 



See code for detailed comments. 


## Versioning

This is version 2.0 The initial release was version 1.3. Version 1.0 was described in: 

Aevermann BD, Novotny M, Bakken T, Miller JA, Diehl AD, Osumi-Sutherland D, Lasken RS, Lein ES, Scheuermann RH.
Cell type discovery using single-cell transcriptomics: implications for ontological representation. 
Hum Mol Genet. 2018 May 1;27(R1):R40-R47. doi: 10.1093/hmg/ddy100.

Aevermann BD, Zhang Y, Novotny M, Bakken TE, Miller JA, Hodge R, Lelieveldt B, Lein ES, Scheuermann RH. “NS-Forest: A machine learning method for the objective identification of minimum marker gene combinations for cell type determination from single cell RNA sequencing” bioRxiv 2020.09.23.308932; doi: https://doi.org/10.1101/2020.09.23.308932


## Authors

* Brian Aevermann baeverma@jcvi.org and Richard Scheuermann RScheuermann@jcvi.org


## License

This project is licensed under the MIT License - see the https://opensource.org/licenses/MIT for details

## Acknowledgments

* Allen Institute of Brain Science
* Chan Zuckerberg Initiative 
* California Institute for Regenerative Medicine 

