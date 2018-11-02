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
cluster versus all. The top ten ranked features from the Random Forest are then tested using f-measure as an objective function. 
For example, during the first step all top ten features are independently evaluated for their discriminatory power at an 
expression value where 75% of the cells have greater than or equal expression. Given that 25% of the cells are lost de facto,
the maximum f-measure for the first step is estimated to be around 0.87 (there will be cases where its higher or lower, such 
as having equal expression across all cells). After the best f-measure is found classifying with one gene than the remaining 
nine genes are tested in combination with the top first gene, again using an expression value where 75% of the cells have expression. 
After the best pair of genes is found, the remaining 8 genes are tested in third position, and onward until the analysis reaches
6 combinations.   

See code for detailed comments. 


## Versioning

The initial release is version 1.3. Version 1.0 was described in: 

Aevermann BD, Novotny M, Bakken T, Miller JA, Diehl AD, Osumi-Sutherland D, Lasken RS, Lein ES, Scheuermann RH.
Cell type discovery using single-cell transcriptomics: implications for ontological representation. 
Hum Mol Genet. 2018 May 1;27(R1):R40-R47. doi: 10.1093/hmg/ddy100.


## Authors

* Brian Aevermann baeverma@jcvi.org and Richard Scheuermann RScheuermann@jcvi.org


## License

This project is licensed under the MIT License - see the LICENSE.txt file for details

## Acknowledgments

* Allen Institute of Brain Science
* Chan Zuckerberg Initiative 
* California Institute for Regenerative Medicine 

