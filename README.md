<img src="NS-Forest-sticker.png" width="110" height="125">

# NS-Forest v4.1

Documentation: https://nsforest.readthedocs.io/en/latest/

BMC Methods Link: https://bmcmethods.biomedcentral.com/articles/10.1186/s44330-024-00015-2

## Download and installation

In terminal: 
```
git clone https://github.com/JCVenterInstitute/NSForest.git

cd NSForest

conda env create -f nsforest.yml

conda activate nsforest

pip install .

```
## Tutorial

Follow the on readthedocs: https://nsforest.readthedocs.io/en/latest/tutorial.html

## Pipeline

<img src="pipeline.PNG">

NS-Forest is an algorithm designed to identify minimum combinations of necessary and sufficient marker genes for a cell type cluster identified in a single cell or single nucleus RNA sequencing experiment that optimizes classification accuracy. NS-Forest proceeds through the following steps (default setting):

1. Data input: An AnnData object (e.g., .h5ad file) with cell type cluster labels. 

2. Binary score calculation: Each gene is assigned a binary score for every cluster. Binary score is a measurement of the binary expression pattern of a gene. A higher binary score means a gene is expressed in one cluster and not others. A lower binary score means a gene is expressed in many clusters and would not be an ideal candidate for a cell type-specific marker gene. 

3. Binary scoring criterion: NS-Forest then filters for genes with high binary scores. Candidate genes are selected if their binary scores are 2 standard deviations above the mean of all genes expressed in the cluster. 

4. Random forest: The top 15 binary score genes are used as input into a random forest classifier, which ranks the genes by Gini Impurity, while producing a classification model for each cluster. 

5. Decision tree evaluation: The top 6 ranked random forest genes are used as input into decision trees where all combinations of input genes are evaluated and the combination with the highest F-beta score is selected. 

6. Output: The NS-Forest algorithm outputs 1-6 marker genes per cluster along with the classification metrics (F-beta, PPV (precision), recall) and the On-Target Fraction expression metric. 

### NS-Forest Marker Gene Evaluation

The final module in the NS-Forest algorithm can also be used to assess the performance of any collection of marker gene combinations identified using any approach.  The marker gene evaluation module includes the following steps (default setting):

1. Data input: 1) An AnnData object (e.g., .h5ad file) with cell type cluster labels. 2) A list of marker genes for every cluster to be evaluated. 

2. Decision tree creation: One-vs-all decision trees are created for each gene in the cluster combination and evaluated for classification accuracy. 

3. Decision tree evaluation: Each gene in the cluster combination is evaluated using these decision trees to determine if the gene gives the correct classification. If even one gene in the cluster combination gives a misclassification, then the prediction is considered incorrect. Note: This strict criteria may lead to PPV = 0 when no true positives (TP) classification are obtained. 

4. Output: The NS-Forest marker gene evaluation outputs the classification metrics (F-beta, PPV (precision), recall) and On-Target Fraction for every cluster combination, which can be used to compare against other marker gene lists.


## Prerequisites
* This is a python script written and tested in python 3.11, scanpy 1.9.6.
* Other required libraries: numpy, pandas, sklearn, plotly, time, tqdm.

## Versions and citations

Earlier versions are managed in [Releases](https://github.com/JCVenterInstitute/NSForest/releases).  

Version 4.0:

Liu A, Peng B, Pankajam A, Duong TE, Pryhuber G, Scheuermann RH, Zhang Y. Discovery of optimal cell type classification marker genes from single cell RNA sequencing data. BMC Methods 1, 15 (2024). https://doi.org/10.1186/s44330-024-00015-2

Version 2:

Aevermann BD, Zhang Y, Novotny M, Keshk M, Bakken TE, Miller JA, Hodge RD, Lelieveldt B, Lein ES, Scheuermann RH. A machine learning method for the discovery of minimum marker gene combinations for cell-type identification from single-cell RNA sequencing. Genome Res. 2021 Jun 4:gr.275569.121. doi: 10.1101/gr.275569.121.

Version 1.3/1.0:

Aevermann BD, Novotny M, Bakken T, Miller JA, Diehl AD, Osumi-Sutherland D, Lasken RS, Lein ES, Scheuermann RH. Cell type discovery using single-cell transcriptomics: implications for ontological representation. Hum Mol Genet. 2018 May 1;27(R1):R40-R47. doi: 10.1093/hmg/ddy100.

## Authors

* Yun (Renee) Zhang zhangy@jcvi.org
* Beverly Peng bpeng@jcvi.org
* Angela Liu aliu@jcvi.org
* Richard Scheuermann richard.scheuermann@nih.gov
* Brian Aevermann baevermann@chanzuckerberg.com

## License

This project is licensed under the [MIT License](https://github.com/JCVenterInstitute/NSForest/blob/master/LICENSE).

## Acknowledgments

* Allen Institute of Brain Science
* Brain Initiative Cell Census Network
* Chan Zuckerberg Initiative
* California Institute for Regenerative Medicine
