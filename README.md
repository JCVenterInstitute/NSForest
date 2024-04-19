# NS-Forest
<img src="NS-Forest-sticker.png" width="110" height="125">
<img src="pipeline.PNG" width="500" height="200">

# NS-Forest v4.0
## Download and installation

### In terminal: 

git clone https://github.com/JCVenterInstitute/NSForest.git

cd NSForest

conda env create -f nsforest.yml

conda activate nsforest

cd tutorials

### How to run: 

Follow the tutorials/tutorial.ipynb. 


## Dev version of NS-Forest v4.0

Follow the [tutorial](https://jcventerinstitute.github.io/celligrate/tutorials/NS-Forest_tutorial.html) to get started.

Download 'NSForest_v4dot0_dev.py' and replace the version in the tutorial. Sample code below.

```
adata_median = preprocessing_medians(adata, cluster_header)
adata_median.varm["medians_" + cluster_header].stack().plot.hist(bins=30, title = 'cluster medians')

adata_median_binary = preprocessing_binary(adata_median, cluster_header, "medians_" + cluster_header)
adata_median_binary.varm["binary_scores_" + cluster_header].stack().plot.hist(bins=30, title='binary scores')

## make a copy of prepared adata
adata_prep = adata_median_binary.copy()

NSForest(adata_prep, cluster_header=cluster_header, n_trees=1000, n_genes_eval=6,
          medians_header = "medians_" + cluster_header, binary_scores_header = "binary_scores_" + cluster_header,
          gene_selection = "BinaryFirst_high", outputfilename="BinaryFirst_high")
```
# NS-Forest v3.9
## Download and installation

NS-Forest can be installed using `pip`:
`sudo pip install nsforest`

If you are using a machine on which you lack administrative access, NS-Forest can be installed locally using `pip`:
`pip install --user nsforest`

NS-Forest can also be installed using `conda`:
`conda install -c ttl074 nsforest`

Will be uploaded to official conda channel soon.

**Prerequisites:**
* This is a python script written and tested in python 3.8, scanpy 1.8.2, anndata 0.8.0.
* Other required libraries: numpy, pandas, sklearn, itertools, time, tqdm.

## Tutorial

Follow the [tutorial](https://jcventerinstitute.github.io/celligrate/tutorials/NS-Forest_tutorial.html) to get started.

If you download 'NSForest_v3dot9_2.py' directly, replace the version to the most updated one in the tutorial.

If you download the `pip` or `conda` package, use the following in the tutorial.

```
import nsforest as ns
ns.NSForest()
```

## Versions and citations

Earlier versions are managed in [Releases](https://github.com/JCVenterInstitute/NSForest/releases).  

Version 2 and beyond:

Aevermann BD, Zhang Y, Novotny M, Keshk M, Bakken TE, Miller JA, Hodge RD, Lelieveldt B, Lein ES, Scheuermann RH. A machine learning method for the discovery of minimum marker gene combinations for cell-type identification from single-cell RNA sequencing. Genome Res. 2021 Jun 4:gr.275569.121. doi: 10.1101/gr.275569.121.

Version 1.3/1.0:

Aevermann BD, Novotny M, Bakken T, Miller JA, Diehl AD, Osumi-Sutherland D, Lasken RS, Lein ES, Scheuermann RH. Cell type discovery using single-cell transcriptomics: implications for ontological representation. Hum Mol Genet. 2018 May 1;27(R1):R40-R47. doi: 10.1093/hmg/ddy100.

## Authors

* Yun (Renee) Zhang zhangy@jcvi.org
* Richard Scheuermann RScheuermann@jcvi.org
* Brian Aevermann baevermann@chanzuckerberg.com
* Angela Liu aliu@jcvi.org
* Beverly Peng bpeng@jcvi.org

## License

This project is licensed under the [MIT License](https://github.com/JCVenterInstitute/NSForest/blob/master/LICENSE).

## Acknowledgments

* BICCN
* Allen Institute of Brain Science
* Chan Zuckerberg Initiative
* California Institute for Regenerative Medicine
