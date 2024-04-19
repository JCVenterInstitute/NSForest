Preprocessing
=============

.. _preprocessing:

To install NSForest via Github: 

.. code-block:: console

  $ cluster_header = "cluster"
  $ medians_header = "medians_" + cluster_header
  $ binary_scores_header = "binary_scores_" + cluster_header

.. code-block:: console

  $ ns.pp.dendrogram(adata, cluster_header, save = "_" + cluster_header + ".png") 

.. code-block:: console

  $ ns.pp.prep_medians(adata, cluster_header, use_mean = False, positive_genes_only = True) 

.. code-block:: console

  $ ns.pp.prep_binary_scores(adata, cluster_header, medians_header)

.. code-block:: console

  $ adata
