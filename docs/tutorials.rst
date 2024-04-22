Tutorials
=========

.. _tutorials:

To run NSForest : 

.. code-block:: console

  $ # Adding this line to actually subset to positive genes if running NSForest
  $ adata = ns.pp.prep_medians(adata, cluster_header, use_mean = True, positive_genes_only = True) 
  $ results = nsforesting.NSForest(adata, cluster_header, medians_header, binary_scores_header, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)
