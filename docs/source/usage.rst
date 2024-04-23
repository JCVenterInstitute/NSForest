Usage
=====

.. _installation:

Installation
------------

To install NSForest via Github: 

.. code-block:: console

   $ git clone https://github.com/JCVenterInstitute/NSForest.git
   $ cd NSForest
   $ conda env create -f nsforest.yml
   $ conda activate nsforest
   $ cd tutorials

To install NSForest with pip: 

.. code-block:: console

   $ pip install nsforest

.. _running:

Preprocessing
-------------

.. autofunction:: preprocessing.get_medians

.. autofunction:: preprocessing.prep_medians

.. autofunction:: preprocessing.prep_binary_scores

.. autofunction:: preprocessing.dendrogram

NSForest
--------

.. autofunction:: nsforesting.NSForest

Evaluating
----------

.. autofunction:: evaluating.DecisionTree

.. autofunction:: evaluating.add_fraction

Plotting
--------

.. autofunction:: plotting.dotplot

.. autofunction:: plotting.stackedviolin

.. autofunction:: plotting.matrixplot

.. autofunction:: plotting.boxplot_fscore

.. autofunction:: plotting.boxplot_ppv

.. autofunction:: plotting.boxplot_ontarget

.. autofunction:: plotting.scatter_w_clusterSize_fscore

.. autofunction:: plotting.scatter_w_clusterSize_ppv

.. autofunction:: plotting.scatter_w_clusterSize_ontarget
