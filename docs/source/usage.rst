Usage
=====

.. _installation:

Installation
------------

To install NSForest with pip: 

.. code-block:: console

   $ pip install nsforest

To install NSForest via Github: 

.. code-block:: console

   $ git clone https://github.com/JCVenterInstitute/NSForest.git
   $ cd NSForest
   $ conda env create -f nsforest.yml
   $ conda activate nsforest
   $ cd tutorials

.. _running:

Preprocessing
-------------

.. autofunction:: preprocessing.get_medians

NSForest
--------

.. autofunction:: nsforesting.NSForest

Evaluating
----------

.. autofunction:: evaluating.DecisionTree

Plotting
--------

.. autofunction:: plotting.dotplot
