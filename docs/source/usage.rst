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

.. _running:

Preprocessing
-------------

.. autofunction:: preprocessing.get_medians

Running NSForest

.. code-block:: console

   (nsforest) $ python3 nsforest -a ${input_folder}/arguments_${prefix}.csv 

Parallelizing NSForest

.. code-block:: console

   (nsforest) $ python3 nsforest -a ${input_folder}/arguments_${prefix}.csv -c "${cluster}"
