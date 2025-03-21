Installation
============

To install NSForest via Github: 

.. code-block:: console

   $ git clone https://github.com/JCVenterInstitute/NSForest.git
   $ cd NSForest
   
Now we create the environment with all the dependencies required to run NSForest, including the ability to run a jupyter lab notebook version of the tutorial. This was successfully run on a Macbook pro, Xtools installed, Apple M3 pro chip, Sonoma 14.6 operating system.

Additionally, emacs is installed, to faciliate editing.   Feel free to inspect the yaml file, nsforest.yml.   Once the environment is created then we activate the environment. 

.. code-block:: console
   $ conda env create -f nsforest.yml
   $ conda activate nsforest
   $ cd tutorials
Once the environment is activated.  We need to install the NSForest package.

.. code-block:: console
   $ pip install .

You should be able to "import nsforest" and run the tutorial.ipynb. 
