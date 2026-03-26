Installation
============

To install NS-Forest via Github: 

.. code-block:: console

   git clone https://github.com/JCVenterInstitute/NSForest.git
   cd NSForest
   
Now we create the environment with all the dependencies required to run NS-Forest, including the ability to run a jupyter lab notebook version of the tutorial. This was successfully run on a Macbook pro, Xtools installed, Apple M3 pro chip, Sonoma 14.6 operating system.

Additionally, emacs is installed to faciliate editing. Once the environment is created then we activate the environment. 

.. code-block:: console

   conda env create -f environment.yml
   conda activate environment

There are 2 ways to run NS-Forest: Jupyter Notebook and command line

Jupyter Notebook: docs/source/tutorial_nsforesting.ipynb

.. code-block:: console
   
   import sys
   import os
   code_folder = "location/of/NSForest/folder"
   sys.path.insert(0, os.path.abspath(code_folder))
   import nsforest as ns
   from nsforest import utils

Command Line: nsforest/__main__.py.

.. code-block:: console

   pip install .
   python3 nsforest

