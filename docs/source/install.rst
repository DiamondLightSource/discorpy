============
Installation
============

Installing from source
----------------------

- Download the source from `github <https://github.com/DiamondLightSource/discorpy>`_ (click-> Code -> Download ZIP).
  Unzip to a local folder.
- Install `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_.
- Open command prompt, navigate to the source folder, run the following
  commands:

  .. code-block:: console

     conda create -n discorpy
     conda activate discorpy
     conda install python
     python setup.py install

Using conda
-----------

- Install Miniconda as instructed above.
- Open terminal or command prompt and run the following commands:

    + If install to an existing environment:

      .. code-block:: console

         conda install -c conda-forge discorpy

    + If install to a new environment:

      .. code-block:: console

         conda create -n discorpy
         conda activate discorpy
         conda install python
         conda install -c conda-forge discorpy

Using pip
---------

- Install Miniconda as instructed above.
- Open terminal or command prompt and run the following commands:
    + If install to an existing environment:

      .. code-block:: console

         pip install discorpy

    + If install to a new environment:

      .. code-block:: console

         conda create -n discorpy
         conda activate discorpy
         conda install python
         pip install discorpy
