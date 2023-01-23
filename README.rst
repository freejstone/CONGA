"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
CONGA: Combined Open and Narrow searches via Group Analysis
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
A tool for discovering peptides with unaccounted for PTMs
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

CONGA is a tool for discovering peptides in mass spectrometry data with rigourous FDR control. It is designed to allow for unexpected post-translational modifications as well as chimeric spectra. Open Group-walk takes three inputs: (1) the set of top-scoring PSMs from a traditional, narrow-window search against a concatenated target-decoy database, (2) the set of top-k PSMs (or fewer if less than k PSMs exist) from an open search, again against a concatenated target-decoy database, and (3) the pairs of target and decoy peptide sequences from the database. Given this, Open Group-walk will return a discovery list of peptides.

Documentation
=============

You can find documentation on how to use CONGA over on `readtheddocs <https://open-groupwalk.readthedocs.io/en/latest/>`_. Alternatively you can find the same documentation under docs/pages in this respository.

Paper
=====

A link to the biorXiv paper will eventually go here.

Installation
============

To install, first create a virtual environment using conda:

``conda create --name conga_env python=3.9``

Then activate this virtual environment:

``conda activate conga_env``

Next download the latest release using `pip`:

``pip install CONGA``

Alternatively you can download the latest release from Github, and install using pip in the same directory as setup.py using:

``pip install .``

Please see the documentation, specifically the tutorial, on how to run CONGA.

Releasing
---------

Releases are published automatically when a tag is pushed to GitHub.

.. code-block:: bash

   # Set next version number
   export RELEASE=x.x.x

   # Create tags
   git commit --allow-empty -m "Release $RELEASE"
   git tag -a $RELEASE -m "Version $RELEASE"

   # Push
   git push upstream --tags