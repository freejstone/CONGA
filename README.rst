"""""""""""""""
Open Group-walk
"""""""""""""""
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
A tool for discovering peptides with unaccounted for PTMs
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Open Group-walk is a tool for discovering peptides in mass spectrometry data with rigourous FDR control. It is designed to allow for unexpected post-translational modifications as well as chimeric spectra. Open Group-walk takes three inputs: (1) the set of top-scoring PSMs from a traditional, narrow-window search against a concatenated target-decoy database, (2) the set of top-k PSMs (or fewer if less than k PSMs exist) from an open search, again against a concatenated target-decoy database, and (3) the pairs of target and decoy peptide sequences from the database. Given this, Open Group-walk will return a discovery list of peptides.

Documentation
=============

You can find documentation on how to use Open Group-walk over on `readtheddocs <https://open-groupwalk.readthedocs.io/en/latest/>`_. Alternatively you can find the same documentation under docs/pages in this respository.

Paper
=====

A link to the biorXiv paper will eventually go here.

Installation
============

To install simply go to Code, Download Zip and extract the Zip file. Installation of the relevant package dependencies are required (see Groupwalk.py). Please see the documentation, specifically the tutorial, on how to run Open Group-walk.
