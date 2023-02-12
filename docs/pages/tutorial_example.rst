""""""""""""""""
Tutorial example
""""""""""""""""

Once CONGA has been installed in your conda-environment, CONGA can be called as a module using ``python3 -m CONGA``. On this page we walk through a simple example of running conga for the three different search programs, Tide, Comet and MSFragger. 



.. toctree::
   :maxdepth: 2
   :caption: Contents:

=========
Crux-Tide
=========

You can download the search files to test the commands below: `Tide-search files <https://github.com/freejstone/open_groupwalk/tree/main/data/tide>`_.

To run CONGA using Tide results with default parameters, use the following command in the terminal.

``python3 -m CONGA path/to/narrow_example1.tide-search.txt path/to/open_example1.tide-search.txt``

CONGA allows for a variety of delta-masses (possibly unaccounted PTMs) and variable modifications to be reported for the same discovered sequence. It does this by clustering together PSMs by considering precursors with similar masses (in m/z units) that are matched to the same discovered "stem" sequence. It then takes the maximal scoring PSM in each cluster and for each modified variant of the "stem" sequence. To cluster the PSMs accurately, it is beneficial to supply a pair of values to ``--isolation-window``, which should be the lower and upper isolation window offsets in m/z units. By default, CONGA uses ``--isolation-window 2,2``. E.g. 

``python3 -m CONGA --isolation-window 0.75,0.75 path/to/narrow_example1.tide-search.txt path/to/open_example1.tide-search.txt``

CONGA's target-decoy competition is conducted over (unmodified) "stem" sequences. In other words, it reports the highest scoring PSM among all PSMs with the same stem sequence. If you want CONGA to report peptides that distinguish between sequences with different modifications, you can turn this feature off by setting ``--account_mods F``. Note that the assumptions of theoretical FDR control do not strictly hold in this case, because target peptides with variable modifications may be scored highly simply because they share similar b- and y-ions in common with the same target peptide but with possibly a different ensemble of variable modifications. In any case, if you wish to proceed, you will need to provide the tab-delimited target-decoy peptide pairs produced by tide-index.

``python3 -m CONGA --account_mods F --isolation-window 0.75,0.75 path/to/narrow_example1.tide-search.txt path/to/open_example1.tide-search.txt path/to/tide-index_example1.peptides.txt``

The default score uses Tailor scores. If you prefer to use Xcorr Scores, you can change the score option ``--score xcorr_score``:

``python3 -m CONGA --score xcorr_score --isolation-window 0.75,0.75 path/to/narrow_example1.tide-search.txt path/to/open_example1.tide-search.txt``

The software will automatically detect if the PSMs in the search files have been searched against a concatenated target-decoy database or separately against each database. If the latter, then you will need to have the name *target* appear in the search file names. In that case, the software will search for the respective decoy-search files by replacing *target* with *decoy* in the search file names in the same directory.

``python3 -m CONGA path/to/narrow_example2.tide-search.target.txt path/to/open_example2.tide-search.target.txt``

=====
Comet
=====

You can download the search files to test the commands below: `Comet-search files <https://github.com/freejstone/open_groupwalk/tree/main/data/comet>`_.

You can use either standalone Comet or the implementation of Comet in Crux. The CONGA software will automatically standardize the difference in terminology between the standalone version and the one in Crux. To run CONGA using Comet-search results with default parameters, use the following command in the terminal.

``python3 -m CONGA --score e-value path/to/narrow_example1.comet.txt path/to/open_example1.comet.txt``

There is no need to provide a list of target-decoy pairs, because Comet always reverses the target peptides to yield decoy peptides. Comet can use either Xcorr scores or E-values. You will need to specify the score, because the default (Tailor score) is specific to Tide.

All other discussion related to CONGA's handle on varaible modifications in :ref:`Crux-Tide` equally applies to Comet search files.

=========
MSFragger
=========

You can download the search files to test the commands below: `MSFragger-search files <https://github.com/freejstone/open_groupwalk/tree/main/data/MS>`_.

Using MSFragger search files is more challenging. MSFragger nor any of the related softwares in Fragpipe produce decoy peptides that *pair* with the target peptides, which is essential for CONGA. Unfortunately not even reversing the entire protein sequence and digesting the protein to produce decoy peptides yields the same result as reversing the target peptides. Instead the user must always provide a tab-delimited target-decoy pairing in the same format as what Tide-index produces. You can do this by creating your own script, with two columns one labelled *target* and the other *decoy*, or more simply by implementing Tide-index. The latter is reasonably achievable by matching the parameters in Tide-index with MSFragger. More specifically we have the following equivalencies between the two softwares.

.. list-table:: Parameters
   :align: center
   :widths: 50 50
   :header-rows: 1

   * - Tide-index
     - MSFragger
   * - --enzyme
     - search_enzyme_name_1
   * - --missed_cleavages
     - allowed_missed_cleavage_1
   * - --max_mods
     - max_variable_mods_combinations
   * - --max_length
     - digest_max_length
   * - --min_length
     - digest_min_length
   * - --min_mass, --max_mass
     - digest_mass_range

To run CONGA using MSFragger-search results with default parameters, use the following command:

``python3 -m CONGA --score hyperscore path/to/narrow.tsv path/to/open.tsv path/to/tide-index.peptides.txt``

Like Comet, the score needs to be specified, because the default is Tailor score. All other discussion related to CONGA's handle on varaible modifications :ref:`Crux-Tide` equally applies to MSFragger-search files.

=========================
Interpreting the log file
=========================

Here we explain what a typical log file looks like from CONGA::

  INFO: CPU: macOS-10.16-x86_64-i386-64bit (1)
  INFO: Version: 1.0.1 (2)
  INFO: 2023-02-09 16:06:16.072313 (3)
  INFO: Command used: /Volumes/Expansion/Documents_from_local_computer/CONGA/CONGA/__main__.py --overwrite T --competition_window 4 data/tide/narrow_example1.tide-search.txt data/tide/open_example1.tide-search.txt (4)
  INFO: Successfully read in arguments
  INFO: Reading in search files.
  INFO: Successfully read in search files.
  INFO: tailor_score successfully found in search files. (5)

  INFO: Tide search files detected. (6)
  INFO: Concatenated search files detected. (7)
  INFO: Filtering for neighbours. (8)
  INFO: Doing dynamic level competition. (9)
  INFO: Constructing groups adaptively. (10)
  INFO:                                    decoys  targets     ratio (11)
  group names:                                                
  narrow                               1706    19495  0.087510
  top 1 PSMs & top (4, 6] mass bins      69      117  0.589744
  top 1 PSMs & top 1 mass bin            16      714  0.022409
  top 1 PSMs & top 3 mass bin            18      176  0.102273
  top 2 or more PSMs                    114      590  0.193220
  left over group                      2226     2497  0.891470
  INFO: Applying group walk. (12)
  INFO: Group walk complete.
  INFO: 19196 peptides discovered at the 1% FDR level. (13)
  INFO: 20076 peptides discovered at the 5% FDR level.
  INFO: Scan multiplicities among the discovered peptides at 1% FDR level: (14)
  INFO:                     Count
  Scan multiplicity:       
  1                   18426
  2                     385
  INFO: Scan multiplicities among the discovered peptides at 5% FDR level:
  INFO:                     Count
  Scan multiplicity:       
  1                   19284
  2                     396
  INFO: Writing peptides at user-specified FDR level to directory.
  INFO: Elapsed time: 49.64 s



#. CPU information.
#. Version number.
#. Date and time.
#. Command used.
#. CONGA checks whether the score specified is found in the search files provided.
#. CONGA detects what type of search engine was used based off the column names of the search files provided.
#. CONGA detects whether a separate search or a concantenated search was conducted based off whether it is able to find any decoys in the search files.
#. CONGA notifies the user that it is about to initiate the first phase of filtering for neighbouring peptides.
#. CONGA notifies the user that it is about to initiate the second phase of undergoing dynamic competition.
#. CONGA notifies the user that it is about to initiate the third phase of partitioning the peptides into groups.
#. CONGA reports the names of the groups it constructed, the number of decoys and targets in each group, and the ratio of decoys to targets in each group.
#. CONGA notifies the user that it is about to initiate the fourth phase of applying the group-walk algorithm.
#. The number of peptides discovered at the 1% FDR level and at the 5% FDR level.
#. A count of the number of times each scan is responsbile for discovering peptide at the 1% and 5% FDR level. As an example, there are 18426 scans each responsible for discovering 1 peptide, and 385 scans responsible for discovering 2 peptides at the 1% level. Hence the total number of discoveries are 18426 + 2*385 = 19196


