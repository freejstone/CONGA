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

You can download the search files to test the commands below: `Tide-search files <https://github.com/freejstone/open_groupwalk/tree/main/data/tide>`_. These files also exist under the CONGA package in the latest release (``CONGA/data/tide``).

  To run CONGA using Tide results with default parameters, use the following command in the terminal.

``python3 -m CONGA path/to/narrow_example1.tide-search.txt path/to/open_example1.tide-search.txt``

CONGA collapses peptides that are equal up to variable modifications by selecting the best-scoring peptide. If you want to run CONGA without this collapsing, and thus allowing the discovery of peptides that differ by some variable modification, you can turn this collapsing off by setting ``--account_mods F``. Note that the assumptions of theoretical FDR control do not strictly hold in this case, because target peptides with variable modifications may be scored highly simply because they share similar b- and y-ions in common with the same target peptide but with possibly a different set of variable modifications. In any case, if you wish to run without this collapsing, you will need to provide the target-decoy peptide pairs produced by tide-index.

``python3 -m CONGA --account_mods F path/to/narrow_example1.tide-search.txt path/to/open_example1.tide-search.txt path/to/tide-index_example1.peptides.txt``

Alternatively, it might be of interest to report all other top-1 matched peptides that differ from a genuinely discovered peptide by some variable modification(s). One can achieve this by setting ``--return_extra_mods T``:

``python3 -m CONGA --return_extra_mods T path/to/narrow_example1.tide-search.txt path/to/open_example1.tide-search.txt``

The default score uses Tailor scores. If you prefer to use Xcorr Scores, you can change the score option ``--score xcorr_score``:

``python3 -m CONGA --score xcorr_score path/to/narrow_example1.tide-search.txt path/to/open_example1.tide-search.txt``

The software will automatically detect if the PSMs in the search files have been searched against a concatenated target-decoy database or separately against each database. If the latter, then you will need to have the name *target* appear in the search file names. In that case, the software will search for the respective decoy-search files by replacing *target* with *decoy* in the search file names in the same directory.

``python3 -m CONGA path/to/narrow_example2.tide-search.target.txt path/to/open_example2.tide-search.target.txt``

=====
Comet
=====

You can download the search files to test the commands below: `Comet-search files <https://github.com/freejstone/open_groupwalk/tree/main/data/comet>`_. These files also exist under the CONGA package in the latest release (``CONGA/data/comet``).

You can use either standalone Comet or the implementation of Comet in Crux. The CONGA software will automatically standardize the difference in terminology between the standalone version and the one in Crux. To run CONGA using Comet-search results with default parameters, use the following command in the terminal.

``python3 -m CONGA --score e-value path/to/narrow_example1.comet.txt path/to/open_example1.comet.txt``

There is no need to provide a list of target-decoy pairs, because Comet always reverses the target peptides to yield decoy peptides. Comet can use either Xcorr scores or E-values. You will need to specify the score, because the default (Tailor score) is specific to Tide.

All other discussion related to CONGA's handle on varaible modifications in :ref:`Crux-Tide` equally applies to Comet search files.

=========
MSFragger
=========

You can download the search files to test the commands below: `MSFragger-search files <https://github.com/freejstone/open_groupwalk/tree/main/data/MS>`_. These files also exist under the CONGA package in the latest release (``CONGA/data/MS``).

Using MSFragger search files is more challenging. MSFragger nor any of the related softwares in Fragpipe produce decoy peptides that *pair* with the target peptides, which is essential for CONGA. Unfortunately not even reversing the entire protein sequence and digesting the protein to produce decoy peptides yields the same result as reversing the target peptides. Instead the user must always provide a target-decoy pairing in the same format as what Tide-index produces. You can do this by creating your own script, with two columns one labelled *target* and the other *decoy*, or more simply by implementing Tide-index. The latter is reasonably achievable by matching the parameters in Tide-index with MSFragger. More specifically we have the following equivalencies between the two softwares.

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

The default precision of the mass-modifications in Tide-index is to two decimal places, while MSFragger uses much higher precision. For this reason, CONGA will round the variable modifications in the MSFragger search files to two decimal places, allowing us to pair the peptides in MSFragger according to the target-decoy pairs produced by Tide-index. To run CONGA using MSFragger-search results with default parameters, use the following command:

``python3 -m CONGA --score hyperscore path/to/narrow.tsv path/to/open.tsv path/to/tide-index.peptides.txt``

Like Comet, the score needs to be specified, because the default is Tailor score. All other discussion related to CONGA's handle on varaible modifications :ref:`Crux-Tide` equally applies to MSFragger-search files.

=========================
Interpreting the log file
=========================

Here we explain what a typical log file looks like from CONGA::

  INFO: CPU: macOS-10.16-x86_64-i386-64bit (1)
  INFO: 2022-12-15 13:00:40.486262 (2)
  INFO: Command used: ../../CONGA.py narrow_example1.tide-search.txt open_example1.tide-search.txt (3)
  INFO: Successfully read in arguments
  INFO: Reading in search files.
  INFO: Successfully read in search files.
  INFO: tailor_score successfully found in search files.  (4)

  INFO: Tide search files detected. (5)
  INFO: Concatenated search files detected. (6)
  INFO: Filtering for neighbours. (7)
  INFO: Doing head to head competition. (8)
  INFO: Constructing groups adaptively. (9)
  INFO:                                    decoys  targets     ratio (10)
  group names:                                                
  narrow                               1161    15646  0.074204
  top 1 PSMs & top (4, 7] mass bins      41       68  0.602941
  top 1 PSMs & top 1 mass bin            27      333  0.081081
  top 2 or more PSMs                     61      430  0.141860
  left over group                      1506     1545  0.974757
  INFO: Applying group walk. (11)
  INFO: Group walk complete.
  INFO: 15183 peptides discovered at the 1% FDR level. (12)
  INFO: 15768 peptides discovered at the 5% FDR level.
  INFO: Scan multiplicities among the discovered peptides at 1% FDR level: (13)
  INFO:                     Count
  Scan multiplicity:       
  1                   14611
  2                     286
  INFO: Scan multiplicities among the discovered peptides at 5% FDR level:
  INFO:                     Count
  Scan multiplicity:       
  1                   15182
  2                     293
  INFO: Writing peptides at user-specified FDR level to directory. 
  INFO: Elapsed time: 23.28 s

#. CPU information.
#. Date and time.
#. Command used.
#. CONGA checks whether the score specified is found in the search files provided.
#. CONGA detects what type of search engine was used based off the column names of the search files provided.
#. CONGA detects whether a separate search or a concantenated search was conducted based off whether it is able to find any decoys in the search files.
#. CONGA notifies the user that it is about to initiate the first phase of filtering for neighbouring peptides.
#. CONGA notifies the user that it is about to initiate the second phase of undergoing peptide-level competition.
#. CONGA notifies the user that it is about to initiate the third phase of partitioning the peptides into groups.
#. CONGA reports the names of the groups it constructed, the number of decoys and targets in each group, and the ratio of decoys to targets in each group.
#. CONGA notifies the user that it is about to initiate the fourth phase of applying the group-walk algorithm.
#. The number of peptides discovered at the 1% FDR level and at the 5% FDR level.
#. A count of the number of times each scan is responsbile for discovering peptide at the 1% and 5% FDR level. As an example, there are 14611 scans each responsible for discovering 1 peptide, and 286 scans responsible for discovering 2 peptides at the 1% level. Hence the total number of discoveries are 14611 + 2*286 = 15183


