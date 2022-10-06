"""""""""""
Parameters
"""""""""""

We provide the following descriptions of the input parameters used by Groupwalk.py.

-------------------
Neighbor filtration
-------------------

* ``-- tops_open <integer>``: The number of top PSMs in the open search file used in the neighbour-filtering process. Default = 5.
* ``--neighbour_remove <T|F>``: If true, for each scan, we successively move down the list of PSMs associated with each scan, ordered in terms of the score from highest to lowest, and with the current peptide is thrown out, and we proceed to the next. Default = T.
* ``--thresh <float>``: The similarity score used as a threshold to filter out neighbouring peptides. Default = 0.25.
* ``--return_filt_search <T|F>``:  Whether or not to return filtered narrow and open search files. Default = F.

-------------------------
Peptide-level competition
-------------------------

* ``--score <string>``: Either 'tailor_score', 'xcorr_xcore', 'e-value' or 'hyperscore'. The score that will be used in the peptide-level competition and subsequent Group construction and Group-walk algorithm. If 'tailor_score', it is assumed the search files are  derived from Tide. If 'xcorr_score', either Tide search or Comet is assumed to be used. If 'e-value', Comet is assumed. If 'hyperscore' it is assumed the search files are derived from MS-Fragger. Default = 'tailor_score'.
* ``--account_mods <T|F>``: To determine whether the group-walk uses the best PSM among the equivalent classes of peptides which are equal up to variable modification, or not. Default = T.

------------------
Group construction
------------------

* ``--precursor_bin_width <value>``: To determine the size of the bins used to create discretize the mass- differences between the sample and theoretical spectra. Default = 1.0005079/4.
* ``--adaptive <T|F>``: To determine whether groups should be chosen using the Kolmogorov-Smirnov test or if a fixed procedure should be used. Default = T.
* ``--min_group_size <integer>``: The number of multiples of K that the used to determine the minimum size of each group. See option ``--K``. Default = 2.
* ``--tops_gw <integer>``: The number of top PSMs for each scan in the open search that will be used by group-walk. Default = 2.
* ``--group_thresh <value>``: The p-value threshold used to determine whether groups are merged or kept separate. Default = 0.01.
* ``--n_top_groups <integer>``: The number of top mass differences used when constructing the groups. Active only if adaptive = F. Default = 4.

--------------------
Group-walk algorithm
--------------------

* ``--K <integer>``: The number of recently observed peptides used to estimate the probability that the next peptide is a target or decoy. Default = 40.
* ``--return_extra_mods <T|F>``: If ``--account_mods T``, all target peptides equal to a peptide up to variable modification used in the group-walk algorithm will also be reported with the same q-value. Default = F.
* ``--return_frontier <T|F>``: The sequence of indices describing the positions of the frontier used by Groupwalk is returned as a .txt file to the output directory. Default = F.

-------------
CPU processes
-------------

* ``--n_processes <integer>``: The number of processes used in the filtering process. Default = 1.

----------------
Input and output
----------------

* ``--output_dir <string>``: The file-path of the output. Default = './'.
* ``--file_root <string>``: The file name of the output. Default = group_walk_results.txt.
* ``--print_chimera <T|F>``: To determine whether we print the number of scans that have more than 1 peptide discovered at the 1% and 5% FDR level to the log file. Default = T.
* ``--print_group_pi0 <T|F>``: To determine whether the group-level proportion of pi_0 is printed to the log file. Default = T.
* ``--return_decoys <T|F>``: Also report decoys. Default = F.
* ``--static_mods <string>``: Of the form X:[+-]A where X is the amino acid, or rather "cterm" or "nterm" if it is a modification on the C-terminal or N-terminal. A is the absolute mass shift in Daltons. [+-] indicates whether the mass shift is positive or negative. C+57.02146 is always included by default. Variable modifications do not need specification (they are accounted for via the search files). List mods in comma separated format, e.g. nterm:10,cterm:-20,L:50. Default = None.
* ``--dcy_prefix <string>``: The prefix used for the decoy proteins. Default = 'decoy\_'.
* ``--seed <int>``: Set random seed. Default = None.