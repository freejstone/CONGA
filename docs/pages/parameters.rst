"""""""""""
Parameters
"""""""""""

CONGA.py supports the following input parameters:

-------------------
Neighbor filtration
-------------------

* ``-- tops_open <integer>``: The number of top PSMs in the open search file used in the neighbour-filtering process. Default = 5.
* ``--neighbour_remove <T|F>``: If true, for each scan, we successively move down the list of PSMs associated with each scan, ordered in terms of the score from highest to lowest,  and compute a similarity score between the current peptide and the previous peptide(s). If one of the similarity score(s) exceeds a certain threshold, then the PSM associated with the current peptide is thrown out, and we proceed to the next.
* ``--thresh <float>``: The similarity score used as a threshold to filter out neighbouring peptides. Default = 0.05.
* ``--return_filt_search <T|F>``:  Whether or not to return filtered narrow and open search files. Default = F.
* ``--mz_error <float>``: Tolerance in mz for deciding whether two peaks are matched. Also used during pyascoring for localization of modifications. Default = 0.05.

-------------------------
Dynamic-level competition
-------------------------

* ``--score <string>``: Either 'tailor_score', 'xcorr_score', 'e-value' or 'hyperscore'. The score that will be used in the peptide-level competition and subsequent Group construction and Group-walk algorithm. If 'tailor_score', it is assumed the search files are derived from Tide. If 'xcorr_score', either Tide or Comet is assumed to be used. If 'e-value', Comet is assumed. If 'hyperscore' it is assumed the search files are derived from MS-Fragger. Default = 'tailor_score'.
* ``--account_mods <T|F>``: To determine whether CONGA uses the best PSM among the equivalent classes of peptides which are equal up to variable modification, or not. Default = T.

------------------
Group construction
------------------

* ``--precursor_bin_width <value>``: To determine the size of the bins used to discretize the mass-differences between the sample and theoretical spectra. Default = 1.0005079/4.
* ``--min_group_size <integer>``: The number of multiples of K that is used to determine the minimum size of each group. See option ``--K``. Default = 2.
* ``--tops_gw <integer>``: The number of top PSMs for each scan in the open search that will be used by group-walk. Default = 2.
* ``--group_thresh <value>``: The p-value threshold used to determine whether groups are merged or kept separate in the KS test. Default = 0.01.

--------------------
Group-walk algorithm
--------------------

* ``--K <integer>``: The number of recently observed peptides used to estimate the probability that the next peptide is a target or decoy. Default = 40.
* ``--FDR_threshold <value>``: The FDR threshold. Default = 0.01.

--------------------
Localization scoring
--------------------

* ``--spectrum_files <string>``: Comma-separated file paths to the mzML spectrum files used during the original MS/MS search. If not specified, localized scoring using pyAscore is not conducted. File path needs to exactly match the file path used during Tide Search/Comet, as these searches produce a file column which is what is used to correctly match the scan number used in CONGA to the scan number if the spectrum files. MS-Fragger does not contain a file column, so it is assumed just a single spectrum file was used. Default = None.
* ``--mods_to_localize <string>``: Of the form X:[+-]A where X is the amino acid. A is the absolute mass shift in Daltons. [+-] indicates whether the mass shift is positive or negative. pyAscore will be used to isolate the most-likely site containing the modification using the user-supplied modifications if the observed delta mass is reasonably close (up to to the isolation window). Else, the observed delta mass is uesd for localization instead or if it produces a better score. Default = None. List mods in comma-separated format, e.g. S:79.966331,T:79.966331,M:15.9949.
* ``--mods_for_correction <string>``: Variable modifications used during MS/MS search. Often the mass of a modification is rounded in the output file from an MS/MS search. This just allows these rounded modifications to be replaced by more accurate values so when it comes to localized scoring via pyAscore, the errors associated to the rounded masses do not compound if several variable modifications exist on a single peptide. Specify the modifications in the same way as the ``--static_mods`` option below. Default = None.

-------------
CPU processes
-------------

* ``--n_processes <integer>``: The number of processes used in the filtering process. This is very valuable as the filtering process is typically the longest part of CONGA. Default = 1.

----------------
Input and output
----------------

* ``--output_dir <string>``: The file-path of the output. Default = './'.
* ``--file_root <string>``: The file prefix of the output files. Default = 'conga'.
* ``--print_chimera <T|F>``: To determine whether we print the number of scans that have more than 1 peptide discovered at the 1% and 5% FDR level to the log file. Default = T.
* ``--isolation_window <value>``: A comma-separated pair of values (in m/z) used to cluster the precursors masses matched to the same (unmodified) stem peptide. For each cluster, and for each modified variant of the stem peptide, the best scoring PSM is selected and subsequently reported in ``file_root.target_mods.txt``. Default = 2,2.
* ``--print_group_pi0 <T|F>``: To determine whether the group-level proportion of pi_0 is printed to the log file. Default = T.
* ``--return_extra_mods <T|F>``: All top 1 PSM-ranked target peptides that are equal to a discovered peptide up to variable modification will be included in the file output. Default = F.
* ``--return_decoys <T|F>``: Also report decoys used to estimate the number of false discoveries. Default = F.
* ``--return_frontier <T|F>``: The sequence of indices describing the positions of the frontier used by the group-walk algorithm is returned as a .txt file to the output directory. Default = F.
* ``--static_mods <string>``: Of the form X:[+-]A where X is the amino acid, or rather "cterm" or "nterm" if it is a modification on the C-terminal or N-terminal of the peptide. A is the absolute mass shift in Daltons. [+-] indicates whether the mass shift is positive or negative. C+57.02146 is always included by default. Variable modifications do not need specification (they are accounted for via the search files). List mods in comma separated format, e.g. nterm:10,cterm:-20,L:50. Default = None.
* ``--return_mass_mod_hist <T|F>``: To determine whether a histogram of unaccounted-for mass-modifications is written to the directory. Default = T.
* ``--dcy_prefix <string>``: The prefix used for the decoy proteins. Default = 'decoy\_'.
* ``--overwrite <T|F>``: If a log file with the same file root is detected in the output directory, CONGA will not overwrite the files unless this option is set to T. Default = F.
* ``--seed <int>``: Set random seed. Default = None.
