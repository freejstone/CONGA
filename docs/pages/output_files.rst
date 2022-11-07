""""""""""""
Output files
""""""""""""

We provide the following descriptions of the (possible) output files. ``file_root.target.txt`` is always returned. 

* ``file_root.target.txt``: A .txt file containing a list of target peptides discovered by Open group-walk at a user-specified FDR threshold.
* ``file_root.decoy.txt``: A .txt file containing a list of decoy peptides used to estimate the number of false discoveries at a user-specified FDR threshold. Returned if ``--return_decoys T``.
* ``file_root.unaccounted-mass-mods.pdf``: A .pdf file containing a histogram of the unaccounted-for mass-modifications in ``file_root.target.txt``. Returned if ``--return_mass_mod_hist T``.
* ``file_root.frontier.txt``: A .txt file containing the complete sequence of positions of the frontier vector used by group-walk. Returned if ``--return_frontier T``.

We provide the following descriptions of the column names found in ``file_root.target.txt``.

* ``peptide``: The sequence of the reported peptide.
* ``scan``: The scan number responsible for identifying the reported peptide.
* ``score``: The score of the PSM between the scan and the peptide.
* ``delta_mass``: The mass difference (in Daltons) between the scan and the peptide.
* ``rank``: The rank of the PSM in the search-file (either the narrow- or open-search file) it was identified from.
* ``search_file``: Indicating whether the PSM was taken from the narrow- or open-search file.
* ``originally_discovered``: This column is made available when the option ``--return_extra_mods T`` is used. In this case, a peptide discovered at an FDR threshold will also have their subsequent variants (other variable modifications) of this peptide reported, so long as these variants are identified by a top 1 peptide-spectrum match (top 1 PSM).
