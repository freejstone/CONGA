""""""""""""
Output files
""""""""""""

CONGA returns files with names of the form ``file_root.content.ext``, where ``file_root`` is provided by the user, ``content`` indicates the contents of the file, and ``ext`` indicates the file type.
The file ``file_root.target.txt`` is always returned, whereas the others are optional.

* ``file_root.target.txt``: A list of target peptides discovered by CONGAs at a user-specified FDR threshold.
* ``file_root.target_mods.txt``: A list of target peptides with distinct delta-masses and variable modifications that are associated to the list of target peptides discovered in ``file_root.target.txt``.
* ``file_root.decoy.txt``: A list of decoy peptides used to estimate the number of false discoveries at a user-specified FDR threshold. Returned if ``--return_decoys T``.
* ``file_root.unaccounted-mass-mods.pdf``: A histogram of the unaccounted-for mass-modifications in ``file_root.target.txt``. Returned if ``--return_mass_mod_hist T``.
* ``file_root.frontier.txt``: The complete sequence of positions of the frontier vector used by group-walk. Returned if ``--return_frontier T``.

The column names found in ``file_root.target.txt`` are as follows.

* ``peptide``: The sequence of the reported peptide.
* ``scan``: The scan number responsible for identifying the reported peptide.
* ``score``: The score of the PSM between the scan and the peptide.
* ``delta_mass``: The mass difference (in Daltons) between the scan and the peptide.
* ``rank``: The rank of the PSM in the search file (either the narrow- or open-search file) it was identified from.
* ``search_file``: Indication of whether the PSM was taken from the narrow- or open-search file.
* ``charge``: Charge of the precursor.
* ``spectrum_neutral_mass``: Neutral mass of the precursor.
* ``modification_info``: Contains the variable modification information of the discovered peptide as a comma-delimited list of "position[mass-modification]".
* ``flag``: Flags whether the peptide discovered has a ``delta_mass`` value that coincides with a loss or gain of an amino acid. (Works only for Tide-search inputs).
* ``originally_discovered``: This column is found in ``file_root.target_mods.txt``. In this case, a peptide discovered at an FDR threshold will also have their subsequent variants (other variable modifications and delta masses) of this peptide reported.
