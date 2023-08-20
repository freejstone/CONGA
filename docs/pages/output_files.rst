""""""""""""
Output files
""""""""""""

CONGA returns files with names of the form ``file_root.content.ext``, where ``file_root`` is provided by the user, ``content`` indicates the contents of the file, and ``ext`` indicates the file type.
The files ``file_root.target.txt`` and ``file_root.target_mods.txt`` are always returned, whereas the others are optional.

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

The extra column names found in ``file_root.target_mods.txt`` are as follows.

* ``originally_discovered``: A peptide discovered at an FDR threshold will also have their subsequent variants (other variable modifications and delta masses) of this peptide reported. This column indicates whether the reported row was originally discovered by CONGA, or if it is one of these subsequent variants.
* ``above_group_threshold``: This column indicates whether the reported peptide exceeds the corresponding group threshold that it belongs to. For ``originally_discovered`` peptides, this is trivially true. For subsequent variants that are reported, this may help with filtering for correct identifications with PTMs for later analysis.
* ``localized_peptide``: The best-scoring localization of the peptide in the ``peptide`` via pyAscore. For PSMs taken from the narrow-search file, this value is the same as the value in ``peptide``.
* ``localized_better``: A boolean indicating whether the peptide in ``localized_peptide`` scored better than peptide in ``peptide``.
*  ``dm_used``: A boolean indicating whether the observed delta-mass (in ``delta_mass``) was used to localize the peptide (in ``localized_peptide``) or whether a user supplied mass-modification was used instead (see ``--mods_for_correction`` in the parameters file).