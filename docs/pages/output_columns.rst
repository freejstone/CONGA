""""""""""""""
Output columns
""""""""""""""

We provide the following descriptions of the input parameters used by Groupwalk.py.

* ``peptide``: The sequence of the reported peptide.
* ``q_value``: The q_value associated to the reported peptide.
* ``scan``: The scan number responsible for identifying the reported peptide.
* ``score``: The score of the PSM between the scan and the peptide.
* ``delta_mass``: The mass difference (in Daltons) between the scan and the peptide.
* ``rank``: The rank of the PSM in the search-file it was identified from.
* ``search_file``: Indicating whether the PSM was taken from the narrow- or open-search file.
* ``originally_discovered``: This column is made available when the option ``--return_extra_mods T`` is used. In this case, a peptide discovered at an FDR threshold will also have subsequent variants (other variable modifications) of this peptide reported satisfying the condition that each variant must exceed the score threshold of the group containing that variant.