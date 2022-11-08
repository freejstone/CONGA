""""""""""""""""
Tutorial example
""""""""""""""""

On this page we walk through a simple example of running Open Group-walk for the three different search programs, Tide, Comet and MSFragger.



.. toctree::
   :maxdepth: 2
   :caption: Contents:

=========
Crux-Tide
=========

You can download the search files to test the commands below: `Tide-search files <https://github.com/freejstone/open_groupwalk/tree/main/docs/pages/files/tide>`_.

  To run Open Group-walk using Tide results with default parameters, use the following command in the terminal.

``python3 Groupwalk.py narrow_example1.tide-search.txt open_example1.tide-search.txt``

Open Group-walk collapses peptides that are equal up to variable modifications by selecting the best-scoring peptide. If you want to run Open Group-walk without this collapsing, and thus allowing the discovery of peptides that differ by some variable modification, you can turn this collapsing off by setting ``--account_mods F``. Note that the assumptions of theoretical FDR control do not strictly hold in this case, because target peptides with variable modifications may be scored highly simply because they share similar b- and y-ions in common with the same target peptide but with possibly a different set of variable modifications. In any case, if you wish to run without this collapsing, you will need to provide the target-decoy peptide pairs produced by tide-index.

``python3 Groupwalk.py --account_mods F narrow_example1.tide-search.txt open_example1.tide-search.txt tide-index_example1.peptides.txt``

Alternatively, it might be of interest to report all other top-1 matched peptides that differ from a genuinely discovered peptide by some variable modification(s). One can achieve this by setting ``--return_extra_mods T``:

``python3 Groupwalk.py --return_extra_mods T narrow_example1.tide-search.txt open_example1.tide-search.txt``

The default score uses Tailor scores. If you prefer to use Xcorr Scores, you can change the score option ``--score xcorr_score``:

``python3 Groupwalk.py --score xcorr_score narrow_example1.tide-search.txt open_example1.tide-search.txt``

The software will automatically detect if the PSMs in the search files have been searched against a concatenated target-decoy database or separately against each database. If the latter, then you will need to have the name *target* appear in the search file names. In that case, the software will search for the respective decoy-search files by replacing *target* with *decoy* in the search file names in the same directory.

``python3 Groupwalk.py narrow_example2.tide-search.target.txt open_example2.tide-search.target.txt``

=====
Comet
=====

You can download the search files to test the commands below: `Comet-search files <https://github.com/freejstone/open_groupwalk/tree/main/docs/pages/files/comet>`_.

You can use either standalone Comet or the implementation of Comet in Crux. The Open Group-walk software will automatically standardize the difference in terminology between the standalone version and the one in Crux. To run Open Group-walk using Comet-search results with default parameters, use the following command in the terminal.

``python3 Groupwalk.py --score e-value narrow_example1.comet.txt open_example1.comet.txt``

There is no need to provide a list of target-decoy pairs, because Comet always reverses the target peptides to yield decoy peptides. Comet can use either Xcorr scores or E-values. You will need to specify the score, because the default (Tailor score) is specific to Tide.

All other discussion in :ref:`Crux-Tide` equally applies to Comet search files.

=========
MSFragger
=========

You can download the search files to test the commands below: `MSFragger-search files <https://github.com/freejstone/open_groupwalk/tree/main/docs/pages/files/MS>`_.

Using MSFragger search files is more challenging. Neither MSFragger nor any of the related software in Fragpipe produce decoy peptides that *pair* with the target peptides, which is essential for Open Group-walk. Unfortunately, not even reversing the entire protein sequence and digesting the protein to produce decoy peptides yields the same result as reversing the target peptides. Instead, the user must always provide a target-decoy pairing in the same format as what tide-Index produces. This is reasonably achievable by matching the parameters in tide-index with MSFragger. More specifically, we have the following equivalencies between the two pieces of software:

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

The default precision of the mass-modifications in tide-index is to two decimal places, while MSFragger uses much higher precision. For this reason, Open Group-walk will round the variable modifications in the MSFragger search files to two decimal places, allowing us to pair the peptides in MSFragger according to the target-decoy pairs produced by tide-index. To run Open Group-walk using MSFragger-search results with default parameters, use the following command:

``python3 Groupwalk.py --score hyperscore narrow.tsv open.tsv tide-index.peptides.txt``

Like Comet, the score needs to be specified, because the default is Tailor score. All other discussion in :ref:`Crux-Tide` equally applies to MSFragger-search files.

