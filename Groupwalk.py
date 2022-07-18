#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 15:27:59 2022

@author: jackfreestone

#Code that executes the group-walk algorithm
"""
import sys
import pandas as pd
import random
import numpy as np
from scipy import stats
import peptides
import multiprocessing
import os 
import logging
#import re


USAGE = """USAGE: Groupwalk.py <narrow> <wide> <matching>

  This program implements the Group-walk algorithm, including the
  creation of groups and the filtering procedure that eliminates
  similar peptide matches per scan, The first input file is the narrow
  search file output, the second input file is the open search file
  output, and the last input file contains the target-decoy peptide
  pairs. Output is a list of PSMs and their corresponding q-values.
  The input search files can be either tab-delimited Tide search files
  or tab-delimited .tsv files from MS-Fragger. The output file is a
  tab-delimited .txt file containing the following information about the PSMs
  used by groupwalk: the peptide sequence, the scan number, the 
  score, the target-decoy label, whether the PSM came from the narrow
  or open search file, and a q-value.

  Options:
      
    --output_dir <string> The file-path of the output
                          Default = './'.
                          
    --file_name <string>  The file-name of the output
                          Default = group_walk_results.txt.

    --K <integer>         The number of recently observed peptides used
                          to estimate the probability that the next
                          peptide is a target or decoy.
                          Default = 40.

    --tops_gw <integer>   The number of top PSMs for each scan in the open 
                          search that will be used by group-walk.
                          Default = 2.
    
    --tops_open <integer> The number of top PSMs in the open search used in
                          the neighbour-filtering process.
                          Default = 5.     
                          

    --score <string>      Either 'tailor' or 'xcorr' or 'hyper'. The score 
                          that will be used by group-walk. If either
                          'tailor' or 'xcorr' is used, it is assumed the
                          search files are derived from Tide. If
                          'hyper' it is assumed the search files are
                          derived from MS-Fragger.
                          Default = 'tailor'.

    --random_h2h <T|F>    To determine whether target-decoy head-to-head
                          competition uses a random selection in case of
                          ties (or deletes both target and decoy from
                          consideration).
                          Default = T.
    
    --account_mods <T|F>  To determine whether the group-walk algorithm
                          selects the best PSM among the equivalent
                          classes of peptides which are equal up to
                          variable modification, or not.
                          Default = T.
                          
    --precursor_bin_width <value>   To determine the size of the bins
                                    used to create discretize the mass-
                                    differences between the sample
                                    and theoretical spectra.
                                    Default = 1.0005079/4.
                                    
    --adaptive <T|F>      To determine whether groups should be chosen
                          using the Kolmogorov-Smirnov test or if a
                          fixed procedure should be used.
                          Default = T.
                          
    --print_chimera <T|F>          To determine whether we print the number
                                   of scans that have more than 1 peptide
                                   discovered at the 1% and 5% FDR level
                                   to the log file.
                                   Default = T.
            
    --group_thresh <value>         The p-value threshold used to determine
                                   whether groups are merged or kept
                                   separate.
                                   Default = 0.01.
                                   
    --print_group_pi0 <T|F>        To determine whether the group-level
                                   proportion of pi_0 is printed to the 
                                   log file.
                                   Default = T.
    
    --min_group_size <integer>     The number of multiples of K that the
                                   used to determine the minimum size of
                                   each group. See option --K.
                                   Default = 2.
                                   
    --n_top_groups <integer>       The number of top mass differences used
                                   when constructing the groups. Active
                                   only if adaptive = False.
                                   Default = 4.
    
    --neighbour_remove <T|F>       If true, for each scan, we successively
                                   move down the list of PSMs associated with
                                   each scan, ordered in terms of the score
                                   from highest to lowest, and compute
                                   a similarity score between the current
                                   peptide and the previous peptide(s).
                                   If one of the similarity score(s) exceeds
                                   a certain threshold, the PSM associated 
                                   with the current peptide is thrown out,
                                   and we proceed to the next.
                                   Default = T.
    
    --thresh <value>      The similarity score used as a threshold to filter
                          out neighbouring peptides.
                          Default = 0.25.
    
    --concat <T|F>        To determine whether the input search files
                          were searched against the concatenated target
                          decoy database, or against the databases
                          separately. If separate the two input search files
                          should be from the target database. Groupwalk
                          will automatically search for the decoy database
                          search files in the same file path by replacing
                          the text "target" with "decoy".
                          Default = T.
                          
    --return_filt_search <T|F>  Whether or not to return filtered narrow
                                and open search files.
                                Default = F.
                                
    --return_frontier <T|F>     The sequence of indices describing the
                                positions of the frontier used by Groupwalk
                                is returned as a .txt file to the output
                                directory.
                                Default = F.
                                    
    --frontier_name <string>    The file-name of the output frontier.
                                Default = frontier_results.txt.
    
    --n_processes <integer>     The number of threads, used in the filtering
                                process.
                                Default = 1.
                                
    --static_mods <string>      Of the form X:[+-]A where X is the amino acid,
                                or rather "cterm" or "nterm" if it is a
                                modification on the C-terminal or N-terminal.
                                A is the absolute mass shift in Daltons.
                                [+-] indicates whether the mass shift is
                                positive or negative. C+57.02146 is always
                                included by default. Variable modifications
                                do not need specification (they are accounted
                                for via the search files). List mods in comma
                                separated format, e.g.
                                nterm:10,cterm:-20,L:50.
                                Default = None.
            
"""

###############################################################################
def parse_static_mods(my_string):
  """
  Parse a static mod string (see USAGE) into a dictinoary.
  Key = amino acid, value = mass offset
  """
  return_value = {}
  
  for one_mod in my_string.split(","):
    words = one_mod.split(":")
    return_value[words[0]] = float(words[1])

  # Add in carbamidomethylation.
  if ("C" not in return_value):
    return_value["C"] = 57.02146
  
  return(return_value)
  
  
###############################################################################
def make_binidx_matchcount_map(mzs, fragment_min_mz, fragment_max_mz, bin_size):
    """
    Utility function for calc_proportion_fragments_incommon.
    Notionally bin the mzs from a list of mzs using the specified bin
    size (don't actually build the binned array), with specified lower
    and upper limits. Return a map from bin indexes to a count of
    fragments
    :param mzs: 
    :param fragment_min_mz:
    :param fragment_max_mz:
    :param bin_size:
    :return:  a map from bin index to a count of fragments in that bin
    """
    binidx_matchcount_map = {}
    for mz in mzs:
        if mz < fragment_min_mz or mz > fragment_max_mz:
            continue
        bin_idx = int((mz - fragment_min_mz) / bin_size)
        if bin_idx not in binidx_matchcount_map:
            binidx_matchcount_map[bin_idx] = 0
        binidx_matchcount_map[bin_idx] += 1
    return binidx_matchcount_map

###############################################################################
def calc_proportion_fragments_incommon(pepseq1, peptide1_mods, nterm1, cterm1,
                                       pepseq2, peptide2_mods, nterm2, cterm2,
                                       binsize, charge1, charge2):
    """
    Determine all the fragment mzs for each peptide. Bin the fragments from
    each peptide.  Calculate the fraction of fragments that fall into a bin
    with a fragment from the other peptide.

    :param pepseq1: First peptide sequence
    :param peptide1_mods: List of mass offsets, same length as peptide.
    :param nterm1: N-terminal modification mass of first peptide.
    :param cterm1: C-terminal modification mass of first peptide.
    :param pepseq2: Second peptide sequence
    :param peptide2_mods: List of mass offsets, same length as peptide.
    :param nterm2: N-terminal modification mass of second peptide.
    :param cterm2: C-terminal modification mass of second peptide.
    :param binsize: Size of m/z bins.
    :return: Fraction of matched peaks.
    """

    # Compute fragments and put them in bins.
    if charge1 in [1, 2]:
        fragment_charge1 = [1]
    elif charge1 >= 3:
        fragment_charge1 = [1, 2]
    
    if charge2 in [1, 2]:
        fragment_charge2 = [1]
    elif charge2 >= 3: #not sure if this is "kosher"
        fragment_charge2 = [1, 2]
    
    
    mzs_1 = peptides.calc_theoretical_peak_mzs(pepseq1, fragment_charge1, peptide1_mods, 
                                               200, 3000, nterm1, cterm1)
    mzs_2 = peptides.calc_theoretical_peak_mzs(pepseq2, fragment_charge2, peptide2_mods,
                                               200, 3000, nterm2, cterm2)
    bin_count_map_1 = make_binidx_matchcount_map(mzs_1, 200, 3000, binsize)
    bin_count_map_2 = make_binidx_matchcount_map(mzs_2, 200, 3000, binsize)

    # Count matched bins.
    n_fragments_in_matched_bins = 0
    for binidx in bin_count_map_1:
        if binidx in bin_count_map_2:
            n_fragments_in_matched_bins += (bin_count_map_1[binidx]
                                            + bin_count_map_2[binidx])

    return float(n_fragments_in_matched_bins) / (len(mzs_1) + len(mzs_2))

###############################################################################
def parse_mods(pepseq_with_mods, static_mods):
    """
    Parse a modified peptide sequence string.

    :param pepseq_with_mods: Peptide string with interpolated bracketed mods.  
    :param static_mods: Dictionary of static mods. 
                        Key = amino acid, value = mass offset.
    :return: A list of amino acids, a list of modification values, and 
             the n-term and c-term terminal modification values.
    """
    aa_seq_list = []
    modifications = []
    nterm_delta = 0.0
    cterm_delta = 0.0

    aaseq_position = -1
    modpepseq_position = 0
    while modpepseq_position < len(pepseq_with_mods):
        my_char = pepseq_with_mods[modpepseq_position]
        if my_char.isalpha():

            # Create an unmodified amino acid and add it to the growing list.
            aa_seq_list.append(my_char)
            modifications.append(0.0)
            aaseq_position += 1
            modpepseq_position += 1
        elif ( my_char == '[' ):
            end_pos = (pepseq_with_mods[modpepseq_position + 1:].index(']') 
                       + modpepseq_position + 1)
            mod_mass = float(pepseq_with_mods[modpepseq_position + 1:end_pos])

            # Store a modification associated with the current position.
            modifications[aaseq_position] = mod_mass
            modpepseq_position = end_pos + 1
        else:
            sys.stderr.write("Invalid character (%s) at position %d.\n"
                             % (my_char, modpepseq_position))
            sys.exit(1)

    # Add in the static mods.
    for index in range(0, len(aa_seq_list)):
        amino = aa_seq_list[index]
        if (amino in static_mods):
          modifications[index] += static_mods[amino]
    if ("nterm" in static_mods):
      nterm_delta = static_mods["nterm"]
    if ("cterm" in static_mods):
      cterm_delta = static_mods["cterm"]


    return(aa_seq_list, modifications, nterm_delta, cterm_delta)
    
###############################################################################
"""
def parse_mods_for_MS_fragger(pepseq_no_mods, mods):
    #split peptide into amino acids
    aa_seq_list = list(pepseq_no_mods)
    modifications = [0.0]*len(aa_seq_list)
    
    #split mod string into list of strings based off commas
    mods_list = mods.split(', ')
    
    for i in range(len(mods_list)):
        #pos = position of the amino acid mass shift is associated to
        #mass_shift = associated mass difference at that amino acid
        (pos, mass_shift) = re.findall(r'\d+', mods_list[i])
        modifications[pos - 1] = mass_shift
    
    return(aa_seq_list, modifications)
"""
###############################################################################
def get_similarity(list1, charges1, list2, charges2, frag_bin_size = 0.05, static_mods = {'C':57.02146}):
  
  peptide_out_1 = []
  peptide_out_2 = []
  similarity_out = []
  start_index = 0

  # Compare each pair of peptides.
  num_pairs = 0
  num_printed = 0
  for index1 in range(0, len(list1)):
    peptide1 = list1[index1]
    (unmodified_peptide1, peptide1_mods, nterm1, cterm1) \
       = parse_mods(peptide1, static_mods)

    for index2 in range(start_index, len(list2)):
      peptide2 = list2[index2]
      num_pairs += 1


      # Don't bother if they are the same peptide.
      if (peptide1.replace("I", "L") == peptide2.replace("I", "L")):
        peptide_out_1.append(peptide1)
        peptide_out_2.append(peptide2)
        similarity_out.append(1.0)
      else:
         
        (unmodified_peptide2, peptide2_mods, nterm2, cterm2) \
         = parse_mods(peptide2, static_mods)
         
        charge1 = charges1[index1]
        charge2 = charges2[index2]

        similarity = calc_proportion_fragments_incommon(
          unmodified_peptide1, peptide1_mods, nterm1, cterm1,
          unmodified_peptide2, peptide2_mods, nterm2, cterm2,
          frag_bin_size, charge1, charge2)
        

        num_printed += 1
        
        peptide_out_1.append(peptide1)
        peptide_out_2.append(peptide2)
        similarity_out.append(similarity)
  return(peptide_out_1, peptide_out_2, similarity_out)
###############################################################################
""" 
def get_similarity_MS(list1, mods1, charges1, list2, mods2, charges2, frag_bin_size = 0.05):
  
  peptide_out_1 = []
  peptide_out_2 = []
  similarity_out = []
  start_index = 0

  # Compare each pair of peptides.
  num_pairs = 0
  num_printed = 0
  for index1 in range(0, len(list1)):
    peptide1 = list1[index1]
    mod1 = mods1[index1]
    (unmodified_peptide1, peptide1_mods) \
       = parse_mods_for_MS_fragger(peptide1, mod1)

    for index2 in range(start_index, len(list2)):
      peptide2 = list2[index2]
      mod2 = mods2[index2]
      num_pairs += 1


      # Don't bother if they are the same peptide.
      if (peptide1.replace("I", "L") == peptide2.replace("I", "L")):
        peptide_out_1.append(peptide1)
        peptide_out_2.append(peptide2)
        similarity_out.append(1.0)
      else:
         
        (unmodified_peptide2, peptide2_mods) \
         = parse_mods_for_MS_fragger(peptide2, mod2)
         
        charge1 = charges1[index1]
        charge2 = charges2[index2]

        similarity = calc_proportion_fragments_incommon(
          unmodified_peptide1, peptide1_mods, 0.0, 0.0,
          unmodified_peptide2, peptide2_mods, 0.0, 0.0,
          frag_bin_size, charge1, charge2)
        

        num_printed += 1
        
        peptide_out_1.append(peptide1)
        peptide_out_2.append(peptide2)
        similarity_out.append(similarity)
  return(peptide_out_1, peptide_out_2, similarity_out)
"""
###############################################################################
def group_walk(winning_scores, labels, all_group_ids, K = 40, return_frontier = False, correction = 1):
    if not (len(winning_scores) == len(labels) == len(all_group_ids)):
        logging.error("winning_scores, labels, all_group_ids are not of the same length.")
        sys.exit("winning_scores, labels, all_group_ids are not of the same length.")
    
    logging.info("Applying group walk.")
    sys.stderr.write("Applying group walk.\n")
        
    ordered_inds = sorted(range(len(winning_scores)), key=lambda k: winning_scores[k]) # record the indices in order of the scores from smallest to largest
    
    #reorder groups and labels according to the winning_scores
    ordered_groups = pd.Series([all_group_ids[i] for i in ordered_inds])
    ordered_labels = pd.Series([labels[i] for i in ordered_inds])
    
    #initailise the q_vals
    q_vals = pd.Series([1.0] * len(winning_scores))
    q_val = 1
    
    #get unique set of group_ids
    group_ids = pd.Series(list(set(sorted(all_group_ids))))
    g_total = len(group_ids)
    if type(all_group_ids) == list:
        totals = np.array([all_group_ids.count(x) for x in list(group_ids)])
    else:
        totals = np.array(all_group_ids.value_counts())
    totals += -1
    
    #create individualised labels, and q_vals for each group_id
    labels_sorted_pd = [0]*g_total
    labels_sorted_ls = [0]*g_total
    q_vals_sorted = [0]*g_total
    for g in range(g_total):
        labels_sorted_pd[g] = ordered_labels[ordered_groups == group_ids[g]] == 1
        labels_sorted_ls[g] = [x == 1 for x in ordered_labels[ordered_groups == group_ids[g]]]
        q_vals_sorted[g] = pd.Series([1.0]*len(labels_sorted_pd[g]))
        #q_vals_sorted[g].set_axis(list(labels_sorted_pd[g].index), inplace = True)
        labels_sorted_pd[g].reset_index(drop = True, inplace = True)
        q_vals_sorted[g].reset_index(drop = True, inplace = True)
        
    #initialise start of frontier and weights
    starts = np.array([0]*g_total)
    weights = np.array([0]*g_total)
    
    #number of decoys + 1 and targets
    decoys_plus_one = sum([x == False for sublist in labels_sorted_pd for x in sublist]) + correction
    rejections = sum(totals + 1) - decoys_plus_one + correction
    
    index_update = 0
    switch = True
    
    frontiers = [[0]*g_total]*(sum(totals + 1) + 1)
    frontiers = pd.DataFrame(frontiers, columns = group_ids)
    counter = 1
    
    candidate_starts_bool = (starts <= totals)
    
    index_row_mapping = pd.Series(range(len(weights)))
    
    check_greater = True
    randind = 0
    
    while any(candidate_starts_bool):
        #print(starts)
        if return_frontier:
            frontiers.loc[counter] = starts
            counter += 1
            
        #estimate FDR
        FDR_e = decoys_plus_one / max(rejections, 1)
        
        #update q_value
        q_val = min(q_val, FDR_e)
        
        #we need at least each group to surpass K peptides
        if (switch and any((starts <= K - 1) & (candidate_starts_bool))):
            candidate_starts = starts[candidate_starts_bool]
            min_val = min(candidate_starts)
            inds = np.where(candidate_starts == min_val)[0]
        else:
            #print('Yes')
            if switch:
                #cacluate all of the group weights for the first time
                for g in index_row_mapping.index:
                    if starts[index_row_mapping[g]] <= totals[index_row_mapping[g]]:
                        weights[g] = K - sum(labels_sorted_ls[g][(starts[g] - K): (starts[g])])
                
                switch = False
            else:
                #now calculate the group weights for the ones that need changing
                if ~check_greater:
                    s = starts[index_update] - 1
                    #the group weight adds the most recent peptide, and deletes the last peptide
                    weights[randind] = weights[randind] - ~(labels_sorted_ls[index_update][s - K]) + ~(labels_sorted_ls[index_update][s])
            inds = np.where(weights == max(weights))[0]
        if len(inds) == 1:
            randind = inds[0]
        else:
            randind = random.choice(inds)
        #if switch:
        #    index_update = randind
        #else:
        #    index_update = index_row_mapping[randind]
        index_update = index_row_mapping[randind]
        #the frontier that needs updating
        s = starts[index_update]
        
        q_vals_sorted[index_update][s] = q_val
        label_at_update = labels_sorted_pd[index_update][s]
        decoys_plus_one = decoys_plus_one - ~label_at_update # updating the number of decoys
        rejections = rejections - label_at_update # updating the number of targets
        starts[index_update] += 1
        
        check_greater = (starts[index_update] > totals[index_update])
        
        if (check_greater):
            candidate_starts_bool = (starts <= totals)
            weights = weights[index_row_mapping.isin([i for i in range(len(starts)) if starts[i] <= totals[i]])]
            index_row_mapping = index_row_mapping[index_row_mapping.isin([i for i in range(len(starts)) if starts[i] <= totals[i]])] 
            index_row_mapping = index_row_mapping.reset_index(drop = True)

    for g in range(g_total):
        q_vals[ordered_groups == group_ids[g]] = list(q_vals_sorted[g])
    
    q_vals[ordered_inds] = list(q_vals)
    
    if return_frontier:
        q_vals = pd.Series(q_vals, name = 'q_vals')
        results = [q_vals, frontiers]
    else:
        results = [pd.Series(q_vals, name = 'q_vals')]
    
    logging.info("Group walk complete.")
    sys.stderr.write("Group walk complete.\n")
    return(results)

###############################################################################
#This appears to work perfectly fine
def filter_scan(search_df, thresh = 0.25, frag_bin_size = 0.05, static_mods = {'C':57.02146}):
    search_df = search_df.reset_index(drop = True)
    n_scans = search_df.shape[0]
    
    peptide_1 = [search_df['sequence'].loc[0]]
    charge_1 = [search_df['charge'].loc[0]]
    drop_scan = [False]*n_scans
    for top in range(1, n_scans):
        peptide_2 = [search_df['sequence'].loc[top]]
        charge_2 = [search_df['charge'].loc[top]]
        results = get_similarity(peptide_1, charge_1, peptide_2, charge_2, frag_bin_size, static_mods)
        if any(sim >= thresh for sim in results[2]):
            drop_scan[top] = True
        else:
            peptide_1 += peptide_2
            charge_1 += charge_2
   
    return drop_scan

###############################################################################
def filter_narrow_open(narrow_target_decoys, open_target_decoys, open_top = 2, thresh = 0.25, n_processes = 1, neighbour_remove = True, tide_used = True, static_mods = {'C':57.02146}):
    '''
    Parameters
    ----------
    narrow_target_decoys : Pandas Dataframe
        Top 1 PSM for concatenated search of target-decoy database
    open_target_decoys : Pandas Dataframe
        Concatenated serch of target-decoy database
    open_top : int
        The number of top PSMs for each scan in the open search
    to be considered. The default is 2.
    thresh : float
        The threshold used to determine neighbouring peptides.
        The default is 0.25.
    n_processes : int
        Number of proccesses to be used. The default is 1.
    neighbour_remove : bool
        Whether we remove neighbours or not. The default is True.
    tide_used : bool
        Whether the input dataframes come from tide or not.
        The default is True.
    static_mods : dict
        Static modifications. The default is {'C':57.02146}.

    Returns
    -------
    target_decoys_all : Pandas Dataframe
        Combined PSMs from both narrow and open search, with neighbours
        for each scan number.
    '''
    open_target_decoys['database'] = 'open'
    narrow_target_decoys['database'] = 'narrow'
    target_decoys_all = pd.concat([narrow_target_decoys, open_target_decoys])
    
    #makes sure the PSMs from the narrow and open search are ordered by scan first, then by their xcorr_score
    #indeed it makes sense to use xcorr_score over tailor score, as tailor_score uses candidate PSMs to 
    #normalise the xcorr_score - this differs from narrow to open. 
    if tide_used:
        target_decoys_all = target_decoys_all.sort_values(by=['scan', 'xcorr_score'], ascending = False)
    else:
        target_decoys_all = target_decoys_all.sort_values(by=['scannum', 'hyperscore'], ascending = False)
    target_decoys_all.reset_index(drop = True, inplace = True)
    
    if not neighbour_remove:
        logging.info("Not removing neighbours.")
        sys.stderr.write("Not removing neighbours.\n")
        return target_decoys_all
    
    logging.info("Filtering for neighbours.")
    sys.stderr.write("Filtering for neighbours.\n")
    
    if n_processes == 1:
        if tide_used:
            results = target_decoys_all.groupby('scan').apply(filter_scan, thresh, 0.05, static_mods)
        else:
            results = target_decoys_all.groupby('scannum').apply(filter_scan, thresh, 0.05, static_mods)
        results = results.sort_index(ascending = False)
        results = results.apply(pd.Series).stack().reset_index()
        
        target_decoys_all['drop_scan'] = results[0]
    else:
        #if more than 1 thread, we partition the dataframe
        if tide_used:
            target_decoys_all['split_col'] = pd.cut(target_decoys_all['scan'], n_processes)
        else:
            target_decoys_all['split_col'] = pd.cut(target_decoys_all['scannum'], n_processes)
        target_decoys_grouped = target_decoys_all.groupby(target_decoys_all.split_col)
        list_of_df = [0]*n_processes
        for i in range(len(target_decoys_all['split_col'].unique())):
            list_of_df[i] = target_decoys_grouped.get_group(target_decoys_all['split_col'].unique()[i])
            list_of_df[i].reset_index(drop = True, inplace = True)
        
        if tide_used:
            def filter_scan_subset(df, task_number, return_dict):
                sys.stderr.write("Starting Process " + str(task_number) + " \n")
                logging.info("Starting Process " + str(task_number))
                results = df.groupby('scan').apply(filter_scan, thresh, 0.05, static_mods)
                results = results.sort_index(ascending = False)
                results = results.apply(pd.Series).stack().reset_index()
                return_dict[task_number] = results[0]
                sys.stderr.write("Process " + str(task_number) + " finished \n")
                logging.info("Starting Process " + str(task_number) + " finished")
        else:
            def filter_scan_subset(df, task_number, return_dict):
                sys.stderr.write("Starting Process " + str(task_number) + " \n")
                logging.info("Starting Process " + str(task_number))
                results = df.groupby('scannum').apply(filter_scan, thresh, 0.05, static_mods)
                results = results.sort_index(ascending = False)
                results = results.apply(pd.Series).stack().reset_index()
                return_dict[task_number] = results[0]
                sys.stderr.write("Process " + str(task_number) + " finished \n")
                logging.info("Starting Process " + str(task_number) + " finished")
        
        processes = [0]*n_processes
        manager = multiprocessing.Manager()
        return_dict = manager.dict()
        for i in range(len(target_decoys_all['split_col'].unique())):
            processes[i] = multiprocessing.Process(target = filter_scan_subset, args = (list_of_df[i], i, return_dict))
            processes[i].start()
        
        for process in processes:
            process.join()
        
        for i in range(len(target_decoys_all['split_col'].unique())):
            target_decoys_all.loc[target_decoys_all['split_col'] == target_decoys_all['split_col'].unique()[i], 'drop_scan'] = list(return_dict[i])
            
    target_decoys_all = target_decoys_all[target_decoys_all['drop_scan'] == False] 
    target_decoys_all.reset_index(drop = True, inplace = True)
    logging.info("Filtering for neighbours is complete.")
    sys.stderr.write("Filtering for neighbours is complete.\n")
        
    return target_decoys_all
###############################################################################

def create_groups(target_decoys, target_decoys_narrow_TDC, peptide_list, K = 40, tops = 2, score = 'tailor', random_h2h = True, account_mods = True, precursor_bin_width = 1.0005079/4, group_thresh = 0.01, adaptive = True, min_group_size = 2, n_top_groups = 4, tide_used = True, print_group_pi0 = True):
    
    #standardize column names
    target_decoys.columns = target_decoys.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '').str.replace('/', '_').str.replace('.', '_')
    target_decoys_narrow_TDC.columns = target_decoys_narrow_TDC.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '').str.replace('/', '_').str.replace('.', '_')
    
    #take only the top 1 from the narrow search
    if tide_used:
        target_decoys_narrow_TDC = target_decoys_narrow_TDC[target_decoys_narrow_TDC['xcorr_rank'] == 1]
    else:
        target_decoys_narrow_TDC = target_decoys_narrow_TDC[target_decoys_narrow_TDC['hit_rank'] == 1]
    target_decoys_narrow_TDC.loc[:, 'database'] = 'narrow'
    
    #rounding in case of any numerical issues from previous filterings
    if score == 'tailor':
        target_decoys['tailor_score'] = round(target_decoys['tailor_score'], 8)
        target_decoys_narrow_TDC['tailor_score'] = round(target_decoys_narrow_TDC['tailor_score'], 8)
    elif score == 'xcorr':
        target_decoys['xcorr_score'] = round(target_decoys['xcorr_score'], 8)
        target_decoys_narrow_TDC['xcorr_score'] = round(target_decoys_narrow_TDC['xcorr_score'], 8)
    elif score == 'hyper':
        target_decoys['hyperscore'] = round(target_decoys['hyperscore'], 8)
        target_decoys_narrow_TDC['hyperscore'] = round(target_decoys_narrow_TDC['hyperscore'], 8)
    
    logging.info("Doing head to head competition.")
    sys.stderr.write("Doing head to head competition.\n")
    
    target_decoys_open = target_decoys[target_decoys['database'] == "open"]
    target_decoys_narrow = target_decoys[target_decoys['database'] == "narrow"]


    if tide_used:
        targets_narrow = target_decoys_narrow[( target_decoys_narrow['target_decoy'] == 'target' )]
        decoys_narrow = target_decoys_narrow[( target_decoys_narrow['target_decoy'] == 'decoy' )]
         #taking the top "tops" PSMs from the open search
        targets_open = target_decoys_open[( target_decoys_open['target_decoy'] == 'target' ) & ( target_decoys_open['xcorr_rank'].isin(range(1, tops + 1)) )]
        decoys_open = target_decoys_open[( target_decoys_open['target_decoy'] == 'decoy' ) & ( target_decoys_open['xcorr_rank'].isin(range(1, tops + 1)) )]    
        
    else:
        targets_narrow = target_decoys_narrow[target_decoys_narrow['sequence'].isin(peptide_list['target'])].copy()
        decoys_narrow = target_decoys_narrow[target_decoys_narrow['sequence'].isin(peptide_list['decoy'])].copy()
        #targets_narrow = target_decoys_narrow[( target_decoys_narrow['peptide'].isin(peptide_list['target'].str.replace("\\[|\\]|\\.|\d+", "")) )]
        #decoys_narrow = target_decoys_narrow[( target_decoys_narrow['peptide'].isin(peptide_list['decoy'].str.replace("\\[|\\]|\\.|\d+", "")) )]
        targets_narrow.loc[:, 'target_decoy'] = 'target'
        decoys_narrow.loc[:, 'target_decoy'] = 'decoy'
        
        targets_open = target_decoys_open[( target_decoys_open['sequence'].isin(peptide_list['target']) ) & ( target_decoys_open['hit_rank'].isin(range(1, tops + 1)) )].copy()
        decoys_open = target_decoys_open[( target_decoys_open['sequence'].isin(peptide_list['decoy']) ) & ( target_decoys_open['hit_rank'].isin(range(1, tops + 1)) )].copy()  
        
        #targets_open = target_decoys_open[( target_decoys_open['peptide'].isin(peptide_list['target'].str.replace("\\[|\\]|\\.|\d+", "")) ) & ( target_decoys_open['hit_rank'].isin(range(1, tops + 1)) )]
        #decoys_open = target_decoys_open[( target_decoys_open['peptide'].isin(peptide_list['decoy'].str.replace("\\[|\\]|\\.|\d+", "")) ) & ( target_decoys_open['hit_rank'].isin(range(1, tops + 1)) )]    
        targets_open.loc[:, 'target_decoy'] = 'target'
        decoys_open.loc[:, 'target_decoy'] = 'decoy'
        
        #decoys_open['sequence_plus_variable_mods'] = decoys_open['peptide'] + decoys_open['modification_info']
        #targets_open['sequence_plus_variable_mods'] = targets_open['peptide'] + targets_open['modification_info']
        #decoys_narrow['sequence_plus_variable_mods'] = decoys_narrow['peptide'] + decoys_narrow['modification_info']
        #targets_narrow['sequence_plus_variable_mods'] = targets_narrow['peptide'] + targets_narrow['modification_info']
        
    #sorting PSMs by their score and deleting duplicate PSMs that share the same sequence (up to variable modification if account_mods == True)
    if score == 'tailor':
        decoys_open = decoys_open.sort_values(by=['tailor_score'], ascending=False)
        targets_open = targets_open.sort_values(by=['tailor_score'], ascending=False)
        decoys_narrow = decoys_narrow.sort_values(by=['tailor_score'], ascending=False)
        targets_narrow = targets_narrow.sort_values(by=['tailor_score'], ascending=False)
    elif score == 'xcorr':
        decoys_open = decoys_open.sort_values(by=['xcorr_score'], ascending=False)
        targets_open = targets_open.sort_values(by=['xcorr_score'], ascending=False)
        decoys_narrow = decoys_narrow.sort_values(by=['xcorr_score'], ascending=False)
        targets_narrow = targets_narrow.sort_values( by=['xcorr_score'], ascending=False)
    elif score == "hyper":
        decoys_open = decoys_open.sort_values(by=['hyperscore'], ascending=False)
        targets_open = targets_open.sort_values(by=['hyperscore'], ascending=False)
        decoys_narrow = decoys_narrow.sort_values(by=['hyperscore'], ascending=False)
        targets_narrow = targets_narrow.sort_values( by=['hyperscore'], ascending=False)
    if account_mods:
        if tide_used:
            decoys_open = decoys_open.drop_duplicates(subset = ['original_target_sequence'])
            targets_open = targets_open.drop_duplicates(subset = ['original_target_sequence'])
            decoys_narrow = decoys_narrow.drop_duplicates(subset = ['original_target_sequence'])
            targets_narrow = targets_narrow.drop_duplicates(subset = ['original_target_sequence'])
        else:
            #MS_fragger's peptide column does not include variable modifications
            decoys_open = decoys_open.drop_duplicates(subset = ['peptide'])
            targets_open = targets_open.drop_duplicates(subset = ['peptide'])
            decoys_narrow = decoys_narrow.drop_duplicates(subset = ['peptide'])
            targets_narrow = targets_narrow.drop_duplicates(subset = ['peptide'])
    else:
        #sequence is both unique for MS_fragger and tide-search results
        decoys_open = decoys_open.drop_duplicates(subset = ['sequence'])
        targets_open = targets_open.drop_duplicates(subset = ['sequence'])
        decoys_narrow = decoys_narrow.drop_duplicates(subset = ['sequence'])
        targets_narrow = targets_narrow.drop_duplicates(subset = ['sequence'])

    #resetting the indices
    decoys_open.reset_index(drop = True, inplace = True)
    targets_open.reset_index(drop = True, inplace = True)
    decoys_narrow.reset_index(drop = True, inplace = True)
    targets_narrow.reset_index(drop = True, inplace = True)
    
    #creating PSM column to uniquely identify the PSMs
    if tide_used:
        decoys_open.loc[:, 'scan_plus_seq'] = decoys_open['scan'].astype(str) + ' ' + decoys_open['sequence']
        targets_open.loc[:, 'scan_plus_seq'] = targets_open['scan'].astype(str) + ' ' + targets_open['sequence']
    else:
        decoys_open.loc[:, 'scan_plus_seq'] = decoys_open['scannum'].astype(str) + ' ' + decoys_open['sequence']
        targets_open.loc[:, 'scan_plus_seq'] = targets_open['scannum'].astype(str) + ' ' + targets_open['sequence']

    #creating PSM column to uniquely identify the PSMs
    target_decoys_narrow_TDC.reset_index(drop = True, inplace = True)
    if tide_used:
        target_decoys_narrow_TDC.loc[:,'scan_plus_seq'] = target_decoys_narrow_TDC['scan'].astype(str) + ' ' + target_decoys_narrow_TDC['sequence']
    else:
        target_decoys_narrow_TDC.loc[:,'scan_plus_seq'] = target_decoys_narrow_TDC['scannum'].astype(str) + ' ' + target_decoys_narrow_TDC['sequence']
    
    #translating into an iterable
    unique_scan_seq_td_narrow = set(target_decoys_narrow_TDC['scan_plus_seq'])
    
    #decoy/target PSMs open geta ssigned "narrow" database if they exist in the narrow search
    decoys_open.loc[decoys_open['scan_plus_seq'].isin(unique_scan_seq_td_narrow), 'database'] = 'narrow'    
    targets_open.loc[targets_open['scan_plus_seq'].isin(unique_scan_seq_td_narrow), 'database'] = 'narrow'
    
    #subsetting the decoy-target PSMs in the open-search that are also found in the narrow search
    scan_plus_seq_decoys_sub = decoys_open[decoys_open['scan_plus_seq'].isin(unique_scan_seq_td_narrow)]['scan_plus_seq']
    scan_plus_seq_targets_sub = targets_open[targets_open['scan_plus_seq'].isin(unique_scan_seq_td_narrow)]['scan_plus_seq']
    
    #translating into iterables
    unique_scan_plus_seq_decoys_sub = set(scan_plus_seq_decoys_sub)
    unique_scan_plus_seq_targets_sub = set(scan_plus_seq_targets_sub)
    
    #take subset of narrow PSMs that contain the decoy PSMs in the open
    #reorder this sset so that they match the decoy PSMs in the open
    target_decoys_narrow_sub = target_decoys_narrow_TDC[target_decoys_narrow_TDC['scan_plus_seq'].isin(unique_scan_plus_seq_decoys_sub)].copy()
    target_decoys_narrow_sub.reset_index(drop = True, inplace = True)
    target_decoys_narrow_sub.scan_plus_seq = target_decoys_narrow_sub.scan_plus_seq.astype("category")
    target_decoys_narrow_sub.scan_plus_seq.cat.set_categories(scan_plus_seq_decoys_sub, inplace = True)
    target_decoys_narrow_sub.reset_index(drop = True, inplace = True)
    
    #assign the narrow decoy PSM scores to the open decoy PSMs that are found in the narrow search
    if score == 'tailor':
        decoys_open.loc[decoys_open['scan_plus_seq'].isin(unique_scan_seq_td_narrow), 'tailor_score'] = list(target_decoys_narrow_sub.sort_values(['scan_plus_seq'])['tailor_score'])
    elif score == 'xcorr':
        decoys_open.loc[decoys_open['scan_plus_seq'].isin(unique_scan_seq_td_narrow), 'xcorr_score'] = list(target_decoys_narrow_sub.sort_values(['scan_plus_seq'])['xcorr_score'])
    elif score == 'hyper':
        decoys_open.loc[decoys_open['scan_plus_seq'].isin(unique_scan_seq_td_narrow), 'hyperscore'] = list(target_decoys_narrow_sub.sort_values(['scan_plus_seq'])['hyperscore'])
    
    #take subset of narrow PSMs that contain the target PSMs in the open
    #reorder this set so that they match the target PSMs in the open
    target_decoys_narrow_sub = target_decoys_narrow_TDC[target_decoys_narrow_TDC['scan_plus_seq'].isin(unique_scan_plus_seq_targets_sub)].copy()
    target_decoys_narrow_sub= target_decoys_narrow_sub.reset_index(drop = True)
    target_decoys_narrow_sub.scan_plus_seq = target_decoys_narrow_sub.scan_plus_seq.astype("category")
    target_decoys_narrow_sub.scan_plus_seq.cat.set_categories(scan_plus_seq_targets_sub, inplace = True)
    target_decoys_narrow_sub.reset_index(drop = True, inplace = True)
    
    #assign the narrow target PSM scores to the open target PSMs that are found in the narrow search 
    if score == 'tailor':
        targets_open.loc[targets_open['scan_plus_seq'].isin(unique_scan_seq_td_narrow), 'tailor_score'] = list(target_decoys_narrow_sub.sort_values(['scan_plus_seq'])['tailor_score'])
    elif score == 'xcorr':
        targets_open.loc[targets_open['scan_plus_seq'].isin(unique_scan_seq_td_narrow), 'xcorr_score'] = list(target_decoys_narrow_sub.sort_values(['scan_plus_seq'])['xcorr_score'])
    elif score == 'hyper':
        targets_open.loc[targets_open['scan_plus_seq'].isin(unique_scan_seq_td_narrow), 'hyperscore'] = list(target_decoys_narrow_sub.sort_values(['scan_plus_seq'])['hyperscore'])
        
        
    #concatenate the open and narrow search PSMs and delete any remaining duplicated sequences (up to variable modification if account_mods == True)
    #narrow takes precedence always, regardless of score
    decoys = pd.concat([decoys_narrow, decoys_open])
    targets = pd.concat([targets_narrow, targets_open])
    
    decoys = decoys.sort_values(by=['database'])
    targets = targets.sort_values(by=['database'])
    
    if account_mods:
        if tide_used:            
            decoys = decoys.drop_duplicates(subset = ['original_target_sequence'])
            targets = targets.drop_duplicates(subset = ['original_target_sequence'])
        else:
            decoys = decoys.drop_duplicates(subset = ['peptide'])
            targets = targets.drop_duplicates(subset = ['peptide'])
    else:
        decoys = decoys.drop_duplicates(subset = ['sequence'])
        targets = targets.drop_duplicates(subset = ['sequence'])

    peptide_list = peptide_list.drop_duplicates()
    peptide_list.reset_index(drop = True, inplace = True)
    
    #give the indices corresponding the index in peptide_list
    #then identify indices that have both target and decoy pairs, indices that have one
    #for indices that have one, produce NA rows
    decoys.sequence = decoys.sequence.astype("category")
    decoys.sequence.cat.set_categories(list(peptide_list['decoy'].dropna()), inplace = True)
    decoys = decoys.sort_values(['sequence'])    
    
    targets.sequence = targets.sequence.astype("category")
    targets.sequence.cat.set_categories(list(peptide_list['target'].dropna()), inplace = True)
    targets = targets.sort_values(['sequence'])
    
    
    decoy_inds = peptide_list.index[peptide_list['decoy'].isin(decoys.sequence)]
    target_inds = peptide_list.index[peptide_list['target'].isin(targets.sequence)]
    pair_inds = list(set(decoy_inds) & set(target_inds))
    decoy_only_inds = list(set(decoy_inds) - set(pair_inds))
    target_only_inds = list(set(target_inds) - set(pair_inds))
    
    decoys.index = decoy_inds
    targets.index = target_inds
    
    additional_decoy_rows = pd.DataFrame(np.nan, index = target_only_inds, columns = decoys.columns)
    additional_target_rows = pd.DataFrame(np.nan, index = decoy_only_inds, columns = targets.columns)
    
    #do head-to-head competition.
    #if one of the sequences in the pair does not exist, it gets -Infinity score
    if score == 'tailor':
        additional_decoy_rows['tailor_score'] = -float('inf') 
        additional_target_rows['tailor_score'] = -float('inf') 
    elif score == "xcorr":
        additional_decoy_rows['xcorr_score'] = -float('inf') 
        additional_target_rows['xcorr_score'] = -float('inf') 
    elif score == "hyper":
        additional_decoy_rows['hyperscore'] = -float('inf') 
        additional_target_rows['hyperscore'] = -float('inf') 
    
    decoys = decoys.append(additional_decoy_rows)
    targets = targets.append(additional_target_rows)
    
    decoys = decoys.sort_index()
    targets = targets.sort_index()
    
    targets.reset_index(drop = True, inplace = True)
    decoys.reset_index(drop = True, inplace = True)
    
    different_databases = ~(decoys['database'] == targets['database']) & ~targets['database'].isnull() & ~decoys['database'].isnull()
    
    if score == 'tailor':
        targets.loc[different_databases & (targets['database'] == "open"),'tailor_score'] = -float('inf')
        decoys.loc[different_databases & (decoys['database'] == "open"), 'tailor_score'] = -float('inf')
    elif score == "xcorr":
        targets.loc[different_databases & (targets['database'] == "open"),'xcorr_score'] = -float('inf')
        decoys.loc[different_databases & (decoys['database'] == "open"), 'xcorr_score'] = -float('inf')
    elif score == "hyper":
        targets.loc[different_databases & (targets['database'] == "open"),'hyperscore'] = -float('inf')
        decoys.loc[different_databases & (decoys['database'] == "open"), 'hyperscore'] = -float('inf')
    
    if score == 'tailor':
        decoys_scores = list(decoys['tailor_score'])
        targets_scores = list(targets['tailor_score'])
    elif score == 'xcorr':
        decoys_scores = list(decoys['xcorr_score'])
        targets_scores = list(targets['xcorr_score'])
    elif score == 'hyper':
        decoys_scores = list(decoys['hyperscore'])
        targets_scores = list(targets['hyperscore'])
    
    
    #Determining indices of winning target peptides, losing target peptides, and ties
    Ws = [i for i in range(len(decoys_scores)) if decoys_scores[i] < targets_scores[i]]
    ties = [i for i in range(len(decoys_scores)) if decoys_scores[i] == targets_scores[i]]
    #break ties randomly or ignore them
    if random_h2h:
        Ws_inds = random.choices([False, True], k = len(ties))
        Ws += [ties[i] for i in range(len(ties)) if Ws_inds[i]]
        Ds = list(set(range(len(decoys_scores))) - set(Ws))
    else:
        Ds = list(set(range(len(decoys_scores))) - set(Ws) - set(ties))
    winning_scores = [max(l1, l2) for l1, l2 in zip(decoys_scores, targets_scores)]
    
    labels = pd.Series([0]*len(winning_scores))
    labels.loc[Ws] = 1
    labels[Ds] = -1
    
    #targets = targets.reset_index(drop = True)
    #decoys = decoys.reset_index(drop = True)
    
    winning_peptides = pd.Series([0]*len(winning_scores))
    winning_peptides[Ws] = targets['sequence'][Ws]
    winning_peptides[Ds] = decoys['sequence'][Ds]
    
    rank = pd.Series([0]*len(winning_scores))
    if tide_used:
        rank[Ws] = targets['xcorr_rank'][Ws]
        rank[Ds] = decoys['xcorr_rank'][Ds]
    else:
        rank[Ws] = targets['hit_rank'][Ws]
        rank[Ds] = decoys['hit_rank'][Ds]
    
    delta_mass = pd.Series([0]*len(winning_scores))
    if tide_used:
        delta_mass[Ws] = targets['spectrum_neutral_mass'][Ws] - targets['peptide_mass'][Ws]
        delta_mass[Ds] = decoys['spectrum_neutral_mass'][Ds] - decoys['peptide_mass'][Ds]
    else:
        delta_mass[Ws] = targets['precursor_neutral_mass'][Ws] - targets['calc_neutral_pep_mass'][Ws]
        delta_mass[Ds] = decoys['precursor_neutral_mass'][Ds] - decoys['calc_neutral_pep_mass'][Ds]
    
    database = pd.Series([0]*len(winning_scores))
    database[Ws] = targets['database'][Ws]
    database[Ds] = decoys['database'][Ds]
    
    scan = pd.Series([0]*len(winning_scores))
    if tide_used:
        scan[Ws] = targets['scan'][Ws]
        scan[Ds] = decoys['scan'][Ds]
    else:
        scan[Ws] = targets['scannum'][Ws]
        scan[Ds] = decoys['scannum'][Ds]
    
    #collect information regarding winning peptide
    df = pd.DataFrame(zip(winning_scores, labels, delta_mass, winning_peptides, rank, database, scan), columns = ['winning_scores', 'labels', 'delta_mass', 'winning_peptides', 'rank', 'database', 'scan'])
    df = df.loc[random.sample(range(df.shape[0]), df.shape[0])]
    
    #creating bins for the mass differences
    breaks_p = np.arange(0, 100 + 2*precursor_bin_width, precursor_bin_width) - precursor_bin_width/2
    breaks_n = list(reversed(-breaks_p))
    breaks = pd.Series(breaks_n[0:(len(breaks_n) - 1)] + list(breaks_p), name = 'bins')
    digitized = np.digitize(df['delta_mass'], breaks)
    df['bins'] = digitized
    df.reset_index(drop = True, inplace = True)
    
    logging.info("Head to head competition is complete.")
    sys.stderr.write("Head to head competition is complete.\n")
    logging.info("Constructing groups.")
    sys.stderr.write("Constructing groups.\n")
    
    #reorder bins according to most frequence
    #bin that includes the narrow is always the top bin (indeed open search has PSMs whose mass difference
    #are close to the narrow bin-width that include a majority of target-PSMs)
    bin_freq =df['bins'].value_counts()
    n_bin_freq = bin_freq.shape[0]
    bin_freq = bin_freq.reset_index()
    bin_freq.columns = ['bin', 'freq']
    
    if adaptive:
        logging.info("Constructing groups adaptively.")
        sys.stderr.write("Constructing groups adaptively.\n")
        #Now to create the adaptive group structure
        all_group_ids = pd.Series([0]*len(winning_scores))
        all_group_ids[df['database'] == 'narrow'] = 'narrow'
        #PSMs with xcorr_rank of 2 or more get "bundled together"
        rank = 1
        rank_val = [1]
        rank_lab = '1 rank'
        i = 0
        #create a candidate group based off the most frequent bins
        #make sure it's at least 2K in size by joining subsequent bins
        #compare this to the left over set of PSMs in the remaining bin
        #if sufficiently different to this left over, along with any existing groups, then create a new group
        #if not different to this left over, terminate
        #if similar to an existing group, join to this group
        largest_bin_considered = 0
        while i < n_bin_freq - 1:
            cand_group = df[(df['bins'] == bin_freq['bin'].loc[i]) & df['rank'].isin(rank_val) & (df['database'] == "open")].index
            while ((len(cand_group) <= min_group_size*K) and (i < n_bin_freq - 1)):
                i += 1
                cand_group = cand_group.append(df[(df['bins'] == bin_freq['bin'].loc[i]) & df['rank'].isin(rank_val) & (df['database'] == "open")].index)
            if (i >= n_bin_freq - 1): break
            
            else_group = df[df['rank'].isin(rank_val) & (all_group_ids == 0)].index.difference(cand_group)
        
            if len(else_group) == 0: break
        
            test_first = stats.ks_2samp(df['winning_scores'].loc[cand_group], df['winning_scores'].loc[else_group], mode = 'asymp')
        
            if test_first.pvalue <= group_thresh:
                unique_gs = all_group_ids.unique()
                unique_gs = set(unique_gs) - set([0, 'narrow'])
                unique_gs = [unique_g for unique_g in unique_gs if rank_lab in unique_g]
                unique_gs = random.sample(unique_gs, len(unique_gs))
                if len(unique_gs) == 0:
                    all_group_ids[cand_group] = rank_lab + ' ' + str(i)
                    largest_bin_considered = i
                    i += 1
                    continue
                else:
                    test_pvals = []
                    for g in unique_gs:
                        other_group_scores = df['winning_scores'][(all_group_ids == g)]
                        test = stats.ks_2samp(df['winning_scores'].loc[cand_group], other_group_scores, mode = 'asymp')
                        test_pvals.append(test.pvalue)
                        if test.pvalue > group_thresh:
                            all_group_ids[cand_group] = g
                            break
                    if all([test_pval <= group_thresh for test_pval in test_pvals]):
                        all_group_ids[cand_group] = rank_lab + ' ' + str(i)
                        largest_bin_considered = i
                    i += 1
            else:
                break
            
        all_group_ids[(all_group_ids == 0) & (df['rank'] == 1)] = "1 rank left over"
        #For PSMs of rank 2 or higher, we use the bins considered determined by rank 1 PSMs and bundle them together in one group
        #If they are not sufficiently many of them, we throw them into the "1 rank left over" group
        if (tops > 1) and (largest_bin_considered > 0):
            inds_higher_rank = ( df['bins'].isin(bin_freq['bin'].loc[0:largest_bin_considered]) ) &  ( df['rank'].isin(range(2, tops + 1)) ) & (df['database'] == "open")
            if sum(inds_higher_rank) <= min_group_size*K:
                all_group_ids[inds_higher_rank] = "1 rank left over"
            else:
                all_group_ids[inds_higher_rank] = "top for rank 2 or more"
            
        
        df['all_group_ids'] = all_group_ids            
    
    else:
        logging.info("Constructing groups non-adaptively using the given top mass-differences.")
        sys.stderr.write("Constructing groups non-adaptively using the given top mass-differences.\n")
        #uses pre-fixed number of top mass bins to use when creating the group
        #n_top_groups = 4 means that we will consider the top 4 bins, and 5th "left over group".
        all_group_ids = pd.Series([0]*len(winning_scores))
        all_group_ids[df['database'] == 'narrow'] = 'narrow'
        
        for rank in range(1, 3):
            if tops == 1 and rank == 2:
                break
            if rank == 1:
                rank_val = [1]
                rank_lab = '1 rank'
            if rank == 2:
                rank_val = range(2, tops + 1)
                rank_lab = '2 rank'
                    
            for i in range(1, min(n_top_groups, bin_freq.shape[0])):
                all_group_ids[(df['bins'] == bin_freq['bin'].loc[i]) & df['rank'].isin(rank_val) & (df['database'] == "open")] = rank_lab + ' ' + str(i)
        
        all_group_ids[(all_group_ids == 0) & df['rank'].isin([1]) & (df['database'] == "open")] = "1 rank left over"
        if tops > 1:
            all_group_ids[(all_group_ids == 0) & df['rank'].isin(range(2, tops + 1)) & (df['database'] == "open")] = "2 rank left over"
        
        df['all_group_ids'] = all_group_ids
    
    df = df[df['all_group_ids'] != 0]
    df.reset_index(drop = True, inplace = True)
    if print_group_pi0:
        label_all_groups = df.groupby(['labels', 'all_group_ids']).size().reset_index(name = "val")
        label_all_groups = label_all_groups.set_index(['all_group_ids', 'labels']).unstack()
        label_all_groups = label_all_groups.iloc[:,0].div(label_all_groups.iloc[:, 1])
        sys.stderr.write(label_all_groups.to_string() + "\n")
        logging.info(label_all_groups.to_string())
    
    logging.info("Constructing groups is complete.")
    sys.stderr.write("Constructing groups is complete.\n")   
    
    return df[['winning_scores', 'labels', 'all_group_ids', 'winning_peptides', 'scan', 'database']]
        

###############################################################################     
def add_modification_to_amino_acid(df, static_mods):
    '''
    Parameters
    ----------
    df : Pandas Dataframe
        Pandas dataframe containing a list of amino acids
        and a list containing modification sites and the
        modification mass.

    Returns
    -------
    A Pandas series containing the amino acids and their
    variable mass modifications. (Static mods ignored).

    '''
    aa_list = df[0]
    mods = df[1]
    
    if type(mods) == float:
        return aa_list
    else:
        for i in np.arange(0, len(mods), 2):
            if (i == 0) & ("nterm" in static_mods):
                continue
            if (i == len(mods) - 2) & ("cterm" in static_mods):
                continue
            
            if mods[i] == "N-term":
                site = 0
            elif mods[i] == "C-term":
                site = len(aa_list) - 1
            else:
                site = int(mods[i]) - 1
            
            mass = mods[i + 1]
            if aa_list[site] in static_mods:
                if abs(static_mods[aa_list[site]] - float(mass)) <= 10**(-4):
                    continue
            aa_list[site] = aa_list[site] + '[' + mass + ']'
            
        return(aa_list)

###############################################################################
def create_peptides_with_mod(unmodified_peptides, modification_info, static_mods):
    '''
    Parameters
    ----------
    unmodified_peptides : Pandas Series
        Unmodified peptides from MS-fragger.
    modification_info : Pandas Series
        Modifications for peptides from MS-fragger.

    Returns
    -------
    Modified peptides in the format of tide-search.

    '''
    modification_info = modification_info.str.findall('-?\d+\.\d+|-?\d+|N\\-term|C\\-term')
    unmodified_peptides = unmodified_peptides.apply(list)
    df = pd.DataFrame({'unmodified_peptides': unmodified_peptides, 'modification_info':modification_info})
    modified_peptides = df.apply(lambda x: add_modification_to_amino_acid(x, static_mods), axis = 1)
    modified_peptides = modified_peptides.apply(''.join)
    return(modified_peptides)
    

###############################################################################
def main():
    global USAGE
    
    
    
    # Set default values for parameters.
    K = 40
    tops_gw = 2
    tops_open = 5
    score = 'tailor'
    random_h2h = True
    account_mods = True
    #print_scan_mult = True
    precursor_bin_width = 1.0005079/4
    adaptive = True
    print_chimera = True
    group_thresh = 0.01
    print_group_pi0 = True
    min_group_size = 2
    n_top_groups = 4
    neighbour_remove = True
    thresh = 0.25
    concat = True
    return_filt_search = False
    correction = 1
    n_processes = 1
    return_frontier = False
    output_dir = './'
    file_name = 'group_walk_results.txt'
    frontier_name = 'frontier_results.txt'
    static_mods = {'C':57.02146}
    
    
    
    
    # Parse the command line.
    sys.argv = sys.argv[1:]
    while (len(sys.argv) > 3):
        next_arg = sys.argv[0]
        sys.argv = sys.argv[1:]
        if (next_arg == "--K"):
            K = int(sys.argv[0])
            sys.argv = sys.argv[1:]
        elif (next_arg == "--tops_open"):
            tops_open = int(sys.argv[0])
            sys.argv = sys.argv[1:]
        elif (next_arg == "--tops_gw"):
            tops_gw = int(sys.argv[0])
            sys.argv = sys.argv[1:]
        elif (next_arg == "--score"):
            score = sys.argv[0]
            sys.argv = sys.argv[1:]
        elif (next_arg == "--random_h2h"):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                random_h2h = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                random_h2h = False
            else:
                sys.stderr.write("Invalid argument for --random_h2h")
                sys.exit(1)
            sys.argv = sys.argv[1:]
        elif (next_arg == "--account_mods"):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                account_mods = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                account_mods = False
            else:
                sys.stderr.write("Invalid argument for --account_mods")
                sys.exit(1)
            sys.argv = sys.argv[1:]
        elif (next_arg == "--precursor_bin_width"):
            precursor_bin_width = float(sys.argv[0])  
            sys.argv = sys.argv[1:]
        elif (next_arg == "--adaptive"):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                adaptive = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                adaptive = False
            else:
                sys.stderr.write("Invalid argument for --adaptive")
                sys.exit(1)
            sys.argv = sys.argv[1:]
        elif (next_arg == "--print_chimera"):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                print_chimera = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                print_chimera = False
            else:
                sys.stderr.write("Invalid argument for --print_chimera")
                sys.exit(1)
            sys.argv = sys.argv[1:]
        elif (next_arg == "--group_thresh"):
            group_thresh = float(sys.argv[0])  
            sys.argv = sys.argv[1:]
        elif (next_arg == "--print_group_pi0"):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                print_group_pi0 = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                print_group_pi0 = False
            else:
                sys.stderr.write("Invalid argument for --print_group_pi0")
                sys.exit(1)
            sys.argv = sys.argv[1:]
        elif (next_arg == "--min_group_size"):
            min_group_size = int(sys.argv[0])  
            sys.argv = sys.argv[1:]
        elif (next_arg == "--n_top_groups"):
            n_top_groups = int(sys.argv[0])  
            sys.argv = sys.argv[1:]
        elif (next_arg == "--neighbour_remove"):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                neighbour_remove = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                neighbour_remove = False
            else:
                sys.stderr.write("Invalid argument for --neighbour_remove")
                sys.exit(1)
            sys.argv = sys.argv[1:]
        elif (next_arg == "--thresh"):
            thresh = float(sys.argv[0])  
            sys.argv = sys.argv[1:]
        elif (next_arg == "--concat"):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                concat = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                concat = False
            else:
                sys.stderr.write("Invalid argument for --concat")
                sys.exit(1)
            sys.argv = sys.argv[1:]
        elif (next_arg == "--return_filt_search"):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                return_filt_search = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                return_filt_search = False
            else:
                sys.stderr.write("Invalid argument for --return_filt_search")
                sys.exit(1)
            sys.argv = sys.argv[1:]
        elif (next_arg == "--n_processes"):
            n_processes = int(sys.argv[0])  
            sys.argv = sys.argv[1:]
        elif (next_arg == "--return_frontier"):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                return_frontier = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                return_frontier = False
            else:
                sys.stderr.write("Invalid argument for --return_frontier")
                sys.exit(1)
            sys.argv = sys.argv[1:]
        elif (next_arg == "--output_dir"):
            output_dir = str(sys.argv[0])  
            sys.argv = sys.argv[1:]
        elif (next_arg == "--file_name"):
            file_name = str(sys.argv[0])  
            sys.argv = sys.argv[1:]
        elif (next_arg == "--frontier_name"):
            frontier_name = str(sys.argv[0])  
            sys.argv = sys.argv[1:]
        elif (next_arg == "--static_mods"):
            static_mods = parse_static_mods(sys.argv[0])
            sys.argv = sys.argv[1:]
        else:
            sys.stderr.write("Invalid option (%s).c" % next_arg)
            sys.exit(1)
    if (len(sys.argv) != 3):
        sys.stderr.write(USAGE)
        sys.exit(1)
    search_file_narrow = sys.argv[0]
    search_file_open = sys.argv[1]
    td_list = sys.argv[2]
    
    logging.basicConfig(filename=output_dir + "/" + file_name.replace('.txt', '.log'), level=logging.DEBUG)
    
    sys.stderr.write("Successfully read in arguments. \n")
    logging.info("Successfully read in arguments")
    
    if score == 'tailor' or score == 'xcorr':
        tide_used = True
    else:
        tide_used = False
        sys.stderr.write("MS-Fragger output files detected. \n")
        logging.info("MS-Fragger output files detected.")

    if concat:
        sys.stderr.write("Reading search files. \n")
        logging.info("Reading search files.")
        narrow_target_decoys = pd.read_table(search_file_narrow)
        open_target_decoys = pd.read_table(search_file_open)
        narrow_target_decoys.columns = narrow_target_decoys.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '').str.replace('/', '_').str.replace('.', '_')
        open_target_decoys.columns = open_target_decoys.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '').str.replace('/', '_').str.replace('.', '_')

        if tide_used:
            if 'hyperscore' in narrow_target_decoys.columns:
                logging.error("Please set --score hyper")
                sys.exit("Please set --score hyper.\n")
            narrow_target_decoys = narrow_target_decoys[narrow_target_decoys['xcorr_rank'] == 1]
            open_target_decoys = open_target_decoys[open_target_decoys['xcorr_rank'].isin(range(1, tops_open + 1))]
        else:
            narrow_target_decoys = narrow_target_decoys[narrow_target_decoys['hit_rank'] == 1]
            open_target_decoys = open_target_decoys[open_target_decoys['hit_rank'].isin(range(1, tops_open + 1))]
        
        narrow_target_decoys.reset_index(drop = True, inplace = True)
        open_target_decoys.reset_index(drop = True, inplace = True)
        
    else:
        narrow_1 = pd.read_table(search_file_narrow)
        open_1 = pd.read_table(search_file_open)
        narrow_2 = pd.read_table(search_file_narrow.replace("target", "decoy"))
        open_2 = pd.read_table(search_file_open.replace("target", "decoy"))
        
        narrow_1.columns = narrow_1.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '').str.replace('/', '_').str.replace('.', '_')
        narrow_2.columns = narrow_2.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '').str.replace('/', '_').str.replace('.', '_')
        open_1.columns = open_1.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '').str.replace('/', '_').str.replace('.', '_')
        open_2.columns = open_2.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '').str.replace('/', '_').str.replace('.', '_')
        
        
        
        #No need for MS-Fragger output as peptide sequence is given without modification
        if tide_used:
            if 'hyperscore' in narrow_target_decoys.columns:
                logging.error("Please set --score hyper")
                sys.exit("Please set --score hyper.\n")
            narrow_1['original_target_sequence'] = narrow_1['sequence'].str.replace("\\[|\\]|\\.|\d+", "")
            open_1['original_target_sequence'] = open_1['sequence'].str.replace("\\[|\\]|\\.|\d+", "")
            
            
            
        if not len(narrow_1.columns) == len(narrow_2.columns):
            logging.error("search files have different number of columns.")
            sys.exit("search files have different number of columns.\n")
        if not len(open_1.columns) == len(open_2.columns):
            logging.error("search files have different number of columns.")
            sys.exit("search files have different number of columns.\n")
        if not len(narrow_1.columns.intersection(narrow_2.columns)) == narrow_1.columns.shape[0]:
            logging.error("search files do not share the same column names.")
            sys.exit("search files do not share the same column names.\n")
        if not len(open_1.columns.intersection(open_2.columns)) == open_1.columns.shape[0]:
            logging.error("search files do not share the same column names.")
            sys.exit("search files do not share the same column names.\n")
        
        sys.stderr.write("concatenating separate searches.\n")
        logging.info("concatenating separate searches.")
        
        narrow_target_decoys = pd.concat([narrow_1, narrow_2])
        open_target_decoys = pd.concat([open_1, open_2])
        
        if tide_used:
            narrow_target_decoys = narrow_target_decoys.sort_values(by=['scan', 'xcorr_score'], ascending=False)
            open_target_decoys = open_target_decoys.sort_values(by=['scan', 'xcorr_score'], ascending=False)
        else:
            narrow_target_decoys = narrow_target_decoys.sort_values(by=['scannum', 'hyperscore'], ascending=False)
            open_target_decoys = open_target_decoys.sort_values(by=['scannum', 'hyperscore'], ascending=False)
            
        narrow_target_decoys = narrow_target_decoys.reset_index(drop = True)
        open_target_decoys = open_target_decoys.reset_index(drop = True)
        
        if tide_used:
            narrow_target_decoys = narrow_target_decoys.assign(xcorr_rank = narrow_target_decoys.groupby('scan').scan.transform(lambda x: range(1, len(x) + 1)))
            open_target_decoys = open_target_decoys.assign(xcorr_rank = open_target_decoys.groupby('scan').scan.transform(lambda x: range(1, len(x) + 1)))
            narrow_target_decoys = narrow_target_decoys[narrow_target_decoys['xcorr_rank'] == 1]
            open_target_decoys = open_target_decoys[open_target_decoys['xcorr_rank'].isin(range(1, tops_open + 1))]
        else:
            narrow_target_decoys = narrow_target_decoys.assign(xcorr_rank = narrow_target_decoys.groupby('scannum').scan.transform(lambda x: range(1, len(x) + 1)))
            open_target_decoys = open_target_decoys.assign(xcorr_rank = open_target_decoys.groupby('scannum').scan.transform(lambda x: range(1, len(x) + 1)))
            narrow_target_decoys = narrow_target_decoys[narrow_target_decoys['hit_rank'] == 1]
            open_target_decoys = open_target_decoys[open_target_decoys['hit_rank'].isin(range(1, tops_open + 1))]
            

        narrow_target_decoys.reset_index(drop = True, inplace = True)
        open_target_decoys.reset_index(drop = True, inplace = True)
        sys.stderr.write("concatenating separate searches complete.\n")
        logging.info("concatenating separate searches complete.")
        
    if not len(narrow_target_decoys.columns) == len(open_target_decoys.columns):
        logging.error("concatenating separate searches.")
        sys.exit("search files have different number of columns.\n")
    if not all([narrow_target_decoys.columns[i] == open_target_decoys.columns[i] for i in range(len(narrow_target_decoys.columns))]):
        logging.error("search files do not share the same column names.")
        sys.exit("search files do not share the same column names.\n")
    
    if not tide_used:
        sys.stderr.write("Rewriting peptide sequence in tide-search format. \n")
        logging.info("Rewriting peptide sequence in tide-search format.")
        narrow_target_decoys['sequence'] = create_peptides_with_mod(narrow_target_decoys['peptide'], narrow_target_decoys['modification_info'], static_mods)
        open_target_decoys['sequence'] = create_peptides_with_mod(open_target_decoys['peptide'], open_target_decoys['modification_info'], static_mods)
        logging.info("Rewriting peptide sequence in tide-search format complete.")
        sys.stderr.write("Rewriting peptide sequence in tide-search format complete. \n")
        
    target_decoys_all = filter_narrow_open(narrow_target_decoys, open_target_decoys, tops_open, thresh, n_processes, neighbour_remove, tide_used, static_mods)
    
    if return_filt_search:
        logging.info("Returning filtered search files in output directory.")
        sys.stderr.write("Returning filtered search files in output directory. \n")
        target_decoys_all_narrow = target_decoys_all[target_decoys_all['database'] == "narrow"]
        target_decoys_all_open = target_decoys_all[target_decoys_all['database'] == "open"]
        target_decoys_all_narrow.to_csv(output_dir + "/" + "filtered_" + search_file_narrow, header=True, index = False, sep = '\t')
        target_decoys_all_open.to_csv(output_dir + "/" + "filtered_" + search_file_open, header=True, index = False, sep = '\t')
        
    logging.info("Reading peptide list.")
    sys.stderr.write("Reading peptide list. \n")
    peptide_list = pd.read_table(td_list)
    df = create_groups(target_decoys_all, narrow_target_decoys, peptide_list, K, tops_gw, score, random_h2h, account_mods, precursor_bin_width, group_thresh, adaptive, min_group_size, n_top_groups, tide_used, print_group_pi0)
    results = group_walk(list(df['winning_scores']), list(df['labels']), list(df['all_group_ids']), K, return_frontier, correction)
    df['q_vals'] = results[0]
    
    power_1 = sum( ( df['q_vals'] <= 0.01 ) & ( df['labels'] == 1) )
    power_5 = sum( ( df['q_vals'] <= 0.05 ) & ( df['labels'] == 1) )
    
    logging.info(str(power_1) + " peptides discovered at the 1% FDR level.")
    logging.info(str(power_5) + " peptides discovered at the 5% FDR level.")
    sys.stderr.write(str(power_1) + " peptides discovered at the 1% FDR level. \n")
    sys.stderr.write(str(power_5) + " peptides discovered at the 5% FDR level. \n")
    
    if print_chimera:
        scan_mult1 = df['scan'][ ( df['q_vals'] <= 0.01 ) & ( df['labels'] == 1) ].value_counts().value_counts()
        logging.info("Printing the scan multiplicity among the discovered peptides at 1% FDR level:")
        logging.info(scan_mult1.to_string())
        sys.stderr.write("Printing the scan multiplicity among the discovered peptides at 1% FDR level: \n")
        sys.stderr.write(scan_mult1.to_string() + "\n")
        
        scan_mult5 = df['scan'][ ( df['q_vals'] <= 0.05 ) & ( df['labels'] == 1) ].value_counts().value_counts()
        logging.info("Printing the scan multiplicity among the discovered peptides at 5% FDR level:")
        logging.info(scan_mult5.to_string())
        sys.stderr.write("Printing the scan multiplicity among the discovered peptides at 5% FDR level: \n")
        sys.stderr.write(scan_mult5.to_string() + "\n")
        
    if output_dir != './':
        if os.path.isdir(output_dir):
            df.to_csv(output_dir + "/" + file_name, header=True, index = False, sep = '\t')
        else:
            os.mkdir(output_dir)
            df.to_csv(output_dir + "/" + file_name, header=True, index = False, sep = '\t')
    if return_frontier:
        results[1].to_csv(output_dir + "/" + frontier_name, header=True, index = False, sep = '\t')

if __name__ == "__main__":
    main()    
    
    
