#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 15:24:39 2023

@author: jackfreestone
"""
import sys
import pandas as pd
import random
import numpy as np
from scipy import stats
from . import peptides
import multiprocessing
import logging
from tqdm import tqdm
import re

###############################################################################
aa_table = pd.DataFrame({
    71.037114: 'A',
   103.009184: 'C',
   115.026943: 'D',
   129.042593: 'E',
   147.068414: 'F',
    57.021464: 'G',
   137.058912: 'H',
   113.084064: 'I, L',
   128.094963: 'K',
   131.040485: 'M',
   114.042927: 'N',
    97.052764: 'P',
   128.058578: 'Q',
   156.101111: 'R',
    87.032028: 'S',
   101.047668: 'T',
    99.068414: 'V',
   186.079313: 'W',
   163.063329: 'Y',
    }.items(), columns = ['mass', 'aa'])
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
    elif charge2 >= 3:
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
            if (modifications[aaseq_position] != 0.0): #To handle when comet has a mass mod due to n- or c-term AND a variable mod on the same site. MSFragger/Tide do not appear to have this doubling up issue!
                modifications[aaseq_position] += mod_mass
            else:
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
def group_walk(winning_scores, labels, all_group_ids, K = 40, return_frontier = False, correction = 1):
    '''
    Parameters
    ----------
    winning_scores : list
        A list of winning scores associated to each PSM.
    labels : list
        A list of winning labels indicating a target or decoy win associated to each PSM.
    all_group_ids : list
        A list of group ids associated to each PSM.
    K : int, optional
        Window size used to estimate the probability the next PSM is a target or decoy in each group. The default is 40.
    return_frontier : bool, optional
        Whether the entire sequence of frontiers is returned. The default is False.
    correction : int, optional
        The additional term used in the numerator to estimate the FDR. The default is 1.

    Returns
    -------
    None.

    '''
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
    group_ids = pd.Series(sorted(list(set(all_group_ids))))
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
def filter_scan(search_df, thresh = 0.05, frag_bin_size = 0.05, static_mods = {'C':57.02146}):
    '''
    Parameters
    ----------
    search_df : Pandas Dataframe
        Dataframe containing just one scan, and multiple PSMs.
    thresh : float, optional
        The similarity score threshold used to filter neighbouring PSMs. The default is 0.05.
    frag_bin_size : float, optional
        Size of m/z bins used to determine the shared number of b- and y-ions. The default is 0.05.
    static_mods : dic, optional
        Dictionary containing all static modifications. The default is {'C':57.02146}.

    Returns
    -------
    drop_scan : list
        A list of booleans indicating which PSMs should be filtered.

    '''
    
    search_df = search_df.reset_index(drop = True)
    n_scans = search_df.shape[0]
    
    peptide_1 = [search_df['sequence'].loc[0]]
    charge_1 = [search_df['charge'].loc[0]]
    drop_scan = [False]*n_scans
    
    #if any(search_df['sequence'][1:].str.contains(peptide_1[0])):
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
def wrapper(df, q, thresh, frag_bin_size, static_mods): #wrapper is used to update the manager queue every time the wrapper is called
    '''
    a wrapper for filter_scan that can be used for multiprocessing
    '''
    result = filter_scan(df, thresh, frag_bin_size, static_mods)
    q.put(1)
    return(result)
###############################################################################
def filter_scan_subset(df, q, task_number, return_dict, tide_used, score, thresh, static_mods):
    '''
    Effectively calls a wrapper for filter_scan that can be used for multiprocessing
    '''
    if tide_used == 'tide':
        sys.stderr.write("Starting Process " + str(task_number) + " \n")
        logging.info("Starting Process " + str(task_number))
        results = df.groupby(['file', 'scan', "charge", "spectrum_neutral_mass"]).apply(lambda x: wrapper(x, q, thresh, 0.05, static_mods))
        results = results.sort_index(ascending = False)
        results = results.apply(pd.Series).stack().reset_index()
        return_dict[task_number] = results[0]
    elif tide_used == 'comet':
        sys.stderr.write("Starting Process " + str(task_number) + " \n")
        logging.info("Starting Process " + str(task_number))
        results = df.groupby(['scan', "charge", "spectrum_neutral_mass"]).apply(lambda x: wrapper(x, q, thresh, 0.05, static_mods))
        results = results.sort_index(ascending = False)
        results = results.apply(pd.Series).stack().reset_index()
        return_dict[task_number] = results[0]
    else:
        sys.stderr.write("Starting Process " + str(task_number) + " \n")
        logging.info("Starting Process " + str(task_number))
        results = df.groupby(['scannum']).apply(lambda x: wrapper(x, q, thresh, 0.05, static_mods))
        results = results.sort_index(ascending = False)
        results = results.apply(pd.Series).stack().reset_index()
        return_dict[task_number] = results[0]
    
###############################################################################
def listener(q, list_of_df, tide_used):
    '''
    constantly checks to see if the queue q has been updated
    '''
    if tide_used == 'tide':
        pbar = tqdm(total = sum([len(j.groupby(['file', 'scan', "charge", "spectrum_neutral_mass"])) for j in list_of_df]))
        for item in iter(q.get, None):
            pbar.update(item)
    elif tide_used == 'comet':
        pbar = tqdm(total = sum([len(j.groupby(['scan', "charge", "spectrum_neutral_mass"])) for j in list_of_df]))
        for item in iter(q.get, None):
            pbar.update(item)
    else:
        pbar = tqdm(total = sum([len(j.groupby(['scannum'])) for j in list_of_df]))
        for item in iter(q.get, None):
            pbar.update(item)
    

###############################################################################
def filter_narrow_open(narrow_target_decoys, open_target_decoys, score, thresh = 0.05, n_processes = 1, neighbour_remove = True, tide_used = 'tide', static_mods = {'C':57.02146}):
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
        The default is 0.05.
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
    
    
    target_decoys_all = pd.concat([narrow_target_decoys, open_target_decoys]) #combine narrow and open PSMs
    
    #drop duplicate PSMs that are subsequently found in the open search
    if tide_used == 'tide':
        target_decoys_all.drop_duplicates(['file', 'scan', 'charge', 'spectrum_neutral_mass', 'sequence'])
    elif tide_used == 'comet':
        target_decoys_all = target_decoys_all.drop_duplicates(['scan', "charge", "spectrum_neutral_mass", 'sequence'])
    else:
        target_decoys_all = target_decoys_all.drop_duplicates(['scannum', 'sequence'])
   
    
    #makes sure the PSMs from the narrow and open search are ordered by scan first, then by their score
    #indeed it makes sense to use xcorr_score over tailor score, as tailor_score uses candidate PSMs to 
    #normalise the xcorr_score - this differs from narrow to open. 
    if tide_used == 'tide':
        target_decoys_all = target_decoys_all.sort_values(by=['file', 'scan', 'charge', 'spectrum_neutral_mass', 'xcorr_score'], ascending = False)
    elif tide_used == 'comet':
        target_decoys_all = target_decoys_all.sort_values(by=['scan', "charge", "spectrum_neutral_mass", 'xcorr_score'], ascending = False)
    else:
        target_decoys_all = target_decoys_all.sort_values(by=['scannum', 'hyperscore'], ascending = False)
    target_decoys_all.reset_index(drop = True, inplace = True)
    
    if not neighbour_remove:
        logging.info("Not removing neighbours.")
        sys.stderr.write("Not removing neighbours.\n")
        return(target_decoys_all)
    
    logging.info("Filtering for neighbours.")
    sys.stderr.write("Filtering for neighbours.\n")
    
    if n_processes == 1:
        tqdm.pandas()
        if tide_used == 'tide':
            results = target_decoys_all.groupby(['file', 'scan', 'charge', 'spectrum_neutral_mass']).progress_apply(lambda x: filter_scan(x, thresh, 0.05, static_mods)) #apply filtering by experimental scan
        elif tide_used == 'comet':
            results = target_decoys_all.groupby(['scan', "charge", "spectrum_neutral_mass"]).progress_apply(lambda x: filter_scan(x, thresh, 0.05, static_mods)) #apply filtering by experimental scan
        else:
            results = target_decoys_all.groupby(['scannum']).progress_apply(lambda x: filter_scan(x, thresh, 0.05, static_mods))  #apply filtering by experimental scan
        
        results = results.sort_index(ascending = False)
        results = results.apply(pd.Series).stack().reset_index()
        
        target_decoys_all['drop_scan'] = results[0]  #create drop_scan column indicating which PSMs are being kept
    else:
        #if more than 1 thread, we partition the dataframe
        if tide_used == "tide" or tide_used == 'comet':
            target_decoys_all['split_col'] = pd.qcut(target_decoys_all['scan'], n_processes)
        else:
            target_decoys_all['split_col'] = pd.qcut(target_decoys_all['scannum'], n_processes)
        target_decoys_grouped = target_decoys_all.groupby(target_decoys_all.split_col)
        list_of_df = [0]*n_processes #create a list of partitioned dataframes
        for i in range(len(target_decoys_all['split_col'].unique())):
            list_of_df[i] = target_decoys_grouped.get_group(target_decoys_all['split_col'].unique()[i])
            list_of_df[i].reset_index(drop = True, inplace = True)
            
        manager = multiprocessing.Manager()
        q = manager.Queue() #creating instance of manager queue
        return_dict = manager.dict() #creating manager dict variable that can be used to store changes to the argument by ALL processes at the same time
        proc = multiprocessing.Process(target=listener, args=(q, list_of_df, tide_used)) 
        proc.start() #start listening for updates to manager queue
        workers = [multiprocessing.Process(target = filter_scan_subset, args=(list_of_df[i], q, i, return_dict, tide_used, score, thresh, static_mods)) for i in range(n_processes)] #now run each of the processes
        for worker in workers:
            worker.start()
        for worker in workers:
            worker.join()
        q.put(None) 
        proc.join()
        
        for i in range(len(target_decoys_all['split_col'].unique())):
            target_decoys_all.loc[target_decoys_all['split_col'] == target_decoys_all['split_col'].unique()[i], 'drop_scan'] = list(return_dict[i]) #update the main dataframe with the results from each process
            
    target_decoys_all = target_decoys_all[target_decoys_all['drop_scan'] == False] #drop the neighbours
    target_decoys_all.reset_index(drop = True, inplace = True)

    return target_decoys_all
###############################################################################
def create_groups(target_decoys, narrow_target_decoys, peptide_list, dcy_prefix = 'decoy_', K = 40, tops = 2, score = 'tailor_score', account_mods = True, any_mods = True, precursor_bin_width = 1.0005079/4, group_thresh = 0.01, adaptive = True, min_group_size = 2, n_top_groups = 4, tide_used = 'tide', print_group_pi0 = True, competition_window = [2,2]):
    '''
    Parameters
    ----------
    target_decoys : Pandas Dataframe
        Dataframe of target-decoy PSMs from both open/narrow searches after neighbour-filtering
    narrow_target_decoys : Pandas Dataframe
        Dataframe of target-decoy PSMs taken directly from the narrow search file.
    peptide_list : Pandas Dataframe or string
        Dataframe of target-decoy pairs, else an empty string.
    dcy_prefix : string, optional
        Decoy prefix used in the protein database. The default is 'decoy_'.
    K : string, optional
        Window size used to estimate the probability the next PSM is a target or decoy in each group. The default is 40.
    tops : int, optional
        The number of top-ranked PSMs in the open search to use. The default is 2.
    score : string, optional
        Which score type is used. The default is 'tailor'.
    account_mods : bool, optional
        Whether winning PSMs are defined up to variable modification or not. The default is True.
    any_mods : bool, optional
        Whether there exists any variable modifications. The default is True.
    precursor_bin_width : float, optional
        The bin width used to bin the mass differences in open search. The default is 1.0005079/4.
    group_thresh : float, optional
        The p-value threshold used to agglomerate groups. The default is 0.01.
    adaptive : bool, optional
        Whether to use an adaptive grouping scheme according to the KS test or a fixed method. The default is True.
    min_group_size : int, optional
        The number of multiples of K that the used to determine the minimum size of each group. The default is 2.
    n_top_groups : int, optional
        If adaptive = False, then n_top_groups is the number of top mass differences used in creating the groups. The default is 4.
    tide_used : bool, optional
        Whether tide-search is used or not. The default is True.
    print_group_pi0 : bool, optional
        Whether to print the within-group decoy-target ratio to console. The default is True.
    competition_window : list, optional
        Window used to cluster the precursor masses (in m/z) for dynamic level competition. The default is [2,2].

    Returns
    -------
    Dataframe containing the winning PSMs after head-to-head competition and their group labels.

    '''
    #take only the top 1 from the narrow search
    if score == 'xcorr_score' or score == 'tailor_score' or score == 'e-value':
        narrow_target_decoys = narrow_target_decoys[narrow_target_decoys['xcorr_rank'] == 1]
    else:
        narrow_target_decoys = narrow_target_decoys[narrow_target_decoys['hit_rank'] == 1]
    narrow_target_decoys.loc[:, 'database'] = 'narrow'
    
    #rounding in case of any numerical issues
    #cannot apply this to e-values, as they go right up to E-12
    if score != 'e-value':
        target_decoys[score] = round(target_decoys[score], 8)
        narrow_target_decoys[score] = round(narrow_target_decoys[score], 8)
    
    logging.info("Doing dynamic level competition.")
    sys.stderr.write("Doing dynamic level competition.\n")

    #get target_decoy column
    target_decoys["target_decoy"] = "target"
    if type(peptide_list) != str:
        target_decoys.loc[target_decoys['sequence'].isin(peptide_list['target']), "target_decoy"] = "target"
        target_decoys.loc[target_decoys['sequence'].isin(peptide_list['decoy']), "target_decoy"] = "decoy"
    else:
        target_decoys.loc[~(target_decoys['protein_id'].str.contains(dcy_prefix)), "target_decoy"] = "target"
        target_decoys.loc[(target_decoys['protein_id'].str.contains(dcy_prefix)), "target_decoy"] = "decoy"
    
    #ensure correct ranks are considered
    if score == 'xcorr_score' or score == 'tailor_score' or score == 'e-value':
        target_decoys = target_decoys[(target_decoys['xcorr_rank'].isin(range(1, tops + 1)) & (target_decoys['database'] == "open")) | \
                                      ((target_decoys['xcorr_rank'] == 1) & (target_decoys['database'] == "narrow"))]
    else:
        target_decoys = target_decoys[(target_decoys['hit_rank'].isin(range(1, tops + 1)) & (target_decoys['database'] == "open")) | \
                                      ((target_decoys['hit_rank'] == 1) & (target_decoys['database'] == "narrow"))]
    
    #get original_target_sequence for ms_fragger
    if tide_used == 'ms_fragger':
        if account_mods:
            pep = 'peptide'
        else:
            pep = 'sequence'
        
        #create equivalent 'original_target_sequence' column for MS fragger results
        target_decoys['original_target_sequence'] = target_decoys[pep]
        peptide_list.rename(columns = {'target':'original_target_sequence','decoy':pep}, inplace = True)
        target_decoys_sub = target_decoys[target_decoys['target_decoy'] == 'decoy'].copy()
        target_decoys_sub.pop('original_target_sequence')
        target_decoys_sub = target_decoys_sub.merge(peptide_list[['original_target_sequence',pep]], how='left', on=pep)
        target_decoys.loc[target_decoys['target_decoy'] == 'decoy', :] = target_decoys_sub.copy()
    if tide_used == 'tide' and (not account_mods) and any_mods:
        #create equivalent 'original_target_sequence' column for MS fragger results but with modifications
        target_decoys['original_target_sequence'] = target_decoys['sequence']
        peptide_list.rename(columns = {'target':'original_target_sequence','decoy':'sequence'}, inplace = True)
        target_decoys_sub = target_decoys[target_decoys['target_decoy'] == 'decoy'].copy()
        target_decoys_sub.pop('original_target_sequence')
        target_decoys_sub = target_decoys_sub.merge(peptide_list[['original_target_sequence','sequence']], how='left', on=pep)
        target_decoys.loc[target_decoys['target_decoy'] == 'decoy', :] = target_decoys_sub.copy()
    
    #create a cluster column
    target_decoys = target_decoys.sample(frac = 1).reset_index(drop = True)
    if tide_used == 'tide' or tide_used == 'comet':
        #original target_sequence already handles whether it is at the sequence level or modified-sequence level
        target_decoys = target_decoys.sort_values(by='spectrum_neutral_mass', ascending=True).reset_index(drop = True)
        target_decoys["mass_plus"] = target_decoys.groupby('original_target_sequence', group_keys = False).apply(lambda x: x.spectrum_neutral_mass.shift(1) + np.maximum(competition_window[0]*x.charge, competition_window[1]*x.charge.shift(1)))
        target_decoys.loc[target_decoys["mass_plus"].isna(), "mass_plus"] = -np.Inf
        target_decoys["condition"] = target_decoys["spectrum_neutral_mass"] > target_decoys["mass_plus"]
        target_decoys["cluster"] = target_decoys.groupby('original_target_sequence', group_keys = False).condition.cumsum()
    else:
        target_decoys = target_decoys.sort_values(by='precursor_neutral_mass', ascending=True)
        target_decoys["mass_plus"] = target_decoys.groupby('original_target_sequence', group_keys = False).apply(lambda x: x.precursor_neutral_mass.shift(1) + np.maximum(competition_window[0]*x.charge, competition_window[1]*x.charge.shift(1)))
        target_decoys.loc[target_decoys["mass_plus"].isna(), "mass_plus"] = -np.Inf
        target_decoys["condition"] = target_decoys["precursor_neutral_mass"] > target_decoys["mass_plus"]
        target_decoys["cluster"] = target_decoys.groupby('original_target_sequence', group_keys = False).condition.cumsum() 
    
    #now doing h2h competition
    if tide_used == 'tide' or tide_used == 'comet':
        target_decoys['scan_plus_seq'] = target_decoys['scan'].astype(str) + ' ' + target_decoys['charge'].astype(str) + ' ' + target_decoys['spectrum_neutral_mass'].astype(str) + ' ' + target_decoys['sequence']
        narrow_target_decoys.reset_index(drop = True, inplace = True)
        narrow_target_decoys['scan_plus_seq'] = narrow_target_decoys['scan'].astype(str) + ' ' + narrow_target_decoys['charge'].astype(str) + ' ' + narrow_target_decoys['spectrum_neutral_mass'].astype(str) + ' ' + narrow_target_decoys['sequence']
        
        #translating into an iterable
        unique_scan_seq_td_narrow = set(narrow_target_decoys['scan_plus_seq'])
        check_1 = all(target_decoys.database[target_decoys['scan_plus_seq'].isin(unique_scan_seq_td_narrow)] == 'narrow')
        
        if not check_1:
            logging.warning("Some PSMs that are found both in the narrow- and open- search file are labelled as 'open' instead of 'narrow'. Likely a numerical issue.")
            sys.stderr.write("WARNING: Some PSMs that are found both in the narrow- and open- search file are labelled as 'open' instead of 'narrow'. Likely a numerical issue.")
        
        #order according to score
        if score != 'e-value':
            target_decoys = target_decoys.sort_values(by=[score], ascending = False)
        else:
            target_decoys = target_decoys.sort_values(by=[score], ascending = True)    
        
        #take best PSM
        target_decoys = target_decoys.drop_duplicates(subset = ['original_target_sequence', "cluster"])
    
        #gather winning information
        winning_scores = target_decoys[score]
            
        labels = target_decoys['target_decoy'].copy()
        labels[labels == 'target'] = 1
        labels[labels == 'decoy'] = -1
        winning_peptides = target_decoys['sequence']
        rank = target_decoys['xcorr_rank']
        delta_mass = target_decoys['spectrum_neutral_mass'] - target_decoys['peptide_mass']
        database = target_decoys['database']
        scan = target_decoys['scan']
        charge = target_decoys['charge']
        spectrum_neutral_mass = target_decoys['spectrum_neutral_mass']
        protein = target_decoys['protein_id']
        if tide_used == 'tide':
            file = target_decoys['file']
            flanking_aa = target_decoys['flanking_aa']
            df = pd.DataFrame(zip(winning_scores, labels, delta_mass, winning_peptides, rank, database, charge, spectrum_neutral_mass, flanking_aa, scan, protein, file), columns = ['winning_scores', 'labels', 'delta_mass', 'winning_peptides', 'rank', 'database', 'charge', 'spectrum_neutral_mass', 'flanking_aa', 'scan', 'protein', 'file'])
        else:
            df = pd.DataFrame(zip(winning_scores, labels, delta_mass, winning_peptides, rank, database, charge, spectrum_neutral_mass, scan, protein), columns = ['winning_scores', 'labels', 'delta_mass', 'winning_peptides', 'rank', 'database', 'charge', 'spectrum_neutral_mass', 'scan', 'protein'])
    else:
        target_decoys['scan_plus_seq'] = target_decoys['scannum'].astype(str) + ' ' + target_decoys['sequence']
        narrow_target_decoys.reset_index(drop = True, inplace = True)
        narrow_target_decoys['scan_plus_seq'] = narrow_target_decoys['scannum'].astype(str) + ' ' + narrow_target_decoys['sequence']
        
        #translating into an iterable
        unique_scan_seq_td_narrow = set(narrow_target_decoys['scan_plus_seq'])
        check_1 = all(target_decoys.database[target_decoys['scan_plus_seq'].isin(unique_scan_seq_td_narrow)] == 'narrow')
        
        if not check_1:
            logging.warning("Some PSMs that are found both in the narrow- and open- search file are labelled as 'open' instead of 'narrow'. Likely a numerical issue.")
            sys.stderr.write("WARNING: Some PSMs that are found both in the narrow- and open- search file are labelled as 'open' instead of 'narrow'. Likely a numerical issue.")
          
        #do h2h competition by taking the best scoring PSM
        target_decoys = target_decoys.sort_values(by=['hyperscore'], ascending = False)
        
        #take best PSM up to variable modification, (narrow takes precedence) if account_mods
        target_decoys = target_decoys.drop_duplicates(subset = ['original_target_sequence', 'cluster'])
        
        winning_scores = target_decoys['hyperscore']
        
        labels = target_decoys['target_decoy'].copy()
        labels[labels == 'target'] = 1
        labels[labels == 'decoy'] = -1
        winning_peptides = target_decoys['sequence']
        rank = target_decoys['hit_rank']
        delta_mass = target_decoys['precursor_neutral_mass'] - target_decoys['calc_neutral_pep_mass']
        database = target_decoys['database']
        scan = target_decoys['scannum']
        charge = target_decoys['charge']
        spectrum_neutral_mass = target_decoys['precursor_neutral_mass']
        protein = target_decoys['protein_id']
        df = pd.DataFrame(zip(winning_scores, labels, delta_mass, winning_peptides, rank, database, charge, spectrum_neutral_mass, scan, protein), columns = ['winning_scores', 'labels', 'delta_mass', 'winning_peptides', 'rank', 'database', 'charge', 'spectrum_neutral_mass', 'scan', 'protein'])
        
    df = df.sample(frac = 1)
    
    #creating bins for the mass differences
    delta_mass_max = max(abs(df['delta_mass']))
    breaks_p = np.arange(0, delta_mass_max + 2*precursor_bin_width, precursor_bin_width) - precursor_bin_width/2
    breaks_n = list(reversed(-breaks_p))
    breaks = pd.Series(breaks_n[0:(len(breaks_n) - 1)] + list(breaks_p), name = 'bins')
    digitized = np.digitize(df['delta_mass'], breaks)
    df['bins'] = digitized
    df.reset_index(drop = True, inplace = True)
    
    #logging.info("Head to head competition is complete.")
    #sys.stderr.write("Head to head competition is complete.\n")
    
    #reorder bins according to most frequent
    bin_freq =df['bins'][(df['database'] != 'narrow') & (df['rank'] == 1)].value_counts()
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
        rank_lab = 'top 1 PSMs'
        i = 0
        #create a candidate group based off the most frequent bins
        #make sure it's at least 2K in size by joining subsequent bins
        #compare this to the left over set of PSMs in the remaining bin
        #if sufficiently different to this left over, along with any existing groups, then create a new group
        #if not different to this left over, terminate
        #if similar to an existing group, join to this group
        largest_bin_considered = 0
        recent_group = 'recent_group'
        while i < n_bin_freq - 1:
            cand_group = df[(df['bins'] == bin_freq['bin'].loc[i]) & df['rank'].isin(rank_val) & (df['database'] == "open")].index
            while ((len(cand_group) < min_group_size*K) and (i < n_bin_freq - 1)):
                i += 1
                cand_group = cand_group.append(df[(df['bins'] == bin_freq['bin'].loc[i]) & df['rank'].isin(rank_val) & (df['database'] == "open")].index)
            if (i >= n_bin_freq - 1) and (len(cand_group) < min_group_size*K):
                if recent_group == 'recent_group':
                    all_group_ids[cand_group] = rank_lab + ' & top '  + str(i + 1) + ' mass bin'
                else:
                    all_group_ids[cand_group] = recent_group
                break
            elif (i >= n_bin_freq - 1) and (len(cand_group) >= min_group_size*K):
                all_group_ids[cand_group] = "left over group"
                break

            else_group = df[df['rank'].isin(rank_val) & (all_group_ids == 0)].index.difference(cand_group)
        
            if len(else_group) == 0: break
        
            test_first = stats.ks_2samp(df['winning_scores'].loc[cand_group], df['winning_scores'].loc[else_group], mode = 'asymp')
        
            if test_first.pvalue <= group_thresh:
                unique_gs = all_group_ids.unique()
                unique_gs = set(unique_gs) - set([0, 'narrow'])
                unique_gs = [unique_g for unique_g in unique_gs if rank_lab in unique_g]
                unique_gs = random.sample(unique_gs, len(unique_gs))
                #unique_gs.sort()
                if len(unique_gs) == 0:
                    if largest_bin_considered == 0 and i == 0:
                        all_group_ids[cand_group] = rank_lab + ' & top '  + str(i + 1) + ' mass bin'
                        recent_group = rank_lab + ' & top '  + str(i + 1) + ' mass bin'
                    elif largest_bin_considered + 1 == i:
                        all_group_ids[cand_group] = rank_lab + ' & top '  + str(i + 1) + ' mass bin'
                        recent_group = rank_lab + ' & top '  + str(i + 1) + ' mass bin'
                    else:
                        all_group_ids[cand_group] = rank_lab + ' & top ('  + str(largest_bin_considered + 1) + ', ' + str(i + 1) + '] mass bins'
                        recent_group = rank_lab + ' & top ('  + str(largest_bin_considered + 1) + ', ' + str(i + 1) + '] mass bins'
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
                            largest_bin_considered = i
                            break
                    if all([test_pval <= group_thresh for test_pval in test_pvals]):
                        if largest_bin_considered == 0 and i == 0:
                            all_group_ids[cand_group] = rank_lab + ' & top '  + str(i + 1) + ' mass bin'
                            recent_group = rank_lab + ' & top '  + str(i + 1) + ' mass bin'
                        elif largest_bin_considered + 1 == i:
                            all_group_ids[cand_group] = rank_lab + ' & top '  + str(i + 1) + ' mass bin'
                            recent_group = rank_lab + ' & top '  + str(i + 1) + ' mass bin'
                        else:
                            all_group_ids[cand_group] = rank_lab + ' & top ('  + str(largest_bin_considered + 1) + ', ' + str(i + 1) + '] mass bins'
                            recent_group = rank_lab + ' & top ('  + str(largest_bin_considered + 1) + ', ' + str(i + 1) + '] mass bins'
                        largest_bin_considered = i
                    i += 1
            else:
                break
            
        all_group_ids[(all_group_ids == 0) & (df['rank'] == 1)] = "left over group"
        #For PSMs of rank 2 or higher, we use the bins considered determined by rank 1 PSMs and bundle them together in one group
        #If they are not sufficiently many of them, we throw them into the "1 rank left over" group
        if (tops > 1) and (largest_bin_considered > 0):
            inds_higher_rank = ( df['bins'].isin(bin_freq['bin'].loc[0:largest_bin_considered]) ) &  ( df['rank'].isin(range(2, tops + 1)) ) & (df['database'] == "open")
            if sum(inds_higher_rank) < min_group_size*K:
                all_group_ids[inds_higher_rank] = "left over group"
            else:
                all_group_ids[inds_higher_rank] = "top 2 or more PSMs"
            
        
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
                rank_lab = 'top PSMs'
            if rank == 2:
                rank_val = range(2, tops + 1)
                rank_lab = 'top 2 or more PSMs'
                    
            for i in range(1, min(n_top_groups, bin_freq.shape[0])):
                all_group_ids[(df['bins'] == bin_freq['bin'].loc[i]) & df['rank'].isin(rank_val) & (df['database'] == "open")] = rank_lab + ' & top '  + str(i + 1) + ' mass bin'
        
        all_group_ids[(all_group_ids == 0) & df['rank'].isin([1]) & (df['database'] == "open")] = "top PSMs left over group"
        if tops > 1:
            all_group_ids[(all_group_ids == 0) & df['rank'].isin(range(2, tops + 1)) & (df['database'] == "open")] = "top 2 or more PSMs left over group"
        
        df['all_group_ids'] = all_group_ids
    
    df = df[df['all_group_ids'] != 0]
    df.reset_index(drop = True, inplace = True)
    if print_group_pi0:
        label_all_groups = df.groupby(['labels', 'all_group_ids']).size().reset_index(name = "val")
        label_all_groups = label_all_groups.set_index(['all_group_ids', 'labels']).unstack()
        label_all_groups.columns = ['decoys', 'targets']
        label_all_groups['ratio'] = label_all_groups.iloc[:,0].div(label_all_groups.iloc[:, 1])
        label_all_groups.index.names = ['group names:']
        if any(df['all_group_ids'] == 'left over group'):
            if len(label_all_groups.index) >= 2:
                label_all_groups = label_all_groups.reindex(list(label_all_groups.index[1:len(label_all_groups.index)]) + [label_all_groups.index[0]])
        sys.stderr.write(label_all_groups.to_string() + "\n")
        logging.info(label_all_groups.to_string())
    
    #logging.info("Constructing groups is complete.")
    #sys.stderr.write("Constructing groups is complete.\n")   
    
    return(df)
###############################################################################     
def add_modification_to_amino_acid(df, static_mods):
    '''
    Parameters
    ----------
    df : Pandas Dataframe
        Pandas dataframe containing a list of amino acids
        and a list containing modification sites and the
        modification mass.
    static_mods : Dict					      
        A dictionary of user specified static modifications. 

    Returns
    -------
    A Pandas series containing the amino acids and their
    variable mass modifications. (Static mods ignored).

    '''
    aa_list = df[0]
    mods = df[1]
    
    if type(mods) == float:
        return aa_list #to throw away when we have NaN appear
    else:
        for i in np.arange(0, len(mods), 2):
            #If the modification determined in MSFragger coin 
            #then we assume it is a static-modification.      
            if (float(mods[i]) == 1) & ("nterm" in static_mods):
                continue
            if (float(mods[i]) == len(aa_list)) & ("cterm" in static_mods):
                continue
            
            if mods[i] == "N-term": #capital N/C refer to terminals of proteins
                site = 0
            elif mods[i] == "C-term":
                site = len(aa_list) - 1
            else:
                site = int(mods[i]) - 1
            
            mass = mods[i + 1]
            if aa_list[site] in static_mods:
                if abs(static_mods[aa_list[site]] - float(mass)) <= 10**(-4):
                    continue
            aa_list[site] = aa_list[site] + '[' + str(round(float(mass), 6)) + ']'
            
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
    static_mods: Dict					      
        Static modifications given by user.    
    
    Returns
    -------
    Modified peptides in the format of tide-search.

    '''
    modification_info = modification_info.str.findall('-?\\d+.\\d+|-?\\d+')
    unmodified_peptides = unmodified_peptides.apply(list)
    df = pd.DataFrame({'unmodified_peptides': unmodified_peptides, 'modification_info':modification_info})
    modified_peptides = df.apply(lambda x: add_modification_to_amino_acid(x, static_mods), axis = 1)
    modified_peptides = modified_peptides.apply(''.join)
    return(modified_peptides)
###############################################################################
def reverse_sequence(sequence, modifications):
    '''
    Parameters
    ----------
    sequence : Pandas Series
        List of sequence with variable modifications in tide-search format.

    Returns
    -------
    Recovers target sequences.
    '''
    results = re.findall('[a-zA-Z]\\[\\d+.\\d+\\]\\[\\d+.\\d+\\]|[a-zA-Z]\\[\\d+.\\d+\\]|[a-zA-Z]\\[\\d+\\]|[a-zA-Z]', sequence)
    if '_N' in modifications:
        n_term_mass_mod_list = re.findall(r"(\d+.\d+)_N", modifications)
        if len(n_term_mass_mod_list) > 0:
            n_term_mass_mod = float(n_term_mass_mod_list[0])
            masses_on_n_term = re.findall(r"(\d+.\d+)", results[0])
            calc_n_term_mass_mod_list = [float(m) for m in masses_on_n_term if abs(float(m) - n_term_mass_mod) <= 10e-3]
            if len(calc_n_term_mass_mod_list) > 0:
                calc_n_term_mass_mod = str(calc_n_term_mass_mod_list[0])
                results[0] = results[0].replace('[' + calc_n_term_mass_mod + ']', '', 1)
    
    results = results[-2::-1] + [results[-1]]
    
    if '_N' in modifications:
        results[0] = results[0] + '[' + calc_n_term_mass_mod + ']'
    
    results = ''.join(results)
    return(results)
###############################################################################
def check_n_term(sequences):
    '''
    Parameters
    ----------
    sequences : Pandas series
        List of sequence with variable modifications.

    Returns
    -------
    List of sequence with n-terminal modification properly relocated.

    '''
    
    n_term_bool = sequences.str.startswith('[')
    results = sequences[n_term_bool].str.split(pat = r'(\])', n = 1, expand = True)
    
    results = results.apply(lambda x: x[2][0] + x[0][::] + x[1][::] + x[2][1:], axis = 1)
    return(results)
###############################################################################
def del_protein(protein, dcy_prefix):
    '''
    Parameters
    ----------
    protein : Pandas series
        List of proteins that contain the peptides in a search file
    dcy_prefix : String
        String indicating the prefix used for decoy peptides

    Returns
    -------
    Boolean series indicating whether proteins contain both target and decoy

    '''
    
    protein_split = protein.split(',')
    check_protein = [dcy_prefix in protein for protein in protein_split]
    
    select_protein = (True in check_protein) and (False in check_protein)
    
    return(select_protein)
###############################################################################
def get_amino_acid_to_warn(df):
    '''
    Parameters
    ----------
    df : Pandas DataFrame
        Data frame containing delta masses, discovered peptides, and the flanking aa

    Returns
    -------
    A string warning the user that some of the discovered mass-modifications coincide with the
    addition or loss of amino acids -- hence potentially biologically irrelevant PTMs

    '''
    mass_diff = df.delta_mass
    peptide = df.peptide
    flanking = df.flanking_aa
    abs_mass_diff = abs(mass_diff)
    diff_diff = abs_mass_diff - aa_table.mass
    get_indices = aa_table.mass[abs(diff_diff) <= 0.2].index
    if len(get_indices) > 0:
        get_index = get_indices[0]
        get_aa = aa_table.aa[get_index]
        if mass_diff > 0:
            if (get_aa != 'I, L') and (get_aa in flanking):
                message = 'Possible addition of ' + get_aa
            elif (get_aa == 'I, L') and ('I' in flanking or 'L' in flanking):
                message = 'Possible addition of ' + get_aa
            else:
                message = ''
        else:
            if (peptide[0] in get_aa) or (peptide[len(peptide) - 1] in get_aa):
                message = 'Possible loss of ' + get_aa
            else:
                message = ''
        return message
    else:
        return ''
###############################################################################
def get_modification_info(peptide):
    '''
    Parameters
    ----------
    peptide : str
        String containing peptide sequence with variable modifications

    Returns
    -------
    A string containing the modification info.

    '''
    peptide_split = re.findall(r"[^\W\d_]\[\d+.\d+\]\[\d+.\d+\]|[^\W\d_]\[\d+.\d+\]|[^\W\d_]", peptide)
    positions = [any(char.isdigit() for char in inputString) for inputString in peptide_split]
    positions = [str(y) for y in (np.where(positions)[0] + 1)]
    mass_mod = re.findall(r"\[\d+.\d+\]\[\d+.\d+\]|\[\d+.\d+\]", peptide)
    modification_info = [''.join(x) for x in zip(positions, mass_mod)]
    modification_info = ','.join(modification_info)
    modification_info = modification_info.replace('][', ',')
    return(modification_info)
###############################################################################