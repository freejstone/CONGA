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
import time
import datetime
from tqdm import tqdm
import re
import platform


USAGE = """USAGE: Groupwalk.py [options] <narrow> <wide> <matching>

  This program implements the Group-walk algorithm, including the
  creation of groups and the filtering procedure that eliminates
  similar peptide matches per scan, The first input file is the narrow
  search file output, the second input file is the open search file
  output, and the last input file contains the target-decoy peptide
  pairs. Output is a list of PSMs and their corresponding q-values.
  The input search files can be either tab-delimited .txt Tide search files,
  tab-delimited .txt Comet search files, or tab-delimited .tsv files from 
  MS-Fragger. Last input file is not required for Comet.
  The output file is a tab-delimited .txt file containing the following
  information about the PSMs used by groupwalk: the peptide sequence, 
  the scan number, the score, the target-decoy label, whether the PSM came 
  from the narrow or open search file, and a q-value.

  Options:
      
    --output_dir <string> The file-path of the output
                          Default = './'.
                          
    --file_root <string>  The file name of the output
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
                          

    --score <string>      Either 'tailor', 'xcorr', 'e-value' or 'hyper'. 
                          The score that will be used by group-walk. If
                          'tailor', it is assumed the search files are 
                          derived from Tide. If 'xcorr', either Tide 
                          search or Comet is assumed to be used. If 
                          'e-value', Comet is assumed. If 'hyper' it 
                          is assumed the search files are derived from
                          MS-Fragger.
                          Default = 'tailor'.
    
    --account_mods <T|F>  To determine whether the group-walk algorithm
                          selects the best PSM among the equivalent
                          classes of peptides which are equal up to
                          variable modification, or not.
                          Default = T.
                          
    --return_all_mods <T|F>         (Not added yet) If --account_mods T and
                                    --return_all_mods T, all target
                                    peptides equal to a peptide up to
                                    variable modification used in the
                                    group-walk algorithm will also
                                    be reported with the same q-value.
                                    Default = F.
                          
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
                                
    --dcy_prefix <string>       The prefix used for the decoy proteins.
                                Default = 'decoy_'
                                
    --seed <int>          Set random seed.
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
    '''
    Parameters
    ----------
    search_df : Pandas Dataframe
        Dataframe containing just one scan, and multiple PSMs.
    thresh : float, optional
        The similarity score threshold used to filter neighbouring PSMs. The default is 0.25.
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
def filter_narrow_open(narrow_target_decoys, open_target_decoys, score, open_top = 2, thresh = 0.25, n_processes = 1, neighbour_remove = True, tide_used = 'tide', static_mods = {'C':57.02146}):
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
    target_decoys_all = pd.concat([narrow_target_decoys, open_target_decoys]) #combine narrow and open PSMs
    
    #makes sure the PSMs from the narrow and open search are ordered by scan first, then by their score
    #indeed it makes sense to use xcorr_score over tailor score, as tailor_score uses candidate PSMs to 
    #normalise the xcorr_score - this differs from narrow to open. 
    if tide_used == 'tide':
        target_decoys_all = target_decoys_all.sort_values(by=['file', 'scan', 'xcorr_score'], ascending = False)
    elif tide_used == 'comet' and score == 'xcorr':
        target_decoys_all = target_decoys_all.sort_values(by=['scan', "charge", "spectrum_neutral_mass", 'xcorr_score'], ascending = False)
    elif tide_used == 'comet' and score == 'e-value':
        target_decoys_all = target_decoys_all.sort_values(by=['scan', "charge", "spectrum_neutral_mass", 'e-value'], ascending = True)
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
            results = target_decoys_all.groupby(['file', 'scan']).progress_apply(lambda x: filter_scan(x, thresh, 0.05, static_mods)) #apply filtering by experimental scan
        elif tide_used == 'comet':
            results = target_decoys_all.groupby(['scan', "charge", "spectrum_neutral_mass"]).progress_apply(lambda x: filter_scan(x, thresh, 0.05, static_mods)) #apply filtering by experimental scan
        else:
            results = target_decoys_all.groupby('scannum').progress_apply(lambda x: filter_scan(x, thresh, 0.05, static_mods))  #apply filtering by experimental scan
        if score == 'e-value': #we sort the results in ascending order since target_decoys_all have been sorted in such a manner
            results = results.sort_index(ascending = True)
        else:
            results = results.sort_index(ascending = False)
        results = results.apply(pd.Series).stack().reset_index()
        
        target_decoys_all['drop_scan'] = results[0]  #create drop_scan column indicating which PSMs are being kept
    else:
        #if more than 1 thread, we partition the dataframe
        if tide_used == "tide":
            target_decoys_all['split_col'] = pd.qcut(target_decoys_all['scan'], n_processes)
        elif tide_used == "comet":
            target_decoys_all['split_col'] = pd.qcut(target_decoys_all['scan'], n_processes)
        else:
            target_decoys_all['split_col'] = pd.qcut(target_decoys_all['scannum'], n_processes)
        target_decoys_grouped = target_decoys_all.groupby(target_decoys_all.split_col)
        list_of_df = [0]*n_processes #create a list of partitioned dataframes
        for i in range(len(target_decoys_all['split_col'].unique())):
            list_of_df[i] = target_decoys_grouped.get_group(target_decoys_all['split_col'].unique()[i])
            list_of_df[i].reset_index(drop = True, inplace = True)
        
        def wrapper(df, q, thresh, frag_bin_size, static_mods): #wrapper is used to update the manager queue every time the wrapper is called
            result = filter_scan(df, thresh, frag_bin_size, static_mods)
            q.put(1)
            return(result)
        
        if tide_used == 'tide':
            def filter_scan_subset(df, q, task_number, return_dict): #function used to apply wrapper to a partitioned dataframe
                sys.stderr.write("Starting Process " + str(task_number) + " \n")
                logging.info("Starting Process " + str(task_number))
                results = df.groupby(['file', 'scan']).apply(lambda x: wrapper(x, q, thresh, 0.05, static_mods))
                results = results.sort_index(ascending = False)
                results = results.apply(pd.Series).stack().reset_index()
                return_dict[task_number] = results[0]
                
            def listener(q): #listener is used to constantly checking manager queue and update the progress bar whenever the manager queue updates
                pbar = tqdm(total = sum([len(j.groupby(['file', 'scan'])) for j in list_of_df]))
                for item in iter(q.get, None):
                    pbar.update(item)
                    
        elif tide_used == 'comet' and score == 'e-value':
            def filter_scan_subset(df, q, task_number, return_dict):
                sys.stderr.write("Starting Process " + str(task_number) + " \n")
                logging.info("Starting Process " + str(task_number))
                results = df.groupby(['scan', "charge", "spectrum_neutral_mass"]).apply(lambda x: wrapper(x, q, thresh, 0.05, static_mods))
                results = results.sort_index(ascending = True) #Note that ascending is true here
                results = results.apply(pd.Series).stack().reset_index()
                return_dict[task_number] = results[0]
                
            def listener(q):
                pbar = tqdm(total = sum([len(j.groupby(['scan', "charge", "spectrum_neutral_mass"])) for j in list_of_df]))
                for item in iter(q.get, None):
                    pbar.update(item)
        elif tide_used == 'comet' and score == 'xcorr':
            def filter_scan_subset(df, q, task_number, return_dict):
                sys.stderr.write("Starting Process " + str(task_number) + " \n")
                logging.info("Starting Process " + str(task_number))
                results = df.groupby(['scan', "charge", "spectrum_neutral_mass"]).apply(lambda x: wrapper(x, q, thresh, 0.05, static_mods))
                results = results.sort_index(ascending = False) #Note that ascending is true here
                results = results.apply(pd.Series).stack().reset_index()
                return_dict[task_number] = results[0]
                
            def listener(q):
                pbar = tqdm(total = sum([len(j.groupby(['scan', "charge", "spectrum_neutral_mass"])) for j in list_of_df]))
                for item in iter(q.get, None):
                    pbar.update(item)
        else:
            def filter_scan_subset(df, q, task_number, return_dict):
                sys.stderr.write("Starting Process " + str(task_number) + " \n")
                logging.info("Starting Process " + str(task_number))
                results = df.groupby('scannum').apply(lambda x: wrapper(x, q, thresh, 0.05, static_mods))
                results = results.sort_index(ascending = False)
                results = results.apply(pd.Series).stack().reset_index()
                return_dict[task_number] = results[0]
        
            def listener(q):
                pbar = tqdm(total = sum([len(j.groupby('scannum')) for j in list_of_df]))
                for item in iter(q.get, None):
                    pbar.update(item)
            
        manager = multiprocessing.Manager()
        q = manager.Queue() #creating instance of manager queue
        return_dict = manager.dict() #creating manager dict variable that can be used to store changes to the argument by ALL processes at the same time
        proc = multiprocessing.Process(target=listener, args=(q,)) 
        proc.start() #start listening for updates to manager queue
        workers = [multiprocessing.Process(target = filter_scan_subset, args=(list_of_df[i], q, i, return_dict)) for i in range(n_processes)] #now run each of the processes
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

def create_groups(target_decoys, narrow_target_decoys, peptide_list, dcy_prefix = 'decoy_', K = 40, tops = 2, score = 'tailor', account_mods = True, any_mods = True, precursor_bin_width = 1.0005079/4, group_thresh = 0.01, adaptive = True, min_group_size = 2, n_top_groups = 4, tide_used = 'tide', print_group_pi0 = True):
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

    Returns
    -------
    Dataframe containing the winning PSMs after head-to-head competition and their group labels.

    '''
    #take only the top 1 from the narrow search
    if tide_used == 'tide':
        narrow_target_decoys = narrow_target_decoys[narrow_target_decoys['xcorr_rank'] == 1]
    elif tide_used == "comet" and score == "xcorr":
        narrow_target_decoys = narrow_target_decoys[narrow_target_decoys['xcorr_rank'] == 1]
    elif tide_used == "comet" and score == 'e-value':
        narrow_target_decoys = narrow_target_decoys[narrow_target_decoys['e_rank'] == 1]
    else:
        narrow_target_decoys = narrow_target_decoys[narrow_target_decoys['hit_rank'] == 1]
    narrow_target_decoys.loc[:, 'database'] = 'narrow'
    
    #rounding in case of any numerical issues
    #cannot apply this to e-values, as they go right up to E-12
    if score == 'tailor':
        target_decoys['tailor_score'] = round(target_decoys['tailor_score'], 8)
        narrow_target_decoys['tailor_score'] = round(narrow_target_decoys['tailor_score'], 8)
    elif score == 'xcorr':
        target_decoys['xcorr_score'] = round(target_decoys['xcorr_score'], 8)
        narrow_target_decoys['xcorr_score'] = round(narrow_target_decoys['xcorr_score'], 8)
    elif score == 'hyper':
        target_decoys['hyperscore'] = round(target_decoys['hyperscore'], 8)
        narrow_target_decoys['hyperscore'] = round(narrow_target_decoys['hyperscore'], 8)
    
    logging.info("Doing head to head competition.")
    sys.stderr.write("Doing head to head competition.\n")
    
    #take the open and narrow PSMs
    target_decoys_open = target_decoys[target_decoys['database'] == "open"]
    target_decoys_narrow = target_decoys[target_decoys['database'] == "narrow"]

    if tide_used == 'tide':
        #get target and decoys separately
        targets_narrow = target_decoys_narrow[( target_decoys_narrow['target_decoy'] == 'target' )]
        decoys_narrow = target_decoys_narrow[( target_decoys_narrow['target_decoy'] == 'decoy' )]
        
        #taking the top ranked PSMs from the open search as specified by user - nominally top 2
        targets_open = target_decoys_open[( target_decoys_open['target_decoy'] == 'target' ) & ( target_decoys_open['xcorr_rank'].isin(range(1, tops + 1)) )]
        decoys_open = target_decoys_open[( target_decoys_open['target_decoy'] == 'decoy' ) & ( target_decoys_open['xcorr_rank'].isin(range(1, tops + 1)) )]    
    elif tide_used == "comet":
        #get target and decoys separately
        targets_narrow = target_decoys_narrow[~(target_decoys_narrow['protein_id'].str.contains(dcy_prefix))].copy()
        decoys_narrow = target_decoys_narrow[(target_decoys_narrow['protein_id'].str.contains(dcy_prefix))].copy()
        targets_narrow['target_decoy'] = 'target'
        decoys_narrow['target_decoy'] = 'decoy'
        
        if score == 'xcorr':
            #taking the top ranked PSMs from the open search as specified by user - nominally top 2
            targets_open = target_decoys_open[(~( target_decoys_open['protein_id'].str.contains(dcy_prefix) )) & ( target_decoys_open['xcorr_rank'].isin(range(1, tops + 1)) )].copy()
            decoys_open = target_decoys_open[( target_decoys_open['protein_id'].str.contains(dcy_prefix) ) & ( target_decoys_open['xcorr_rank'].isin(range(1, tops + 1)) )].copy()  
            targets_open['target_decoy'] = 'target'
            decoys_open['target_decoy'] = 'decoy'
        else:
            #taking the top ranked PSMs from the open search as specified by user - nominally top 2
            targets_open = target_decoys_open[(~( target_decoys_open['protein_id'].str.contains(dcy_prefix) )) & ( target_decoys_open['e_rank'].isin(range(1, tops + 1)) )].copy()
            decoys_open = target_decoys_open[( target_decoys_open['protein_id'].str.contains(dcy_prefix) ) & ( target_decoys_open['e_rank'].isin(range(1, tops + 1)) )].copy()  
            targets_open['target_decoy'] = 'target'
            decoys_open['target_decoy'] = 'decoy'
    else:
        #get target and decoys separately
        targets_narrow = target_decoys_narrow[~(target_decoys_narrow['protein'].str.contains(dcy_prefix))].copy()
        decoys_narrow = target_decoys_narrow[(target_decoys_narrow['protein'].str.contains(dcy_prefix))].copy()
        targets_narrow['target_decoy'] = 'target'
        decoys_narrow['target_decoy'] = 'decoy'
        
        #taking the top ranked PSMs from the open search as specified by user - nominally top 2
        targets_open = target_decoys_open[(~( target_decoys_open[target_decoys_open['protein'].str.contains(dcy_prefix)] )) & ( target_decoys_open['hit_rank'].isin(range(1, tops + 1)) )].copy()
        decoys_open = target_decoys_open[( target_decoys_open[target_decoys_open['protein'].str.contains(dcy_prefix)] ) & ( target_decoys_open['hit_rank'].isin(range(1, tops + 1)) )].copy()  
        targets_open['target_decoy'] = 'target'
        decoys_open['target_decoy'] = 'decoy'
            
    #sorting PSMs by their score
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
    elif score == 'e-value':
        decoys_open = decoys_open.sort_values(by=['e-value'], ascending=True)
        targets_open = targets_open.sort_values(by=['e-value'], ascending=True)
        decoys_narrow = decoys_narrow.sort_values(by=['e-value'], ascending=True)
        targets_narrow = targets_narrow.sort_values( by=['e-value'], ascending=True)
    elif score == "hyper":
        decoys_open = decoys_open.sort_values(by=['hyperscore'], ascending=False)
        targets_open = targets_open.sort_values(by=['hyperscore'], ascending=False)
        decoys_narrow = decoys_narrow.sort_values(by=['hyperscore'], ascending=False)
        targets_narrow = targets_narrow.sort_values( by=['hyperscore'], ascending=False)
    if account_mods:
        if tide_used == 'tide' or tide_used == 'comet':
            #taking best PSM score up to variable modification in each of the database
            decoys_open = decoys_open.drop_duplicates(subset = ['original_target_sequence'])
            targets_open = targets_open.drop_duplicates(subset = ['original_target_sequence'])
            decoys_narrow = decoys_narrow.drop_duplicates(subset = ['original_target_sequence'])
            targets_narrow = targets_narrow.drop_duplicates(subset = ['original_target_sequence'])
        else:
            #MS_fragger's peptide column does not include variable modifications, so we can use that
            decoys_open = decoys_open.drop_duplicates(subset = ['peptide'])
            targets_open = targets_open.drop_duplicates(subset = ['peptide'])
            decoys_narrow = decoys_narrow.drop_duplicates(subset = ['peptide'])
            targets_narrow = targets_narrow.drop_duplicates(subset = ['peptide'])
    else:
        #taking best PSM score up to sequence in each of the database
        decoys_open = decoys_open.drop_duplicates(subset = ['sequence'])
        targets_open = targets_open.drop_duplicates(subset = ['sequence'])
        decoys_narrow = decoys_narrow.drop_duplicates(subset = ['sequence'])
        targets_narrow = targets_narrow.drop_duplicates(subset = ['sequence'])

    #resetting the indices
    decoys_open.reset_index(drop = True, inplace = True)
    targets_open.reset_index(drop = True, inplace = True)
    decoys_narrow.reset_index(drop = True, inplace = True)
    targets_narrow.reset_index(drop = True, inplace = True)
    
    #creating PSM column called scan_plus_seq to uniquely identify the PSMs in the open search
    if tide_used == 'tide':
        decoys_open['scan_plus_seq'] = decoys_open['file'].astype(str) + decoys_open['scan'].astype(str) + ' ' + decoys_open['sequence']
        targets_open['scan_plus_seq'] = targets_open['file'].astype(str) + targets_open['scan'].astype(str) + ' ' + targets_open['sequence']
    elif tide_used == 'comet':
        decoys_open['scan_plus_seq'] = decoys_open['scan'].astype(str) + ' ' + decoys_open['charge'].astype(str) + ' ' + decoys_open['spectrum_neutral_mass'].astype(str) + ' ' + decoys_open['sequence']
        targets_open['scan_plus_seq'] = targets_open['scan'].astype(str) + ' ' + targets_open['charge'].astype(str) + ' ' + targets_open['spectrum_neutral_mass'].astype(str) + ' ' + targets_open['sequence']
    else:
        decoys_open['scan_plus_seq'] = decoys_open['scannum'].astype(str) + ' ' + decoys_open['sequence']
        targets_open['scan_plus_seq'] = targets_open['scannum'].astype(str) + ' ' + targets_open['sequence']

    #creating PSM column to uniquely identify the PSMs in the narrow search
    narrow_target_decoys.reset_index(drop = True, inplace = True)
    if tide_used == 'tide':
        narrow_target_decoys['scan_plus_seq'] = narrow_target_decoys['file'].astype(str) + narrow_target_decoys['scan'].astype(str) + ' ' + narrow_target_decoys['sequence']
    elif tide_used == 'comet':
        narrow_target_decoys['scan_plus_seq'] = narrow_target_decoys['scan'].astype(str) + ' ' + narrow_target_decoys['charge'].astype(str) + ' ' + narrow_target_decoys['spectrum_neutral_mass'].astype(str) + ' ' + narrow_target_decoys['sequence']
    else:
        narrow_target_decoys['scan_plus_seq'] = narrow_target_decoys['scannum'].astype(str) + ' ' + narrow_target_decoys['sequence']
    
    #translating into an iterable
    unique_scan_seq_td_narrow = set(narrow_target_decoys['scan_plus_seq'])
    
    #decoy/target PSMs open get assigned "narrow" database if they exist in the narrow search
    decoys_open.loc[decoys_open['scan_plus_seq'].isin(unique_scan_seq_td_narrow), 'database'] = 'narrow'    
    targets_open.loc[targets_open['scan_plus_seq'].isin(unique_scan_seq_td_narrow), 'database'] = 'narrow'
    
    #subsetting the decoy-target PSMs in the open-search that are also found in the narrow search
    scan_plus_seq_decoys_sub = decoys_open[decoys_open['scan_plus_seq'].isin(unique_scan_seq_td_narrow)]['scan_plus_seq']
    scan_plus_seq_targets_sub = targets_open[targets_open['scan_plus_seq'].isin(unique_scan_seq_td_narrow)]['scan_plus_seq']
    
    #translating into iterables
    unique_scan_plus_seq_decoys_sub = set(scan_plus_seq_decoys_sub)
    unique_scan_plus_seq_targets_sub = set(scan_plus_seq_targets_sub)
    
    #take subset of narrow PSMs that contain the decoy PSMs in the open
    #reorder this set so that they match the decoy PSMs in the open
    target_decoys_narrow_sub = narrow_target_decoys[narrow_target_decoys['scan_plus_seq'].isin(unique_scan_plus_seq_decoys_sub)].copy()
    target_decoys_narrow_sub.reset_index(drop = True, inplace = True)
    target_decoys_narrow_sub.scan_plus_seq = target_decoys_narrow_sub.scan_plus_seq.astype("category")
    target_decoys_narrow_sub.scan_plus_seq.cat.set_categories(scan_plus_seq_decoys_sub, inplace = True)
    target_decoys_narrow_sub.reset_index(drop = True, inplace = True)
    
    #assign the narrow decoy PSM scores to the open decoy PSMs that are found in the narrow search
    if score == 'tailor':
        decoys_open.loc[decoys_open['scan_plus_seq'].isin(unique_scan_seq_td_narrow), 'tailor_score'] = list(target_decoys_narrow_sub.sort_values(['scan_plus_seq'])['tailor_score'])
    elif score == 'xcorr':
        decoys_open.loc[decoys_open['scan_plus_seq'].isin(unique_scan_seq_td_narrow), 'xcorr_score'] = list(target_decoys_narrow_sub.sort_values(['scan_plus_seq'])['xcorr_score'])
    elif score == 'e-value':
        decoys_open.loc[decoys_open['scan_plus_seq'].isin(unique_scan_seq_td_narrow), 'e-value'] = list(target_decoys_narrow_sub.sort_values(['scan_plus_seq'])['e-value'])
    elif score == 'hyper':
        decoys_open.loc[decoys_open['scan_plus_seq'].isin(unique_scan_seq_td_narrow), 'hyperscore'] = list(target_decoys_narrow_sub.sort_values(['scan_plus_seq'])['hyperscore'])
    
    #take subset of narrow PSMs that contain the target PSMs in the open
    #reorder this set so that they match the target PSMs in the open
    target_decoys_narrow_sub = narrow_target_decoys[narrow_target_decoys['scan_plus_seq'].isin(unique_scan_plus_seq_targets_sub)].copy()
    target_decoys_narrow_sub = target_decoys_narrow_sub.reset_index(drop = True)
    target_decoys_narrow_sub.scan_plus_seq = target_decoys_narrow_sub.scan_plus_seq.astype("category")
    target_decoys_narrow_sub.scan_plus_seq.cat.set_categories(scan_plus_seq_targets_sub, inplace = True)
    target_decoys_narrow_sub.reset_index(drop = True, inplace = True)
    
    #assign the narrow target PSM scores to the open target PSMs that are found in the narrow search 
    if score == 'tailor':
        targets_open.loc[targets_open['scan_plus_seq'].isin(unique_scan_seq_td_narrow), 'tailor_score'] = list(target_decoys_narrow_sub.sort_values(['scan_plus_seq'])['tailor_score'])
    elif score == 'xcorr':
        targets_open.loc[targets_open['scan_plus_seq'].isin(unique_scan_seq_td_narrow), 'xcorr_score'] = list(target_decoys_narrow_sub.sort_values(['scan_plus_seq'])['xcorr_score'])
    elif score == 'e-value':
        targets_open.loc[targets_open['scan_plus_seq'].isin(unique_scan_seq_td_narrow), 'e-value'] = list(target_decoys_narrow_sub.sort_values(['scan_plus_seq'])['e-value'])
    elif score == 'hyper':
        targets_open.loc[targets_open['scan_plus_seq'].isin(unique_scan_seq_td_narrow), 'hyperscore'] = list(target_decoys_narrow_sub.sort_values(['scan_plus_seq'])['hyperscore'])
        
    
    #rejoining the decoys and targets from the two databases
    #note if account_mods = T, this will still mean that PSMs are not unique up to variable mod from the two databases
    #note if account_mods = F, this will still mean that PSMs are not unique up to sequence from the two databases.
    decoys = pd.concat([decoys_narrow, decoys_open])
    targets = pd.concat([targets_narrow, targets_open])
    
    if tide_used == 'tide':
        if account_mods or (not any_mods):
            #combine targets/decoys
            target_decoys_final = pd.concat([targets, decoys])
            #randomly shuffle order
            target_decoys_final = target_decoys_final.sample(frac = 1).reset_index(drop=True)
            
        else:
            #At this point, there are possibly two PSMs to the same peptide in 'decoys' and 'targets', they 
            #are the narrow and open version PSMs
            #Take the narrow PSM
            decoys = decoys.drop_duplicates(subset = ['sequence'])
            targets = targets.drop_duplicates(subset = ['sequence'])
            #we reset index since we refer to the peptide_list to match the corresponding target pairs to
            #decoys
            decoys.reset_index(drop = True, inplace = True)
            targets.reset_index(drop = True, inplace = True)
            #Changes original_target_sequence to include variable mods
            targets['original_target_sequence'] = targets['sequence']
            #peptide_list is required here to pair target-decoys
            peptide_list_sub = peptide_list[peptide_list['decoy'].isin(decoys['sequence'])].copy()
            peptide_list_sub.decoy = peptide_list_sub.decoy.astype('category')
            peptide_list_sub.decoy.cat.set_categories(decoys['sequence'], inplace = True)
            peptide_list_sub = peptide_list_sub.sort_values(['decoy'])
            peptide_list_sub.reset_index(drop = True, inplace = True)
            #find the associated pair for the decoy sequences
            decoys['original_target_sequence'] = peptide_list_sub['target']
            
            #combine targets and decoys
            target_decoys_final = pd.concat([targets, decoys])
            #randomly shuffle order, so that ties between targets and decoys are broken
            target_decoys_final = target_decoys_final.sample(frac = 1).reset_index(drop=True)
            
        #order according to database, so narrow always take precedence in head-to-head competitions
        if score == 'tailor':
            target_decoys_final = target_decoys_final.sort_values(by=['database', 'tailor_score'], ascending = [True, False])
        elif score == 'xcorr':
            target_decoys_final = target_decoys_final.sort_values(by=['database', 'xcorr_score'], ascending = [True, False])
        
        #take best PSM up to variable modification, (narrow takes precedence) if account_mods
        #else take best PSM up to sequence (narrow takes precedence)
        target_decoys_final = target_decoys_final.drop_duplicates(subset = ['original_target_sequence'])
        
        #gather winning information
        if score == 'tailor':
            winning_scores = target_decoys_final['tailor_score']
        elif score == 'xcorr':
            winning_scores = target_decoys_final['xcorr_score']
            
        labels = target_decoys_final['target_decoy'].copy()
        labels[labels == 'target'] = 1
        labels[labels == 'decoy'] = -1
        winning_peptides = target_decoys_final['sequence']
        rank = target_decoys_final['xcorr_rank']
        delta_mass = target_decoys_final['spectrum_neutral_mass'] - target_decoys_final['peptide_mass']
        database = target_decoys_final['database']
        scan = target_decoys_final['scan']
        file = target_decoys_final['file']
        df = pd.DataFrame(zip(winning_scores, labels, delta_mass, winning_peptides, rank, database, scan, file), columns = ['winning_scores', 'labels', 'delta_mass', 'winning_peptides', 'rank', 'database', 'scan', 'file'])
    elif tide_used == 'comet':
        #combine targets/decoys
        target_decoys_final = pd.concat([targets, decoys])
        #randomly shuffle order
        target_decoys_final = target_decoys_final.sample(frac = 1).reset_index(drop=True)
        
        #order according to database, so narrow always take precedence
        if score == 'e-value':
            target_decoys_final = target_decoys_final.sort_values(by=['database', 'e-value'], ascending = [True, True])
        elif score == 'xcorr':
            target_decoys_final = target_decoys_final.sort_values(by=['database', 'xcorr_score'], ascending = [True, False])
        
        #take best PSM up to variable modification, (narrow takes precedence) if account_mods
        #else take best PSM up to sequence (narrow takes precedence)
        if account_mods:
            target_decoys_final = target_decoys_final.drop_duplicates(subset = ['original_target_sequence'])
        else:
            target_decoys_final = target_decoys_final.drop_duplicates(subset = ['target_sequence'])
        
        #gather winning information
        if score == 'e-value':
            winning_scores = target_decoys_final['e-value']
        elif score == 'xcorr':
            winning_scores = target_decoys_final['xcorr_score']
            
        labels = target_decoys_final['target_decoy'].copy()
        labels[labels == 'target'] = 1
        labels[labels == 'decoy'] = -1
        winning_peptides = target_decoys_final['sequence']
        if score == 'e-value':
            rank = target_decoys_final['e_rank']
        if score == 'xcorr':
            rank = target_decoys_final['xcorr_rank']
        delta_mass = target_decoys_final['spectrum_neutral_mass'] - target_decoys_final['peptide_mass']
        database = target_decoys_final['database']
        scan = target_decoys_final['scan'].astype(str) + ' ' + target_decoys_final['charge'].astype(str) + ' ' + target_decoys_final['spectrum_neutral_mass'].astype(str)
        df = pd.DataFrame(zip(winning_scores, labels, delta_mass, winning_peptides, rank, database, scan), columns = ['winning_scores', 'labels', 'delta_mass', 'winning_peptides', 'rank', 'database', 'scan_charge_sp_neutral_mass'])
 
    
    else:
        #unnecessary, but in case
        peptide_list = peptide_list.drop_duplicates()
        peptide_list.reset_index(drop = True, inplace = True)
        
        if account_mods:
            #Take the narrow PSM between narrow and open up to variable modification.
            decoys = decoys.drop_duplicates(subset = ['peptide'])
            targets = targets.drop_duplicates(subset = ['peptide'])
            
            #we reset index since we refer to the peptide_list to match the corresponding target pairs to
            #decoys
            decoys.reset_index(drop = True, inplace = True)
            targets.reset_index(drop = True, inplace = True)
            
            #create equivalent 'original_target_sequence' column for MS fragger results
            targets['original_target_sequence'] = targets['peptide']
            peptide_list_sub = peptide_list[peptide_list['decoy'].isin(decoys['peptide'])]
            peptide_list_sub.decoy = peptide_list_sub.decoy.astype('category')
            peptide_list_sub.decoy.cat.set_categories(decoys['peptide'], inplace = True)
            peptide_list_sub = peptide_list_sub.sort_values(['decoy'])
            peptide_list_sub.reset_index(drop = True, inplace = True)
            #find the associated pair for the decoy sequences
            decoys['original_target_sequence'] = peptide_list_sub['target']
        else:
            #Take the narrow PSM between narrow and open equal to the same sequence.
            decoys = decoys.drop_duplicates(subset = ['sequence'])
            targets = targets.drop_duplicates(subset = ['sequence'])
            
            #we reset index since we refer to the peptide_list to match the corresponding target pairs to
            #decoys
            decoys.reset_index(drop = True, inplace = True)
            targets.reset_index(drop = True, inplace = True)
            
            #create 'original_target_sequence' column for MS fragger results on the individual sequence-level
            targets['original_target_sequence'] = targets['sequence']
            
            peptide_list_sub = peptide_list[peptide_list['decoy'].isin(decoys['sequence'])]
            peptide_list_sub.decoy = peptide_list_sub.decoy.astype('category')
            peptide_list_sub.decoy.cat.set_categories(decoys['sequence'], inplace = True)
            peptide_list_sub = peptide_list_sub.sort_values(['decoy'])
            peptide_list_sub.reset_index(drop = True, inplace = True)
            #find the associated pair for the decoy sequences
            decoys['original_target_sequence'] = peptide_list_sub['target']
            
        #combine targets/decoys
        target_decoys_final = pd.concat([targets, decoys])
        #randomly shuffle order
        target_decoys_final = target_decoys_final.sample(frac = 1)
        #do h2h competition by taking the best scoring PSM
        target_decoys_final = target_decoys_final.sort_values(by=['database', 'hyperscore'], ascending = [True, False])
        
        #take best PSM up to variable modification, (narrow takes precedence) if account_mods
        #else take best PSM up to sequence (narrow takes precedence)
        target_decoys_final = target_decoys_final.drop_duplicates(subset = ['original_target_sequence'])
        
        winning_scores = target_decoys_final['hyperscore']
        
        labels = target_decoys_final['target_decoy']
        labels[labels == 'target'] = 1
        labels[labels == 'decoy'] = -1
        winning_peptides = target_decoys_final['sequence']
        rank = target_decoys_final['hit_rank']
        delta_mass = target_decoys_final['precursor_neutral_mass'] - target_decoys_final['calc_neutral_pep_mass']
        database = target_decoys_final['database']
        scan = target_decoys_final['scannum']
        df = pd.DataFrame(zip(winning_scores, labels, delta_mass, winning_peptides, rank, database, scan), columns = ['winning_scores', 'labels', 'delta_mass', 'winning_peptides', 'rank', 'database', 'scan'])
        
    df = df.sample(frac = 1)
    
    #creating bins for the mass differences
    breaks_p = np.arange(0, 100 + 2*precursor_bin_width, precursor_bin_width) - precursor_bin_width/2
    breaks_n = list(reversed(-breaks_p))
    breaks = pd.Series(breaks_n[0:(len(breaks_n) - 1)] + list(breaks_p), name = 'bins')
    digitized = np.digitize(df['delta_mass'], breaks)
    df['bins'] = digitized
    df.reset_index(drop = True, inplace = True)
    
    #logging.info("Head to head competition is complete.")
    #sys.stderr.write("Head to head competition is complete.\n")
    
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
        rank_lab = 'top 1 PSMs'
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
                    if largest_bin_considered == 0 and i == 0:
                        all_group_ids[cand_group] = rank_lab + ' & top '  + str(i + 1) + ' mass bin'
                    elif largest_bin_considered + 1 == i:
                        all_group_ids[cand_group] = rank_lab + ' & top '  + str(i + 1) + ' mass bin'
                    else:
                        all_group_ids[cand_group] = rank_lab + ' & top ('  + str(largest_bin_considered + 1) + ', ' + str(i + 1) + '] mass bins'
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
                        if largest_bin_considered == 0 and i == 0:
                            all_group_ids[cand_group] = rank_lab + ' & top '  + str(i + 1) + ' mass bin'
                        elif largest_bin_considered + 1 == i:
                            all_group_ids[cand_group] = rank_lab + ' & top '  + str(i + 1) + ' mass bin'
                        else:
                            all_group_ids[cand_group] = rank_lab + ' & top ('  + str(largest_bin_considered + 1) + ', ' + str(i + 1) + '] mass bins'
                        largest_bin_considered = i
                    i += 1
            else:
                break
            
        all_group_ids[(all_group_ids == 0) & (df['rank'] == 1)] = "left over group"
        #For PSMs of rank 2 or higher, we use the bins considered determined by rank 1 PSMs and bundle them together in one group
        #If they are not sufficiently many of them, we throw them into the "1 rank left over" group
        if (tops > 1) and (largest_bin_considered > 0):
            inds_higher_rank = ( df['bins'].isin(bin_freq['bin'].loc[0:largest_bin_considered]) ) &  ( df['rank'].isin(range(2, tops + 1)) ) & (df['database'] == "open")
            if sum(inds_higher_rank) <= min_group_size*K:
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
def reverse_sequence(sequences):
    '''
    Parameters
    ----------
    sequence : Pandas Series
        Recovers target sequences from decoy sequences for Comet.

    Returns
    -------
    Recovers target sequences.
    '''
    results = sequences.str.findall('[a-zA-Z]\\[\d+\.\d+\\]|[a-zA-Z]\\[\d+\\]|[a-zA-Z]')
    results = results.apply(lambda x: x[-2::-1] + [x[-1]])
    results = results.apply(''.join)
    return(results)

###############################################################################
def main():
    global USAGE
    
    start_time = time.time()
    
    # Set default values for parameters.
    K = 40
    tops_gw = 2
    tops_open = 5
    score = 'tailor'
    #random_h2h = True
    account_mods = True
    precursor_bin_width = 1.0005079/4
    adaptive = True
    print_chimera = True
    group_thresh = 0.01
    print_group_pi0 = True
    min_group_size = 2
    n_top_groups = 4
    neighbour_remove = True
    thresh = 0.25
    #concat = True
    return_filt_search = False
    correction = 1
    n_processes = 1
    return_frontier = False
    output_dir = './'
    file_root = 'group_walk_results.txt'
    frontier_name = 'frontier_results.txt'
    static_mods = {'C':57.02146}
    dcy_prefix = 'decoy_'
    seed = None
    
    command_line = ' '.join(sys.argv)
    
    
    # Parse the command line.
    sys.argv = sys.argv[1:]
    while (any('--' in string for string in sys.argv)):
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
        elif (next_arg == "--file_root"):
            file_root = str(sys.argv[0])  
            sys.argv = sys.argv[1:]
        elif (next_arg == "--frontier_name"):
            frontier_name = str(sys.argv[0])  
            sys.argv = sys.argv[1:]
        elif (next_arg == "--static_mods"):
            static_mods = parse_static_mods(sys.argv[0])
            sys.argv = sys.argv[1:]
        elif (next_arg == "--dcy_prefix"):
            dcy_prefix = str(sys.argv[0])
            sys.argv = sys.argv[1:]
        elif (next_arg == '--seed'):
            seed = int(sys.argv[0])
            sys.argv = sys.argv[1:]
        else:
            sys.stderr.write("Invalid option (%s).c" % next_arg)
            sys.exit(1)
    if (len(sys.argv) == 2):
        search_file_narrow = sys.argv[0]
        search_file_open = sys.argv[1]
        td_list = ''
    elif (len(sys.argv) == 3):
        search_file_narrow = sys.argv[0]
        search_file_open = sys.argv[1]
        td_list = sys.argv[2]
    else:
        sys.stderr.write(USAGE)
        sys.exit(1)
    
    #setting seed for reproducibility
    if type(seed) == int:
        random.seed(seed)
    
    if '.txt' not in file_root:
        file_root += '.txt'
    if '.txt' not in frontier_name:
        frontier_name += '.txt'
    
    #check if output directory exists, if not create and store log file there.
    if output_dir != './':
        if os.path.isdir(output_dir):
            logging.basicConfig(filename=output_dir + "/" + file_root.replace('.txt', '.log.txt'), level=logging.DEBUG, format = '%(levelname)s: %(message)s')
        else:
            os.mkdir(output_dir)
            logging.basicConfig(filename=output_dir + "/" + file_root.replace('.txt', '.log.txt'), level=logging.DEBUG, format = '%(levelname)s: %(message)s')

    #print CPU info
    logging.info('CPU: ' + str(platform.platform()))
    sys.stderr.write('CPU: ' + str(platform.platform()) + " \n")
    
    #print date time info
    logging.info(str(datetime.datetime.now()))
    sys.stderr.write(str(datetime.datetime.now()) + " \n")
    
    #print command used
    logging.info('Command used: ' + command_line)
    sys.stderr.write('Command used: ' + command_line + "\n")
    
    sys.stderr.write("Successfully read in arguments. \n")
    logging.info("Successfully read in arguments")
    
    if not account_mods:
        sys.stderr.write("Warning: no longer accounting for variable modification. FDR control not guaranteed. \n")
        logging.warning("No longer accounting for variable modification. FDR control not guaranteed.")
    
    sys.stderr.write("Reading in search files. \n")
    logging.info("Reading in search files.")
    
    #read in files and standardize column names
    narrow_target_decoys = pd.read_table(search_file_narrow)
    open_target_decoys = pd.read_table(search_file_open)
    narrow_target_decoys.columns = narrow_target_decoys.columns.str.strip().str.lower().str.replace(' ', '_', regex = False).str.replace('(', '', regex = False).str.replace(')', '', regex = False).str.replace('/', '_', regex = False).str.replace('.', '_', regex = False)
    open_target_decoys.columns = open_target_decoys.columns.str.strip().str.lower().str.replace(' ', '_', regex = False).str.replace('(', '', regex = False).str.replace(')', '', regex = False).str.replace('/', '_', regex = False).str.replace('.', '_', regex = False)
    
    sys.stderr.write("Successfully read in search files. \n")
    logging.info("Successfully read in search files.")
    
    #check if score matches with search file type given by the user
    #and establishes score and search file type
    if score == 'tailor' or score == 'xcorr' or score == 'e-value':
        check_1 = ('hyperscore' in narrow_target_decoys.columns)
        check_2 = ('hyperscore' in open_target_decoys.columns)
        if check_1 and check_2:
            tide_used = 'ms_fragger'
            sys.stderr.write("MS-Fragger output files detected. \n")
            logging.info("MS-Fragger output files detected.")
            score = 'hyper'
            sys.stderr.write("Using hyperscore. \n")
            logging.info("Using hyperscore.")
        elif (not check_1) and (not check_2):
            check_again_1 = ('e-value' in narrow_target_decoys.columns)
            check_again_2 = ('e-value' in open_target_decoys.columns)
            if check_again_1 and check_again_2:
                tide_used = 'comet'
                sys.stderr.write("Comet output files detected. \n")
                logging.info("Comet output files detected.")
                if score == 'tailor':
                    logging.info("--score tailor was specified yet Comet was used. Please use either e-value or xcorr.")
                    sys.exit("--score tailor was specified yet Comet was used. Please use either e-value or xcorr. \n")
            elif (not check_again_1) and (not check_again_2):
                tide_used = 'tide'
                sys.stderr.write("Tide search output files detected. \n")
                logging.info("Tide search output files detected.")
                if score == 'e-value':
                    logging.info("--score e-value was specified yet Tide search was used. Please use either tailor or xcorr.")
                    sys.exit("--score e-value was specified yet Tide search was used. Please use either tailor or xcorr. \n")
            else:
                sys.stderr.write("e-value was detected in one file but not the other. \n")
                logging.info("e-value was detected in one file but not the other.")
                logging.info("Both search files should contain e-value.")
                sys.exit("Both search files should contain e-value. \n")
        else:
            sys.stderr.write("Hyperscore was detected in one file but not the other. \n")
            logging.info("Hyperscore was detected in one file but not the other.")
            logging.info("Both search files should contain hyperscore.")
            sys.exit("Both search files should contain hyperscore. \n")
    else:
        check_1 = ('hyperscore' in narrow_target_decoys.columns)
        check_2 = ('hyperscore' in open_target_decoys.columns)
        if check_1 and check_2:
            tide_used = 'ms_fragger'
            sys.stderr.write("MS-Fragger output files detected. \n")
            logging.info("MS-Fragger output files detected.")
        else:
            logging.info("--score hyper was specified yet hyperscore was not detected in both search files.")
            sys.exit("--score hyper was specified yet hyperscore was not detected in both search files. \n")
        
   
    #check to see if concantenated search file or separate search file used
    if tide_used == 'tide' or tide_used == 'comet':
        check_1 = any(narrow_target_decoys['protein_id'].str.contains(dcy_prefix))
        check_2 = any(open_target_decoys['protein_id'].str.contains(dcy_prefix))
        if check_1 and check_2:
            logging.info('Concatenated search files detected.')
            sys.stderr.write('Concatenated search files detected. \n')
            concat = True
        elif (not check_1) and (not check_2):
            logging.info('Separate search files detected.')
            sys.stderr.write('Separate search files detected. \n')
            concat = False
        else:
            logging.info('Decoy peptides detected in one file, but not detected in the other.')
            logging.info('Please make both search files contain the PSMs against the target-decoy database, or just the target database.')
            sys.stderr.write('Decoy peptides detected in one file, but not detected in the other. \n')
            sys.exit('Please make both search files contain the PSMs against the target-decoy database, or just the target database. \n')

    else:
        check_1 = any(narrow_target_decoys['protein'].str.contains(dcy_prefix))
        check_2 = any(open_target_decoys['protein'].str.contains(dcy_prefix))
        if check_1 and check_2:
            logging.info('Concatenated search files detected.')
            sys.stderr.write('Concatenated search files detected. \n')
            concat = True
        elif (not check_1) and (not check_2):
            logging.info('Separate search files detected.')
            sys.stderr.write('Separate search files detected. \n')
            concat = False
        else:
            logging.info('Decoy peptides detected in one file, but not detected in the other.')
            logging.info('Please make both search files contain the PSMs against the target-decoy database, or just the target database.')
            sys.stderr.write('Decoy peptides detected in one file, but not detected in the other. \n')
            sys.exit('Please make both search files contain the PSMs against the target-decoy database, or just the target database. \n')
    
    if concat:
        #take the top 1 PSM for narrow search and top 'tops_open' PSM for open search
        if tide_used == 'tide' or (tide_used == 'comet' and score == 'xcorr'):
            narrow_target_decoys = narrow_target_decoys[narrow_target_decoys['xcorr_rank'] == 1]
            open_target_decoys = open_target_decoys[open_target_decoys['xcorr_rank'].isin(range(1, tops_open + 1))]
        elif tide_used == 'comet' and score == 'e-value':
            narrow_target_decoys['e_rank'] = narrow_target_decoys.groupby(["scan", "charge", "spectrum_neutral_mass"])["e-value"].rank("first", ascending=True)
            narrow_target_decoys = narrow_target_decoys[narrow_target_decoys['e_rank'] == 1]
            open_target_decoys['e_rank'] = open_target_decoys.groupby(["scan", "charge", "spectrum_neutral_mass"])["e-value"].rank("first", ascending=True)
            open_target_decoys = open_target_decoys[open_target_decoys['e_rank'].isin(range(1, tops_open + 1))]
        else:
            narrow_target_decoys = narrow_target_decoys[narrow_target_decoys['hit_rank'] == 1]
            open_target_decoys = open_target_decoys[open_target_decoys['hit_rank'].isin(range(1, tops_open + 1))]
        
        if tide_used == 'tide':
            #This is redundant currently, as original_target_sequence is being printed WITHOUT modification (but this is a current bug in tide-search, so we run this anyways)
            narrow_target_decoys['original_target_sequence'] = narrow_target_decoys['original_target_sequence'].str.replace("\\[|\\]|\\.|\d+", "")
            open_target_decoys['original_target_sequence'] = open_target_decoys['original_target_sequence'].str.replace("\\[|\\]|\\.|\d+", "")
            
        narrow_target_decoys.reset_index(drop = True, inplace = True)
        open_target_decoys.reset_index(drop = True, inplace = True)
        
    else:
        narrow_1 = narrow_target_decoys
        open_1 = open_target_decoys
        
        logging.info('Searching for decoy files in same directory.')
        sys.stderr.write('Searching for decoy files in same directory. \n')
        
        #search for corresponding decoy search files
        narrow_2 = pd.read_table(search_file_narrow.replace("target", "decoy"))
        open_2 = pd.read_table(search_file_open.replace("target", "decoy"))
        #standardizing column names
        narrow_2.columns = narrow_2.columns.str.strip().str.lower().str.replace(' ', '_', regex = False).str.replace('(', '', regex = False).str.replace(')', '', regex = False).str.replace('/', '_', regex = False).str.replace('.', '_', regex = False)
        open_2.columns = open_2.columns.str.strip().str.lower().str.replace(' ', '_', regex = False).str.replace('(', '', regex = False).str.replace(')', '', regex = False).str.replace('/', '_', regex = False).str.replace('.', '_', regex = False)
        
        logging.info('Successfully read in decoy search files.')
        sys.stderr.write('Successfully read in decoy search files. \n')
        
        if tide_used == 'tide' or tide_used == 'comet':
            #check to see if decoy search files indeed contain decoy peptides
            check_1 = any(narrow_target_decoys['protein_id'].str.contains(dcy_prefix))
            check_2 = any(open_target_decoys['protein_id'].str.contains(dcy_prefix))
            if (not check_1) or (not check_2):
                logging.info('No decoys were found in one of the decoy search files.')
                sys.exit('No decoys were found in one of the decoy search files. \n')
        else:
            check_1 = any(narrow_target_decoys['protein'].str.contains(dcy_prefix))
            check_2 = any(open_target_decoys['protein'].str.contains(dcy_prefix))
            if (not check_1) or (not check_2):
                logging.info('No decoys were found in one of the decoy search files.')
                sys.exit('No decoys were found in one of the decoy search files. \n')
        
        #Creating original_target_sequence column for target files too, so that they exist when we concatenate our target and decoy search files
        if tide_used == 'tide':
            narrow_1['original_target_sequence'] = narrow_1['sequence'].str.replace("\\[|\\]|\\.|\d+", "")
            open_1['original_target_sequence'] = open_1['sequence'].str.replace("\\[|\\]|\\.|\d+", "")
            narrow_2['original_target_sequence'] = narrow_2['original_target_sequence'].str.replace("\\[|\\]|\\.|\d+", "")
            open_2['original_target_sequence'] = open_2['original_target_sequence'].str.replace("\\[|\\]|\\.|\d+", "")

        if not len(narrow_1.columns) == len(narrow_2.columns):
            logging.error("Search files have different number of columns.")
            sys.exit("Search files have different number of columns.\n")
        if not len(open_1.columns) == len(open_2.columns):
            logging.error("Search files have different number of columns.")
            sys.exit("Search files have different number of columns.\n")
        if not len(narrow_1.columns.intersection(narrow_2.columns)) == narrow_1.columns.shape[0]:
            logging.error("Search files do not share the same column names.")
            sys.exit("Search files do not share the same column names.\n")
        if not len(open_1.columns.intersection(open_2.columns)) == open_1.columns.shape[0]:
            logging.error("Search files do not share the same column names.")
            sys.exit("Search files do not share the same column names.\n")
        
        sys.stderr.write("Concatenating separate searches.\n")
        logging.info("Concatenating separate searches.")
        
        #concatenating separate searches
        narrow_target_decoys = pd.concat([narrow_1, narrow_2])
        open_target_decoys = pd.concat([open_1, open_2])
        
        narrow_target_decoys = narrow_target_decoys.reset_index(drop = True)
        open_target_decoys = open_target_decoys.reset_index(drop = True)
        
        #randomly ordering
        narrow_target_decoys = narrow_target_decoys.sample(frac = 1)
        open_target_decoys = open_target_decoys.sample(frac = 1)
        
        #ordering the PSMs in each scan from best score to worst score
        if tide_used == 'tide':
            narrow_target_decoys = narrow_target_decoys.sort_values(by=['file', 'scan', 'xcorr_score'], ascending=False)
            open_target_decoys = open_target_decoys.sort_values(by=['file', 'scan', 'xcorr_score'], ascending=False)
        elif tide_used == 'comet' and score == 'xcorr':
            narrow_target_decoys = narrow_target_decoys.sort_values(by=['scan', 'charge', 'spectrum_neutral_mass', 'xcorr_score'], ascending=False)
            open_target_decoys = open_target_decoys.sort_values(by=['scan', 'charge', 'spectrum_neutral_mass', 'xcorr_score'], ascending=False)
        elif tide_used == 'comet' and score == 'e-value':
            narrow_target_decoys = narrow_target_decoys.sort_values(by=['scan', 'e-value'], ascending=True)
            open_target_decoys = open_target_decoys.sort_values(by=['scan', 'e-value'], ascending=True)
        else:
            narrow_target_decoys = narrow_target_decoys.sort_values(by=['scannum', 'hyperscore'], ascending=False)
            open_target_decoys = open_target_decoys.sort_values(by=['scannum', 'hyperscore'], ascending=False)
            
        narrow_target_decoys = narrow_target_decoys.reset_index(drop = True)
        open_target_decoys = open_target_decoys.reset_index(drop = True)
        
        #create a new ranking in the concantenated search and take the top 1 in narrow and top 'tops_open' in open.
        if tide_used == 'tide':
            narrow_target_decoys['xcorr_rank'] = narrow_target_decoys.groupby(["file", "scan"])["xcorr_score"].rank("first", ascending=False)
            narrow_target_decoys = narrow_target_decoys[narrow_target_decoys['xcorr_rank'] == 1]
            open_target_decoys['xcorr_rank'] = open_target_decoys.groupby(["file", "scan"])["xcorr_score"].rank("first", ascending=False)
            open_target_decoys = open_target_decoys[open_target_decoys['xcorr_rank'].isin(range(1, tops_open + 1))]
        elif tide_used == 'comet' and score == 'xcorr':
            narrow_target_decoys['xcorr_rank'] = narrow_target_decoys.groupby(["scan", "charge", "spectrum_neutral_mass"])["xcorr_score"].rank("first", ascending=False)
            narrow_target_decoys = narrow_target_decoys[narrow_target_decoys['xcorr_rank'] == 1]
            open_target_decoys['xcorr_rank'] = open_target_decoys.groupby(["scan", "charge", "spectrum_neutral_mass"])["xcorr_score"].rank("first", ascending=False)
            open_target_decoys = open_target_decoys[open_target_decoys['xcorr_rank'].isin(range(1, tops_open + 1))]
        elif tide_used == 'comet' and score == 'e-value':
            narrow_target_decoys['e_rank'] = narrow_target_decoys.groupby(["scan", "charge", "spectrum_neutral_mass"])["e-value"].rank("first", ascending=True)
            narrow_target_decoys = narrow_target_decoys[narrow_target_decoys['e_rank'] == 1]
            open_target_decoys['e_rank'] = open_target_decoys.groupby(["scan", "charge", "spectrum_neutral_mass"])["e-value"].rank("first", ascending=True)
            open_target_decoys = open_target_decoys[open_target_decoys['e_rank'].isin(range(1, tops_open + 1))]
        else:
            narrow_target_decoys['hit_rank'] = narrow_target_decoys.groupby('scan')["hyperscore"].rank("first", ascending=False)
            open_target_decoys['hit_rank'] = open_target_decoys.groupby('scan')["hyperscore"].rank("first", ascending=False)
            narrow_target_decoys = narrow_target_decoys[narrow_target_decoys['hit_rank'] == 1]
            open_target_decoys = open_target_decoys[open_target_decoys['hit_rank'].isin(range(1, tops_open + 1))]
            

        narrow_target_decoys.reset_index(drop = True, inplace = True)
        open_target_decoys.reset_index(drop = True, inplace = True)

    if not len(narrow_target_decoys.columns) == len(open_target_decoys.columns):
        logging.error("Search files have different number of columns.")
        sys.exit("Search files have different number of columns.\n")
    if not all([narrow_target_decoys.columns[i] == open_target_decoys.columns[i] for i in range(len(narrow_target_decoys.columns))]):
        logging.error("Search files do not share the same column names.")
        sys.exit("Search files do not share the same column names.\n")
    
    if tide_used == 'ms_fragger':
        sys.stderr.write("Rewriting peptide sequence in tide-search format. \n")
        logging.info("Rewriting peptide sequence in tide-search format.")
        #create the same "sequence" column as in tide-search
        narrow_target_decoys['sequence'] = create_peptides_with_mod(narrow_target_decoys['peptide'], narrow_target_decoys['modification_info'], static_mods)
        open_target_decoys['sequence'] = create_peptides_with_mod(open_target_decoys['peptide'], open_target_decoys['modification_info'], static_mods)
    elif tide_used == "comet":
        #Comet uses a reverse sequence for the peptides
        if account_mods:
            logging.info('Creating original_target_sequence column as in tide-search.')
            sys.stderr.write("Creating original_target_sequence column as in tide-search. \n")
            #creating original_target_sequence as in tide-search
            narrow_target_decoys['original_target_sequence'] = narrow_target_decoys['sequence']
            narrow_target_decoys.loc[narrow_target_decoys['protein_id'].str.contains(dcy_prefix), 'original_target_sequence'] = narrow_target_decoys[narrow_target_decoys['protein_id'].str.contains(dcy_prefix)].original_target_sequence.apply(lambda x: x[-2::-1] + x[-1])         
            open_target_decoys['original_target_sequence'] = open_target_decoys['sequence']
            open_target_decoys.loc[open_target_decoys['protein_id'].str.contains(dcy_prefix), 'original_target_sequence'] = open_target_decoys[open_target_decoys['protein_id'].str.contains(dcy_prefix)].original_target_sequence.apply(lambda x: x[-2::-1] + x[-1])
            
            sys.stderr.write("Rewriting peptide sequence in tide-search format. \n")
            logging.info("Rewriting peptide sequence in tide-search format.")
            #create the same "sequence" column as in tide-search
            narrow_target_decoys['sequence'] = narrow_target_decoys['modified_sequence'].apply(lambda x: re.search('\\.(.*)\\.', x).group(1))
            open_target_decoys['sequence'] = open_target_decoys['modified_sequence'].apply(lambda x: re.search('\\.(.*)\\.', x).group(1))
        else:
            sys.stderr.write("Rewriting peptide sequence in tide-search format. \n")
            logging.info("Rewriting peptide sequence in tide-search format.")
            #create the same "sequence" column as in tide-search
            narrow_target_decoys['sequence'] = narrow_target_decoys['modified_sequence'].apply(lambda x: re.search('\\.(.*)\\.', x).group(1))
            open_target_decoys['sequence'] = open_target_decoys['modified_sequence'].apply(lambda x: re.search('\\.(.*)\\.', x).group(1))
            
            sys.stderr.write("Pairing targets and decoys. \n")
            logging.info("Pairing targets and decoys.")
            #create a target_sequence column used to pair the targets and decoys together at the sequence-level
            narrow_target_decoys['target_sequence'] = narrow_target_decoys['sequence']
            narrow_target_decoys.loc[narrow_target_decoys['protein_id'].str.contains(dcy_prefix), 'target_sequence'] = reverse_sequence(narrow_target_decoys['sequence'][narrow_target_decoys['protein_id'].str.contains(dcy_prefix)])
            open_target_decoys['target_sequence'] = open_target_decoys['sequence']
            open_target_decoys.loc[open_target_decoys['protein_id'].str.contains(dcy_prefix), 'target_sequence'] = reverse_sequence(open_target_decoys['sequence'][open_target_decoys['protein_id'].str.contains(dcy_prefix)])

        
        
    
    target_decoys_all = filter_narrow_open(narrow_target_decoys, open_target_decoys, score, tops_open, thresh, n_processes, neighbour_remove, tide_used, static_mods)
    
    #check if there are any variable modifications
    any_mods = any('[' in pep for pep in target_decoys_all['sequence'])
    
    #import the peptide list only if required
    if account_mods or (not any_mods) or tide_used == "comet":
        peptide_list = ''
    else:
        logging.info("Reading peptide list.")
        sys.stderr.write("Reading peptide list. \n")
        if td_list == '':
            logging.info("Peptide list not found.")
            sys.exit("Peptide list not found. \n")
        else:
            peptide_list = pd.read_table(td_list)
    
    #print filtered search results if requested
    if return_filt_search:
        logging.info("Returning filtered search files in output directory.")
        sys.stderr.write("Returning filtered search files in output directory. \n")
        target_decoys_all_narrow = target_decoys_all[target_decoys_all['database'] == "narrow"]
        target_decoys_all_open = target_decoys_all[target_decoys_all['database'] == "open"]
        target_decoys_all_narrow.to_csv(output_dir + "/" + "filtered_" + search_file_narrow, header=True, index = False, sep = '\t')
        target_decoys_all_open.to_csv(output_dir + "/" + "filtered_" + search_file_open, header=True, index = False, sep = '\t')
    
    #create groups
    df = create_groups(target_decoys_all, narrow_target_decoys, peptide_list, dcy_prefix, K, tops_gw, score, account_mods, any_mods, precursor_bin_width, group_thresh, adaptive, min_group_size, n_top_groups, tide_used, print_group_pi0)
    
    #e-values score things in reverse
    if tide_used == 'comet' and score == 'e-value':
        df['winning_scores'] = -df['winning_scores']
    #apply group-walk
    results = group_walk(list(df['winning_scores']), list(df['labels']), list(df['all_group_ids']), K, return_frontier, correction)
    
    #report power for 1 and 5% FDR
    df['q_vals'] = results[0]
    power_1 = sum( ( df['q_vals'] <= 0.01 ) & ( df['labels'] == 1) )
    power_5 = sum( ( df['q_vals'] <= 0.05 ) & ( df['labels'] == 1) )
    
    logging.info(str(power_1) + " peptides discovered at the 1% FDR level.")
    logging.info(str(power_5) + " peptides discovered at the 5% FDR level.")
    sys.stderr.write(str(power_1) + " peptides discovered at the 1% FDR level. \n")
    sys.stderr.write(str(power_5) + " peptides discovered at the 5% FDR level. \n")
    
    #print Scan multiplicity: 
    if print_chimera:
        if tide_used == 'tide':
            scan_mult1 = df['scan'][ ( df['q_vals'] <= 0.01 ) & ( df['labels'] == 1) ].groupby([df['file'], df['scan']]).value_counts().value_counts()
            scan_mult1 = pd.DataFrame(scan_mult1)
            scan_mult1.columns = ['Count']
            scan_mult1.index.names = ['Scan multiplicity:']
            
            logging.info("Scan multiplicities among the discovered peptides at 1% FDR level:")
            logging.info(scan_mult1.to_string())
            sys.stderr.write("Scan multiplicities among the discovered peptides at 1% FDR level: \n")
            sys.stderr.write(scan_mult1.to_string() + "\n")
            
            scan_mult5 = df['scan'][ ( df['q_vals'] <= 0.05 ) & ( df['labels'] == 1) ].groupby([df['file'], df['scan']]).value_counts().value_counts()
            scan_mult5 = pd.DataFrame(scan_mult5)
            scan_mult5.columns = ['Count']
            scan_mult5.index.names = ['Scan multiplicity:']
            
            logging.info("Scan multiplicities among the discovered peptides at 5% FDR level:")
            logging.info(scan_mult5.to_string())
            sys.stderr.write("Scan multiplicities among the discovered peptides at 5% FDR level: \n")
            sys.stderr.write(scan_mult5.to_string() + "\n")
        elif tide_used == 'comet':
            scan_mult1 = df['scan_charge_sp_neutral_mass'][ ( df['q_vals'] <= 0.01 ) & ( df['labels'] == 1) ].value_counts().value_counts()
            scan_mult1 = pd.DataFrame(scan_mult1)
            scan_mult1.columns = ['Count']
            scan_mult1.index.names = ['Scan multiplicity:']

            logging.info("Scan multiplicities among the discovered peptides at 1% FDR level:")
            logging.info(scan_mult1.to_string())
            sys.stderr.write("Scan multiplicities among the discovered peptides at 1% FDR level: \n")
            sys.stderr.write(scan_mult1.to_string() + "\n")
            
            scan_mult5 = df['scan_charge_sp_neutral_mass'][ ( df['q_vals'] <= 0.05 ) & ( df['labels'] == 1) ].value_counts().value_counts()
            scan_mult5 = pd.DataFrame(scan_mult5)
            scan_mult5.columns = ['Count']
            scan_mult5.index.names = ['Scan multiplicity:']

            
            logging.info("Scan multiplicities among the discovered peptides at 5% FDR level:")
            logging.info(scan_mult5.to_string())
            sys.stderr.write("Scan multiplicities among the discovered peptides at 5% FDR level: \n")
            sys.stderr.write(scan_mult5.to_string() + "\n")
        else:
            scan_mult1 = df['scan'][ ( df['q_vals'] <= 0.01 ) & ( df['labels'] == 1) ].value_counts().value_counts()
            scan_mult1 = pd.DataFrame(scan_mult1)
            scan_mult1.columns = ['Count']
            scan_mult1.index.names = ['Scan multiplicity:']
            
            logging.info("Scan multiplicities among the discovered peptides at 1% FDR level:")
            logging.info(scan_mult1.to_string())
            sys.stderr.write("Scan multiplicities among the discovered peptides at 1% FDR level: \n")
            sys.stderr.write(scan_mult1.to_string() + "\n")
            
            scan_mult5 = df['scan'][ ( df['q_vals'] <= 0.05 ) & ( df['labels'] == 1) ].value_counts().value_counts()
            scan_mult5 = pd.DataFrame(scan_mult5)
            scan_mult5.columns = ['Count']
            scan_mult5.index.names = ['Scan multiplicity:']
            
            logging.info("Scan multiplicities among the discovered peptides at 5% FDR level:")
            logging.info(scan_mult5.to_string())
            sys.stderr.write("Scan multiplicities among the discovered peptides at 5% FDR level: \n")
            sys.stderr.write(scan_mult5.to_string() + "\n")
    
    
    if tide_used == 'comet' and score == 'e-value':
        df['winning_scores'] = -df['winning_scores'] 
        
    #output the q-values for the remaining winning PSMs
    if output_dir != './':
        if os.path.isdir(output_dir):
            df.to_csv(output_dir + "/" + file_root, header=True, index = False, sep = '\t')
        else:
            os.mkdir(output_dir)
            df.to_csv(output_dir + "/" + file_root, header=True, index = False, sep = '\t')
    if return_frontier:
        results[1].to_csv(output_dir + "/" + frontier_name, header=True, index = False, sep = '\t')
    
    end_time = time.time()
    
    logging.info("Elapsed time: " + str(round(end_time - start_time, 2)) + " s")
    sys.stderr.write("Elapsed time: " + str(round(end_time - start_time, 2)) + " s \n")
    
if __name__ == "__main__":
    main()    
    
    
