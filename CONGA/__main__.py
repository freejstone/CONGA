#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 15:27:59 2022

@author: jackfreestone

#Code that executes the group-walk algorithm
"""
import os 
import time
import datetime
import platform
import sys
import random
import logging
import re
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from CONGA.utils import CONGA_functions as cg
from . import version
__version__ = version.get_versions()['version']


USAGE = """USAGE: python3 -m CONGA [options] <narrow> <wide> <matching>

  This script implements the CONGA algorithm, including the
  creation of groups and the filtering procedure that eliminates
  similar peptide matches per scan. The first input file is the narrow
  search file output, the second input file is the open search file
  output, and the last input file contains the target-decoy peptide
  pairs. Output is a list of peptides discovered at a user-specified FDR-level.
  The input search files can be either tab-delimited .txt Tide-search files,
  tab-delimited .txt Comet search files, or tab-delimited .tsv files from 
  MSFragger. The target-decoy pairs should be given in a format similar
  to Tide-index's peptide list. Target-decoy peptide pairs is not required for
  Comet since Comet reverses target peptide sequences. The target-decoy peptide
  pairs is also not required for Tide since we can use the original.target.sequence
  column in the output search file (unless --account_mods F is used).

  Options:
      
    --output_dir <string> The file-path of the output
                          Default = './'.
                          
    --file_root <string>  The file prefix of the output files
                          Default = 'conga'.

    --FDR_threshold <value>     A user-specified threshold for which the reported
                                peptides will have FDR below this threshold.
                                Default = 0.01.

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
                          
    --score <string>      Either 'tailor_score', 'xcorr_score', 'e-value' or
                          'hyperscore'. The score that will be used in the 
                          peptide-level competition and subsequent Group 
                          construction and Group-walk algorithm. If 
                          'tailor_score', it is assumed the search files are 
                          derived from Tide. If 'xcorr_score', either Tide 
                          search or Comet is assumed to be used. If 'e-value',
                          Comet is assumed. If 'hyperscore' it is assumed the
                          search files are derived from MS-Fragger.
                          Default = 'tailor_score'.
    
    --account_mods <T|F>  To determine whether the group-walk algorithm
                          selects the best PSM among the equivalent
                          classes of peptides which are equal up to
                          variable modification, or not.
                          Default = T.
                          
    --competition_window <value>    A value (in m/z) used to agglomerate precursors
                                    according to their mass for h2h competition.
                                    The largest of the left- and right- isolation window
                                    offset should go here.
                                    Default = 2.
                          
    --return_extra_mods <T|F>       All top 1 PSM-ranked target peptides that
                                    are equal to a discovered peptide
                                    up to variable modification will
                                    be included in the file output.
                                    Default = F.
                          
    --precursor_bin_width <value>   To determine the size of the bins
                                    used to discretize the mass-
                                    differences between the sample
                                    and theoretical spectra.
                                    Default = 1.0005079/4.
                          
    --print_chimera <T|F>          To determine whether we print the number
                                   of scans that have more than 1 peptide
                                   discovered at the 1% and 5% FDR level
                                   to the log file.
                                   Default = T.
            
    --group_thresh <value>         The p-value threshold used to determine
                                   whether groups are merged or kept
                                   separate in the KS test.
                                   Default = 0.01.
                                   
    --print_group_pi0 <T|F>        To determine whether the group-level
                                   proportion of pi_0 is printed to the 
                                   log file.
                                   Default = T.
    
    --min_group_size <integer>     The number of multiples of K that is
                                   used to determine the minimum size of
                                   each group. See option --K.
                                   Default = 2.
    
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
                          Default = 0.05.
                          
    --return_filt_search <T|F>  Whether or not to return filtered narrow
                                and open search files.
                                Default = F.
                                
    --return_frontier <T|F>     The sequence of indices describing the
                                positions of the frontier used by Groupwalk
                                is returned as a .txt file to the output
                                directory.
                                Default = F.
    
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
                                
    --return_mass_mod_hist <T|F>   A matplotlib histogram of mass 
                                   modifications is returned to the user's
                                   directory.
                                   Default = F.
                                
    --dcy_prefix <string>       The prefix used for the decoy proteins.
                                Default = 'decoy_'

    --return_decoys <T|F>       Also report decoys used to estimate the 
                                number of false discoveries. 
                                Default = F. 
                                
    --overwrite <T|F>     Gives option to overwrite existing files in
                          directory.
                          Default = F.
                                
    --seed <int>          Set random seed.
                          Default = None.
            
"""
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
def print_info(command_line, output_dir, file_root, overwrite, account_mods, search_file_narrow, search_file_open):
    #check if output directory exists, if not create and store log file there.
    if os.path.isdir(output_dir):
        if os.path.exists(path = output_dir + "/" + file_root + ".log.txt") and overwrite:
            os.remove(output_dir + "/" + file_root + ".log.txt")
            logging.basicConfig(filename=output_dir + "/" + file_root + ".log.txt", level=logging.DEBUG, format = '%(levelname)s: %(message)s')
        elif os.path.exists(path = output_dir + "/" + file_root + ".log.txt") and not overwrite:
            log_file = output_dir + "/" + file_root + ".log.txt"
            sys.exit("The file %s already exists and cannot be overwritten. Use --overwrite T to replace or choose a different output file name. \n" %(log_file))
        else:
            logging.basicConfig(filename=output_dir + "/" + file_root + ".log.txt", level=logging.DEBUG, format = '%(levelname)s: %(message)s')
    else:
        os.mkdir(output_dir)
        logging.basicConfig(filename=output_dir + "/" + file_root  + ".log.txt", level=logging.DEBUG, format = '%(levelname)s: %(message)s')
    
    #print CPU info
    logging.info('CPU: ' + str(platform.platform()))
    sys.stderr.write('CPU: ' + str(platform.platform()) + " \n")
    
     #print version
    logging.info('Version: ' + str(__version__))
    sys.stderr.write('Version: ' + str(__version__) + " \n")
    
    #print date time info
    logging.info(str(datetime.datetime.now()))
    sys.stderr.write(str(datetime.datetime.now()) + " \n")
    
    #print command used
    logging.info('Command used: ' + command_line)
    sys.stderr.write('Command used: ' + command_line + "\n")
    
    sys.stderr.write("Successfully read in arguments. \n")
    logging.info("Successfully read in arguments")
    
    if not account_mods:
        logging.warning("No longer accounting for variable modification. FDR control not guaranteed if variable modifications exist.")
        sys.stderr.write("No longer accounting for variable modification. FDR control not guaranteed if variable modifications exist. \n")
    
    if (not os.path.isfile(search_file_narrow)) or (not os.path.isfile(search_file_open)):
        logging.info("One of the search files does not exist.")
        sys.exit("One of the search files does not exist. \n")
###############################################################################
def read_search_files(search_file_narrow, search_file_open):
    #read in files and standardize column names
    #Comet has an initial comment starting with 'C' -- all others do not
    with open(search_file_narrow) as f:
        first_line_narrow = f.readline()
    with open(search_file_open) as f:
        first_line_open = f.readline()
    
    if 'Comet' in first_line_narrow:
        narrow_target_decoys = pd.read_table(search_file_narrow, skiprows = 1)
    else:
        narrow_target_decoys = pd.read_table(search_file_narrow)
    
    if 'Comet' in first_line_open:
        open_target_decoys = pd.read_table(search_file_open, skiprows = 1)
    else:
        open_target_decoys = pd.read_table(search_file_open)
    
    narrow_target_decoys.columns = narrow_target_decoys.columns.str.strip().str.lower().str.replace(' ', '_', regex = False).str.replace('(', '', regex = False).str.replace(')', '', regex = False).str.replace('/', '_', regex = False).str.replace('.', '_', regex = False)
    open_target_decoys.columns = open_target_decoys.columns.str.strip().str.lower().str.replace(' ', '_', regex = False).str.replace('(', '', regex = False).str.replace(')', '', regex = False).str.replace('/', '_', regex = False).str.replace('.', '_', regex = False)
    
    if ('xcorr' in narrow_target_decoys.columns):
        narrow_target_decoys.rename(columns = {'xcorr':'xcorr_score'}, inplace = True)
        open_target_decoys.rename(columns = {'xcorr':'xcorr_score'}, inplace = True)
    return(narrow_target_decoys, open_target_decoys)
###############################################################################
def rename_features(narrow_target_decoys, open_target_decoys, tide_used):
    #rename protein id for consistency
    if tide_used == 'ms_fragger':
        narrow_target_decoys.rename(columns = {'protein':'protein_id'}, inplace = True)
        open_target_decoys.rename(columns = {'protein':'protein_id'}, inplace = True)
        if ('alternative_proteins' in narrow_target_decoys.columns):
            narrow_target_decoys['protein_id'] = narrow_target_decoys['protein_id'].astype(str) + ',' + narrow_target_decoys['alternative_proteins'].astype(str)
        if ('alternative_proteins' in open_target_decoys.columns):
            open_target_decoys['protein_id'] = open_target_decoys['protein_id'].astype(str) + ',' + open_target_decoys['alternative_proteins'].astype(str)
        
    #making standalone comet agree with crux comet 
    if tide_used == 'comet':
        if ('protein' in narrow_target_decoys.columns):
            narrow_target_decoys.rename(columns = {'protein':'protein_id'}, inplace = True)
            open_target_decoys.rename(columns = {'protein':'protein_id'}, inplace = True)
        if ('exp_neutral_mass' in narrow_target_decoys.columns):
            narrow_target_decoys.rename(columns = {'exp_neutral_mass':'spectrum_neutral_mass'}, inplace = True)
            open_target_decoys.rename(columns = {'exp_neutral_mass':'spectrum_neutral_mass'}, inplace = True)
        if ('plain_peptide' in narrow_target_decoys.columns):
            narrow_target_decoys.rename(columns = {'plain_peptide':'sequence'}, inplace = True)
            open_target_decoys.rename(columns = {'plain_peptide':'sequence'}, inplace = True)
        if ('modified_peptide' in narrow_target_decoys.columns):
            narrow_target_decoys.rename(columns = {'modified_peptide':'modified_sequence'}, inplace = True)
            open_target_decoys.rename(columns = {'modified_peptide':'modified_sequence'}, inplace = True)
        if ('calc_neutral_mass' in narrow_target_decoys.columns):
            narrow_target_decoys.rename(columns = {'calc_neutral_mass':'peptide_mass'}, inplace = True)
            open_target_decoys.rename(columns = {'calc_neutral_mass':'peptide_mass'}, inplace = True)
    
    return(narrow_target_decoys, open_target_decoys)
###############################################################################
def main():
    global USAGE
    
    start_time = time.time()
    
    # Set default values for parameters.
    K = 40
    FDR_threshold = 0.01
    tops_gw = 2
    tops_open = 5
    score = 'tailor_score'
    account_mods = True
    competition_window = 2
    return_extra_mods = False
    precursor_bin_width = 1.0005079/4
    adaptive = True
    print_chimera = True
    group_thresh = 0.01
    print_group_pi0 = True
    min_group_size = 2
    n_top_groups = 4
    neighbour_remove = True
    thresh = 0.05
    return_filt_search = False
    correction = 1
    n_processes = 1
    return_frontier = False
    output_dir = '.'
    file_root = 'conga'
    static_mods = {'C':57.02146}
    return_mass_mod_hist = False
    dcy_prefix = 'decoy_'
    return_decoys = False
    overwrite = False
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
        elif (next_arg == "--FDR_threshold"):
            FDR_threshold = float(sys.argv[0])
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
        elif (next_arg == "--competition_window"):
            competition_window = float(sys.argv[0])
            sys.argv = sys.argv[1:]
        elif (next_arg == "--return_extra_mods"):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                return_extra_mods = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                return_extra_mods = False
            else:
                sys.stderr.write("Invalid argument for --return_extra_mods")
                sys.exit(1)
            sys.argv = sys.argv[1:]
        elif (next_arg == "--precursor_bin_width"):
            precursor_bin_width = float(sys.argv[0])  
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
        elif (next_arg == "--static_mods"):
            static_mods = cg.parse_static_mods(sys.argv[0])
            sys.argv = sys.argv[1:]
        elif (next_arg == "--return_mass_mod_hist"):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                return_mass_mod_hist = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                return_mass_mod_hist = False
            else:
                sys.stderr.write("Invalid argument for --return_mass_mod_hist")
                sys.exit(1)
            sys.argv = sys.argv[1:]
        elif (next_arg == "--dcy_prefix"):
            dcy_prefix = str(sys.argv[0])
            sys.argv = sys.argv[1:]
        elif (next_arg == "--return_decoys"):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                return_decoys = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                return_decoys = False
            else:
                sys.stderr.write("Invalid argument for --return_decoys")
                sys.exit(1)
            sys.argv = sys.argv[1:]
        elif (next_arg == "--overwrite"):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                overwrite = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                overwrite = False
            else:
                sys.stderr.write("Invalid argument for --overwrite")
                sys.exit(1)
            sys.argv = sys.argv[1:]
        elif (next_arg == '--seed'):
            seed = int(sys.argv[0])
            sys.argv = sys.argv[1:]
        else:
            sys.stderr.write("Invalid option (%s)" % next_arg)
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
        sys.stderr.write('Version: ' + str(__version__) + " \n")
        sys.stderr.write(USAGE)
        sys.exit(1)
    
    #setting seed for reproducibility
    if type(seed) == int:
        random.seed(seed)
        np.random.seed(seed)
    
    #print meta information, checking directory and printing warnings
    print_info(command_line, output_dir, file_root, overwrite, account_mods, search_file_narrow, search_file_open)
    
    sys.stderr.write("Reading in search files. \n")
    logging.info("Reading in search files.")
    
    #reading in search files
    narrow_target_decoys, open_target_decoys = read_search_files(search_file_narrow, search_file_open)
    
    sys.stderr.write("Successfully read in search files. \n")
    logging.info("Successfully read in search files.")
    
    #check if score matches with search file type given by the user
    #and establishes score and search file type
    check_1 = (score in narrow_target_decoys.columns)
    check_2 = (score in open_target_decoys.columns)
    if check_1 and check_2:
        sys.stderr.write(score + " successfully found in search files. \n")
        logging.info(score + " successfully found in search files. \n")
    else:
        logging.info(score + " was not found in at least one of the search files. Please specify the score using --score option. \n")
        sys.exit(score + " was not found in at least one of the search files. Please specify the score using --score option. \n")
    if ('hyperscore' in narrow_target_decoys.columns) and ('hyperscore' in open_target_decoys.columns):
        tide_used = 'ms_fragger'
        sys.stderr.write("MS-Fragger search files detected. \n")
        logging.info("MS-Fragger search files detected.")
    elif ('e-value' in narrow_target_decoys.columns) and ('e-value' in open_target_decoys.columns):
        tide_used = 'comet'
        sys.stderr.write("Comet search files detected. \n")
        logging.info("Comet search files detected.")
    else:
        tide_used = 'tide'
        sys.stderr.write("Tide search files detected. \n")
        logging.info("Tide search files detected.")     
    
    #standardising feature names
    narrow_target_decoys, open_target_decoys = rename_features(narrow_target_decoys, open_target_decoys, tide_used)
    
    #check to see if concantenated search file or separate search file used
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
        
    if concat:
        #take the top 1 PSM for narrow search and top 'tops_open' PSM for open search
        if tide_used == 'tide':
            narrow_target_decoys = narrow_target_decoys[narrow_target_decoys['xcorr_rank'] == 1]
            open_target_decoys = open_target_decoys[open_target_decoys['xcorr_rank'].isin(range(1, tops_open + 1))]
        elif (tide_used == 'comet'):
            #native comet does not have xcorr_rank, so do this regardless
            narrow_target_decoys['xcorr_rank'] = narrow_target_decoys.groupby(["scan", "charge", "spectrum_neutral_mass"])["xcorr_score"].rank("first", ascending=False)
            narrow_target_decoys = narrow_target_decoys[narrow_target_decoys['xcorr_rank'] == 1]
            open_target_decoys['xcorr_rank'] = open_target_decoys.groupby(["scan", "charge", "spectrum_neutral_mass"])["xcorr_score"].rank("first", ascending=False)
            open_target_decoys = open_target_decoys[open_target_decoys['xcorr_rank'].isin(range(1, tops_open + 1))]
        else:
            narrow_target_decoys = narrow_target_decoys[narrow_target_decoys['hit_rank'] == 1]
            open_target_decoys = open_target_decoys[open_target_decoys['hit_rank'].isin(range(1, tops_open + 1))]
        
        if tide_used == 'tide':
            #This is redundant currently, as original_target_sequence is being printed WITHOUT modification (but this is may change in tide-search, so we run this anyways)
            narrow_target_decoys['original_target_sequence'] = narrow_target_decoys['original_target_sequence'].str.replace("\\[|\\]|\\.|\\d+", "", regex = True)
            open_target_decoys['original_target_sequence'] = open_target_decoys['original_target_sequence'].str.replace("\\[|\\]|\\.|\\d+", "", regex = True)
            
        narrow_target_decoys.reset_index(drop = True, inplace = True)
        open_target_decoys.reset_index(drop = True, inplace = True)
        
    else:
        narrow_1 = narrow_target_decoys
        open_1 = open_target_decoys
        
        logging.info('Searching for decoy files in same directory.')
        sys.stderr.write('Searching for decoy files in same directory. \n')
        
        if (not os.path.isfile(search_file_narrow.replace("target", "decoy"))) or (not os.path.isfile(search_file_open.replace("target", "decoy"))):
            logging.info("One of the decoy files does not exist.")
            sys.exit("One of the decoy files does not exist. \n")
        
        #search for corresponding decoy search files
        search_file_narrow_2 = search_file_narrow.replace("target", "decoy") 
        search_file_open_2 = search_file_open.replace("target", "decoy") 
        narrow_2, open_2 = read_search_files(search_file_narrow_2, search_file_open_2)
        
        logging.info('Successfully read in decoy search files.')
        sys.stderr.write('Successfully read in decoy search files. \n')
        
        #standardising feature names
        narrow_2, open_2 = rename_features(narrow_2, open_2, tide_used)
            
        check_1 = any(narrow_2['protein_id'].str.contains(dcy_prefix))
        check_2 = any(open_2['protein_id'].str.contains(dcy_prefix))
        if (not check_1):
            logging.info('No decoys were found in %s.' % search_file_narrow_2)
            sys.exit('No decoys were found in %s. \n' % search_file_narrow_2)
  
        #Creating original_target_sequence column for target files too, so that they exist when we concatenate our target and decoy search files
        if tide_used == 'tide':
            narrow_1['original_target_sequence'] = narrow_1['sequence'].str.replace("\\[|\\]|\\.|\\d+", "", regex = True)
            open_1['original_target_sequence'] = open_1['sequence'].str.replace("\\[|\\]|\\.|\\d+", "", regex = True)
            narrow_2['original_target_sequence'] = narrow_2['original_target_sequence'].str.replace("\\[|\\]|\\.|\\d+", "", regex = True)
            open_2['original_target_sequence'] = open_2['original_target_sequence'].str.replace("\\[|\\]|\\.|\\d+", "", regex = True)

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
        elif tide_used == 'comet':
            narrow_target_decoys = narrow_target_decoys.sort_values(by=['scan', 'xcorr_score'], ascending=False)
            open_target_decoys = open_target_decoys.sort_values(by=['scan', 'xcorr_score'], ascending=False)
        else:
            narrow_target_decoys = narrow_target_decoys.sort_values(by=['scannum', 'hyperscore'], ascending=False)
            open_target_decoys = open_target_decoys.sort_values(by=['scannum', 'hyperscore'], ascending=False)
            
        narrow_target_decoys = narrow_target_decoys.reset_index(drop = True)
        open_target_decoys = open_target_decoys.reset_index(drop = True)
        
        #create a new ranking in the concantenated search and take the top 1 in narrow and top 'tops_open' in open.
        if tide_used == 'tide':
            narrow_target_decoys['xcorr_rank'] = narrow_target_decoys.groupby(["file", "scan", "charge", "spectrum_neutral_mass"])["xcorr_score"].rank("first", ascending=False)
            narrow_target_decoys = narrow_target_decoys[narrow_target_decoys['xcorr_rank'] == 1]
            open_target_decoys['xcorr_rank'] = open_target_decoys.groupby(["file", "scan", "charge", "spectrum_neutral_mass"])["xcorr_score"].rank("first", ascending=False)
            open_target_decoys = open_target_decoys[open_target_decoys['xcorr_rank'].isin(range(1, tops_open + 1))]
        elif tide_used == 'comet':
            narrow_target_decoys['xcorr_rank'] = narrow_target_decoys.groupby(["scan", "charge", "spectrum_neutral_mass"])["xcorr_score"].rank("first", ascending=False)
            narrow_target_decoys = narrow_target_decoys[narrow_target_decoys['xcorr_rank'] == 1]
            open_target_decoys['xcorr_rank'] = open_target_decoys.groupby(["scan", "charge", "spectrum_neutral_mass"])["xcorr_score"].rank("first", ascending=False)
            open_target_decoys = open_target_decoys[open_target_decoys['xcorr_rank'].isin(range(1, tops_open + 1))]
        else:
            narrow_target_decoys['hit_rank'] = narrow_target_decoys.groupby(['scannum'])["hyperscore"].rank("first", ascending=False)
            open_target_decoys['hit_rank'] = open_target_decoys.groupby(['scannum'])["hyperscore"].rank("first", ascending=False)
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
        narrow_target_decoys['sequence'] = cg.create_peptides_with_mod(narrow_target_decoys['peptide'], narrow_target_decoys['modification_info'], static_mods)
        open_target_decoys['sequence'] = cg.create_peptides_with_mod(open_target_decoys['peptide'], open_target_decoys['modification_info'], static_mods)
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
        narrow_target_decoys['sequence'] = narrow_target_decoys['sequence'].str.replace("n", "").str.replace("c", "")
        open_target_decoys['sequence'] = open_target_decoys['sequence'].str.replace("n", "").str.replace("c", "")
        narrow_target_decoys.loc[narrow_target_decoys['sequence'].str.startswith('['), 'sequence'] = cg.check_n_term(narrow_target_decoys['sequence'][narrow_target_decoys['sequence'].str.startswith('[')])
        open_target_decoys.loc[open_target_decoys['sequence'].str.startswith('['), 'sequence'] = cg.check_n_term(open_target_decoys['sequence'][open_target_decoys['sequence'].str.startswith('[')])
            
        if not account_mods:
            sys.stderr.write("Pairing targets and decoys. \n")
            logging.info("Pairing targets and decoys.")
            #create a target_sequence column used to pair the targets and decoys together at the sequence-level
            narrow_target_decoys['original_target_sequence'] = narrow_target_decoys['sequence']
            narrow_target_decoys.loc[narrow_target_decoys['protein_id'].str.contains(dcy_prefix), 'original_target_sequence'] = narrow_target_decoys[narrow_target_decoys['protein_id'].str.contains(dcy_prefix)].apply(lambda x: cg.reverse_sequence(x.sequence, x.modifications), axis = 1)
            open_target_decoys['original_target_sequence'] = open_target_decoys['sequence']
            open_target_decoys.loc[open_target_decoys['protein_id'].str.contains(dcy_prefix), 'original_target_sequence'] = open_target_decoys[open_target_decoys['protein_id'].str.contains(dcy_prefix)].apply(lambda x: cg.reverse_sequence(x.sequence, x.modifications), axis = 1)
    
    target_decoys_all = cg.filter_narrow_open(narrow_target_decoys, open_target_decoys, score, thresh, n_processes, neighbour_remove, tide_used, static_mods)
      
    if tide_used == 'comet':
        target_decoys_all['check_protein'] = target_decoys_all['protein_id'].apply(lambda x: cg.del_protein(x, dcy_prefix))
        target_decoys_all = target_decoys_all[target_decoys_all['check_protein'] == False]
        target_decoys_all = target_decoys_all.drop('check_protein', axis = 1)
    
    #check if there are any variable modifications
    any_mods = any('[' in pep for pep in target_decoys_all['sequence'])
    
    #import the peptide list only if required
    if ((tide_used == 'tide') and (account_mods or not any_mods)) or tide_used == "comet":
        peptide_list = ''
    else:
        logging.info("Reading peptide list.")
        sys.stderr.write("Reading peptide list. \n")
        if td_list == '':
            logging.info("Peptide list not provided.")
            sys.exit("Peptide list not provided. \n")
        else:
            peptide_list = pd.read_table(td_list)
            if 'decoy(s)' in peptide_list.columns: #too handle newest version of tide-index
                peptide_list.rename(columns = {'decoy(s)':'decoy'}, inplace = True)
    
    #print filtered search results if requested
    if return_filt_search:
        logging.info("Returning filtered search files in output directory.")
        sys.stderr.write("Returning filtered search files in output directory. \n")
        target_decoys_all_narrow = target_decoys_all[target_decoys_all['database'] == "narrow"]
        target_decoys_all_open = target_decoys_all[target_decoys_all['database'] == "open"]
        target_decoys_all_narrow.to_csv(output_dir + "/" + file_root + ".narrow.filtered.txt", header=True, index = False, sep = '\t')
        target_decoys_all_open.to_csv(output_dir + "/" + file_root + ".open.filtered.txt", header=True, index = False, sep = '\t')
    
    if type(peptide_list) != str:
        peptide_list = peptide_list[~(peptide_list['target'] == peptide_list['decoy'])] #where the two peptides agree
        peptide_list = peptide_list.drop_duplicates(['target'])
    
    #create groups
    df = cg.create_groups(target_decoys_all, narrow_target_decoys, peptide_list, dcy_prefix, K, tops_gw, score, account_mods, any_mods, precursor_bin_width, group_thresh, adaptive, min_group_size, n_top_groups, tide_used, print_group_pi0, competition_window)
    
    #e-values score things in reverse so reverse this for groupwalk
    if tide_used == 'comet' and score == 'e-value':
        df['winning_scores'] = -df['winning_scores']
        
    #apply group-walk
    results = cg.group_walk(list(df['winning_scores']), list(df['labels']), list(df['all_group_ids']), K, return_frontier, correction)
    
    df['q_vals'] = results[0]
    df = df.sort_values(by = ['q_vals', 'winning_scores'], ascending = [True, False])
    df = df.drop('bins', axis = 1)
    
    #report power for 1 and 5% FDR (we do this before the reporting of extra variable modifications)
    df.reset_index(drop = True, inplace = True)
    
    power_1 = sum( ( df['q_vals'] <= 0.01 ) & ( df['labels'] == 1) )
    power_5 = sum( ( df['q_vals'] <= 0.05 ) & ( df['labels'] == 1) )
    
    logging.info(str(power_1) + " peptides discovered at the 1% FDR level.")
    logging.info(str(power_5) + " peptides discovered at the 5% FDR level.")
    sys.stderr.write(str(power_1) + " peptides discovered at the 1% FDR level. \n")
    sys.stderr.write(str(power_5) + " peptides discovered at the 5% FDR level. \n")
    
    #print Scan multiplicity: 
    if print_chimera:
            if power_1 > 0:
                if tide_used == 'tide':
                    scan_mult1 = df['scan'][ ( df['q_vals'] <= 0.01 ) & ( df['labels'] == 1) ].groupby([df['file'], df['scan'], df['charge'], df['spectrum_neutral_mass']]).value_counts().value_counts()
                elif tide_used == 'comet' or tide_used == 'ms_fragger':
                    scan_mult1 = df['scan'][ ( df['q_vals'] <= 0.01 ) & ( df['labels'] == 1) ].groupby([df['scan'], df['charge'], df['spectrum_neutral_mass']]).value_counts().value_counts()
                scan_mult1 = pd.DataFrame(scan_mult1)
                scan_mult1.columns = ['Count']
                scan_mult1.index.names = ['Scan multiplicity:']
                
                logging.info("Scan multiplicities among the discovered peptides at 1% FDR level:")
                logging.info(scan_mult1.to_string())
                sys.stderr.write("Scan multiplicities among the discovered peptides at 1% FDR level: \n")
                sys.stderr.write(scan_mult1.to_string() + "\n")
            if power_5 > 0:
                if tide_used == 'tide':
                    scan_mult5 = df['scan'][ ( df['q_vals'] <= 0.05 ) & ( df['labels'] == 1) ].groupby([df['file'], df['scan'], df['charge'], df['spectrum_neutral_mass']]).value_counts().value_counts()
                elif tide_used == 'comet' or tide_used == 'ms_fragger':
                    scan_mult5 = df['scan'][ ( df['q_vals'] <= 0.05 ) & ( df['labels'] == 1) ].groupby([df['scan'], df['charge'], df['spectrum_neutral_mass']]).value_counts().value_counts()
                scan_mult5 = pd.DataFrame(scan_mult5)
                scan_mult5.columns = ['Count']
                scan_mult5.index.names = ['Scan multiplicity:']
                
                logging.info("Scan multiplicities among the discovered peptides at 5% FDR level:")
                logging.info(scan_mult5.to_string())
                sys.stderr.write("Scan multiplicities among the discovered peptides at 5% FDR level: \n")
                sys.stderr.write(scan_mult5.to_string() + "\n")
       
    df = df[df['q_vals'] <= FDR_threshold]
    
    if return_extra_mods and (not any_mods):
        logging.info("--return_extra_mods was set to T however no variable modifications were found.")
        sys.stderr.write("--return_extra_mods was set to T however no variable modifications were found. \n")
    elif return_extra_mods and any_mods:
        logging.info("Also reporting sequences with other variable modifications (--return_extra_mods T).")
        sys.stderr.write("Also reporting sequences with other variable modifications (--return_extra_mods T). \n")
        
        #gets the STEM sequence
        df['sequence_without_mods_td'] = df['winning_peptides'].str.replace("\\[|\\]|\\.|\\d+", "", regex = True) + "_" + df['labels'].astype(str)

        if type(peptide_list) == str:
            target_decoys_all['labels'] = (~(target_decoys_all['protein_id'].str.contains(dcy_prefix))).astype(int)
            target_decoys_all.loc[target_decoys_all['labels'] == 0, 'labels'] = -1
        else:
            target_decoys_all['labels'] = target_decoys_all['sequence'].isin(peptide_list['target'])
            target_decoys_all['labels'].replace({True:1, False:-1}, inplace = True)

        target_decoys_all['sequence_without_mods_td'] = target_decoys_all['sequence'].str.replace("\\[|\\]|\\.|\\d+", "", regex = True) + "_" + target_decoys_all['labels'].astype(str)
        
        #get all target PSMs whose sequence equals to one in df up to modification but isn't exactly equal to the same sequence in df.
        target_decoys_all_sub = target_decoys_all[( target_decoys_all['sequence_without_mods_td'].isin(df['sequence_without_mods_td']) ) & ( target_decoys_all['labels'] == 1 ) & ~( target_decoys_all['sequence'].isin(df['winning_peptides']) )].copy()
        
        if score != 'e-value':
            target_decoys_all_sub = target_decoys_all_sub.sort_values(by=[score], ascending=False)
        else:
            target_decoys_all_sub = target_decoys_all_sub.sort_values(by=[score], ascending=True)
        
        #print all top 1 PSMs that differ by variable mods, even if there are double ups
        #target_decoys_all_sub = target_decoys_all_sub.drop_duplicates(subset = ['sequence'])
        
        if (tide_used == 'tide') or (tide_used == 'comet'):
            rank = 'xcorr_rank'
        else:
            rank = 'hit_rank'
    
        target_decoys_all_sub = target_decoys_all_sub[target_decoys_all_sub[rank] == 1]
        
        if len(target_decoys_all_sub) > 0:
        
            if tide_used == 'tide':
                winning_scores = target_decoys_all_sub[score]
                labels = target_decoys_all_sub['labels']
                winning_peptides = target_decoys_all_sub['sequence']
                rank = target_decoys_all_sub['xcorr_rank']
                delta_mass = target_decoys_all_sub['spectrum_neutral_mass'] - target_decoys_all_sub['peptide_mass']
                database = target_decoys_all_sub['database']
                scan = target_decoys_all_sub['scan']
                file = target_decoys_all_sub['file']
                charge = target_decoys_all_sub['charge']
                spectrum_neutral_mass = target_decoys_all_sub['spectrum_neutral_mass']
                flanking_aa = target_decoys_all_sub['flanking_aa']
                protein = target_decoys_all_sub['protein_id']
                df_extra = pd.DataFrame(zip(winning_scores, labels, delta_mass, winning_peptides, rank, database, charge, spectrum_neutral_mass, flanking_aa, scan, file), columns = ['winning_scores', 'labels', 'delta_mass', 'winning_peptides', 'rank', 'database', 'charge', 'spectrum_neutral_mass', 'flanking_aa', 'scan', 'file'])
                
            elif tide_used == 'comet':
                winning_scores = target_decoys_all_sub[score]
                labels = target_decoys_all_sub['labels']
                winning_peptides = target_decoys_all_sub['sequence']
                rank = target_decoys_all_sub['xcorr_rank']
                delta_mass = target_decoys_all_sub['spectrum_neutral_mass'] - target_decoys_all_sub['peptide_mass']
                database = target_decoys_all_sub['database']
                scan = target_decoys_all_sub['scan']
                charge = target_decoys_all_sub['charge']
                spectrum_neutral_mass = target_decoys_all_sub['spectrum_neutral_mass']
                protein = target_decoys_all_sub['protein_id']
                df_extra = pd.DataFrame(zip(winning_scores, labels, delta_mass, winning_peptides, rank, database, charge, spectrum_neutral_mass, scan), columns = ['winning_scores', 'labels', 'delta_mass', 'winning_peptides', 'rank', 'database', 'charge', 'spectrum_neutral_mass', 'scan'])
            
            elif tide_used == 'ms_fragger':
                winning_scores = target_decoys_all_sub['hyperscore']
                labels = target_decoys_all_sub['labels']
                winning_peptides = target_decoys_all_sub['sequence']
                rank = target_decoys_all_sub['hit_rank']
                delta_mass = target_decoys_all_sub['precursor_neutral_mass'] - target_decoys_all_sub['calc_neutral_pep_mass']
                database = target_decoys_all_sub['database']
                scan = target_decoys_all_sub['scannum']
                spectrum_neutral_mass = target_decoys_all_sub['precursor_neutral_mass']
                protein = target_decoys_all_sub['protein_id']
                df_extra = pd.DataFrame(zip(winning_scores, labels, delta_mass, winning_peptides, rank, database, spectrum_neutral_mass, scan, protein), columns = ['winning_scores', 'labels', 'delta_mass', 'winning_peptides', 'rank', 'database', 'spectrum_neutral_mass', 'scan', 'protein'])
            
            df['originally_discovered'] = True
            df_extra['originally_discovered'] = False
            df = pd.concat([df, df_extra])
        
        else:
            logging.info("There are no other sequences with variable modifications that are the top match for an experimentral spectrum.")
            sys.stderr.write("There are no other sequences with variable modifications that are the top match for an experimentral spectrum. \n")
        
        
    if tide_used == 'comet' and score == 'e-value':
        df['winning_scores'] = -df['winning_scores']
        
    #changing the name of some of the features for user-readability
    #doing this here since I still want to use more detailed names for creating discrepancies when reviewing code e.g winning_scores vs. score or winning_peptides vs. peptide
    df.pop('all_group_ids') #I don't believe this will be needed in the output for users honestly
    df_decoys = df[df['labels'] == -1].copy()
    df = df[df['labels'] == 1]
    df.pop('labels')
    df.rename(columns = {'winning_scores':'score', 'winning_peptides':'peptide', 'q_vals':'q_value', 'database':'search_file'}, inplace = True)
    if return_decoys:
        df_decoys.pop('labels')
        df_decoys.rename(columns = {'winning_scores':'score', 'winning_peptides':'peptide', 'q_vals':'q_value', 'database':'search_file'}, inplace = True)
        df_decoys = df_decoys.round({'delta_mass':4})

    #ordering by q-value then peptides second
    if 'originally_discovered' in df.columns:
        df['sequence_without_mods_td'] = df['peptide'].str.replace("\\[|\\]|\\.|\\d+", "", regex = True)
        df['min_q_vals'] = df.groupby(['sequence_without_mods_td']).q_value.transform('min')
        df = df.sort_values(by = ['min_q_vals', 'originally_discovered'], ascending = [True, False])
        df.pop('sequence_without_mods_td')
        df.pop('min_q_vals')
    
    #rounding the mass differences
    df = df.round({'delta_mass':4})

    #dropping q_values
    df.pop('q_value')
        
    #output the discovered peptides
    if return_extra_mods and any_mods:
        logging.info("Writing peptides at user-specified FDR level along with extrato directory.")
        sys.stderr.write("Writing peptides at user-specified FDR level to directory. \n")
    else:
        logging.info("Writing peptides at user-specified FDR level to directory.")
        sys.stderr.write("Writing peptides at user-specified FDR level to directory. \n")
    
    for aa in static_mods.keys():
        if aa == 'I' or aa == 'L': #Problem if someone uses two different static mods for I and L
            get_indices = aa_table.mass[aa_table.aa == 'I, L'].index
            if len(get_indices) > 0:
                aa_table.mass.at[get_indices[0]] += static_mods[aa]
        else:
            get_indices = aa_table.mass[aa_table.aa == aa].index
            if len(get_indices) > 0:
                aa_table.mass.at[get_indices[0]] += static_mods[aa]

    #dropping flanking_aa
    if tide_used == 'tide':
        if df.shape[0] > 0:
            df['flag'] = df.apply(cg.get_amino_acid_to_warn, axis = 1)
        else:
            df['flag'] = pd.Series([], dtype = 'object')
        df.pop('flanking_aa')
    
    df['modification_info'] = df['peptide'].apply(cg.get_modification_info)
    
    if output_dir != './':
        if os.path.isdir(output_dir):
            df.to_csv(output_dir + "/" + file_root + ".target.txt", header=True, index = False, sep = '\t')
            if return_decoys:
                df_decoys.to_csv(output_dir + "/" + file_root + ".decoy.txt", header=True, index = False, sep = '\t')
        else:
            os.mkdir(output_dir)
            df.to_csv(output_dir + "/" + file_root + ".target.txt", header=True, index = False, sep = '\t')
            if return_decoys:
                df_decoys.to_csv(output_dir + "/" + file_root + ".decoy.txt", header=True, index = False, sep = '\t')
    else:
        df.to_csv(output_dir + "/" + file_root + ".target.txt", header=True, index = False, sep = '\t')
        if return_decoys:
                df_decoys.to_csv(output_dir + "/" + file_root + ".decoy.txt", header=True, index = False, sep = '\t')
        
    if return_frontier:
        results[1].to_csv(output_dir + "/" + file_root + ".frontier.txt", header=True, index = False, sep = '\t')
    
    if return_mass_mod_hist:
        #Producing histogram for unaccounted for modifications
        df_mass_mods = df[df['search_file'] == 'open'].copy()
        if len(df_mass_mods) > 0:
            logging.info("Creating histograms for unaccounted mass-modifications.")
            sys.stderr.write("Creating histograms for unaccounted mass-modifications. \n")
            delta_mass_min = np.floor(min(df_mass_mods['delta_mass']))
            delta_mass_max = np.ceil(max(df_mass_mods['delta_mass']))
            plt.hist(df_mass_mods['delta_mass'], bins = np.arange(delta_mass_min, delta_mass_max + 2, 1) - 0.5, )
            plt.title('Frequency of unaccounted mass-modifications')
            plt.xlabel('Mass-modifications (Da)')
            plt.ylabel('Frequency')
            plt.savefig(output_dir + "/" + file_root + '.unaccounted-mass-mods.pdf')
        else:
            logging.info("No unaccounted mass-modifications.")
            sys.stderr.write("No unaccounted mass-modifications. \n")
        
        #Producing histogram for accounted for modifications
        variable_mods = df['peptide'].apply(lambda x: re.findall('\\[(.+?)\\]', x))
        variable_mods_inds = variable_mods.apply(lambda x: len(x) > 0)
        variable_mods = variable_mods[variable_mods_inds]
        variable_mods = variable_mods.sort_values().apply(lambda x: sorted(x))
        if len(variable_mods) > 0:
            logging.info("Counting the number of accounted-for variable mass-modifications.")
            sys.stderr.write("Counting the number of accounted-for variable mass-modifications. \n")
            variable_mods_table = pd.DataFrame(variable_mods.value_counts())
            variable_mods_table.columns = ['Count']
            variable_mods_table.index.names = ['Modifications:']
            logging.info(variable_mods_table.to_string())
            sys.stderr.write(variable_mods_table.to_string() + "\n")
        else:
            logging.info("No variable mass-modifications.")
            sys.stderr.write("No variable mass-modifications. \n")
        
    end_time = time.time()
    
    logging.info("Elapsed time: " + str(round(end_time - start_time, 2)) + " s")
    sys.stderr.write("Elapsed time: " + str(round(end_time - start_time, 2)) + " s \n")
    
if __name__ == "__main__":
    main()    

