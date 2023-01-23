#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 21 22:38:47 2023

@author: jackfreestone

Testing primarily utility/helper functions used in CONGA.py. Overall testing is accomplished by ensuring
that the CONDA modules completes successfully using a few example data sets via Github actions.

"""

from CONGA.utils import CONGA_functions as cg
import random
import pandas as pd
import numpy as np

def test_reverse_sequence():
    assert cg.reverse_sequence('ABCD', '') == 'CBAD'
    assert cg.reverse_sequence('ABCD[10.1]', '4_V_10.1') == 'CBAD[10.1]'
    assert cg.reverse_sequence('A[9.2]BCD[10.1]', '1_V_9.2, 4_V_10.1') == 'CBA[9.2]D[10.1]'
    assert cg.reverse_sequence('A[9.2][8.3]BCD[10.1]', '1_V_9.2, 1_V_8.3_N, 4_V_10.1') == 'C[8.3]BA[9.2]D[10.1]'
    
def test_parse_static_mods():
    assert cg.parse_static_mods('S:+79.345') == {'S': 79.345, 'C': 57.02146}
    assert cg.parse_static_mods('S:+79.345,nterm:30.987,cterm:-53.135') == {'S': 79.345, 'nterm': 30.987, 'cterm': -53.135, 'C': 57.02146}

def test_parse_mods():
    assert cg.parse_mods('A[1.2345]BCD', {'S': 79.345}) == (['A', 'B', 'C', 'D'], [1.2345, 0.0, 0.0, 0.0], 0.0, 0.0)
    assert cg.parse_mods('A[1.2345][1.2345]BCD', {'S': 79.345}) == (['A', 'B', 'C', 'D'], [2.469, 0.0, 0.0, 0.0], 0.0, 0.0)
    assert cg.parse_mods('A[1.2345][1.2345]BCD', {'S': 79.345, 'A': 1.2345}) ==  (['A', 'B', 'C', 'D'], [3.7035, 0.0, 0.0, 0.0], 0.0, 0.0)
    
def test_get_similarity():
    assert cg.get_similarity(['PEPTIDE'], [2], ['PEPTIDE'], [2]) == (['PEPTIDE'], ['PEPTIDE'], [1.0])
    assert cg.get_similarity(['PEPTI'], [2], ['LEPTI'], [2]) == (['PEPTI'], ['LEPTI'], [0.5])
    
def TDC_flex_c(winning_scores, labels):
    ordered_inds = sorted(range(len(winning_scores)), key=lambda k: winning_scores[k], reverse = True)
    ordered_labels = [labels[i] for i in ordered_inds]
    nTD = np.cumsum(ordered_labels)
    nDD = np.cumsum([(not l) for l in ordered_labels])
    nTD = [max(1, t) for t in nTD]
    fdps = (1 + nDD)/nTD
    fdps = pd.Series([float(min(1, fdp)) for fdp in fdps])
    fdps = fdps[::-1]
    qvals = fdps.cummin() 
    qvals = qvals[::-1]
    original_inds = sorted(range(len(winning_scores)), key=lambda k: ordered_inds[k], reverse = False)
    qvals_original_order = [qvals[i] for i in original_inds]
    return (qvals_original_order)

winning_scores = list(np.random.normal(size = 20))
labels = random.choices([True, False], k = 20)

def test_group_walk():
    results = cg.group_walk(winning_scores, labels, [1]*len(winning_scores))
    list_qval = [float(result) for result in results[0]]
    assert list_qval == TDC_flex_c(winning_scores, labels)

def test_create_peptides_with_mod():
    unmod_peptides = pd.Series(['NQVNMWEPCK', 'ASDFBDCDC'])
    mods = pd.Series(['1N(229.163), 9C(57.021465), 10C(100.12345)', ''])
    static_mods = {'C': 57.02146}
    assert cg.create_peptides_with_mod(unmod_peptides, mods, static_mods).to_list() == ['N[229.163]QVNMWEPCK[100.12345]', 'ASDFBDCDC']

def test_reverse_sequence():
    assert cg.reverse_sequence('P[10.01][1.2345]EPTIDE[11.1]', '1_V_10.011_N, 1_V_1.23450, 7_V_11.1002_C') == 'D[10.01]ITPEP[1.2345]E[11.1]'
    
def test_check_n_term():
    sequences = pd.Series(['[1.2345]ABCDEFG', '[1.2345]BCDEFG'])
    assert cg.check_n_term(sequences).to_list() == ['A[1.2345]BCDEFG', 'B[1.2345]CDEFG']
    
def test_del_protein():
    dcy_prefix = 'decoy_'
    cg.del_protein('decoy_PROTEIN1, PROTEIN2', dcy_prefix) == True
    cg.del_protein('decoy_PROTEIN1', dcy_prefix) == False
    cg.del_protein('PROTEIN3, PROTEIN4', dcy_prefix) == False
    cg.del_protein('PROTEIN5, decoy_PROTEIN6', dcy_prefix) == True

def test_get_amino_acid_to_warn():
    delta_mass = pd.Series([128.094963, -114.03])
    peptide = pd.Series(['ABCD', 'ABCD[13.1]N'])
    flanking_aa = pd.Series(['KT', 'RA'])
    df = pd.DataFrame(zip(delta_mass, peptide, flanking_aa), columns = ['delta_mass', 'peptide', 'flanking_aa'])
    df.apply(cg.get_amino_acid_to_warn, axis = 1) == pd.Series(['Possible addition of K', 'Possible loss of N'])
    
def test_get_modification_info():
    cg.get_modification_info('A[10.1]BCD[15.333]EF') == '1[10.1],4[15.333]'
    cg.get_modification_info('A[10.1][1.2345]BCD[15.333]EF') == '1[10.1,1.2345],4[15.333]'