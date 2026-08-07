#! /usr/bin/env python

import sys
import os
import logging

_LOG = logging.getLogger(__name__)

def get_overlap_and_diff(seq1, seq2, missing_symbols = ("?", "-", "N", "n")):
    align_length = len(seq1)
    if align_length != len(seq2):
        raise Exception("Sequences are not aligned")
    if align_length < 1:
        return 0.0, None
    num_shared_sites = 0
    num_shared_site_diffs = 0
    for i in range(align_length):
        if (seq1[i] in missing_symbols) or (seq2[i] in missing_symbols):
            continue
        num_shared_sites += 1
        if seq1[i].upper() != seq2[i].upper():
            num_shared_site_diffs += 1
    prop_shared = num_shared_sites / float(align_length)
    prop_shared_diffs = None
    if num_shared_sites > 0:
        prop_shared_diffs = num_shared_site_diffs / float(num_shared_sites)
    return prop_shared, prop_shared_diffs

def get_shared_indices(seqs, symbol_set):
    gap_chars = {'?', '-'}
    try:
        indices = [i for i, col in enumerate(zip(*seqs, strict = True)) if all(c in symbol_set for c in col)]
    except ValueError as e:
        sys.stderr.write("get_shared_indices requires aligned sequences")
        raise e
    return indices

def get_missing_column_indices(seqs, missing_symbols = {'?', '-'}):
    return get_shared_indices(seqs, missing_symbols)

def remove_missing_columns(labeled_seqs, missing_symbols = {'?', '-'}):
    seq_len = len(labeled_seqs[0][1])
    seqs = (seq for lab, seq in labeled_seqs)
    remove_set = set(get_missing_column_indices(seqs, missing_symbols))
    indices_to_keep = [i for i in range(seq_len) if i not in remove_set]
    new_seqs = [[lab, [seq[i] for i in indices_to_keep]] for lab, seq in labeled_seqs]
    return new_seqs
