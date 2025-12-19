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
